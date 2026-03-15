package overlay

import (
	"context"
	"errors"
	"fmt"
	"io"
	"os"
	"os/exec"
	"path/filepath"
	"strings"
	"syscall"

	"github.com/Justype/condatainer/internal/utils"
)

// ---------------------------------------------------------
// 1. Filesystem Profiles and Options
// ---------------------------------------------------------

// Profile defines the tuning parameters for the filesystem.
type Profile struct {
	InodeRatio   int // -i: bytes per inode
	ReservedPerc int // -m: percentage (0-100)
}

// CreateOptions holds all configuration options for creating an overlay image.
type CreateOptions struct {
	Path           string  // Path to the overlay image file
	SizeMB         int     // Size in megabytes
	UID            int     // Owner user ID
	GID            int     // Owner group ID
	Profile        Profile // Filesystem tuning profile
	Sparse         bool    // Create sparse file (vs allocated)
	FilesystemType string  // Filesystem type: "ext3" or "ext4" (default: "ext3")
	Quiet          bool    // Suppress detailed specs output
}

var (
	// ProfileSmall is tuned for environments with many small files (Conda packages, Python/site-packages).
	// 1 inode per 4KB keeps metadata lean so you don't hit inode exhaustion inside Conda envs.
	ProfileSmall = Profile{InodeRatio: 4096, ReservedPerc: 3}

	// ProfileDefault is a general-purpose profile (Linux default 16K).
	// 1 inode per 16KB balances metadata and data for typical workloads.
	ProfileDefault = Profile{InodeRatio: 16384, ReservedPerc: 3}
	// ProfileLarge is optimized for large files (genomes, FASTA, databases).
	// 1 inode per 1MB drastically lowers metadata pressure for huge single files.
	ProfileLarge = Profile{InodeRatio: 1048576, ReservedPerc: 3}
)

// GetProfile returns the Profile struct for a simple string name.
func GetProfile(name string) Profile {
	switch strings.ToLower(name) {
	case "small", "conda", "python":
		return ProfileSmall
	case "large", "data", "genome":
		return ProfileLarge
	default:
		return ProfileDefault
	}
}

// ---------------------------------------------------------
// 2. Internal helpers
// ---------------------------------------------------------

// validateOpts normalizes and validates CreateOptions in-place.
func validateOpts(opts *CreateOptions) error {
	if opts.FilesystemType == "" {
		opts.FilesystemType = "ext3"
	}
	fsType := strings.ToLower(opts.FilesystemType)
	if fsType != "ext3" && fsType != "ext4" {
		return fmt.Errorf("unsupported filesystem type '%s': must be ext3 or ext4", opts.FilesystemType)
	}
	opts.FilesystemType = fsType
	return nil
}

// typeLabel returns a human-readable label for the overlay type (e.g. "Allocated Fakeroot ").
func typeLabel(opts *CreateOptions) string {
	label := "Allocated "
	if opts.Sparse {
		label = "Sparse "
	}
	if opts.UID == 0 && opts.GID == 0 {
		label += "Fakeroot "
	}
	return label
}

// createOverlayFile runs dd + mke2fs + debugfs to build a raw overlay at filePath.
// sparse controls whether dd creates a sparse file (fast, saves local space) or a fully-allocated one.
func createOverlayFile(ctx context.Context, opts *CreateOptions, filePath string, sparse bool) error {
	cleanup := func() { os.Remove(filePath) }

	// 1. Create raw file (dd)
	ddArgs := []string{"if=/dev/zero", "of=" + filePath, "bs=1M"}
	if sparse {
		ddArgs = append(ddArgs, "count=0", fmt.Sprintf("seek=%d", opts.SizeMB))
	} else {
		ddArgs = append(ddArgs, fmt.Sprintf("count=%d", opts.SizeMB), "status=progress")
	}
	if err := runCommand(ctx, "create_file", opts.Path, "dd", ddArgs...); err != nil {
		cleanup()
		return err
	}

	// 2. Format filesystem (mke2fs)
	if err := runCommand(ctx, "format", opts.Path, "mke2fs",
		"-t", opts.FilesystemType,
		"-i", fmt.Sprintf("%d", opts.Profile.InodeRatio),
		"-m", fmt.Sprintf("%d", opts.Profile.ReservedPerc),
		"-F", filePath); err != nil {
		cleanup()
		return err
	}

	// 3. Inject overlay structure (debugfs)
	var script strings.Builder
	script.WriteString("mkdir upper\n")
	script.WriteString("mkdir work\n")
	fmt.Fprintf(&script, "set_inode_field upper uid %d\n", opts.UID)
	fmt.Fprintf(&script, "set_inode_field upper gid %d\n", opts.GID)
	fmt.Fprintf(&script, "set_inode_field work uid %d\n", opts.UID)
	fmt.Fprintf(&script, "set_inode_field work gid %d\n", opts.GID)
	script.WriteString("quit\n")

	cmd := exec.CommandContext(ctx, "debugfs", "-w", filePath)
	cmd.Stdin = strings.NewReader(script.String())
	output, execErr := cmd.CombinedOutput()
	if execErr != nil {
		cleanup()
		return &Error{
			Op:      "inject structure",
			Path:    opts.Path,
			Tool:    "debugfs",
			Output:  string(output),
			BaseErr: execErr,
		}
	}

	// 4. Set permissions
	if err := os.Chmod(filePath, utils.PermFile); err != nil {
		utils.PrintDebug("Failed to set permissions on overlay: %v", err)
	}

	return nil
}

// Linux lseek whence values for sparse file support (available since kernel 3.1).
const (
	seekData = 3 // SEEK_DATA: seek to next data region
	seekHole = 4 // SEEK_HOLE: seek to next hole
)

// sparseAwareCopy copies src → dst preserving sparse holes.
// It uses SEEK_DATA/SEEK_HOLE to skip over holes, writing only data segments.
// The destination file is truncated to match the source size so holes are implicit.
func sparseAwareCopy(src, dst *os.File) error {
	srcFd := int(src.Fd())

	// Get total size for the final truncate.
	info, err := src.Stat()
	if err != nil {
		return err
	}
	size := info.Size()

	var offset int64
	buf := make([]byte, 1<<20) // 1 MiB copy buffer
	for {
		// Find next data segment starting at offset.
		dataStart, err := syscall.Seek(srcFd, offset, seekData)
		if err != nil {
			break // ENXIO means no more data segments
		}
		holeStart, err := syscall.Seek(srcFd, dataStart, seekHole)
		if err != nil {
			holeStart = size // data runs to EOF
		}

		// Seek both files to the data region.
		if _, err := src.Seek(dataStart, io.SeekStart); err != nil {
			return err
		}
		if _, err := dst.Seek(dataStart, io.SeekStart); err != nil {
			return err
		}

		// Copy the data segment.
		remaining := holeStart - dataStart
		for remaining > 0 {
			n := min(int64(len(buf)), remaining)
			nr, err := src.Read(buf[:n])
			if nr > 0 {
				if _, werr := dst.Write(buf[:nr]); werr != nil {
					return werr
				}
				remaining -= int64(nr)
			}
			if err == io.EOF {
				break
			} else if err != nil {
				return fmt.Errorf("read from tmp overlay: %w", err)
			}
		}
		offset = holeStart
	}

	// Truncate dst to exact source size so trailing holes are implicit.
	return dst.Truncate(size)
}

// moveFile moves src to dst. It tries os.Rename first (same-filesystem, instant).
// If that fails due to a cross-device link (e.g. local /tmp → LustreFS), it falls
// back to a copy followed by removal of src.
// When sparse=true it uses a hole-aware copy so the destination stays sparse.
// When sparse=false it uses io.Copy which materialises all holes (fully allocated).
// Returns (copied=true) when a copy was performed, (copied=false) when renamed.
func moveFile(src, dst string, sparse bool) (copied bool, err error) {
	if err := os.Rename(src, dst); err == nil {
		return false, nil
	} else if !errors.Is(err, syscall.EXDEV) {
		return false, fmt.Errorf("rename %s → %s: %w", src, dst, err)
	}

	// Cross-filesystem: copy then remove src.
	in, err := os.Open(src)
	if err != nil {
		return false, fmt.Errorf("open tmp overlay: %w", err)
	}
	defer in.Close()

	out, err := os.OpenFile(dst, os.O_WRONLY|os.O_CREATE|os.O_TRUNC, utils.PermFile)
	if err != nil {
		return false, fmt.Errorf("create destination overlay: %w", err)
	}

	var copyErr error
	if sparse {
		copyErr = sparseAwareCopy(in, out)
	} else {
		_, copyErr = io.Copy(out, in)
	}
	if copyErr != nil {
		out.Close()
		os.Remove(dst)
		return false, fmt.Errorf("copy overlay to destination: %w", copyErr)
	}
	if err := out.Close(); err != nil {
		os.Remove(dst)
		return false, fmt.Errorf("flush destination overlay: %w", err)
	}

	os.Remove(src)
	return true, nil
}

// ---------------------------------------------------------
// 3. Public API
// ---------------------------------------------------------

// MoveOverlayCopied moves src to dst and reports whether a copy was performed.
// sparse=true preserves holes (space-efficient, slightly slower);
// sparse=false materialises all zeros via io.Copy (fully allocated, faster sequential write).
// Returns copied=true when a cross-filesystem copy was done, copied=false when os.Rename was used.
func MoveOverlayCopied(src, dst string, sparse bool) (copied bool, err error) {
	return moveFile(src, dst, sparse)
}

// AllocateOverlay pre-allocates disk blocks for a sparse overlay file using fallocate.
// After a cross-filesystem move the file is already allocated (io.Copy fills holes);
// after a same-filesystem rename it may still be sparse — fallocate handles both cases efficiently.
// A warning is printed if fallocate is unavailable; the overlay remains functional but sparse.
func AllocateOverlay(ctx context.Context, path string, sizeMB int) {
	utils.PrintMessage("Allocating %s MB at %s", utils.StyleNumber(sizeMB), utils.StylePath(path))
	cmd := exec.CommandContext(ctx, "fallocate", "-l", fmt.Sprintf("%dM", sizeMB), path)
	if output, err := cmd.CombinedOutput(); err != nil {
		utils.PrintWarning("fallocate failed, overlay may be sparse: %s", strings.TrimSpace(string(output)))
	}
}

// CreateInTmp builds the overlay at a local temp path and returns that path.
// The overlay is always created sparse to minimize local disk usage.
// The caller is responsible for:
//   - Running any additional initialization (e.g. conda init) on the returned tmp path.
//   - Calling MoveOverlay to move the file to its final destination.
//   - Calling AllocateOverlay if the final overlay should be fully allocated (not sparse).
//   - Removing the temp file on error.
func CreateInTmp(ctx context.Context, opts *CreateOptions) (tmpPath string, err error) {
	tmpPath, err = createAtTmp(ctx, opts)
	return
}

// createAtTmp builds the overlay at a local temp path under utils.GetTmpDir() and returns
// that path. The overlay is always created sparse to minimize local disk usage.
// The caller owns the file and must move it or remove it.
func createAtTmp(ctx context.Context, opts *CreateOptions) (tmpPath string, err error) {
	if err := validateOpts(opts); err != nil {
		return "", err
	}

	if err := checkDependencies([]string{"dd", "mke2fs", "debugfs"}); err != nil {
		return "", err
	}

	label := typeLabel(opts)

	// Prepare local tmp path
	tmpDir := utils.GetTmpDir()
	if err := os.MkdirAll(tmpDir, utils.PermDir); err != nil {
		return "", fmt.Errorf("failed to create tmp dir %s: %w", tmpDir, err)
	}
	tmpPath = filepath.Join(tmpDir, filepath.Base(opts.Path))

	if !opts.Quiet {
		utils.PrintMessage("Creating %soverlay %s", utils.StyleInfo(label), utils.StylePath(opts.Path))
		utils.PrintMessage("Size %s MB | Filesystem %s | Profile Inode Ratio %s Reserved %s%%",
			utils.StyleNumber(opts.SizeMB), utils.StyleInfo(opts.FilesystemType),
			utils.StyleNumber(opts.Profile.InodeRatio),
			utils.StyleNumber(opts.Profile.ReservedPerc))
		utils.PrintMessage("Building at local tmp: %s", utils.StylePath(tmpPath))
	}

	// Always create sparse in tmp to save local disk space.
	// If the final overlay must be allocated, call AllocateOverlay after MoveOverlay.
	if err := createOverlayFile(ctx, opts, tmpPath, true); err != nil {
		return "", err
	}

	return tmpPath, nil
}

// CreateDirectly builds the overlay at opts.Path without using a local tmp directory.
// The overlay is created with opts.Sparse as-is (no deferred allocation step).
// Use this when the target filesystem is fast (e.g. local SSD via --no-tmp).
func CreateDirectly(ctx context.Context, opts *CreateOptions) error {
	if err := validateOpts(opts); err != nil {
		return err
	}

	if err := checkDependencies([]string{"dd", "mke2fs", "debugfs"}); err != nil {
		return err
	}

	label := typeLabel(opts)

	if !opts.Quiet {
		utils.PrintMessage("Creating %soverlay %s", utils.StyleInfo(label), utils.StylePath(opts.Path))
		utils.PrintMessage("Size %s MB | Filesystem %s | Profile Inode Ratio %s Reserved %s%%",
			utils.StyleNumber(opts.SizeMB), utils.StyleInfo(opts.FilesystemType),
			utils.StyleNumber(opts.Profile.InodeRatio),
			utils.StyleNumber(opts.Profile.ReservedPerc))
	}

	if err := createOverlayFile(ctx, opts, opts.Path, opts.Sparse); err != nil {
		return err
	}

	utils.PrintSuccess("Created %soverlay %s", utils.StyleInfo(label), utils.StylePath(opts.Path))
	return nil
}

// CreateWithOptions builds the overlay at a local tmp path and moves it to opts.Path.
// If opts.Sparse is false, AllocateOverlay is called after the move to fill sparse holes.
// For cases where additional work must happen before the final move (e.g. conda init),
// use CreateInTmp + MoveOverlay + AllocateOverlay instead.
func CreateWithOptions(ctx context.Context, opts *CreateOptions) error {
	tmpPath, err := createAtTmp(ctx, opts)
	if err != nil {
		return err
	}
	defer os.Remove(tmpPath)

	label := typeLabel(opts)

	if !opts.Quiet {
		utils.PrintMessage("Moving overlay to %s", utils.StylePath(opts.Path))
	}
	copied, err := moveFile(tmpPath, opts.Path, opts.Sparse)
	if err != nil {
		return fmt.Errorf("failed to install overlay: %w", err)
	}

	// Allocate space only if not sparse and io.Copy was not used.
	// When io.Copy was used (cross-filesystem), all zero-bytes were already written physically
	// so the destination is fully allocated — fallocate is unnecessary and may not be supported.
	if !opts.Sparse && !copied {
		AllocateOverlay(ctx, opts.Path, opts.SizeMB)
	}

	utils.PrintSuccess("Created %soverlay %s", utils.StyleInfo(label), utils.StylePath(opts.Path))
	return nil
}
