package overlay

import (
	"bytes"
	"context"
	"errors"
	"fmt"
	"os"
	"os/exec"
	"path/filepath"
	"strings"
	"syscall"

	"github.com/Justype/condatainer/internal/logging"
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
	ProfileSmall = Profile{InodeRatio: 4096, ReservedPerc: 3}

	// ProfileDefault is a general-purpose profile (Linux default 16K).
	ProfileDefault = Profile{InodeRatio: 16384, ReservedPerc: 3}

	// ProfileLarge is optimized for large files (genomes, FASTA, databases).
	ProfileLarge = Profile{InodeRatio: 1048576, ReservedPerc: 3}
)

// GetProfile returns the Profile for a name (small/balanced/large, or an alias),
// or nil if the name is not found.
func GetProfile(name string) *Profile {
	switch strings.ToLower(name) {
	case "small", "conda", "python":
		return &ProfileSmall
	case "large", "data", "genome":
		return &ProfileLarge
	case "balanced", "default", "":
		return &ProfileDefault
	default:
		return nil
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
	log := logging.FromContext(ctx)
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
	var debugBuf bytes.Buffer
	if w := logging.WriterFromCtx(ctx); w != nil {
		cmd.Stdout = w
		cmd.Stderr = w
	} else {
		cmd.Stdout = &debugBuf
		cmd.Stderr = &debugBuf
	}
	if execErr := cmd.Run(); execErr != nil {
		cleanup()
		return &Error{
			Op:      "inject structure",
			Path:    opts.Path,
			Tool:    "debugfs",
			Output:  debugBuf.String(),
			BaseErr: execErr,
		}
	}

	// 4. Set permissions
	if err := os.Chmod(filePath, utils.PermFile); err != nil {
		log.Debug(fmt.Sprintf("failed to set permissions on overlay: %v", err))
	}

	return nil
}

// crossFsCopy copies src to dst using the system cp command, which:
//   - Is a subprocess → context-cancellable (Ctrl+C works during large copies)
//   - sparse=true:  cp --sparse=always — detects holes, destination stays sparse
//   - sparse=false: cp --sparse=never  — writes all zeros, destination is fully allocated
//
// Uses cmd.Start + a Wait goroutine instead of CombinedOutput so that cancellation
// returns immediately. CombinedOutput blocks in cmd.Wait even after SIGKILL because
// on network filesystems (Lustre/NFS) the cp process can enter D-state
// (uninterruptible I/O sleep) and not honour SIGKILL until the I/O resolves.
// The background Wait goroutine will eventually clean up once the process exits.
func crossFsCopy(ctx context.Context, src, dst string, sparse bool) error {
	sparseFlag := "--sparse=never"
	if sparse {
		sparseFlag = "--sparse=always"
	}
	var outBuf bytes.Buffer
	cmd := exec.Command("cp", sparseFlag, src, dst)
	cmd.Stdout = &outBuf
	cmd.Stderr = &outBuf
	if err := cmd.Start(); err != nil {
		return fmt.Errorf("cp %s → %s: %w", src, dst, err)
	}

	waitErr := make(chan error, 1)
	go func() { waitErr <- cmd.Wait() }()

	select {
	case err := <-waitErr:
		if err != nil {
			os.Remove(dst) //nolint:errcheck
			return fmt.Errorf("cp %s → %s: %w\n%s", src, dst, err, strings.TrimSpace(outBuf.String()))
		}
	case <-ctx.Done():
		cmd.Process.Kill() //nolint:errcheck
		os.Remove(dst)     //nolint:errcheck
		return ctx.Err()
	}

	if err := os.Chmod(dst, utils.PermFile); err != nil {
		logging.FromContext(ctx).Debug(fmt.Sprintf("failed to set permissions on overlay: %v", err))
	}
	return nil
}

// moveFile moves src to dst. It tries os.Rename first (same-filesystem, instant).
// If that fails due to a cross-device link (e.g. local /tmp → LustreFS), it falls
// back to crossFsCopy which uses cp and is context-cancellable (Ctrl+C works).
// Returns (copied=true) when a copy was performed, (copied=false) when renamed.
func moveFile(ctx context.Context, src, dst string, sparse bool) (copied bool, err error) {
	if err := os.MkdirAll(filepath.Dir(dst), utils.PermDir); err != nil {
		return false, fmt.Errorf("create destination directory: %w", err)
	}

	if err := os.Rename(src, dst); err == nil {
		return false, nil
	} else if !errors.Is(err, syscall.EXDEV) {
		return false, fmt.Errorf("rename %s → %s: %w", src, dst, err)
	}

	if err := crossFsCopy(ctx, src, dst, sparse); err != nil {
		return false, err
	}
	os.Remove(src)
	return true, nil
}

// ---------------------------------------------------------
// 3. Public API
// ---------------------------------------------------------

// MoveOverlayCopied moves src to dst and reports whether a copy was performed.
func MoveOverlayCopied(ctx context.Context, src, dst string, sparse bool) (copied bool, err error) {
	return moveFile(ctx, src, dst, sparse)
}

// AllocateOverlay pre-allocates disk blocks for a sparse overlay file using fallocate.
func AllocateOverlay(ctx context.Context, path string, sizeMB int) {
	log := logging.FromContext(ctx)
	log.Info(fmt.Sprintf("allocating %d MiB at %s", sizeMB, path))
	cmd := exec.CommandContext(ctx, "fallocate", "-l", fmt.Sprintf("%dM", sizeMB), path)
	if output, err := cmd.CombinedOutput(); err != nil {
		log.Warn(fmt.Sprintf("fallocate failed, overlay may be sparse: %s", strings.TrimSpace(string(output))))
	}
}

// CreateInTmp builds the overlay at a local temp path and returns that path.
// The overlay is always created sparse to minimize local disk usage.
// The caller is responsible for running additional init, moving to final destination,
// and calling AllocateOverlay if needed. The caller must remove the temp file on error.
func CreateInTmp(ctx context.Context, opts *CreateOptions) (tmpPath string, err error) {
	tmpPath, err = createAtTmp(ctx, opts)
	return
}

// createAtTmp builds the overlay at a local temp path under utils.GetTmpDir() and returns
// that path. Always created sparse to minimize local disk usage.
func createAtTmp(ctx context.Context, opts *CreateOptions) (tmpPath string, err error) {
	if err := validateOpts(opts); err != nil {
		return "", err
	}

	if err := checkDependencies([]string{"dd", "mke2fs", "debugfs"}); err != nil {
		return "", err
	}

	label := typeLabel(opts)
	log := logging.FromContext(ctx)

	tmpDir := utils.GetTmpDir()
	if err := utils.EnsureTmpSubdir(tmpDir); err != nil {
		return "", fmt.Errorf("failed to create tmp dir %s: %w", tmpDir, err)
	}
	tmpPath = filepath.Join(tmpDir, filepath.Base(opts.Path))

	if !opts.Quiet {
		log.Info(fmt.Sprintf("creating %soverlay %s", label, opts.Path))
		log.Info(fmt.Sprintf("size %d MiB | filesystem %s | inode ratio %d | reserved %d%%",
			opts.SizeMB, opts.FilesystemType, opts.Profile.InodeRatio, opts.Profile.ReservedPerc))
		log.Info(fmt.Sprintf("building at local tmp: %s", tmpPath))
	}

	if err := createOverlayFile(ctx, opts, tmpPath, true); err != nil {
		utils.RemoveDirIfEmpty(tmpDir)
		return "", err
	}

	return tmpPath, nil
}

// CreateDirectly builds the overlay at opts.Path without using a local tmp directory.
func CreateDirectly(ctx context.Context, opts *CreateOptions) error {
	if err := validateOpts(opts); err != nil {
		return err
	}

	if err := checkDependencies([]string{"dd", "mke2fs", "debugfs"}); err != nil {
		return err
	}

	if err := os.MkdirAll(filepath.Dir(opts.Path), utils.PermDir); err != nil {
		return fmt.Errorf("create destination directory: %w", err)
	}

	label := typeLabel(opts)
	log := logging.FromContext(ctx)

	if !opts.Quiet {
		log.Info(fmt.Sprintf("creating %soverlay %s", label, opts.Path))
		log.Info(fmt.Sprintf("size %d MiB | filesystem %s | inode ratio %d | reserved %d%%",
			opts.SizeMB, opts.FilesystemType, opts.Profile.InodeRatio, opts.Profile.ReservedPerc))
	}

	if err := createOverlayFile(ctx, opts, opts.Path, opts.Sparse); err != nil {
		return err
	}

	log.Info(fmt.Sprintf("created %soverlay %s", label, opts.Path), "kind", "success")
	return nil
}

// CreateWithOptions builds the overlay at a local tmp path and moves it to opts.Path.
// If opts.Sparse is false, AllocateOverlay is called after the move to fill sparse holes.
// For cases where additional work must happen before the final move (e.g. conda init),
// use CreateInTmp + MoveOverlayCopied + AllocateOverlay instead.
func CreateWithOptions(ctx context.Context, opts *CreateOptions) error {
	tmpPath, err := createAtTmp(ctx, opts)
	if err != nil {
		return err
	}

	label := typeLabel(opts)
	log := logging.FromContext(ctx)

	if tmpPath != opts.Path {
		defer func() {
			os.Remove(tmpPath)
			utils.RemoveDirIfEmpty(filepath.Dir(tmpPath))
		}()

		if !opts.Quiet {
			log.Info(fmt.Sprintf("moving overlay to %s", opts.Path))
		}
		copied, err := moveFile(ctx, tmpPath, opts.Path, opts.Sparse)
		if err != nil {
			return fmt.Errorf("failed to install overlay: %w", err)
		}

		if !opts.Sparse && !copied {
			AllocateOverlay(ctx, opts.Path, opts.SizeMB)
		}
	}

	log.Info(fmt.Sprintf("created %soverlay %s", label, opts.Path), "kind", "success")
	return nil
}
