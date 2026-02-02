package overlay

import (
	"context"
	"fmt"
	"os"
	"os/exec"
	"strings"

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

// getProfile returns the struct based on a simple string name.
func getProfile(name string) Profile {
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
// 2. Public API
// ---------------------------------------------------------

// CreateForCurrentUser generates an overlay owned by the current host user.
func CreateForCurrentUser(ctx context.Context, path string, sizeMB int, profileType string, sparse bool, fsType string, quiet bool) error {
	opts := &CreateOptions{
		Path:           path,
		SizeMB:         sizeMB,
		UID:            os.Getuid(),
		GID:            os.Getgid(),
		Profile:        getProfile(profileType),
		Sparse:         sparse,
		FilesystemType: fsType,
		Quiet:          quiet,
	}
	return CreateWithOptions(ctx, opts)
}

// CreateForRoot generates an overlay owned by root (UID=0, GID=0).
func CreateForRoot(ctx context.Context, path string, sizeMB int, profileType string, sparse bool, fsType string, quiet bool) error {
	opts := &CreateOptions{
		Path:           path,
		SizeMB:         sizeMB,
		UID:            0,
		GID:            0,
		Profile:        getProfile(profileType),
		Sparse:         sparse,
		FilesystemType: fsType,
		Quiet:          quiet,
	}
	return CreateWithOptions(ctx, opts)
}

// CreateWithOptions generates an overlay image using the provided options.
// Supports ext3 and ext4 filesystem types.
func CreateWithOptions(ctx context.Context, opts *CreateOptions) error {
	// 0. Validate and normalize filesystem type
	if opts.FilesystemType == "" {
		opts.FilesystemType = "ext3" // Default to ext3 for compatibility
	}
	fsType := strings.ToLower(opts.FilesystemType)
	if fsType != "ext3" && fsType != "ext4" {
		return fmt.Errorf("unsupported filesystem type '%s': must be ext3 or ext4", opts.FilesystemType)
	}
	opts.FilesystemType = fsType
	// 1. Check for required system tools
	if err := checkDependencies([]string{"dd", "mke2fs", "debugfs"}); err != nil {
		return err
	}

	typeStr := "Allocated "
	if opts.Sparse {
		typeStr = "Sparse "
	}
	if opts.UID == 0 && opts.GID == 0 {
		typeStr += "Fakeroot "
	}
	if !opts.Quiet {
		styleType := utils.StyleInfo(typeStr)
		stylePath := utils.StylePath(opts.Path)

		utils.PrintMessage("Creating %soverlay %s",
			styleType, stylePath)
		utils.PrintMessage("Size %s MB | Filesystem %s | Profile Inode Ratio %s Reserved %s%%",
			utils.StyleNumber(opts.SizeMB), utils.StyleInfo(opts.FilesystemType),
			utils.StyleNumber(opts.Profile.InodeRatio),
			utils.StyleNumber(opts.Profile.ReservedPerc))
	}

	// 2. Create File (dd)
	// Strategy switch based on 'sparse' flag
	ddArgs := []string{"if=/dev/zero", "of=" + opts.Path, "bs=1M"}

	if opts.Sparse {
		// Sparse: Seek to end, write nothing. Instant.
		ddArgs = append(ddArgs, "count=0", fmt.Sprintf("seek=%d", opts.SizeMB))
	} else {
		// Allocated: Write zeros. Slow but guarantees space.
		ddArgs = append(ddArgs, fmt.Sprintf("count=%d", opts.SizeMB), "status=progress")
	}

	if err := runCommand(ctx, "create_file", opts.Path, "dd", ddArgs...); err != nil {
		// Cleanup the file if dd was interrupted or failed (especially for non-sparse)
		_ = os.Remove(opts.Path)
		return err
	}

	// 3. Format filesystem with Profile (mke2fs)
	if err := runCommand(ctx, "format", opts.Path, "mke2fs",
		"-t", opts.FilesystemType,
		"-i", fmt.Sprintf("%d", opts.Profile.InodeRatio),
		"-m", fmt.Sprintf("%d", opts.Profile.ReservedPerc),
		"-F", opts.Path); err != nil {
		_ = os.Remove(opts.Path) // Cleanup on fail
		return err
	}

	// 4. Inject Structure (debugfs)
	var script strings.Builder
	script.WriteString("mkdir upper\n")
	script.WriteString("mkdir work\n")
	script.WriteString(fmt.Sprintf("set_inode_field upper uid %d\n", opts.UID))
	script.WriteString(fmt.Sprintf("set_inode_field upper gid %d\n", opts.GID))
	script.WriteString(fmt.Sprintf("set_inode_field work uid %d\n", opts.UID))
	script.WriteString(fmt.Sprintf("set_inode_field work gid %d\n", opts.GID))
	script.WriteString("quit\n")

	cmd := exec.CommandContext(ctx, "debugfs", "-w", opts.Path)
	cmd.Stdin = strings.NewReader(script.String())

	output, err := cmd.CombinedOutput()
	if err != nil {
		// If structure injection fails (e.g. cancelled), the image is corrupt/incomplete.
		_ = os.Remove(opts.Path)

		return &Error{
			Op:      "inject structure",
			Path:    opts.Path,
			Tool:    "debugfs",
			Output:  string(output),
			BaseErr: err,
		}
	}

	// Set permissions on overlay image (664 for group-writable)
	if err := os.Chmod(opts.Path, 0o664); err != nil {
		utils.PrintDebug("Failed to set permissions on %s: %v", opts.Path, err)
	}

	utils.PrintSuccess("Created %soverlay %s", utils.StyleInfo(typeStr), utils.StylePath(opts.Path))
	return nil
}
