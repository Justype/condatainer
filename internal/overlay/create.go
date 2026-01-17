package overlay

import (
	"fmt"
	"os"
	"os/exec"
	"strings"

	"github.com/Justype/condatainer/internal/utils"
)

// ---------------------------------------------------------
// 1. Filesystem Profiles (Tuning)
// ---------------------------------------------------------

// Profile defines the tuning parameters for the ext3 filesystem.
type Profile struct {
	BlockSize    int // -b: usually 4096
	InodeRatio   int // -i: bytes per inode
	ReservedPerc int // -m: percentage (0-100)
}

var (
	// ProfileConda is optimized for many small files (Python, Node, build artifacts).
	// Ratio 4096 = 1 inode per 4KB. Ensures you don't run out of inodes.
	ProfileConda = Profile{BlockSize: 4096, InodeRatio: 4096, ReservedPerc: 0}

	// ProfileDefault is the standard Linux balance (1 inode per 16KB).
	// Good for general purpose usage.
	ProfileDefault = Profile{BlockSize: 4096, InodeRatio: 16384, ReservedPerc: 0}

	// ProfileData is optimized for large files (Genomics, Video, Database).
	// Ratio 1MB = 1 inode per 1MB. Saves significant metadata space on large disks.
	ProfileData = Profile{BlockSize: 4096, InodeRatio: 1048576, ReservedPerc: 0}
)

// getProfile returns the struct based on a simple string name.
func getProfile(name string) Profile {
	switch strings.ToLower(name) {
	case "conda", "small", "python":
		return ProfileConda
	case "data", "genome", "large":
		return ProfileData
	default:
		return ProfileDefault
	}
}

// ---------------------------------------------------------
// 2. Public API
// ---------------------------------------------------------
// CreateForCurrentUser generates an overlay owned by the current host user.
func CreateForCurrentUser(path string, sizeMB int, profileType string, sparse bool) error {
	uid := os.Getuid()
	gid := os.Getgid()
	return Create(path, sizeMB, uid, gid, getProfile(profileType), sparse)
}

// CreateForRoot generates an overlay owned by root (UID=0, GID=0).
func CreateForRoot(path string, sizeMB int, profileType string, sparse bool) error {
	return Create(path, sizeMB, 0, 0, getProfile(profileType), sparse)
}

// Create generates an ext3 overlay image using specific tuning parameters.
func Create(path string, sizeMB int, uid, gid int, profile Profile, sparse bool) error {
	// 0. Check for required system tools
	if err := checkDependencies([]string{"dd", "mke2fs", "debugfs"}); err != nil {
		return err
	}

	typeStr := "Allocated"
	if sparse {
		typeStr = "Sparse"
	}

	utils.PrintDebug("Creating %s overlay: %s (%s MB) [Profile: -b %d -i %d -m %d]",
		typeStr,
		utils.StylePath(path), utils.StyleNumber(sizeMB),
		profile.BlockSize, profile.InodeRatio, profile.ReservedPerc)

	// 1. Create File (dd)
	// Strategy switch based on 'sparse' flag
	ddArgs := []string{"if=/dev/zero", "of=" + path, "bs=1M"}

	if sparse {
		// Sparse: Seek to end, write nothing. Instant.
		ddArgs = append(ddArgs, "count=0", fmt.Sprintf("seek=%d", sizeMB))
	} else {
		// Allocated: Write zeros. Slow but guarantees space.
		ddArgs = append(ddArgs, fmt.Sprintf("count=%d", sizeMB), "status=progress")
	}

	if err := runCommand("create_file", path, "dd", ddArgs...); err != nil {
		return err
	}

	// 2. Format Ext3 with Profile (mke2fs)
	utils.PrintDebug("Formatting ext3 with optimized settings...")

	if err := runCommand("format", path, "mke2fs",
		"-t", "ext3",
		"-b", fmt.Sprintf("%d", profile.BlockSize),
		"-i", fmt.Sprintf("%d", profile.InodeRatio),
		"-m", fmt.Sprintf("%d", profile.ReservedPerc),
		"-F", path); err != nil {
		_ = os.Remove(path) // Cleanup on fail
		return err
	}

	// 3. Inject Structure (debugfs)
	var script strings.Builder
	script.WriteString("mkdir upper\n")
	script.WriteString("mkdir work\n")
	script.WriteString(fmt.Sprintf("set_inode_field upper uid %d\n", uid))
	script.WriteString(fmt.Sprintf("set_inode_field upper gid %d\n", gid))
	script.WriteString(fmt.Sprintf("set_inode_field work uid %d\n", uid))
	script.WriteString(fmt.Sprintf("set_inode_field work gid %d\n", gid))
	script.WriteString("quit\n")

	utils.PrintDebug("Injecting overlay structure via debugfs...")

	cmd := exec.Command("debugfs", "-w", path)
	cmd.Stdin = strings.NewReader(script.String())

	output, err := cmd.CombinedOutput()
	if err != nil {
		return &Error{
			Op:      "inject structure",
			Path:    path,
			Tool:    "debugfs",
			Output:  string(output),
			BaseErr: err,
		}
	}

	utils.PrintSuccess("Created %s overlay: %s", typeStr, utils.StylePath(path))
	return nil
}
