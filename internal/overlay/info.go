package overlay

import (
	"os"
	"os/exec"
	"strconv"
	"strings"
	"syscall"
)

// Stats holds metadata about the overlay filesystem.
type Stats struct {
	// Filesystem metadata
	FilesystemType  string // "ext3", "ext4"
	BlockSize       int64
	TotalBlocks     int64
	FreeBlocks      int64
	ReservedBlocks  int64
	TotalInodes     int64
	FreeInodes      int64
	FilesystemState string // "clean", "dirty", etc.
	LastMounted     string
	CreatedTime     string
	FilesystemUUID  string
	JournalSize     int64 // in blocks
	InodeSize       int64

	// File metadata
	FileSizeBytes  int64 // Apparent size of the file
	FileBlocksUsed int64 // Actual blocks allocated on disk (for sparse detection)
	IsSparse       bool  // True if file is sparse

	// Overlay ownership
	UpperUID int
	UpperGID int
}

// GetStats reads the filesystem superblock to calculate usage without mounting.
func GetStats(path string) (*Stats, error) {
	if err := checkDependencies([]string{"tune2fs"}); err != nil {
		return nil, err
	}

	// tune2fs -l lists the superblock info
	cmd := exec.Command("tune2fs", "-l", path)
	cmd.Env = append(os.Environ(), "LC_ALL=C", "LC_TIME=C")
	out, err := cmd.Output()
	if err != nil {
		return nil, &Error{Op: "read stats", Path: path, Tool: "tune2fs", BaseErr: err}
	}

	stats := &Stats{
		UpperUID: -1,
		UpperGID: -1,
	}
	lines := strings.Split(string(out), "\n")

	for _, line := range lines {
		parts := strings.SplitN(line, ":", 2)
		if len(parts) != 2 {
			continue
		}
		key := strings.TrimSpace(parts[0])
		val := strings.TrimSpace(parts[1])

		switch key {
		case "Filesystem magic number":
			if val == "0xEF53" {
				stats.FilesystemType = "ext3/ext4"
			}
		case "Filesystem features":
			if strings.Contains(val, "has_journal") {
				if strings.Contains(val, "extent") {
					stats.FilesystemType = "ext4"
				} else {
					stats.FilesystemType = "ext3"
				}
			}
		case "Block size":
			stats.BlockSize, _ = strconv.ParseInt(val, 10, 64)
		case "Block count":
			stats.TotalBlocks, _ = strconv.ParseInt(val, 10, 64)
		case "Reserved block count":
			stats.ReservedBlocks, _ = strconv.ParseInt(val, 10, 64)
		case "Free blocks":
			stats.FreeBlocks, _ = strconv.ParseInt(val, 10, 64)
		case "Inode count":
			stats.TotalInodes, _ = strconv.ParseInt(val, 10, 64)
		case "Free inodes":
			stats.FreeInodes, _ = strconv.ParseInt(val, 10, 64)
		case "Inode size":
			stats.InodeSize, _ = strconv.ParseInt(val, 10, 64)
		case "Filesystem state":
			stats.FilesystemState = val
		case "Last mounted on":
			stats.LastMounted = val
		case "Filesystem created":
			stats.CreatedTime = val
		case "Filesystem UUID":
			stats.FilesystemUUID = val
		case "Journal size":
			// Parse journal size (e.g., "32M")
			val = strings.TrimSuffix(val, "M")
			if size, err := strconv.ParseInt(val, 10, 64); err == nil {
				stats.JournalSize = size * 1024 * 1024 / stats.BlockSize // Convert to blocks
			}
		}
	}

	// Get file stats for sparse detection
	if fileInfo, err := os.Stat(path); err == nil {
		stats.FileSizeBytes = fileInfo.Size()

		// Get actual blocks allocated (for sparse detection)
		if stat, ok := fileInfo.Sys().(*syscall.Stat_t); ok {
			stats.FileBlocksUsed = stat.Blocks * 512 // Blocks are in 512-byte units
			// File is sparse if allocated size is significantly less than apparent size
			stats.IsSparse = stats.FileBlocksUsed < stats.FileSizeBytes*9/10
		}
	}

	// Get upper directory ownership using debugfs
	stats.UpperUID, stats.UpperGID = getUpperOwnership(path)

	return stats, nil
}

// getUpperOwnership extracts UID/GID of the upper directory using debugfs.
func getUpperOwnership(path string) (uid, gid int) {
	uid, gid = -1, -1

	debugfsPath, err := exec.LookPath("debugfs")
	if err != nil {
		return
	}

	cmd := exec.Command(debugfsPath, "-R", "stat upper", path)
	output, err := cmd.CombinedOutput()
	if err != nil {
		return
	}

	lines := strings.Split(string(output), "\n")
	for _, line := range lines {
		if strings.Contains(line, "User:") {
			parts := strings.Fields(line)
			if len(parts) >= 4 {
				uid, _ = strconv.Atoi(parts[1])
				gid, _ = strconv.Atoi(parts[3])
			}
			break
		}
	}
	return
}

// Usage returns the calculated used space in bytes and percentage.
func (s *Stats) Usage() (usedBytes int64, percent float64) {
	if s.TotalBlocks == 0 {
		return 0, 0
	}
	usedBlocks := s.TotalBlocks - s.FreeBlocks
	usedBytes = usedBlocks * s.BlockSize
	percent = (float64(usedBlocks) / float64(s.TotalBlocks)) * 100.0
	return
}

// InodeUsage returns the calculated inode usage percentage.
func (s *Stats) InodeUsage() (percent float64) {
	if s.TotalInodes == 0 {
		return 0
	}
	usedInodes := s.TotalInodes - s.FreeInodes
	return (float64(usedInodes) / float64(s.TotalInodes)) * 100.0
}

// UIDStatus represents the result of InspectImageUIDStatus.
type UIDStatus int

const (
	UIDStatusUnknown UIDStatus = iota
	// UIDStatusRoot indicates the 'upper' directory is owned by root (UID 0).
	UIDStatusRoot
	// UIDStatusCurrentUser indicates it's owned by the current user.
	UIDStatusCurrentUser
	// UIDStatusDifferentUser indicates it's owned by a different user.
	UIDStatusDifferentUser
)

func (s UIDStatus) String() string {
	switch s {
	case UIDStatusRoot:
		return "root"
	case UIDStatusCurrentUser:
		return "current user"
	case UIDStatusDifferentUser:
		return "different user"
	default:
		return "unknown"
	}
}

// InspectImageUIDStatus checks the ownership of the 'upper' directory inside an ext3/4 image.
// Returns:
//   - UIDStatusRoot: owned by root (UID 0)
//   - UIDStatusCurrentUser: owned by current user
//   - UIDStatusDifferentUser: owned by a different user
//   - UIDStatusUnknown: failed to check or not an image file
//
// This is used to determine if --fakeroot should be automatically enabled.
func InspectImageUIDStatus(imgPath string) UIDStatus {
	// Check if debugfs is available
	debugfsPath, err := exec.LookPath("debugfs")
	if err != nil {
		return UIDStatusUnknown // debugfs not found
	}

	// Use debugfs to check the UID of the 'upper' directory
	cmd := exec.Command(debugfsPath, "-R", "stat upper", imgPath)
	output, err := cmd.CombinedOutput()
	if err != nil {
		// Failed to execute debugfs (might not be an ext3/4 image)
		return UIDStatusUnknown
	}

	// Parse the output to extract UID
	// Expected format: "User:  1000   Group:  1000"
	lines := strings.Split(string(output), "\n")
	for _, line := range lines {
		if strings.Contains(line, "User:") {
			parts := strings.Fields(line)
			// Format is typically: "User:  <uid>  Group:  <gid>"
			if len(parts) >= 2 {
				uidStr := parts[1]
				uid, err := strconv.Atoi(uidStr)
				if err != nil {
					return UIDStatusUnknown // Failed to parse UID
				}

				// Compare with current user's UID
				currentUID := os.Getuid()
				switch uid {
				case 0:
					return UIDStatusRoot
				case currentUID:
					return UIDStatusCurrentUser
				default:
					return UIDStatusDifferentUser
				}
			}
		}
	}

	return UIDStatusUnknown // UID information not found in output
}
