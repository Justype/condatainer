package overlay

import (
	"os"
	"os/exec"
	"strconv"
	"strings"
)

// Stats holds metadata about the overlay filesystem.
type Stats struct {
	BlockSize       int64
	TotalBlocks     int64
	FreeBlocks      int64
	TotalInodes     int64
	FreeInodes      int64
	FilesystemState string // "clean", "dirty", etc.
	LastMounted     string
}

// GetStats reads the filesystem superblock to calculate usage without mounting.
func GetStats(path string) (*Stats, error) {
	if err := checkDependencies([]string{"tune2fs"}); err != nil {
		return nil, err
	}

	// tune2fs -l lists the superblock info
	cmd := exec.Command("tune2fs", "-l", path)
	out, err := cmd.Output()
	if err != nil {
		return nil, &Error{Op: "read stats", Path: path, Tool: "tune2fs", BaseErr: err}
	}

	stats := &Stats{}
	lines := strings.Split(string(out), "\n")

	for _, line := range lines {
		parts := strings.SplitN(line, ":", 2)
		if len(parts) != 2 {
			continue
		}
		key := strings.TrimSpace(parts[0])
		val := strings.TrimSpace(parts[1])

		switch key {
		case "Block size":
			stats.BlockSize, _ = strconv.ParseInt(val, 10, 64)
		case "Block count":
			stats.TotalBlocks, _ = strconv.ParseInt(val, 10, 64)
		case "Free blocks":
			stats.FreeBlocks, _ = strconv.ParseInt(val, 10, 64)
		case "Inode count":
			stats.TotalInodes, _ = strconv.ParseInt(val, 10, 64)
		case "Free inodes":
			stats.FreeInodes, _ = strconv.ParseInt(val, 10, 64)
		case "Filesystem state":
			stats.FilesystemState = val
		case "Last mounted on":
			stats.LastMounted = val
		}
	}

	return stats, nil
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
