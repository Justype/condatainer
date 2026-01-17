package overlay

import (
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
