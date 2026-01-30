package cmd

import (
	"fmt"
	"os"
	"os/exec"
	"path/filepath"
	"regexp"
	"strings"
	"syscall"
	"time"

	"github.com/Justype/condatainer/internal/utils"
	"github.com/spf13/cobra"
)

var infoOverlayCmd = &cobra.Command{
	Use:   "info [overlay]",
	Short: "Show information about a specific overlay",
	Long: `Display detailed information about an installed overlay or external overlay file.

Shows file size, compression type, timestamps, mount path, and environment variables.`,
	Example: `  condatainer info samtools/1.22    # Show info for installed overlay
  condatainer info env.img          # Show info for local overlay file`,
	Args:         cobra.ExactArgs(1),
	SilenceUsage: true, // Runtime errors should not show usage
	RunE:         runInfoOverlay,
}

func init() {
	rootCmd.AddCommand(infoOverlayCmd)
}

func runInfoOverlay(cmd *cobra.Command, args []string) error {
	overlayArg := args[0]
	normalized := utils.NormalizeNameVersion(overlayArg)

	// Try to find as installed overlay first
	installedOverlays, err := getInstalledOverlaysMap()
	if err != nil {
		return err
	}

	var overlayPath string
	if path, ok := installedOverlays[normalized]; ok {
		overlayPath = path
	} else {
		// Try as external file path
		overlayPath, _ = filepath.Abs(overlayArg)
		if !utils.FileExists(overlayPath) {
			return fmt.Errorf("overlay file %s not found", utils.StylePath(overlayPath))
		}
	}

	// Display information
	fmt.Printf("Information for %s:\n", utils.StyleName(filepath.Base(overlayPath)))

	// Writable status
	if utils.IsImg(overlayPath) {
		fmt.Println("Writable: Yes")
	} else {
		fmt.Println("Writable: No")
	}

	// File size
	fileInfo, err := os.Stat(overlayPath)
	if err != nil {
		return fmt.Errorf("failed to stat file: %w", err)
	}
	size := fileInfo.Size()
	fmt.Printf("Size: %s\n", formatSize(size))

	// Compression type
	compression := getCompressionType(overlayPath)
	fmt.Printf("Compression: %s\n", compression)

	// Timestamps
	if stat, ok := fileInfo.Sys().(*syscall.Stat_t); ok {
		// Birth time (creation time) - only on some systems
		if stat.Ctim.Sec > 0 {
			ctime := time.Unix(stat.Ctim.Sec, stat.Ctim.Nsec)
			fmt.Printf("Status change: %s\n", ctime.Format("2006-01-02 15:04:05"))
		}
	}

	// Modification time
	mtime := fileInfo.ModTime()
	fmt.Printf("Last modified: %s\n", mtime.Format("2006-01-02 15:04:05"))

	// Access time
	if stat, ok := fileInfo.Sys().(*syscall.Stat_t); ok {
		atime := time.Unix(stat.Atim.Sec, stat.Atim.Nsec)
		fmt.Printf("Last accessed: %s\n", atime.Format("2006-01-02 15:04:05"))
	}

	// Mount path
	if utils.IsSqf(overlayPath) {
		name := strings.TrimSuffix(filepath.Base(overlayPath), filepath.Ext(overlayPath))
		name = strings.ReplaceAll(name, "--", "/")
		name = strings.ReplaceAll(name, "=", "/")
		fmt.Printf("Mount path: /cnt/%s\n", name)
	} else if utils.IsImg(overlayPath) {
		// Check ownership status
		status := inspectImageUIDStatus(overlayPath)
		switch status {
		case 0:
			fmt.Println("Files owner: root owned (--fakeroot required)")
		case 1:
			fmt.Println("Files owner: current UID")
		case 2:
			fmt.Println("Files owner: different UID (chown required)")
		default:
			fmt.Println("Files owner: unknown")
		}
		fmt.Println("Mount path: /ext3/env")
	}

	// Environment variables
	envPath := overlayPath + ".env"
	if utils.FileExists(envPath) {
		fmt.Println("Environment variables:")
		data, err := os.ReadFile(envPath)
		if err == nil {
			lines := strings.Split(string(data), "\n")
			for _, line := range lines {
				line = strings.TrimSpace(line)
				if line != "" {
					fmt.Printf("  - %s\n", line)
				}
			}
		}
	}

	return nil
}

// formatSize formats file size in human-readable format
func formatSize(size int64) string {
	const (
		KB = 1024
		MB = 1024 * KB
		GB = 1024 * MB
	)

	switch {
	case size < KB:
		return fmt.Sprintf("%d Bytes", size)
	case size < MB:
		return fmt.Sprintf("%.2f KiB", float64(size)/KB)
	case size < GB:
		return fmt.Sprintf("%.2f MiB", float64(size)/MB)
	default:
		return fmt.Sprintf("%.2f GiB", float64(size)/GB)
	}
}

// getCompressionType detects compression type using the file command
func getCompressionType(path string) string {
	cmd := exec.Command("file", path)
	output, err := cmd.Output()
	if err != nil {
		return "not available"
	}

	re := regexp.MustCompile(`,\s+(\w+)\s+compressed`)
	matches := re.FindStringSubmatch(string(output))
	if len(matches) > 1 {
		return matches[1]
	}

	return "not available"
}

// inspectImageUIDStatus checks the ownership status of files in an img overlay
// Returns: 0 = root owned, 1 = current UID, 2 = different UID, -1 = unknown
func inspectImageUIDStatus(imagePath string) int {
	// Use debugfs to inspect the image without mounting
	cmd := exec.Command("debugfs", "-R", "ls -l /ext3/env", imagePath)
	output, err := cmd.Output()
	if err != nil {
		return -1
	}

	currentUID := os.Getuid()
	lines := strings.Split(string(output), "\n")

	for _, line := range lines {
		// Parse debugfs output to extract UID
		// Example line format varies, but typically contains UID info
		fields := strings.Fields(line)
		if len(fields) > 2 {
			// Check if files are owned by root (UID 0)
			if strings.Contains(line, "(0)") || strings.HasPrefix(fields[2], "0") {
				return 0
			}
			// Check if owned by current user
			if strings.Contains(line, fmt.Sprintf("(%d)", currentUID)) {
				return 1
			}
			// Different UID
			return 2
		}
	}

	return -1
}
