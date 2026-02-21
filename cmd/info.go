package cmd

import (
	"fmt"
	"os"
	"os/exec"
	"path/filepath"
	"strconv"
	"strings"
	"time"

	"github.com/Justype/condatainer/internal/overlay"
	"github.com/Justype/condatainer/internal/utils"
	"github.com/spf13/cobra"
)

var infoOverlayCmd = &cobra.Command{
	Use:   "info [overlay]",
	Short: "Show information about a specific overlay",
	Long: `Display detailed information about an installed overlay or external overlay file.

For SquashFS (.sqf) overlays: shows compression details, inode count, block size, and mount path.
For ext3 (.img) overlays: shows filesystem stats, disk/inode usage, block size, and ownership.`,
	Example: `  condatainer info samtools/1.22    # Show info for installed overlay
  condatainer info env.img          # Show info for local overlay file`,
	Args:              cobra.ExactArgs(1),
	SilenceUsage:      true,
	ValidArgsFunction: completeInfoArgs,
	RunE:              runInfoOverlay,
}

func init() {
	rootCmd.AddCommand(infoOverlayCmd)
}

// completeInfoArgs restricts file completion to .sqf and .img files with folder navigation.
func completeInfoArgs(cmd *cobra.Command, args []string, toComplete string) ([]string, cobra.ShellCompDirective) {
	if len(args) > 0 {
		return nil, cobra.ShellCompDirectiveNoFileComp
	}
	return []string{"sqf", "img"}, cobra.ShellCompDirectiveFilterFileExt
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

	if utils.IsSqf(overlayPath) {
		return displaySqfInfo(overlayPath)
	} else if utils.IsImg(overlayPath) {
		return displayImgInfo(overlayPath)
	}

	return fmt.Errorf("unsupported file type: %s", filepath.Ext(overlayPath))
}

// SquashFSStats holds metadata parsed from `unsquashfs -stat`.
type SquashFSStats struct {
	CreatedTime         string
	FilesystemSizeBytes int64
	Compression         string
	CompressionLevel    int // 0 means not specified by the file
	BlockSize           int64
	NumFragments        int64
	NumInodes           int64
	NumIDs              int64
	DuplicatesRemoved   bool
	ExportableNFS       bool
}

// getSquashFSStats runs `unsquashfs -stat` on the given path and parses its output.
func getSquashFSStats(path string) (*SquashFSStats, error) {
	cmd := exec.Command("unsquashfs", "-stat", path)
	cmd.Env = append(os.Environ(), "LC_ALL=C", "LC_TIME=C")
	out, err := cmd.Output()
	if err != nil {
		return nil, fmt.Errorf("unsquashfs -stat failed: %w", err)
	}

	stats := &SquashFSStats{}
	lines := strings.Split(string(out), "\n")

	for i, raw := range lines {
		line := strings.TrimSpace(raw)
		switch {
		case strings.HasPrefix(line, "Creation or last append time "):
			stats.CreatedTime = strings.TrimPrefix(line, "Creation or last append time ")

		case strings.HasPrefix(line, "Filesystem size "):
			// "Filesystem size 24511225 bytes (23936.74 Kbytes / 23.38 Mbytes)"
			parts := strings.Fields(line)
			if len(parts) >= 3 {
				stats.FilesystemSizeBytes, _ = strconv.ParseInt(parts[2], 10, 64)
			}

		case strings.HasPrefix(line, "Compression "):
			stats.Compression = strings.TrimPrefix(line, "Compression ")
			// Next line may be indented "compression-level N"
			if i+1 < len(lines) {
				next := strings.TrimSpace(lines[i+1])
				if strings.HasPrefix(next, "compression-level ") {
					lvl, _ := strconv.Atoi(strings.TrimPrefix(next, "compression-level "))
					stats.CompressionLevel = lvl
				}
			}

		case strings.HasPrefix(line, "Block size "):
			stats.BlockSize, _ = strconv.ParseInt(strings.TrimPrefix(line, "Block size "), 10, 64)

		case strings.HasPrefix(line, "Number of fragments "):
			stats.NumFragments, _ = strconv.ParseInt(strings.TrimPrefix(line, "Number of fragments "), 10, 64)

		case strings.HasPrefix(line, "Number of inodes "):
			stats.NumInodes, _ = strconv.ParseInt(strings.TrimPrefix(line, "Number of inodes "), 10, 64)

		case strings.HasPrefix(line, "Number of ids "):
			stats.NumIDs, _ = strconv.ParseInt(strings.TrimPrefix(line, "Number of ids "), 10, 64)

		case strings.Contains(line, "Duplicates are removed"):
			stats.DuplicatesRemoved = true

		case strings.HasPrefix(line, "Filesystem is exportable via NFS"):
			stats.ExportableNFS = true
		}
	}

	return stats, nil
}

// normalizeTime parses ctime()-style timestamps produced by tune2fs and unsquashfs
// (both run with LC_ALL=C) and reformats them as "2006-01-02 15:04:05".
// Falls back to the original string if parsing fails.
func normalizeTime(s string) string {
	// Already ISO — pass through immediately.
	if len(s) >= 19 && s[4] == '-' && s[7] == '-' {
		return s[:19]
	}

	// C-locale ctime() format: "Mon Jan _2 15:04:05 2006"
	// The day field is space-padded for single-digit days.
	formats := []string{
		"Mon Jan _2 15:04:05 2006", // space-padded day (canonical C locale)
		"Mon Jan  2 15:04:05 2006", // two-space pad (explicit match)
		"Mon Jan 02 15:04:05 2006", // zero-padded day
		"Mon Jan 2 15:04:05 2006",  // bare single digit
	}
	for _, f := range formats {
		if t, err := time.Parse(f, s); err == nil {
			return t.Format("2006-01-02 15:04:05")
		}
	}
	return s
}

// getSqfOverlayType classifies a SquashFS overlay based on its top-level directory
// entries and filename:
//   - OS Overlay     : contains .singularity.d
//   - Module Overlay : contains cnt and filename uses -- separators (name--version)
//   - Bundle Overlay : contains cnt but filename has no --
func getSqfOverlayType(overlayPath string) string {
	if sqfPathExists(overlayPath, ".singularity.d") {
		return "OS Overlay"
	}
	if sqfPathExists(overlayPath, "cnt") {
		base := strings.TrimSuffix(filepath.Base(overlayPath), filepath.Ext(overlayPath))
		if strings.Contains(base, "--") {
			return "Module Overlay"
		}
		return "Bundle Overlay"
	}
	return ""
}

// displaySqfInfo prints rich info for a SquashFS overlay.
func displaySqfInfo(overlayPath string) error {
	fileInfo, err := os.Stat(overlayPath)
	if err != nil {
		return fmt.Errorf("failed to stat file: %w", err)
	}

	sqStats, unsquashfsErr := getSquashFSStats(overlayPath)
	overlayType := getSqfOverlayType(overlayPath)

	// File section
	fmt.Println(utils.StyleTitle("File"))
	fmt.Printf("  %-14s %s\n", "Name:", utils.StyleName(filepath.Base(overlayPath)))
	fmt.Printf("  %-14s %s\n", "Path:", utils.StylePath(overlayPath))
	fmt.Printf("  %-14s %s\n", "Size:", utils.FormatBytes(fileInfo.Size()))
	if overlayType != "" {
		fmt.Printf("  %-14s %s (Read-Only)\n", "Type:", utils.StyleInfo(overlayType))
	}
	if sqStats != nil && sqStats.CreatedTime != "" {
		fmt.Printf("  %-14s %s\n", "Created:", normalizeTime(sqStats.CreatedTime))
	} else {
		fmt.Printf("  %-14s %s\n", "Modified:", fileInfo.ModTime().Format("2006-01-02 15:04:05"))
	}

	// SquashFS section
	fmt.Println(utils.StyleTitle("SquashFS"))
	if unsquashfsErr == nil && sqStats != nil {
		comprStr := sqStats.Compression
		if sqStats.CompressionLevel > 0 {
			comprStr = fmt.Sprintf("%s (level %d)", sqStats.Compression, sqStats.CompressionLevel)
		}
		fmt.Printf("  %-14s %s\n", "Compression:", utils.StyleInfo(comprStr))
		fmt.Printf("  %-14s %s\n", "Block Size:", utils.FormatBytes(sqStats.BlockSize))
		fmt.Printf("  %-14s %d\n", "Inodes:", sqStats.NumInodes)
		fmt.Printf("  %-14s %d\n", "Fragments:", sqStats.NumFragments)
		if sqStats.DuplicatesRemoved {
			fmt.Printf("  %-14s %s\n", "Deduplication:", utils.StyleInfo("enabled"))
		}
	} else {
		fmt.Printf("  %-14s unsquashfs not available\n", "Details:")
	}

	// Mount path section
	name := strings.TrimSuffix(filepath.Base(overlayPath), filepath.Ext(overlayPath))
	name = strings.ReplaceAll(name, "--", "/")
	name = strings.ReplaceAll(name, "=", "/")
	if overlayType == "Module Overlay" || overlayType == "Bundle Overlay" {
		fmt.Println(utils.StyleTitle("Mount"))
		fmt.Printf("  %-14s /cnt/%s\n", "Path:", name)
	}

	// Environment variables section
	displayEnvVars(overlayPath)

	return nil
}

// displayImgInfo prints rich info for an ext3 overlay image,
// identical to what `condatainer overlay info` shows.
func displayImgInfo(overlayPath string) error {
	stats, err := overlay.GetStats(overlayPath)
	if err != nil {
		return fmt.Errorf("failed to get stats: %w", err)
	}

	usedBytes, diskPct := stats.Usage()
	totalBytes := stats.TotalBlocks * stats.BlockSize
	reservedBytes := stats.ReservedBlocks * stats.BlockSize
	freeBytes := stats.FreeBlocks*stats.BlockSize - reservedBytes
	usedInodes := stats.TotalInodes - stats.FreeInodes
	inodePct := stats.InodeUsage()
	fileInfo, statErr := os.Stat(overlayPath)

	// File section
	fmt.Println(utils.StyleTitle("File"))
	fmt.Printf("  %-14s %s\n", "Name:", utils.StyleName(filepath.Base(overlayPath)))
	fmt.Printf("  %-14s %s\n", "Path:", utils.StylePath(overlayPath))
	fmt.Printf("  %-14s %s\n", "Size:", utils.FormatBytes(stats.FileSizeBytes))
	if stats.IsSparse {
		fmt.Printf("  %-14s %s (Writable; %s on disk)\n", "Type:", utils.StyleInfo("Workspace Overlay"), utils.FormatBytes(stats.FileBlocksUsed))
	} else {
		fmt.Printf("  %-14s %s (Writable)\n", "Type:", utils.StyleInfo("Workspace Overlay"))
	}
	// Filesystem section
	fmt.Println(utils.StyleTitle("Filesystem"))
	fmt.Printf("  %-14s %s\n", "Format:", utils.StyleInfo(stats.FilesystemType))
	fmt.Printf("  %-14s %s\n", "State:", stats.FilesystemState)
	fmt.Printf("  %-14s %d bytes\n", "Block Size:", stats.BlockSize)
	// fmt.Printf("  %-14s %s\n", "UUID:", stats.FilesystemUUID)
	fmt.Printf("  %-14s %s\n", "Created:", normalizeTime(stats.CreatedTime))
	if statErr == nil {
		fmt.Printf("  %-14s %s\n", "Modified:", fileInfo.ModTime().Format("2006-01-02 15:04:05"))
	}
	if stats.LastMounted != "" && stats.LastMounted != "<not available>" {
		fmt.Printf("  %-14s %s\n", "Last Mounted:", normalizeTime(stats.LastMounted))
	}

	// Ownership section
	fmt.Println(utils.StyleTitle("Ownership"))
	if stats.UpperUID == 0 && stats.UpperGID == 0 {
		fmt.Printf("  %-14s %s (use with --fakeroot)\n", "Inner Owner:", utils.StyleInfo("root"))
	} else if stats.UpperUID >= 0 {
		fmt.Printf("  %-14s UID=%d GID=%d\n", "Inner Owner:", stats.UpperUID, stats.UpperGID)
	} else {
		fmt.Printf("  %-14s unknown\n", "Inner Owner:")
	}

	// Disk Usage section
	fmt.Println(utils.StyleTitle("Disk Usage"))
	fmt.Printf("  %-14s %s / %s (%.2f%%)\n", "Used:",
		utils.FormatBytes(usedBytes),
		utils.FormatBytes(totalBytes),
		diskPct)
	if stats.ReservedBlocks > 0 {
		reservedPct := float64(stats.ReservedBlocks) / float64(stats.TotalBlocks) * 100
		fmt.Printf("  %-14s %s (%.1f%% for root)\n", "Reserved:", utils.FormatBytes(reservedBytes), reservedPct)
	}
	fmt.Printf("  %-14s %s\n", "Free:", utils.FormatBytes(freeBytes))

	// Inode Usage section
	fmt.Println(utils.StyleTitle("Inode Usage"))
	fmt.Printf("  %-14s %d / %d (%.2f%%)\n", "Used:",
		usedInodes,
		stats.TotalInodes,
		inodePct)
	fmt.Printf("  %-14s %d\n", "Free:", stats.FreeInodes)

	// Mount path section
	fmt.Println(utils.StyleTitle("Mount"))
	fmt.Printf("  %-14s %s\n", "Path:", "/ext3/env")

	// Environment variables section
	displayEnvVars(overlayPath)

	return nil
}

// displayEnvVars prints the Environment section from the .env sidecar file, if present.
// Lines of the form #ENVNOTE:VAR=description are shown as inline annotations on the
// corresponding VAR=... line instead of being printed as raw entries.
func displayEnvVars(overlayPath string) {
	envPath := overlayPath + ".env"
	if !utils.FileExists(envPath) {
		return
	}
	data, err := os.ReadFile(envPath)
	if err != nil {
		return
	}

	// First pass: collect notes keyed by variable name.
	notes := map[string]string{}
	var varLines []string
	for _, line := range strings.Split(string(data), "\n") {
		line = strings.TrimSpace(line)
		if line == "" {
			continue
		}
		if after, ok := strings.CutPrefix(line, "#ENVNOTE:"); ok {
			// Format: #ENVNOTE:VARNAME=description
			if k, v, found := strings.Cut(after, "="); found {
				notes[strings.TrimSpace(k)] = strings.TrimSpace(v)
			}
		} else {
			varLines = append(varLines, line)
		}
	}

	if len(varLines) == 0 {
		return
	}

	fmt.Println(utils.StyleTitle("Environment"))
	for _, line := range varLines {
		fmt.Printf("  - %s\n", line)
		varName, _, _ := strings.Cut(line, "=")
		if note, ok := notes[varName]; ok {
			fmt.Printf("    %s\n", utils.StyleNote("# "+note))
		}
	}
}
