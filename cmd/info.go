package cmd

import (
	"fmt"
	"os"
	"path/filepath"
	"strings"
	"time"

	"github.com/Justype/condatainer/internal/config"
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
	Example: `  condatainer info samtools/1.22 # Show info for installed overlay
  condatainer info env.img       # Show info for local overlay file`,
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
	} else if !strings.Contains(normalized, "/") && config.Global.DefaultDistro != "" {
		// Bare name not found: try <default_distro>/<name> (e.g. "igv" → "ubuntu24/igv", "base_image" → "ubuntu24/base_image")
		if path, ok := installedOverlays[config.Global.DefaultDistro+"/"+normalized]; ok {
			overlayPath = path
		}
	}

	if overlayPath == "" {
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

// displaySqfInfo prints rich info for a SquashFS overlay.
func displaySqfInfo(overlayPath string) error {
	fileInfo, err := os.Stat(overlayPath)
	if err != nil {
		return fmt.Errorf("failed to stat file: %w", err)
	}

	sqStats, unsquashfsErr := overlay.GetSquashFSStats(overlayPath)
	overlayType := overlay.GetOverlayType(overlayPath)

	// File section
	fmt.Println(utils.StyleTitle("File"))
	fmt.Printf("  %-14s %s\n", "Name:", utils.StyleName(filepath.Base(overlayPath)))
	displayWhatis(overlayPath)
	fmt.Printf("  %-14s %s\n", "Path:", utils.StylePath(overlayPath))
	fmt.Printf("  %-14s %s\n", "Size:", utils.FormatBytes(fileInfo.Size()))
	if overlayType != "" {
		fmt.Printf("  %-14s %s (Read-Only)\n", "Type:", utils.StyleInfo(overlayType))
	}
	if overlayType == "OS Overlay" {
		if osInfo := overlay.GetOSInfo(overlayPath); osInfo != nil {
			fmt.Printf("  %-14s %s\n", "Distro:", osInfo.String())
		}
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

	// Conda environment section (sqf overlays use /cnt/<name/version>)
	displayCondaEnv(overlayPath, "/cnt/"+name)

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
	displayWhatis(overlayPath)
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

	// Conda environment section (img overlays always use /ext3/env)
	displayCondaEnv(overlayPath, "/ext3/env")

	// Environment variables section
	displayEnvVars(overlayPath)

	return nil
}

// readEnvFile parses the .env sidecar and returns whatis, notes, and var lines.
func readEnvFile(overlayPath string) (whatis string, notes map[string]string, varLines []string) {
	notes = map[string]string{}
	data, err := os.ReadFile(overlayPath + ".env")
	if err != nil {
		return
	}
	for _, line := range strings.Split(string(data), "\n") {
		line = strings.TrimSpace(line)
		if line == "" {
			continue
		}
		if after, ok := strings.CutPrefix(line, "#WHATIS:"); ok {
			whatis = strings.TrimSpace(after)
		} else if after, ok := strings.CutPrefix(line, "#ENVNOTE:"); ok {
			if k, v, found := strings.Cut(after, "="); found {
				notes[strings.TrimSpace(k)] = strings.TrimSpace(v)
			}
		} else {
			varLines = append(varLines, line)
		}
	}
	return
}

// displayWhatis prints the Whatis line from the .env sidecar in the File section.
func displayWhatis(overlayPath string) {
	whatis, _, _ := readEnvFile(overlayPath)
	if whatis != "" {
		fmt.Printf("  %-14s %s\n", "Whatis:", whatis)
	}
}

// displayCondaEnv reads conda-meta/history from inside the overlay and prints
// a "Conda Env" section with the channels and explicitly-installed packages.
// Silently does nothing if the history file is absent or unreadable.
func displayCondaEnv(overlayPath, envPrefix string) {
	info := overlay.ReadCondaInfo(overlayPath, envPrefix)
	if info == nil {
		return
	}

	fmt.Println(utils.StyleTitle("Conda Env"))
	if len(info.Channels) > 0 {
		printWrapped("Channels:", strings.Join(info.Channels, " "), 4)
	}
	printWrapped("Packages:", strings.Join(info.Specs, " "), 4)
}

// printWrapped prints label + content (comma-separated tokens) wrapping at word
// boundaries to fit the terminal width, across at most maxLines lines.
// Any tokens that don't fit are summarised as "... and N more".
func printWrapped(label, content string, maxLines int) {
	prefix := "  " + label + " "
	indent := "    "
	tw := terminalWidth()
	availFirst := max(tw-len(prefix), 20)
	availRest := max(tw-len(indent), 20)

	tokens := strings.Fields(content)
	var lines []string
	var cur strings.Builder
	remaining := 0

	for _, tok := range tokens {
		if maxLines > 0 && len(lines) >= maxLines {
			remaining++
			continue
		}
		avail := availRest
		if len(lines) == 0 {
			avail = availFirst
		}
		sep := ""
		if cur.Len() > 0 {
			sep = "   "
		}
		if cur.Len()+len(sep)+len(tok) <= avail {
			cur.WriteString(sep)
			cur.WriteString(tok)
		} else {
			if cur.Len() > 0 {
				lines = append(lines, cur.String())
				cur.Reset()
			}
			if maxLines > 0 && len(lines) >= maxLines {
				remaining++
				continue
			}
			cur.WriteString(tok)
		}
	}
	if cur.Len() > 0 {
		lines = append(lines, cur.String())
	}

	// If there are remaining tokens, ensure the suffix fits on the last line.
	// Bleed tokens back from the last line into remaining until it fits.
	if remaining > 0 && len(lines) > 0 {
		suffixAvail := availRest
		if len(lines) == 1 {
			suffixAvail = availFirst
		}
		for {
			suffix := fmt.Sprintf(" ... and %d more", remaining)
			last := lines[len(lines)-1]
			if len(last)+len(suffix) <= suffixAvail {
				break
			}
			idx := strings.LastIndex(last, "   ")
			if idx < 0 {
				break // single token on last line; let it overflow
			}
			lines[len(lines)-1] = last[:idx]
			remaining++
		}
	}

	for i, line := range lines {
		suffix := ""
		if i == len(lines)-1 && remaining > 0 {
			suffix = fmt.Sprintf(" ... and %d more", remaining)
		}
		if i == 0 {
			fmt.Printf("%s%s%s\n", prefix, line, suffix)
		} else {
			fmt.Printf("%s%s%s\n", indent, line, suffix)
		}
	}
}

// displayEnvVars prints the Environment section (env vars only) from the .env sidecar.
func displayEnvVars(overlayPath string) {
	_, notes, varLines := readEnvFile(overlayPath)
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
