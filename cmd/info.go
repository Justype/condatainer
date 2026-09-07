package cmd

import (
	"fmt"
	"os"
	"path/filepath"
	"sort"
	"strings"

	"github.com/Justype/condatainer/catalog"
	"github.com/Justype/condatainer/internal/artifact/compare"
	"github.com/Justype/condatainer/internal/artifact/meta"
	"github.com/Justype/condatainer/internal/config"
	"github.com/Justype/condatainer/internal/image/ext3"
	"github.com/Justype/condatainer/internal/image/sif"
	"github.com/Justype/condatainer/internal/image/squashfs"
	"github.com/Justype/condatainer/internal/runtime/container"
	"github.com/Justype/condatainer/internal/store"
	"github.com/Justype/condatainer/internal/utils"
	"github.com/spf13/cobra"
)

// overlayInfoHelp is shared by 'info' and 'overlay info', which run the same command.
const overlayInfoHelp = `Show details about an installed overlay or a local overlay file.

- SquashFS (.sqf): compression details, inode count, block size, mount path
- ext3 (.img): filesystem stats, disk/inode usage, block size, ownership`

// overlayInfoExample is shared by 'info' and 'overlay info'.
const overlayInfoExample = `  condatainer info samtools/1.22 # Show info for installed overlay
  condatainer info env.img       # Show info for local overlay file`

var infoOverlayCmd = &cobra.Command{
	Use:               "info [flags] <overlay>",
	Short:             "Show details about an overlay",
	Long:              overlayInfoHelp,
	Example:           overlayInfoExample,
	Args:              cobra.ExactArgs(1),
	SilenceUsage:      true,
	ValidArgsFunction: completeOverlayArg,
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
	return []string{"sqf", "img", "ext3"}, cobra.ShellCompDirectiveFilterFileExt
}

func runInfoOverlay(cmd *cobra.Command, args []string) error {
	overlayArg := args[0]
	normalized := catalog.Normalize(overlayArg)

	// Try to find as installed overlay first
	installedOverlays, err := getInstalledOverlaysMap()
	if err != nil {
		return err
	}

	var overlayPath string
	if path, ok := installedOverlays[normalized]; ok {
		overlayPath = path
	} else if !strings.Contains(normalized, "/") && config.ResolvedDefaultDistro() != "" {
		// Bare name not found: try <base>/<name> (e.g. "build-essential" → "ubuntu24/build-essential", "base_image" → "ubuntu24/base_image")
		if path, ok := installedOverlays[config.ResolvedDefaultDistro()+"/"+normalized]; ok {
			overlayPath = path
		}
	}

	// The base is not in the installed-overlay map — that map drives `remove`, and
	// the container root is not something to delete by name. Resolve it here only.
	if overlayPath == "" && normalized == config.BaseRecipeName() {
		overlayPath = config.FindBaseImage()
	}

	if overlayPath == "" {
		// Try as external file path
		overlayPath, _ = filepath.Abs(overlayArg)
		if !utils.FileExists(overlayPath) {
			return fmt.Errorf("overlay file %s not found", utils.StylePath(overlayPath))
		}
	}

	switch {
	case utils.IsSqf(overlayPath), utils.IsSif(overlayPath):
		return displayImageInfo(overlayPath)
	case utils.IsImg(overlayPath):
		return displayImgInfo(overlayPath)
	}

	return fmt.Errorf("unsupported file type: %s", filepath.Ext(overlayPath))
}

// imageType reports the recipe type the image was built as — app, base, os or
// data — taken from its embedded metadata. An image built before that metadata
// has none, and its layout cannot tell app from data, so that reports "" for
// unknown.
func imageType(imagePath string) catalog.Type {
	if rt, err := meta.ReadRuntime(imagePath); err == nil {
		switch rt.Type {
		case catalog.TypeBase, catalog.TypeOS, catalog.TypeApp, catalog.TypeData, catalog.TypeEnv:
			return rt.Type
		}
	}
	return ""
}

// typeLabel is how a type is written for a reader. Only env differs from its
// stored value: it is spelled out so it cannot be read as the environment
// variables the rest of the tool calls env, and so a .sqf and the .img it came
// from say the same word.
func typeLabel(typ catalog.Type) string {
	switch typ {
	case catalog.TypeEnv:
		return "environment"
	case "":
		return "unknown"
	}
	return string(typ)
}

// displayPayload prints where the image's payload sits inside the container and
// what it deletes from the layers below.
//
// An environment shows no prefix of its own: every one of them records the same
// EnvPrefix, so printing it says nothing that the type has not already said.
func displayPayload(overlayPath string, typ catalog.Type, prefix string) {
	whiteouts, convention := 0, ""
	if typ == catalog.TypeEnv {
		if m, err := meta.ReadManifest(overlayPath); err == nil && m.Snapshot != nil {
			whiteouts, convention = m.Snapshot.Whiteouts, m.Snapshot.Convention
		}
	}
	showPrefix := prefix != "" && (typ == catalog.TypeApp || typ == catalog.TypeData)
	if !showPrefix && whiteouts == 0 {
		return
	}

	fmt.Println(utils.StyleTitle("Payload"))
	if showPrefix {
		fmt.Printf("  %-14s %s\n", "Prefix:", prefix)
	}
	if whiteouts > 0 {
		fmt.Printf("  %-14s %d (from %s markers)\n", "Deletions:", whiteouts, convention)
	}
}

// imagePrefix reports where the payload sits inside the container, or "" when the
// artifact does not say. Only the artifact is asked: a filename is an address,
// never a naming claim, and decoding one changes the answer when a file is
// renamed.
func imagePrefix(imagePath string) string {
	rt, err := meta.ReadRuntime(imagePath)
	if err != nil {
		return ""
	}
	return rt.Prefix
}

// normalizeTime parses ctime()-style timestamps produced by tune2fs and unsquashfs
// (both run with LC_ALL=C) and reformats them as "2006-01-02 15:04:05".
// Falls back to the original string if parsing fails.
func normalizeTime(s string) string {
	if t, ok := squashfs.ParseStatTime(s); ok {
		return t.Format("2006-01-02 15:04:05")
	}
	return s
}

// displayImageInfo prints rich info for an immutable image. A foreign
// SIF's payload starts partway into the file, so its archive reads take that
// offset; a .sqf is the same reads at 0.
func displayImageInfo(overlayPath string) error {
	fileInfo, err := os.Stat(overlayPath)
	if err != nil {
		return fmt.Errorf("failed to stat file: %w", err)
	}

	var offset int64
	if utils.IsSif(overlayPath) {
		part, err := sif.PrimarySystemPartition(overlayPath)
		if err != nil {
			return fmt.Errorf("not a readable SIF: %w", err)
		}
		offset = part.Offset
	}

	sqStats, unsquashfsErr := squashfs.GetSquashFSStatsAt(overlayPath, offset)
	typ := imageType(overlayPath)
	description, notes, varLines, diagnostics := readEnvFile(overlayPath)
	printDiagnostics(diagnostics)

	// File section
	fmt.Println(utils.StyleTitle("File"))
	fmt.Printf("  %-14s %s\n", "Name:", utils.StyleName(filepath.Base(overlayPath)))
	displayDescription(description)
	fmt.Printf("  %-14s %s\n", "Path:", utils.StylePath(overlayPath))
	fmt.Printf("  %-14s %s\n", "Size:", utils.FormatSize(fileInfo.Size()))
	fmt.Printf("  %-14s %s (Read-Only)\n", "Type:", utils.StyleInfo(typeLabel(typ)))
	if typ == catalog.TypeOS || typ == catalog.TypeBase {
		if osInfo := squashfs.GetOSInfoAt(overlayPath, offset); osInfo != nil {
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
		fmt.Printf("  %-14s %s\n", "Block Size:", utils.FormatSize(sqStats.BlockSize))
		fmt.Printf("  %-14s %d\n", "Inodes:", sqStats.NumInodes)
		fmt.Printf("  %-14s %d\n", "Fragments:", sqStats.NumFragments)
		if sqStats.DuplicatesRemoved {
			fmt.Printf("  %-14s %s\n", "Deduplication:", utils.StyleInfo("enabled"))
		}
	} else {
		fmt.Printf("  %-14s unsquashfs not available\n", "Details:")
	}

	// Where the payload sits inside the container. Nothing is mounted there: the
	// whole image is mounted as one overlay layer, and this is the prefix its
	// contents appear under.
	prefix := imagePrefix(overlayPath)
	displayPayload(overlayPath, typ, prefix)

	// Conda environment section, read from under that same prefix. Without one
	// there is nowhere to look, and guessing would report another image's layout.
	if prefix != "" {
		displayCondaEnv(overlayPath, prefix)
	}

	// Environment variables section
	displayEnvVars(notes, varLines)

	// Identity addresses one exact build; equivalence says what may stand in for
	// it. Both are what `store` and a project lock resolve by.
	displayKeys(overlayPath)

	// Diagnostic provenance is off the mount path and optional for images built
	// before manifests recorded build tools.
	displayBuildTools(overlayPath)

	return nil
}

// displayKeys prints the artifact's identity and equivalence keys.
//
// The keys are regenerated from the image rather than read from the manifest,
// so what is shown is what `store` and `project` will resolve by. An image
// without readable metadata prints nothing, like the Build section.
func displayKeys(imagePath string) {
	artifact, err := compare.Read(imagePath)
	if err != nil {
		return
	}
	// FormatKeyRef, not Empty: a half-recorded key renders as "" and is not a key.
	identity, equiv := store.FormatKeyRef(artifact.IdentityRef()), store.FormatKeyRef(artifact.EquivRef())
	if identity == "" && equiv == "" {
		return
	}

	fmt.Println(utils.StyleTitle("Keys"))
	if artifact.Name != "" {
		fmt.Printf("  %-14s %s\n", "Artifact:", utils.StyleName(artifact.Name))
	}
	if identity != "" {
		fmt.Printf("  %-14s %s\n", "Identity:", utils.StyleInfo(identity))
	}
	if equiv != "" {
		fmt.Printf("  %-14s %s\n", "Equivalence:", utils.StyleInfo(equiv))
	}
}

func displayBuildTools(imagePath string) {
	m, err := meta.ReadManifest(imagePath)
	if err != nil || m.Build.Tools.Empty() {
		return
	}

	fmt.Println(utils.StyleTitle("Build"))
	tools := m.Build.Tools
	if !tools.Condatainer.Empty() {
		fmt.Printf("  %-14s %s\n", "Condatainer:", utils.StyleInfo(tools.Condatainer.Version))
	}
	if !tools.Apptainer.Empty() {
		value := tools.Apptainer.Version
		if tools.Apptainer.Name != "" && tools.Apptainer.Name != "apptainer" {
			name := tools.Apptainer.Name
			if name == "singularity" {
				name = "Singularity"
			}
			value = name + " " + value
		}
		fmt.Printf("  %-14s %s\n", "Apptainer:", utils.StyleInfo(value))
	}
	if !tools.Micromamba.Empty() {
		fmt.Printf("  %-14s %s\n", "Micromamba:", utils.StyleInfo(tools.Micromamba.Version))
	}
}

// displayImgInfo prints rich info for an ext3 overlay image,
// identical to what `condatainer overlay info` shows.
func displayImgInfo(overlayPath string) error {
	stats, err := ext3.GetStats(overlayPath)
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
	description, notes, varLines, diagnostics := readEnvFile(overlayPath)
	printDiagnostics(diagnostics)

	// File section
	fmt.Println(utils.StyleTitle("File"))
	fmt.Printf("  %-14s %s\n", "Name:", utils.StyleName(filepath.Base(overlayPath)))
	displayDescription(description)
	fmt.Printf("  %-14s %s\n", "Path:", utils.StylePath(overlayPath))
	fmt.Printf("  %-14s %s\n", "Size:", utils.FormatSize(stats.FileSizeBytes))
	// Not a recipe type: an .img is a mutable working overlay, never built from a
	// recipe, so it carries no embedded metadata to take a type from.
	if stats.IsSparse {
		fmt.Printf("  %-14s %s (Writable; %s on disk)\n", "Type:", utils.StyleInfo("environment"), utils.FormatSize(stats.FileBlocksUsed))
	} else {
		fmt.Printf("  %-14s %s (Writable)\n", "Type:", utils.StyleInfo("environment"))
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
		utils.FormatSize(usedBytes),
		utils.FormatSize(totalBytes),
		diskPct)
	if stats.ReservedBlocks > 0 {
		reservedPct := float64(stats.ReservedBlocks) / float64(stats.TotalBlocks) * 100
		fmt.Printf("  %-14s %s (%.1f%% for root)\n", "Reserved:", utils.FormatSize(reservedBytes), reservedPct)
	}
	fmt.Printf("  %-14s %s\n", "Free:", utils.FormatSize(freeBytes))

	// Inode Usage section
	fmt.Println(utils.StyleTitle("Inode Usage"))
	fmt.Printf("  %-14s %d / %d (%.2f%%)\n", "Used:",
		usedInodes,
		stats.TotalInodes,
		inodePct)
	fmt.Printf("  %-14s %d\n", "Free:", stats.FreeInodes)

	// Where a writable overlay's conda prefix sits. Nothing is mounted there —
	// the image is attached as one overlay layer — and the value is fixed rather
	// than read, because an .img carries no manifest to read it from.
	fmt.Println(utils.StyleTitle("Payload"))
	fmt.Printf("  %-14s %s\n", "Prefix:", meta.EnvPrefix)

	// Conda environment section, under that same prefix
	displayCondaEnv(overlayPath, meta.EnvPrefix)

	// Environment variables section
	displayEnvVars(notes, varLines)

	return nil
}

// readEnvFile resolves an image's description, notes, var lines, and non-fatal
// diagnostics from its runtime metadata, with {prefix} resolved to the install
// prefix. A writable .img reads its .env sidecar instead, merged with a paired
// snapshot's own variables when one is autoloaded. Var lines are sorted
// KEY=VALUE. Called once per `info` invocation — its caller prints the
// diagnostics once, rather than once per section that wants a piece of this.
func readEnvFile(overlayPath string) (description string, notes map[string]string, varLines []string, diagnostics []container.Diagnostic) {
	description, configs, notes, diagnostics := container.ResolveOverlayEnv(overlayPath)
	keys := make([]string, 0, len(configs))
	for k := range configs {
		keys = append(keys, k)
	}
	sort.Strings(keys)
	for _, k := range keys {
		varLines = append(varLines, k+"="+configs[k])
	}
	return description, notes, varLines, diagnostics
}

// printDiagnostics prints non-fatal container-setup messages the way
// cmd/internal/ui.RenderExecPlan does for exec/run.
func printDiagnostics(diagnostics []container.Diagnostic) {
	for _, d := range diagnostics {
		switch d.Level {
		case "warn", "warning":
			utils.PrintWarning("%s", d.Message)
		case "note":
			utils.PrintNote("%s", d.Message)
		case "error":
			utils.PrintError("%s", d.Message)
		default:
			utils.PrintMessage("%s", d.Message)
		}
	}
}

// displayDescription prints the Description line in the File section.
func displayDescription(description string) {
	if description != "" {
		fmt.Printf("  %-14s %s\n", "Description:", description)
	}
}

// displayCondaEnv reads conda-meta/history from inside the overlay and prints
// a "Conda Env" section with the channels and explicitly-installed packages.
// Silently does nothing if the history file is absent or unreadable.
//
// container.PairedInfo reads both histories as one continuous log when
// overlayPath is a .img paired with an autoloaded snapshot, rather than the
// .img alone, which would show a nearly-empty history for an environment
// that is actually full.
func displayCondaEnv(overlayPath, envPrefix string) {
	info := container.PairedInfo(overlayPath, envPrefix)
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
func displayEnvVars(notes map[string]string, varLines []string) {
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
