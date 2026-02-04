package cmd

import (
	"context"
	"errors"
	"os"
	"path/filepath"
	"sort"
	"strings"

	"github.com/Justype/condatainer/internal/apptainer"
	"github.com/Justype/condatainer/internal/build"
	"github.com/Justype/condatainer/internal/config"
	"github.com/Justype/condatainer/internal/exec"
	"github.com/Justype/condatainer/internal/utils"
	"github.com/spf13/cobra"
)

type execCommand struct {
	cobra.Command
}

func (c *execCommand) InitDefaultHelpFlag() {}

var (
	execOverlays    []string
	execKeep        bool
	execWritableImg bool
	execEnvSettings []string
	execBaseImage   string
	execBindPaths   []string
	execFakeroot    bool
	eReadOnly       bool
	eNoAutoload     bool
)

var execCmd = &execCommand{
	Command: cobra.Command{
		Use:          "exec [flags] [overlays...] [--] [command...]",
		Aliases:      []string{"e"},
		Short:        "Execute a command using overlays",
		SilenceUsage: true,
		Long: `Execute a command inside a container with specified overlays.

The exec command can intelligently parse arguments to separate overlays from commands.
Use -o/--overlay to explicitly specify overlays, or let the command parse them automatically.

Examples:
    # Run bash with samtools overlay
    condatainer exec samtools/1.22

    # Explicit overlay specification for standalone tools
    condatainer exec -o samtools/1.22 samtools

    # Use writable .img overlay
    condatainer exec -w env.img bash

    # Set environment variables
    condatainer exec --env MYVAR=value samtools/1.22 bash

If you omit the command, it falls back to bash. When no -o/--overlay flags are provided,
positional arguments are treated as overlays.`,
		RunE: runExec,
	},
}

func init() {
	rootCmd.AddCommand(&execCmd.Command)

	execCmd.Flags().StringSliceVarP(&execOverlays, "overlay", "o", nil, "overlay file to mount (can be used multiple times)")
	execCmd.Flags().BoolVarP(&execKeep, "keep", "k", false, "disable automatic command parsing to installed overlays")
	execCmd.Flags().BoolVarP(&execWritableImg, "writable", "w", false, "mount .img overlays as writable (default: read-only)")
	execCmd.Flags().StringSliceVar(&execEnvSettings, "env", nil, "set environment variable 'KEY=VALUE' (can be used multiple times)")
	execCmd.Flags().StringVarP(&execBaseImage, "base-image", "b", "", "base image to use instead of default")
	execCmd.Flags().StringSliceVar(&execBindPaths, "bind", nil, "bind path 'HOST:CONTAINER' (can be used multiple times)")
	execCmd.Flags().BoolVar(&execFakeroot, "fakeroot", false, "run container with fakeroot privileges")
	execCmd.Flags().BoolVarP(&eReadOnly, "read-only", "r", false, "mount .img overlays as read-only (only applies when using the 'e' shortcut)")
	execCmd.Flags().BoolVarP(&eNoAutoload, "no-autoload", "n", false, "disable autoloading 'env.img' from current directory (only applies to the 'e' shortcut)")

	// --- NEW: Dynamic Completion Logic ---
	execCmd.RegisterFlagCompletionFunc("overlay", overlayFlagCompletion(true, true))
	execCmd.RegisterFlagCompletionFunc("base-image", baseImageFlagCompletion())
	// Positional args completion depends on which command is used
	execCmd.ValidArgsFunction = func(cmd *cobra.Command, args []string, toComplete string) ([]string, cobra.ShellCompDirective) {
		// Check if called via 'e' alias
		// During completion, os.Args contains the __complete command
		// We need to check if 'e' appears in the args or command context
		isShortcut := false

		// Check os.Args for the actual command used
		for i, arg := range os.Args {
			if arg == "__complete" && i+1 < len(os.Args) && os.Args[i+1] == "e" {
				isShortcut = true
				break
			}
		}

		if isShortcut {
			// For 'e': complete overlays before '--', default after
			if containsDashDash(args) {
				return nil, cobra.ShellCompDirectiveDefault
			}
			return overlaySuggestions(true, true, toComplete)
		}
		// For 'exec': use default file completion for positional args
		return nil, cobra.ShellCompDirectiveDefault
	}
	// -------------------------------------

	// Stop flag parsing after the first positional argument so tool arguments (like --help) bypass exec.
	execCmd.Flags().SetInterspersed(false)
}

func runExec(cmd *cobra.Command, args []string) error {
	if execHelpRequested(args) {
		return cmd.Help()
	}

	if err := ensureBaseImage(cmd.Context()); err != nil {
		return err
	}

	installedOverlays, err := exec.InstalledOverlays()
	if err != nil {
		return err
	}

	overlayFinal := []string{}
	commandFinal := []string{}
	isShortcut := cmd.CalledAs() == "e"

	if !isShortcut {
		if !execKeep && len(execOverlays) == 0 && len(args) > 0 {
			utils.PrintDebug("[EXEC] Parsing commands to separate overlays and command...")
			commandStarted := false
			for _, arg := range args {
				if commandStarted {
					commandFinal = append(commandFinal, arg)
					continue
				}
				if arg == "--" {
					commandStarted = true
					continue
				}
				if utils.IsOverlay(arg) {
					overlayFinal = append(overlayFinal, arg)
					continue
				}
				if _, ok := installedOverlays[utils.NormalizeNameVersion(arg)]; ok {
					utils.PrintWarning("Convert command %s to overlay", utils.StyleName(arg))
					overlayFinal = append(overlayFinal, arg)
					continue
				}
				commandStarted = true
				commandFinal = append(commandFinal, arg)
			}
		} else {
			commandFinal = args
		}
	} else {
		overlayFinal, commandFinal = parseShortcutArgs(args)
	}

	overlayFinal = append(overlayFinal, execOverlays...)

	if isShortcut {
		hasImgOverlay := false
		for _, overlay := range overlayFinal {
			if utils.IsImg(overlay) {
				hasImgOverlay = true
				break
			}
		}
		if !hasImgOverlay && !eNoAutoload {
			if pwd, err := os.Getwd(); err == nil {
				localEnvPath := filepath.Join(pwd, "env.img")
				if utils.FileExists(localEnvPath) {
					utils.PrintMessage("Autoload env.img at %s", utils.StylePath(localEnvPath))
					overlayFinal = append(overlayFinal, localEnvPath)
				}
			}
		}
	}

	if len(commandFinal) == 0 {
		commandFinal = []string{"bash"}
	}

	// Resolve user-specified base image if provided
	var baseImageResolved string
	if execBaseImage != "" {
		// If it's an existing file path, use it directly
		// Otherwise, try to resolve it as an overlay name
		if utils.FileExists(execBaseImage) {
			baseImageResolved = execBaseImage
		} else {
			resolvedBase, err := exec.ResolveOverlayPaths([]string{execBaseImage})
			if err == nil && len(resolvedBase) > 0 {
				baseImageResolved = resolvedBase[0]
			} else {
				// Keep the original value (might be an absolute path)
				baseImageResolved = execBaseImage
			}
		}
	}

	resolvedOverlays, err := exec.ResolveOverlayPaths(overlayFinal)
	if err != nil {
		return err
	}

	writableImg := execWritableImg
	if isShortcut {
		writableImg = !eReadOnly
	}
	options := exec.Options{
		Overlays:     resolvedOverlays,
		Command:      commandFinal,
		WritableImg:  writableImg,
		EnvSettings:  execEnvSettings,
		BindPaths:    execBindPaths,
		Fakeroot:     execFakeroot,
		BaseImage:    baseImageResolved, // Empty string triggers GetBaseImage() in ensureDefaults()
		ApptainerBin: config.Global.ApptainerBin,
	}

	if err := exec.Run(cmd.Context(), options); err != nil {
		if errors.Is(err, context.Canceled) || errors.Is(cmd.Context().Err(), context.Canceled) {
			return nil
		}
		// Propagate exit code from container command
		if appErr, ok := err.(*apptainer.ApptainerError); ok {
			if code := appErr.ExitCode(); code >= 0 {
				os.Exit(code)
			}
		}
		return err
	}
	return nil
}

func execHelpRequested(args []string) bool {
	if len(args) == 0 {
		return false
	}
	switch args[0] {
	case "--help", "-h":
		return true
	default:
		return false
	}
}

func ensureBaseImage(ctx context.Context) error {
	return apptainer.EnsureBaseImage(ctx)
}

func overlayFlagCompletion(includeData bool, includeImg bool) func(*cobra.Command, []string, string) ([]string, cobra.ShellCompDirective) {
	return func(cmd *cobra.Command, args []string, toComplete string) ([]string, cobra.ShellCompDirective) {
		return overlaySuggestions(includeData, includeImg, toComplete)
	}
}

func baseImageFlagCompletion() func(*cobra.Command, []string, string) ([]string, cobra.ShellCompDirective) {
	return func(cmd *cobra.Command, args []string, toComplete string) ([]string, cobra.ShellCompDirective) {
		// -b shows ONLY system overlays (from whitelist) + local .sqf/.sif (no .img)
		return systemOverlaySuggestions(toComplete)
	}
}

// systemOverlaySuggestions returns only system overlays (def-built) and local overlay files
func systemOverlaySuggestions(toComplete string) ([]string, cobra.ShellCompDirective) {
	installed, err := exec.InstalledOverlays()
	if err != nil {
		return nil, cobra.ShellCompDirectiveError
	}

	choices := map[string]struct{}{}

	// Only include overlays from the whitelist (system/def-built overlays)
	whitelist := build.GetDefBuiltWhitelist()
	for name := range installed {
		normalized := utils.NormalizeNameVersion(name)
		if whitelist[normalized] {
			if toComplete == "" || strings.HasPrefix(name, toComplete) {
				choices[name] = struct{}{}
			}
		}
	}

	// Add local .sqf/.sif files (not .img) - for -b flag
	for _, candidate := range localImageSuggestions(toComplete) {
		choices[candidate] = struct{}{}
	}

	suggestions := make([]string, 0, len(choices))
	for choice := range choices {
		suggestions = append(suggestions, choice)
	}
	sort.Strings(suggestions)

	return suggestions, cobra.ShellCompDirectiveNoFileComp
}

// localImageSuggestions returns local .sqf and .sif files (for -b flag)
func localImageSuggestions(toComplete string) []string {
	pathDir, _ := filepath.Split(toComplete)
	dirForRead := pathDir
	if dirForRead == "" {
		dirForRead = "."
	}

	entries, err := os.ReadDir(dirForRead)
	if err != nil {
		return nil
	}

	suggestions := []string{}
	for _, entry := range entries {
		if entry.IsDir() {
			continue
		}
		name := entry.Name()

		// Include .sqf and .sif files only (not .img)
		if !utils.IsSqf(name) && !utils.IsSif(name) {
			continue
		}

		candidate := name
		if pathDir != "" {
			candidate = pathDir + name
		}
		if toComplete != "" && !strings.HasPrefix(candidate, toComplete) {
			continue
		}
		suggestions = append(suggestions, candidate)
	}
	return suggestions
}

func overlaySuggestions(includeData bool, includeImg bool, toComplete string) ([]string, cobra.ShellCompDirective) {
	installed, err := exec.InstalledOverlays()
	if err != nil {
		return nil, cobra.ShellCompDirectiveError
	}

	choices := map[string]struct{}{}
	for name, path := range installed {
		if includeData || isAppOverlay(path) {
			if toComplete == "" || strings.HasPrefix(name, toComplete) {
				choices[name] = struct{}{}
			}
		}
	}

	for _, candidate := range localOverlaySuggestions(toComplete, includeImg) {
		choices[candidate] = struct{}{}
	}

	suggestions := make([]string, 0, len(choices))
	for choice := range choices {
		suggestions = append(suggestions, choice)
	}
	sort.Strings(suggestions)

	// Always return NoFileComp to prevent default file listing
	return suggestions, cobra.ShellCompDirectiveNoFileComp
}

func parseShortcutArgs(args []string) ([]string, []string) {
	overlays := []string{}
	command := []string{}
	commandStarted := false
	for i := 0; i < len(args); i++ {
		arg := args[i]
		if arg == "--" {
			commandStarted = true
			continue
		}
		if commandStarted {
			command = append(command, arg)
			continue
		}
		if arg == "-r" || arg == "--read-only" {
			eReadOnly = true
			continue
		}
		if arg == "-n" || arg == "--no-autoload" {
			eNoAutoload = true
			continue
		}
		if arg == "--fakeroot" {
			execFakeroot = true
			continue
		}
		if arg == "-b" || arg == "--base-image" {
			if i+1 < len(args) {
				execBaseImage = args[i+1]
				i++
			}
			continue
		}
		if strings.HasPrefix(arg, "-b=") {
			execBaseImage = strings.TrimPrefix(arg, "-b=")
			continue
		}
		if strings.HasPrefix(arg, "--base-image=") {
			execBaseImage = strings.TrimPrefix(arg, "--base-image=")
			continue
		}
		overlays = append(overlays, arg)
	}
	return overlays, command
}

func containsDashDash(args []string) bool {
	for _, arg := range args {
		if arg == "--" {
			return true
		}
	}
	return false
}

func localOverlaySuggestions(toComplete string, includeImg bool) []string {
	pathDir, _ := filepath.Split(toComplete)
	dirForRead := pathDir
	if dirForRead == "" {
		dirForRead = "."
	}

	entries, err := os.ReadDir(dirForRead)
	if err != nil {
		return nil
	}

	suggestions := []string{}
	for _, entry := range entries {
		if entry.IsDir() {
			continue
		}
		name := entry.Name()

		// Filter based on file type
		if !utils.IsOverlay(name) {
			continue
		}

		// Filter .sif and .img based on includeImg parameter
		if utils.IsSif(name) && includeImg {
			// Skip .sif files when includeImg is true (for -o flag)
			continue
		}
		if utils.IsImg(name) && !includeImg {
			// Skip .img files when includeImg is false (for -b flag)
			continue
		}

		candidate := name
		if pathDir != "" {
			candidate = pathDir + name
		}
		if toComplete != "" && !strings.HasPrefix(candidate, toComplete) {
			continue
		}
		suggestions = append(suggestions, candidate)
	}
	return suggestions
}

func isAppOverlay(path string) bool {
	// Check if overlay name is in the .def-built whitelist
	// .def-built overlays (created via 'sif dump') are considered "app" overlays
	// regardless of their location
	if isDefBuiltOverlay(path) {
		return true
	}

	// Check if path is in any of the image search directories
	// Regular overlays (created via mksquashfs) are only considered "app" if
	// they're in the images directory
	for _, imagesDir := range config.GetImageSearchPaths() {
		if pathWithinDir(path, imagesDir) {
			return true
		}
	}
	return false
}

// isDefBuiltOverlay checks if an overlay was built from a .def file
// by checking against the build-scripts metadata whitelist
func isDefBuiltOverlay(overlayPath string) bool {
	baseName := filepath.Base(overlayPath)
	// Remove extension (e.g., build-essential.sqf -> build-essential)
	nameWithoutExt := strings.TrimSuffix(baseName, filepath.Ext(baseName))
	// Convert -- to / (e.g., build-essential -> build-essential)
	normalized := strings.ReplaceAll(nameWithoutExt, "--", "/")

	// Check whitelist from build package
	whitelist := build.GetDefBuiltWhitelist()
	return whitelist[normalized]
}

func pathWithinDir(path, dir string) bool {
	if dir == "" {
		return false
	}
	rel, err := filepath.Rel(dir, path)
	if err != nil {
		return false
	}
	if rel == "." {
		return true
	}
	return !(rel == ".." || strings.HasPrefix(rel, ".."+string(os.PathSeparator)))
}
