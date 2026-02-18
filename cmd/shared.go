package cmd

import (
	"context"
	"os"
	"path/filepath"
	"sort"
	"strings"

	"github.com/Justype/condatainer/internal/apptainer"
	"github.com/Justype/condatainer/internal/build"
	"github.com/Justype/condatainer/internal/config"
	"github.com/Justype/condatainer/internal/container"
	"github.com/Justype/condatainer/internal/utils"
	"github.com/spf13/cobra"
)

// CommonFlags holds the common flags used by exec, instance start, and instance exec
type CommonFlags struct {
	Overlays    []string
	WritableImg bool
	EnvSettings []string
	BaseImage   string
	BindPaths   []string
	Fakeroot    bool
}

// RegisterCommonFlags registers common flags on a cobra command
func RegisterCommonFlags(cmd *cobra.Command, flags *CommonFlags) {
	cmd.Flags().StringSliceVarP(&flags.Overlays, "overlay", "o", nil, "overlay file to mount (can be used multiple times)")
	cmd.Flags().BoolVarP(&flags.WritableImg, "writable", "w", false, "mount .img overlays as writable (default: read-only)")
	cmd.Flags().BoolVar(&flags.WritableImg, "writable-img", false, "Alias for --writable")
	cmd.Flags().StringSliceVar(&flags.EnvSettings, "env", nil, "set environment variable 'KEY=VALUE' (can be used multiple times)")
	cmd.Flags().StringVarP(&flags.BaseImage, "base-image", "b", "", "base image to use instead of default")
	cmd.Flags().StringSliceVar(&flags.BindPaths, "bind", nil, "bind path 'HOST:CONTAINER' (can be used multiple times)")
	cmd.Flags().BoolVarP(&flags.Fakeroot, "fakeroot", "f", false, "run container with fakeroot privileges")

	// Register completions
	cmd.RegisterFlagCompletionFunc("overlay", overlayFlagCompletion(true, true))
	cmd.RegisterFlagCompletionFunc("base-image", baseImageFlagCompletion())

	// Stop flag parsing after the first positional argument
	cmd.Flags().SetInterspersed(false)

	// Allow unknown flags to pass through (for apptainer)
	cmd.FParseErrWhitelist.UnknownFlags = true
}

// KnownFlags returns a map of all known condatainer flags
func KnownFlags() map[string]bool {
	return map[string]bool{
		"--overlay": true, "-o": true,
		"--writable": true, "--writable-img": true, "-w": true,
		"--env":        true,
		"--bind":       true,
		"--base-image": true, "-b": true,
		"--fakeroot": true, "-f": true,
		"--debug": true,
		"--local": true,
		"--quiet": true, "-q": true,
		"--yes": true, "-y": true,
	}
}

// Exit codes used by various commands
const (
	// Generic error code
	ExitCodeError = 1
	// Returned when jobs were submitted to a scheduler (overlays will be created asynchronously)
	ExitCodeJobsSubmitted = 3
)

// ExitIfJobsSubmitted exits the process with ExitCodeJobsSubmitted if a scheduler
// submission was requested and the graph shows job IDs were created.
func ExitIfJobsSubmitted(graph *build.BuildGraph) {
	if config.Global.SubmitJob && graph != nil && len(graph.GetJobIDs()) > 0 {
		utils.PrintMessage("%d scheduler job(s) submitted. exiting with code %d",
			len(graph.GetJobIDs()), ExitCodeJobsSubmitted)
		os.Exit(ExitCodeJobsSubmitted)
	}
}

// ExitWithError prints an error and exits with ExitCodeError
func ExitWithError(format string, a ...interface{}) {
	utils.PrintError(format, a...)
	os.Exit(ExitCodeError)
}

// ParseCommandArgs parses arguments from os.Args after a given subcommand
// Returns: commands (positional args), apptainerFlags (unknown flags for apptainer)
func ParseCommandArgs(subcommand string) ([]string, []string) {
	var commands, apptainerFlags []string

	// Find where subcommand appears in os.Args
	cmdIdx := -1
	for i, arg := range os.Args {
		if arg == subcommand {
			cmdIdx = i
			break
		}
	}

	if cmdIdx == -1 {
		// Fallback: no args found
		return commands, apptainerFlags
	}

	knownFlags := KnownFlags()

	// Parse from after subcommand in os.Args
	commandStarted := false
	for i := cmdIdx + 1; i < len(os.Args); i++ {
		arg := os.Args[i]

		// Once we hit the command, everything after is part of the command
		if commandStarted {
			commands = append(commands, arg)
			continue
		}

		// Skip known flags (already handled by cobra)
		if knownFlags[arg] || isKnownFlagWithEquals(knownFlags, arg) {
			// Value flags (space-separated)
			if knownFlags[arg] && needsValue(arg) && i+1 < len(os.Args) {
				i++ // Skip value
			}
			continue
		}

		// Unknown flag → must use --flag=value format for apptainer pass-through
		if strings.HasPrefix(arg, "-") {
			apptainerFlags = append(apptainerFlags, arg)
			continue
		}

		// First non-flag argument starts the command
		commandStarted = true
		commands = append(commands, arg)
	}

	return commands, apptainerFlags
}

// ParseInstanceArgs parses arguments for instance start/run commands
// Returns: instanceName, apptainerFlags
func ParseInstanceArgs(subcommand string) (string, []string) {
	var instanceName string
	var apptainerFlags []string

	// Find where "instance <subcommand>" appears in os.Args
	cmdIdx := -1
	for i := 0; i < len(os.Args)-1; i++ {
		if os.Args[i] == "instance" && os.Args[i+1] == subcommand {
			cmdIdx = i + 1 // Point to the subcommand
			break
		}
	}

	if cmdIdx == -1 {
		// Fallback: no args found
		return instanceName, apptainerFlags
	}

	knownFlags := KnownFlags()

	// Parse from after subcommand in os.Args
	for i := cmdIdx + 1; i < len(os.Args); i++ {
		arg := os.Args[i]

		// Skip known flags (already handled by cobra)
		if knownFlags[arg] || isKnownFlagWithEquals(knownFlags, arg) {
			// Value flags (space-separated)
			if knownFlags[arg] && needsValue(arg) && i+1 < len(os.Args) {
				i++ // Skip value
			}
			continue
		}

		// Unknown flag → pass through to apptainer
		if strings.HasPrefix(arg, "-") {
			apptainerFlags = append(apptainerFlags, arg)
			continue
		}

		// First (and only) non-flag argument is the instance name
		if instanceName == "" {
			instanceName = arg
		}
		// Ignore any additional positional arguments - overlays must use -o flag
	}

	return instanceName, apptainerFlags
}

func isKnownFlagWithEquals(knownFlags map[string]bool, arg string) bool {
	if !strings.Contains(arg, "=") {
		return false
	}
	// Check if prefix matches a known flag
	parts := strings.SplitN(arg, "=", 2)
	return knownFlags[parts[0]]
}

func needsValue(flag string) bool {
	valueFlags := map[string]bool{
		"-o": true, "--overlay": true,
		"--env": true, "--bind": true,
		"-b": true, "--base-image": true,
	}
	return valueFlags[flag]
}

// ensureBaseImage ensures the base image exists
func ensureBaseImage(ctx context.Context) error {
	return apptainer.EnsureBaseImage(ctx, false, false)
}

// ResolveBaseImage resolves a base image path to an absolute path
// Returns empty string if baseImage is empty
func ResolveBaseImage(baseImage string) string {
	if baseImage == "" {
		return ""
	}

	if utils.FileExists(baseImage) {
		return baseImage
	}

	resolvedBase, err := container.ResolveOverlayPaths([]string{baseImage})
	if err == nil && len(resolvedBase) > 0 {
		return resolvedBase[0]
	}

	return baseImage
}

// PrepareCommandAndHidePrompt prepares the command array and determines if prompt should be hidden
// If commands is empty, defaults to ["bash"] with hidePrompt=false
// If commands has 1 element, returns hidePrompt=false
// Otherwise returns hidePrompt=true
func PrepareCommandAndHidePrompt(commands []string) ([]string, bool) {
	hidePrompt := true
	if len(commands) == 0 {
		commands = []string{"bash"}
		hidePrompt = false
	} else if len(commands) == 1 {
		hidePrompt = false
	}
	return commands, hidePrompt
}

// ============================================================================
// Shell Completion Functions
// ============================================================================

// overlayFlagCompletion returns completion function for -o/--overlay flag
func overlayFlagCompletion(includeData bool, includeImg bool) func(*cobra.Command, []string, string) ([]string, cobra.ShellCompDirective) {
	return func(cmd *cobra.Command, args []string, toComplete string) ([]string, cobra.ShellCompDirective) {
		return overlaySuggestions(includeData, includeImg, toComplete)
	}
}

// baseImageFlagCompletion returns completion function for -b/--base-image flag
func baseImageFlagCompletion() func(*cobra.Command, []string, string) ([]string, cobra.ShellCompDirective) {
	return func(cmd *cobra.Command, args []string, toComplete string) ([]string, cobra.ShellCompDirective) {
		return systemOverlaySuggestions(toComplete)
	}
}

// systemOverlaySuggestions returns only system overlays (def-built) and local overlay files
func systemOverlaySuggestions(toComplete string) ([]string, cobra.ShellCompDirective) {
	installed, err := container.InstalledOverlays()
	if err != nil {
		return nil, cobra.ShellCompDirectiveError
	}

	choices := map[string]struct{}{}

	// Only include overlays from the def list (system/def-built overlays)
	defList := build.GetDefBuiltList()
	for name := range installed {
		normalized := utils.NormalizeNameVersion(name)
		if defList[normalized] {
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

// overlaySuggestions returns overlay suggestions including installed overlays and local files
func overlaySuggestions(includeData bool, includeImg bool, toComplete string) ([]string, cobra.ShellCompDirective) {
	installed, err := container.InstalledOverlays()
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

	// Add space after completions (all suggestions are files now, no folders)
	return suggestions, cobra.ShellCompDirectiveNoFileComp
}

// findLocalFilesWithFilter is a shared helper for finding local files recursively
// Includes directories for navigation and recursively finds files up to maxDepth
func findLocalFilesWithFilter(toComplete string, maxDepth int, fileFilter func(name string) bool) []string {
	pathDir, _ := filepath.Split(toComplete)
	dirForRead := pathDir
	if dirForRead == "" {
		dirForRead = "."
	}

	suggestions := []string{}

	// Recursively find files up to maxDepth
	var findFiles func(dir string, prefix string, currentDepth int)
	findFiles = func(dir string, prefix string, currentDepth int) {
		entries, err := os.ReadDir(dir)
		if err != nil {
			return
		}

		for _, entry := range entries {
			name := entry.Name()

			// Skip hidden files/directories (starting with .)
			if strings.HasPrefix(name, ".") {
				continue
			}

			candidate := name
			if prefix != "" {
				candidate = prefix + name
			}

			if entry.IsDir() {
				// Recurse into subdirectories if within depth limit (don't add directories to suggestions)
				if currentDepth < maxDepth {
					findFiles(filepath.Join(dir, name), candidate+"/", currentDepth+1)
				}
			} else {
				// Apply file filter
				if fileFilter(name) {
					if toComplete == "" || strings.HasPrefix(candidate, toComplete) {
						suggestions = append(suggestions, candidate)
					}
				}
			}
		}
	}

	findFiles(dirForRead, pathDir, 0)
	return suggestions
}

// localImageSuggestions returns local .sqf and .sif files (for -b flag)
// Includes directories for navigation and recursively finds files up to 1 level deep
func localImageSuggestions(toComplete string) []string {
	return findLocalFilesWithFilter(toComplete, 1, func(name string) bool {
		// Include .sqf and .sif files only (not .img)
		return utils.IsSqf(name) || utils.IsSif(name)
	})
}

// localOverlaySuggestions returns local overlay files (for -o flag)
// Includes directories for navigation and recursively finds files up to 1 level deep
func localOverlaySuggestions(toComplete string, includeImg bool) []string {
	return findLocalFilesWithFilter(toComplete, 1, func(name string) bool {
		// Filter based on file type
		if !utils.IsOverlay(name) {
			return false
		}

		// Filter .sif and .img based on includeImg parameter
		if utils.IsSif(name) && includeImg {
			// Skip .sif files when includeImg is true (for -o flag)
			return false
		}
		if utils.IsImg(name) && !includeImg {
			// Skip .img files when includeImg is false (for -b flag)
			return false
		}

		return true
	})
}

// isAppOverlay checks if an overlay path is considered an "app" overlay
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

	// Check def list from build package
	defList := build.GetDefBuiltList()
	return defList[normalized]
}

// pathWithinDir checks if a path is within a directory
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
