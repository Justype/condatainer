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
	"github.com/Justype/condatainer/internal/container"
	"github.com/Justype/condatainer/internal/exec"
	"github.com/Justype/condatainer/internal/utils"
	"github.com/spf13/cobra"
)

type execCommand struct {
	cobra.Command
}

func (c *execCommand) InitDefaultHelpFlag() {}

var execFlags CommonFlags

var execCmd = &execCommand{
	Command: cobra.Command{
		Use:          "exec [flags] [overlays...] [command...]",
		Short:        "Execute a command using overlays",
		SilenceUsage: true,
		Long: `Execute a command inside a container with specified overlays.

- Use -o/--overlay to explicitly specify overlays. All positional arguments are treated as commands.
- If command is omitted, it defaults to bash.
- Additional Apptainer flags can be passed using --flag=value format (no space)`,
		Example: `  # Run samtools command with overlay
  condatainer exec -o samtools/1.22 samtools view file.bam

  # Use writable .img overlay
  condatainer exec -w -o env.img

  # Set environment variables
  condatainer exec --env MYVAR=value -o samtools/1.22 bash

  # Pass apptainer flags (use --flag=value format)
  condatainer exec --nv --home=/custom -o samtools/1.22 python gpu_script.py`,
		RunE: runExec,
	},
}

func init() {
	rootCmd.AddCommand(&execCmd.Command)

	// Register common flags
	RegisterCommonFlags(&execCmd.Command, &execFlags)

	// For 'exec': use default file completion for positional args
	execCmd.ValidArgsFunction = func(cmd *cobra.Command, args []string, toComplete string) ([]string, cobra.ShellCompDirective) {
		return nil, cobra.ShellCompDirectiveDefault
	}
}

func runExec(cmd *cobra.Command, args []string) error {
	if execHelpRequested(args) {
		return cmd.Help()
	}

	if err := ensureBaseImage(cmd.Context()); err != nil {
		return err
	}

	// Parse arguments - treat all positional args as commands
	commandFinal, apptainerFlags := ParseCommandArgs("exec")

	// Use overlays from -o flag
	overlayFinal := execFlags.Overlays

	hidePrompt := true
	if len(commandFinal) == 0 {
		commandFinal = []string{"bash"}
		hidePrompt = false
	} else if len(commandFinal) == 1 {
		hidePrompt = false
	}

	// Resolve user-specified base image if provided
	var baseImageResolved string
	if execFlags.BaseImage != "" {
		if utils.FileExists(execFlags.BaseImage) {
			baseImageResolved = execFlags.BaseImage
		} else {
			resolvedBase, err := container.ResolveOverlayPaths([]string{execFlags.BaseImage})
			if err == nil && len(resolvedBase) > 0 {
				baseImageResolved = resolvedBase[0]
			} else {
				baseImageResolved = execFlags.BaseImage
			}
		}
	}

	resolvedOverlays, err := container.ResolveOverlayPaths(overlayFinal)
	if err != nil {
		return err
	}

	options := exec.Options{
		Overlays:       resolvedOverlays,
		Command:        commandFinal,
		WritableImg:    execFlags.WritableImg, // Default false, unless -w specified
		EnvSettings:    execFlags.EnvSettings,
		BindPaths:      execFlags.BindPaths,
		ApptainerFlags: apptainerFlags, // Pass through unknown flags
		Fakeroot:       execFlags.Fakeroot,
		BaseImage:      baseImageResolved,
		ApptainerBin:   config.Global.ApptainerBin,
		HidePrompt:     hidePrompt,
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

	// Always return NoFileComp to prevent default file listing
	return suggestions, cobra.ShellCompDirectiveNoFileComp
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

	// Check def list from build package
	defList := build.GetDefBuiltList()
	return defList[normalized]
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
