package cmd

import (
	"os"
	"path/filepath"
	"sort"
	"strings"

	"github.com/Justype/condatainer/internal/apptainer"
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
	execCmd.Flags().BoolVarP(&execWritableImg, "writable-img", "w", false, "mount .img overlays as writable (default: read-only)")
	execCmd.Flags().StringSliceVar(&execEnvSettings, "env", nil, "set environment variable 'KEY=VALUE' (can be used multiple times)")
	execCmd.Flags().StringVarP(&execBaseImage, "base-image", "b", "", "base image to use instead of default")
	execCmd.Flags().StringSliceVar(&execBindPaths, "bind", nil, "bind path 'HOST:CONTAINER' (can be used multiple times)")
	execCmd.Flags().BoolVar(&execFakeroot, "fakeroot", false, "run container with fakeroot privileges")
	execCmd.Flags().BoolVarP(&eReadOnly, "read-only", "r", false, "mount .img overlays as read-only (only applies when using the 'e' shortcut)")
	execCmd.Flags().BoolVarP(&eNoAutoload, "no-autoload", "n", false, "disable autoloading 'env.img' from current directory (only applies to the 'e' shortcut)")

	// --- NEW: Dynamic Completion Logic ---
	execCmd.RegisterFlagCompletionFunc("overlay", overlayFlagCompletion(true))
	execCmd.RegisterFlagCompletionFunc("base-image", overlayFlagCompletion(false))
	execCmd.ValidArgsFunction = positionalOverlayCompletion(true)
	// -------------------------------------

	// Stop flag parsing after the first positional argument so tool arguments (like --help) bypass exec.
	execCmd.Flags().SetInterspersed(false)
}

func runExec(cmd *cobra.Command, args []string) error {
	if execHelpRequested(args) {
		return cmd.Help()
	}

	if err := ensureBaseImage(); err != nil {
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

	if execBaseImage != "" {
		config.Global.BaseImage = execBaseImage
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
		BaseImage:    config.Global.BaseImage,
		ApptainerBin: config.Global.ApptainerBin,
	}

	return exec.Run(options)
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

func ensureBaseImage() error {
	return apptainer.EnsureBaseImage()
}

func overlayFlagCompletion(includeData bool) func(*cobra.Command, []string, string) ([]string, cobra.ShellCompDirective) {
	return func(cmd *cobra.Command, args []string, toComplete string) ([]string, cobra.ShellCompDirective) {
		return overlaySuggestions(includeData, toComplete)
	}
}

func positionalOverlayCompletion(includeData bool) func(*cobra.Command, []string, string) ([]string, cobra.ShellCompDirective) {
	return func(cmd *cobra.Command, args []string, toComplete string) ([]string, cobra.ShellCompDirective) {
		if containsDashDash(args) {
			return nil, cobra.ShellCompDirectiveDefault
		}
		return overlaySuggestions(includeData, toComplete)
	}
}

func overlaySuggestions(includeData bool, toComplete string) ([]string, cobra.ShellCompDirective) {
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

	for _, candidate := range localOverlaySuggestions(toComplete) {
		choices[candidate] = struct{}{}
	}

	if len(choices) == 0 {
		return nil, cobra.ShellCompDirectiveDefault
	}

	suggestions := make([]string, 0, len(choices))
	for choice := range choices {
		suggestions = append(suggestions, choice)
	}
	sort.Strings(suggestions)
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

func localOverlaySuggestions(toComplete string) []string {
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
		if entry.IsDir() || !utils.IsOverlay(entry.Name()) {
			continue
		}
		candidate := entry.Name()
		if pathDir != "" {
			candidate = pathDir + entry.Name()
		}
		if toComplete != "" && !strings.HasPrefix(candidate, toComplete) {
			continue
		}
		suggestions = append(suggestions, candidate)
	}
	return suggestions
}

func isAppOverlay(path string) bool {
	return pathWithinDir(path, config.Global.ImagesDir)
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
