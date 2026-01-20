package cmd

import (
	"fmt"
	"os"
	"path/filepath"
	"sort"
	"strings" // <--- ADDED: Required for string manipulation

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
	execCmd.RegisterFlagCompletionFunc("overlay", func(cmd *cobra.Command, args []string, toComplete string) ([]string, cobra.ShellCompDirective) {
		// Use InstalledOverlays so we include both images and ref-images and reuse existing normalization.
		installed, err := exec.InstalledOverlays()
		if err != nil {
			// If we can't list overlays, indicate an error to the shell completion system
			return nil, cobra.ShellCompDirectiveError
		}

		suggestions := []string{}
		for name := range installed {
			if strings.HasPrefix(name, toComplete) {
				suggestions = append(suggestions, name)
			}
		}
		sort.Strings(suggestions)
		if len(suggestions) == 0 {
			// No matches — allow shell to perform default file completion
			return nil, cobra.ShellCompDirectiveDefault
		}

		// Return suggestions and tell Shell NOT to fall back to standard file completion
		return suggestions, cobra.ShellCompDirectiveNoFileComp
	})
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

	if execBaseImage != "" {
		config.Global.BaseImage = execBaseImage
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
		overlayFinal = append(overlayFinal, args...)
		commandFinal = []string{"bash"}
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
	if config.Global.BaseImage == "" {
		return fmt.Errorf("base image not configured")
	}
	if !utils.FileExists(config.Global.BaseImage) {
		return fmt.Errorf("base image not found at %s", config.Global.BaseImage)
	}
	return nil
}
