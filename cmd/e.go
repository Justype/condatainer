package cmd

import (
	"context"
	"errors"
	"os"
	"path/filepath"
	"strings"

	"github.com/Justype/condatainer/internal/apptainer"
	"github.com/Justype/condatainer/internal/config"
	"github.com/Justype/condatainer/internal/container"
	"github.com/Justype/condatainer/internal/exec"
	"github.com/Justype/condatainer/internal/utils"
	"github.com/spf13/cobra"
)

var (
	eReadOnly    bool
	eNoAutoload  bool
	eBaseImage   string
	eFakeroot    bool
	eEnvSettings []string
	eBindPaths   []string
)

// eCmd is a quick shortcut for executing commands with overlays
var eCmd = &cobra.Command{
	Use:   "e [flags] [overlays...] [-- command...]",
	Short: "Shortcut for exec with overlays, writable by default",
	Long: `Quick shortcut for executing commands with overlays.

Overlays and flags go before --, commands go after --.
If no -- is provided, all positional arguments are treated as overlays.
- Writable by default (use -r for read-only)
- Auto-loads env.img (unless a different .img is provided or -n is specified)
- Defaults to bash if no command specified

Note: Additional Apptainer flags must use --flag=value format (no space)`,
	Example: `  condatainer e               # Auto-load env.img if present, run bash
  condatainer e testing.img   # Use a specific writable overlay

  # Multiple overlays, run bash (auto loads env.img)
  condatainer e samtools/1.22 bcftools/1.20

  # Run specific command (requires --)
  condatainer e samtools/1.22 -- samtools view file.bam

  # Pass apptainer flags
  condatainer e --home=/custom samtools/1.22`,
	SilenceUsage: true,
	RunE:         runE,
}

func init() {
	rootCmd.AddCommand(eCmd)

	eCmd.Flags().BoolVarP(&eReadOnly, "read-only", "r", false, "Mount .img overlays as read-only (default: writable)")
	eCmd.Flags().BoolVarP(&eNoAutoload, "no-autoload", "n", false, "Disable auto-loading env.img from current directory")
	eCmd.Flags().StringVarP(&eBaseImage, "base-image", "b", "", "Base image to use instead of default")
	eCmd.Flags().BoolVarP(&eFakeroot, "fakeroot", "f", false, "Run container with fakeroot privileges")
	eCmd.Flags().StringSliceVar(&eEnvSettings, "env", nil, "Set environment variable 'KEY=VALUE' (can be used multiple times)")
	eCmd.Flags().StringSliceVar(&eBindPaths, "bind", nil, "Bind path 'HOST:CONTAINER' (can be used multiple times)")

	// Allow flags to be interspersed with overlays
	eCmd.Flags().SetInterspersed(true)

	// Allow unknown flags to pass through (for apptainer)
	eCmd.FParseErrWhitelist.UnknownFlags = true

	// Completion
	eCmd.RegisterFlagCompletionFunc("base-image", baseImageFlagCompletion())
	eCmd.ValidArgsFunction = func(cmd *cobra.Command, args []string, toComplete string) ([]string, cobra.ShellCompDirective) {
		// During completion, check if -- appears in os.Args
		// This is more reliable than checking the processed args array
		for _, arg := range os.Args {
			if arg == "--" {
				// After --, use default file completion
				return nil, cobra.ShellCompDirectiveDefault
			}
		}
		// Before --, complete overlays
		return overlaySuggestions(true, true, toComplete)
	}
}

func runE(cmd *cobra.Command, args []string) error {
	if err := ensureBaseImage(cmd.Context()); err != nil {
		return err
	}

	// STRICT parsing: everything before -- is overlay/flag, after is command
	overlays, commands, apptainerFlags, err := parseEArgs(args)
	if err != nil {
		return err
	}

	// Auto-load env.img if not disabled and no .img in overlays
	if !eNoAutoload {
		hasImgOverlay := false
		for _, overlay := range overlays {
			if utils.IsImg(overlay) {
				hasImgOverlay = true
				break
			}
		}
		if !hasImgOverlay {
			if pwd, err := os.Getwd(); err == nil {
				localEnvPath := filepath.Join(pwd, "env.img")
				if utils.FileExists(localEnvPath) {
					utils.PrintNote("Autoload env.img at %s", utils.StylePath(localEnvPath))
					overlays = append(overlays, localEnvPath)
				}
			}
		}
	}

	// Prepare command and determine if prompt should be hidden
	commands, hidePrompt := PrepareCommandAndHidePrompt(commands)

	// Resolve base image if provided
	baseImageResolved := ResolveBaseImage(eBaseImage)

	// Resolve overlays
	resolvedOverlays, err := container.ResolveOverlayPaths(overlays)
	if err != nil {
		return err
	}

	// Build options
	options := exec.Options{
		Overlays:       resolvedOverlays,
		Command:        commands,
		WritableImg:    !eReadOnly,                                            // Default writable unless -r specified
		EnvSettings:    append(liveJobResourceEnvSettings(), eEnvSettings...), // Inject live job resources, then user env vars
		BindPaths:      eBindPaths,
		ApptainerFlags: apptainerFlags,
		Fakeroot:       eFakeroot,
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

func parseEArgs(args []string) (overlays, commands, apptainerFlags []string, err error) {
	// Parse os.Args directly to catch unknown flags that Cobra filtered out
	eIdx := -1
	for i, arg := range os.Args {
		if arg == "e" {
			eIdx = i
			break
		}
	}

	knownFlags := map[string]bool{
		"--read-only": true, "-r": true,
		"--no-autoload": true, "-n": true,
		"--env":        true,
		"--bind":       true,
		"--base-image": true, "-b": true,
		"--fakeroot": true, "-f": true,
		"--debug": true,
		"--local": true,
		"--quiet": true, "-q": true,
		"--yes": true, "-y": true,
	}

	commandStarted := false

	if eIdx == -1 {
		// Fallback: no os.Args parsing, just return cobra's positional args as overlays
		for _, arg := range args {
			if arg == "--" {
				commandStarted = true
				continue
			}
			if commandStarted {
				commands = append(commands, arg)
				continue
			}
			overlays = append(overlays, arg)
		}
		return overlays, commands, apptainerFlags, nil
	}

	for i := eIdx + 1; i < len(os.Args); i++ {
		arg := os.Args[i]

		if arg == "--" {
			commandStarted = true
			continue
		}

		if commandStarted {
			// Everything after -- is command
			commands = append(commands, arg)
			continue
		}

		// Before -- (overlay/flag section)

		// Skip known flags (already handled by cobra)
		if knownFlags[arg] || isKnownFlagWithEquals(knownFlags, arg) {
			// Check if flag needs value
			if (arg == "-b" || arg == "--base-image" || arg == "--env" || arg == "--bind") && i+1 < len(os.Args) {
				i++ // Skip value
			}
			continue
		}

		// Unknown flag → must use --flag=value format
		if strings.HasPrefix(arg, "-") {
			if strings.Contains(arg, "=") {
				// Has value: --home=/path
				apptainerFlags = append(apptainerFlags, arg)
			} else {
				// No value: boolean flag like --nv
				apptainerFlags = append(apptainerFlags, arg)
			}
			continue
		}

		// Not a flag, before -- → overlay
		overlays = append(overlays, arg)
	}

	return overlays, commands, apptainerFlags, nil
}
