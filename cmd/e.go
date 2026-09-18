package cmd

import (
	"context"
	"errors"
	"os"
	"path/filepath"
	"slices"
	"strings"

	"github.com/Justype/condatainer/cmd/internal/ui"
	"github.com/Justype/condatainer/internal/helper"
	"github.com/Justype/condatainer/internal/image"
	"github.com/Justype/condatainer/internal/runtime/apptainer"
	"github.com/Justype/condatainer/internal/runtime/container"
	"github.com/Justype/condatainer/internal/runtime/exec"
	"github.com/Justype/condatainer/internal/utils"
	"github.com/spf13/cobra"
)

var (
	eReadOnly    bool
	eNoAutoload  bool
	eFakeroot    bool
	eActivation  string
	eEnvSettings []string
	eBindPaths   []string
	eProjectDir  string
)

// eCmd is a quick shortcut for executing commands with overlays
var eCmd = &cobra.Command{
	Use:   "e [flags] [overlays...] [-- command...]",
	Short: "Shortcut for exec with overlays, writable by default",
	Long: `Quick shortcut for executing commands with overlays.

Overlays and flags go before --, commands go after --.
- Writable by default
- Auto-loads env.img from the current directory (unless a .img is provided or -n).
  If only its frozen snapshot (env.sqf) exists, that is loaded instead, read-only.

Note: Additional Apptainer flags must use --flag=value format (no space)`,
	Example: `  condatainer e                # Auto-load env.img if present, run bash
  condatainer e testing.img    # Use a specific writable overlay
  condatainer e samtools/1.22 bcftools/1.20   # Multiple overlays
  condatainer e samtools/1.22 -- samtools     # Run command
  condatainer e --home=/custom samtools/1.22  # Pass apptainer flags`,
	SilenceUsage: true,
	RunE:         runE,
}

func init() {
	rootCmd.AddCommand(eCmd)

	eCmd.Flags().BoolVarP(&eReadOnly, "read-only", "r", false, "Mount .img overlays as read-only (default: writable)")
	eCmd.Flags().BoolVarP(&eNoAutoload, "no-autoload", "n", false, "Disable auto-loading env.img from current directory")
	eCmd.Flags().BoolVarP(&eFakeroot, "fakeroot", "f", false, "Run container with fakeroot privileges")
	eCmd.Flags().StringVar(&eActivation, "activation", "all", "Which activate.d scripts to source before the command: all, env, or none")
	eCmd.Flags().StringSliceVar(&eEnvSettings, "env", nil, "Set environment variable 'KEY=VALUE' (repeatable)")
	eCmd.Flags().StringSliceVar(&eBindPaths, "bind", nil, "Bind path 'HOST:CONTAINER' (repeatable)")
	RegisterProjectFlags(eCmd, &eProjectDir)

	// Allow flags to be interspersed with overlays
	eCmd.Flags().SetInterspersed(true)

	// Allow unknown flags to pass through (for apptainer)
	eCmd.FParseErrWhitelist.UnknownFlags = true

	// Completion
	eCmd.ValidArgsFunction = func(cmd *cobra.Command, args []string, toComplete string) ([]string, cobra.ShellCompDirective) {
		// During completion, check if -- appears in os.Args
		// This is more reliable than checking the processed args array
		if slices.Contains(os.Args, "--") {
			// After --, use default file completion
			return nil, cobra.ShellCompDirectiveDefault
		}
		// Before --, complete overlays
		return overlaySuggestions(true, true, toComplete)
	}
}

func runE(cmd *cobra.Command, args []string) error {
	if err := applyProjectRelocation(eProjectDir); err != nil {
		return err
	}

	// STRICT parsing: everything before -- is overlay/flag, after is command
	overlays, commands, apptainerFlags, err := parseEArgs(args)
	if err != nil {
		return err
	}

	// Auto-load env.img if not disabled and no .img in overlays
	if !eNoAutoload {
		hasImgOverlay := slices.ContainsFunc(overlays, utils.IsImg)
		if !hasImgOverlay {
			if pwd, err := os.Getwd(); err == nil {
				switch candidate := utils.FindEnvOverlay("", pwd); {
				case candidate != "":
					// The probe asks for the lock the mount will take, so a
					// writable session is not told an overlay is free that it
					// then cannot have.
					err := image.CheckAvailable(candidate, !eReadOnly)
					switch {
					case errors.Is(err, image.ErrProtected):
						utils.PrintWarning("%s is write-protected, running without it; -r mounts it read-only",
							filepath.Base(candidate))
					case err != nil:
						utils.PrintWarning("%s is in use, running without it", filepath.Base(candidate))
					default:
						utils.PrintNote("Autoload environment overlay at %s", utils.StylePath(candidate))
						overlays = append(overlays, candidate)
					}
				default:
					// No .img: fall back to its read-only snapshot, if any.
					if snapshot := helper.FindEnvSnapshot(pwd); snapshot != "" {
						utils.PrintNote("No env.img here; auto-loading its snapshot %s, read-only. `overlay create env.img` for a writable copy.",
							utils.StylePath(snapshot))
						overlays = append(overlays, snapshot)
					}
				}
			}
		}
	}

	// Prepare command and determine if prompt should be hidden
	commands, hidePrompt := PrepareCommandAndHidePrompt(commands)

	// Resolve overlays
	overlays, err = projectOverlays(cmd.Context(), overlays)
	if err != nil {
		return err
	}
	resolvedOverlays, err := container.ResolveOverlayPaths(overlays)
	if err != nil {
		return err
	}

	baseImageResolved, err := ensureRootBaseImage(cmd.Context(), resolvedOverlays)
	if err != nil {
		return err
	}

	activation, err := parseActivation(eActivation)
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
		HidePrompt:     hidePrompt,
		Activation:     activation,
	}

	plan, err := exec.Prepare(cmd.Context(), options)
	if err != nil {
		return err
	}
	ui.RenderExecPlan(plan)
	if err := exec.RunPrepared(cmd.Context(), plan, exec.IO{Stdin: os.Stdin, Stdout: os.Stdout, Stderr: os.Stderr}); err != nil {
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
		"--env":      true,
		"--bind":     true,
		"--fakeroot": true, "-f": true,
		"--activation": true,
		"--debug":      true,
		"--quiet":      true, "-q": true,
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
			if (arg == "--env" || arg == "--bind" || arg == "--activation") && i+1 < len(os.Args) {
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
