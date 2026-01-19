package cmd

import (
	"fmt"
	"os"
	"path/filepath"

	"github.com/spf13/cobra"

	"github.com/Justype/condatainer/internal/config"
	"github.com/Justype/condatainer/internal/exec"
	"github.com/Justype/condatainer/internal/utils"
)

var (
	// exec command flags
	execOverlays    []string
	execKeep        bool
	execWritableImg bool
	execEnvSettings []string
	execBaseImage   string
	execBindPaths   []string
	execFakeroot    bool

	// e command flags
	eReadOnly    bool
	eBaseImage   string
	eNoAutoload  bool
	eFakeroot    bool
)

// execCmd represents the exec command
var execCmd = &cobra.Command{
	Use:   "exec [flags] [overlays...] [--] [command...]",
	Short: "Execute a command using overlays",
	Long: `Execute a command inside a container with specified overlays.

The exec command can intelligently parse arguments to separate overlays from commands.
Use -o/--overlay to explicitly specify overlays, or let the command parse them automatically.

Examples:
  # Run bash with numpy overlay
  condatainer exec numpy/1.24

  # Run python script with multiple overlays
  condatainer exec numpy/1.24 scipy/1.10 -- python script.py

  # Explicit overlay specification
  condatainer exec -o numpy/1.24 -o scipy/1.10 -- python script.py

  # Use writable .img overlay
  condatainer exec -w env.img bash

  # Set environment variables
  condatainer exec --env MYVAR=value numpy/1.24 bash`,
	RunE: runExec,
}

// eCmd represents the quick exec (e) command
var eCmd = &cobra.Command{
	Use:   "e [flags] [overlays...]",
	Short: "Run bash using writable overlays (quick exec)",
	Long: `Quick bash execution with writable .img overlays.

The 'e' command is a shortcut for interactive shell sessions with overlays.
By default, .img overlays are mounted as writable, and 'env.img' in the current
directory is automatically loaded if present.

Examples:
  # Run bash with writable env.img (auto-loaded from current dir)
  condatainer e

  # Run bash with specific overlays
  condatainer e myenv.img numpy/1.24

  # Disable autoloading
  condatainer e -n numpy/1.24

  # Read-only mode
  condatainer e -r env.img`,
	RunE: runE,
}

func init() {
	rootCmd.AddCommand(execCmd)
	rootCmd.AddCommand(eCmd)

	// exec command flags
	execCmd.Flags().StringSliceVarP(&execOverlays, "overlay", "o", nil, "overlay file to mount (can be used multiple times)")
	execCmd.Flags().BoolVarP(&execKeep, "keep", "k", false, "disable automatic command parsing to installed overlays")
	execCmd.Flags().BoolVarP(&execWritableImg, "writable-img", "w", false, "mount .img overlays as writable (default: read-only)")
	execCmd.Flags().StringSliceVar(&execEnvSettings, "env", nil, "set environment variable 'KEY=VALUE' (can be used multiple times)")
	execCmd.Flags().StringVarP(&execBaseImage, "base-image", "b", "", "base image to use instead of default")
	execCmd.Flags().StringSliceVar(&execBindPaths, "bind", nil, "bind path 'HOST:CONTAINER' (can be used multiple times)")
	execCmd.Flags().BoolVar(&execFakeroot, "fakeroot", false, "run container with fakeroot privileges")

	// e command flags
	eCmd.Flags().BoolVarP(&eReadOnly, "read-only", "r", false, "mount .img overlays as read-only (default: writable)")
	eCmd.Flags().StringVarP(&eBaseImage, "base-image", "b", "", "base image to use instead of default")
	eCmd.Flags().BoolVarP(&eNoAutoload, "no-autoload", "n", false, "disable autoloading 'env.img' from current directory")
	eCmd.Flags().BoolVar(&eFakeroot, "fakeroot", false, "run container with fakeroot privileges")
}

func runExec(cmd *cobra.Command, args []string) error {
	// Ensure base image exists
	if err := ensureBaseImage(); err != nil {
		return err
	}

	// Change base image if specified
	if execBaseImage != "" {
		config.Global.BaseImage = execBaseImage
	}

	// Get installed overlays for name resolution
	installedOverlays, err := exec.InstalledOverlays()
	if err != nil {
		return err
	}

	overlayFinal := []string{}
	commandFinal := []string{}

	// Parse arguments to separate overlays from commands
	// If --keep is not set and no explicit overlays, try to parse args
	if !execKeep && len(execOverlays) == 0 && len(args) > 0 {
		utils.PrintDebug("[EXEC] Parsing commands to separate overlays and command...")
		for _, arg := range args {
			if utils.IsOverlay(arg) {
				// It's an overlay file
				overlayFinal = append(overlayFinal, arg)
			} else if _, ok := installedOverlays[utils.NormalizeNameVersion(arg)]; ok {
				// It's an installed overlay name
				utils.PrintWarning("Convert command %s to overlay", utils.StyleName(arg))
				overlayFinal = append(overlayFinal, arg)
			} else {
				// It's part of the command
				commandFinal = append(commandFinal, arg)
			}
		}
	} else {
		// Keep mode or explicit overlays: all args are command
		commandFinal = args
	}

	// Add explicit overlays
	overlayFinal = append(overlayFinal, execOverlays...)

	// Default command is bash
	if len(commandFinal) == 0 {
		commandFinal = []string{"bash"}
	}

	// Resolve overlay names to paths
	resolvedOverlays, err := exec.ResolveOverlayPaths(overlayFinal)
	if err != nil {
		return err
	}

	// Run the command
	options := exec.Options{
		Overlays:      resolvedOverlays,
		Command:       commandFinal,
		WritableImg:   execWritableImg,
		EnvSettings:   execEnvSettings,
		BindPaths:     execBindPaths,
		Fakeroot:      execFakeroot,
		BaseImage:     config.Global.BaseImage,
		ApptainerBin:  config.Global.ApptainerBin,
	}

	return exec.Run(options)
}

func runE(cmd *cobra.Command, args []string) error {
	// Ensure base image exists
	if err := ensureBaseImage(); err != nil {
		return err
	}

	// Change base image if specified
	if eBaseImage != "" {
		config.Global.BaseImage = eBaseImage
	}

	overlayFinal := args

	// Autoload env.img if no .img overlay is specified
	hasImgOverlay := false
	for _, overlay := range overlayFinal {
		if utils.IsImg(overlay) {
			hasImgOverlay = true
			break
		}
	}

	if !hasImgOverlay && !eNoAutoload {
		pwd, err := os.Getwd()
		if err == nil {
			localEnvPath := filepath.Join(pwd, "env.img")
			if utils.FileExists(localEnvPath) {
				utils.PrintMessage("Autoload env.img at %s", utils.StylePath(localEnvPath))
				overlayFinal = append(overlayFinal, localEnvPath)
			}
		}
	}

	// Resolve overlay names to paths
	resolvedOverlays, err := exec.ResolveOverlayPaths(overlayFinal)
	if err != nil {
		return err
	}

	// Run bash
	options := exec.Options{
		Overlays:     resolvedOverlays,
		Command:      []string{"bash"},
		WritableImg:  !eReadOnly, // e command defaults to writable
		Fakeroot:     eFakeroot,
		BaseImage:    config.Global.BaseImage,
		ApptainerBin: config.Global.ApptainerBin,
	}

	return exec.Run(options)
}

// ensureBaseImage checks if the base image exists
func ensureBaseImage() error {
	if config.Global.BaseImage == "" {
		return fmt.Errorf("base image not configured")
	}
	if !utils.FileExists(config.Global.BaseImage) {
		return fmt.Errorf("base image not found at %s", config.Global.BaseImage)
	}
	return nil
}
