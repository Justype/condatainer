package cmd

import (
	"context"
	"errors"
	"os"
	"time"

	"github.com/Justype/condatainer/cmd/internal/ui"
	"github.com/Justype/condatainer/internal/runtime/apptainer"
	"github.com/Justype/condatainer/internal/runtime/container"
	"github.com/Justype/condatainer/internal/runtime/exec"
	"github.com/spf13/cobra"
)

type execCommand struct {
	cobra.Command
}

var execFlags CommonFlags
var execGpuRequested bool
var execActivation string
var execProjectDir string
var execStopGrace int

var execCmd = &execCommand{
	Command: cobra.Command{
		Use:          "exec [flags] [command...]",
		Short:        "Execute a command with overlays",
		SilenceUsage: true,
		Long: `Execute a command inside a container with specified overlays.

- Overlays are given with -o/--overlay
- All positional arguments form the command (default: bash)

Note: Additional Apptainer flags must use --flag=value format (no space)`,
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
	execCmd.Flags().BoolVar(&execGpuRequested, "gpu", false, "Force GPU flags (--nv/--rocm) even if autoload_gpu is disabled")
	execCmd.Flags().StringVar(&execActivation, "activation", "all", "Which activate.d scripts to source before the command: all, env, or none")
	execCmd.Flags().IntVar(&execStopGrace, "stop-grace", 0, "Seconds the container gets to exit after this process is told to stop")
	execCmd.Flags().MarkHidden("stop-grace") //nolint:errcheck
	RegisterProjectFlags(&execCmd.Command, &execProjectDir)

	// For 'exec': use default file completion for positional args
	execCmd.ValidArgsFunction = func(cmd *cobra.Command, args []string, toComplete string) ([]string, cobra.ShellCompDirective) {
		return nil, cobra.ShellCompDirectiveDefault
	}
}

func runExec(cmd *cobra.Command, args []string) error {
	if execHelpRequested(args) {
		return cmd.Help()
	}
	if err := applyProjectRelocation(execProjectDir); err != nil {
		return err
	}

	ResolveFlagAlias(cmd, "writable", "writable-img")

	// Parse arguments - treat all positional args as commands
	commandFinal, apptainerFlags := ParseCommandArgs("exec")

	// Use overlays from -o flag
	overlayFinal := execFlags.Overlays

	// Prepare command and determine if prompt should be hidden
	commandFinal, hidePrompt := PrepareCommandAndHidePrompt(commandFinal)

	overlayFinal, err := projectOverlays(cmd.Context(), overlayFinal)
	if err != nil {
		return err
	}
	resolvedOverlays, err := container.ResolveOverlayPaths(overlayFinal)
	if err != nil {
		return err
	}

	baseImageResolved, err := ensureRootBaseImage(cmd.Context(), resolvedOverlays)
	if err != nil {
		return err
	}

	resolvedOverlays, bindLibexec, err := nestedRun(cmd.Context(), resolvedOverlays)
	if err != nil {
		return err
	}

	activation, err := parseActivation(execActivation)
	if err != nil {
		return err
	}

	options := exec.Options{
		Overlays:       resolvedOverlays,
		BindLibexec:    bindLibexec,
		Command:        commandFinal,
		WritableImg:    execFlags.WritableImg,                                          // Default false, unless -w specified
		EnvSettings:    append(liveJobResourceEnvSettings(), execFlags.EnvSettings...), // Inject live job resources, then user env vars
		BindPaths:      execFlags.BindPaths,
		ApptainerFlags: apptainerFlags, // Pass through unknown flags
		Fakeroot:       execFlags.Fakeroot,
		BaseImage:      baseImageResolved,
		HidePrompt:     hidePrompt,
		GpuRequested:   execGpuRequested,
		Activation:     activation,
		StopGrace:      time.Duration(execStopGrace) * time.Second,
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
