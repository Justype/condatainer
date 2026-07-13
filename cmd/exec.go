package cmd

import (
	"context"
	"errors"
	"os"

	"github.com/Justype/condatainer/cmd/internal/ui"
	"github.com/Justype/condatainer/internal/apptainer"
	"github.com/Justype/condatainer/internal/config"
	"github.com/Justype/condatainer/internal/container"
	"github.com/Justype/condatainer/internal/exec"
	"github.com/spf13/cobra"
)

type execCommand struct {
	cobra.Command
}

var execFlags CommonFlags

var execCmd = &execCommand{
	Command: cobra.Command{
		Use:          "exec [flags] [command...]",
		Short:        "Execute a command using overlays",
		SilenceUsage: true,
		Long: `Execute a command inside a container with specified overlays.

Overlays are given with -o/--overlay; all positional arguments form
the command (default: bash).

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

	// Prepare command and determine if prompt should be hidden
	commandFinal, hidePrompt := PrepareCommandAndHidePrompt(commandFinal)

	// Resolve user-specified base image if provided
	baseImageResolved := ResolveBaseImage(execFlags.BaseImage)

	resolvedOverlays, err := container.ResolveOverlayPaths(overlayFinal)
	if err != nil {
		return err
	}

	options := exec.Options{
		Overlays:       resolvedOverlays,
		Command:        commandFinal,
		WritableImg:    execFlags.WritableImg,                                          // Default false, unless -w specified
		EnvSettings:    append(liveJobResourceEnvSettings(), execFlags.EnvSettings...), // Inject live job resources, then user env vars
		BindPaths:      execFlags.BindPaths,
		ApptainerFlags: apptainerFlags, // Pass through unknown flags
		Fakeroot:       execFlags.Fakeroot,
		BaseImage:      baseImageResolved,
		ApptainerBin:   config.Global.ApptainerBin,
		HidePrompt:     hidePrompt,
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
