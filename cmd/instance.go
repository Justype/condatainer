package cmd

import (
	"context"
	"errors"
	"fmt"
	"strings"

	"github.com/Justype/condatainer/internal/apptainer"
	"github.com/Justype/condatainer/internal/config"
	"github.com/Justype/condatainer/internal/container"
	"github.com/Justype/condatainer/internal/instance"
	"github.com/Justype/condatainer/internal/utils"
	"github.com/spf13/cobra"
)

// Instance command variables for start/run
var instanceFlags CommonFlags

// Instance exec specific flags
var instanceExecEnvSettings []string

// instanceCmd is the parent command for instance management
var instanceCmd = &cobra.Command{
	Use:   "instance",
	Short: "Manage Apptainer instances",
	Long: `Manage named Apptainer instances.

Instances are persistent containers that run in the background.
You can start, stop, list, and execute commands in instances.`,
}

// instanceListCmd lists all running instances
var instanceListCmd = &cobra.Command{
	Use:          "list",
	Short:        "List all running and named Apptainer instances",
	SilenceUsage: true,
	RunE: func(cmd *cobra.Command, args []string) error {
		if err := apptainer.SetBin(config.Global.ApptainerBin); err != nil {
			return err
		}
		return apptainer.InstanceList(cmd.Context())
	},
}

// instanceStartCmd starts a named instance
var instanceStartCmd = &cobra.Command{
	Use:   "start [flags] <name>",
	Short: "Start a named instance of the given container image",
	Long: `Start a named instance with specified overlays.

The instance runs in the background and can be accessed later via 'instance exec'.

Use -o/--overlay to explicitly specify overlays.

Examples:
    # Start with multiple overlays
    condatainer instance start -o samtools/1.22 -o bcftools/1.22 myinstance

    # Start with writable .img overlay
    condatainer instance start -w -o env.img myinstance

    # Pass apptainer flags
    condatainer instance start --home=/ext3/home -o samtools/1.22 myinstance

Note: 
- Unknown flags must use --flag=value format (no space) for unambiguous parsing`,
	SilenceUsage: true,
	RunE:         runInstanceStart,
}

// Instance stop specific flags
var instanceStopAll bool
var instanceStopForce bool
var instanceStopSignal string
var instanceStopTimeout int

// instanceStopCmd stops a running instance
var instanceStopCmd = &cobra.Command{
	Use:   "stop [flags] [name]",
	Short: "Stop a named instance of a given container image",
	Long: `Stop a running instance.

Examples:
    # Stop a specific instance
    condatainer instance stop myinstance
    condatainer instance stop mysql*
    condatainer instance stop --all

    # Force stop
    condatainer instance stop --force myinstance

    # Send SIGTERM (SIGTERM, TERM, and 15 are all valid)
    condatainer instance stop --signal TERM myinstance

    # Custom timeout
    condatainer instance stop --timeout 30 myinstance`,
	SilenceUsage: true,
	RunE: func(cmd *cobra.Command, args []string) error {
		// Validate arguments
		if !instanceStopAll && len(args) == 0 {
			return fmt.Errorf("instance name is required (or use --all to stop all instances)")
		}

		if instanceStopAll && len(args) > 0 {
			return fmt.Errorf("cannot specify instance name with --all flag")
		}

		// Build apptainer flags
		var apptainerFlags []string
		if instanceStopAll {
			apptainerFlags = append(apptainerFlags, "--all")
		}
		if instanceStopForce {
			apptainerFlags = append(apptainerFlags, "--force")
		}
		if instanceStopSignal != "" {
			apptainerFlags = append(apptainerFlags, "--signal", instanceStopSignal)
		}
		if instanceStopTimeout > 0 {
			apptainerFlags = append(apptainerFlags, "--timeout", fmt.Sprintf("%d", instanceStopTimeout))
		}

		// Add instance names/patterns
		apptainerFlags = append(apptainerFlags, args...)

		opts := instance.StopOptions{
			ApptainerBin:   config.Global.ApptainerBin,
			ApptainerFlags: apptainerFlags,
			All:            instanceStopAll,
			InstanceNames:  args,
		}

		if err := instance.Stop(cmd.Context(), opts); err != nil {
			return err
		}

		utils.PrintSuccess("Instance(s) stopped successfully")
		return nil
	},
}

// instanceStatsCmd gets stats for a named instance
var instanceStatsCmd = &cobra.Command{
	Use:          "stats <name>",
	Short:        "Get stats for a named instance",
	Args:         cobra.ExactArgs(1),
	SilenceUsage: true,
	RunE: func(cmd *cobra.Command, args []string) error {
		instanceName := args[0]
		if err := apptainer.SetBin(config.Global.ApptainerBin); err != nil {
			return err
		}
		return apptainer.InstanceStats(cmd.Context(), instanceName)
	},
}

// instanceExecCmd executes a command in a running instance
var instanceExecCmd = &cobra.Command{
	Use:   "exec [flags] <name> <command> [args...]",
	Short: "Execute commands in a named instance",
	Long: `Execute a command inside a running instance.

The instance must be already running (started with 'instance start').

Examples:
    # Run a command in an instance
    condatainer instance exec myinstance samtools view file.bam

    # Set additional environment variables
    condatainer instance exec --env VAR1=val1 --env VAR2=val2 myinstance bash

    # Pass apptainer flags (use --flag=value format)
    condatainer instance exec --home=/custom myinstance bash

Note: Unknown flags must use --flag=value format (no space) for unambiguous parsing.`,
	SilenceUsage: true,
	RunE:         runInstanceExec,
}

func init() {
	rootCmd.AddCommand(instanceCmd)

	// Add subcommands
	instanceCmd.AddCommand(instanceListCmd)
	instanceCmd.AddCommand(instanceStartCmd)
	instanceCmd.AddCommand(instanceStopCmd)
	instanceCmd.AddCommand(instanceStatsCmd)
	instanceCmd.AddCommand(instanceExecCmd)

	// Configure start command flags using shared helper
	RegisterCommonFlags(instanceStartCmd, &instanceFlags)

	// Configure stop command flags
	instanceStopCmd.Flags().BoolVarP(&instanceStopAll, "all", "a", false, "stop all user's instances")
	instanceStopCmd.Flags().BoolVarP(&instanceStopForce, "force", "F", false, "force kill instance (may corrupt data)")
	instanceStopCmd.Flags().StringVarP(&instanceStopSignal, "signal", "s", "", "signal to send to the instance")
	instanceStopCmd.Flags().IntVarP(&instanceStopTimeout, "timeout", "t", 0, "timeout before force kill (seconds)")
	instanceStartCmd.ValidArgsFunction = func(cmd *cobra.Command, args []string, toComplete string) ([]string, cobra.ShellCompDirective) {
		return nil, cobra.ShellCompDirectiveDefault
	}

	// Configure instance exec command flags
	instanceExecCmd.Flags().StringSliceVar(&instanceExecEnvSettings, "env", nil, "set environment variable 'KEY=VALUE' (can be used multiple times)")
	// Stop flag parsing after the first positional argument
	instanceExecCmd.Flags().SetInterspersed(false)
	// Allow unknown flags to pass through (for apptainer)
	instanceExecCmd.FParseErrWhitelist.UnknownFlags = true
	// Enable default file completion for positional args
	instanceExecCmd.ValidArgsFunction = func(cmd *cobra.Command, args []string, toComplete string) ([]string, cobra.ShellCompDirective) {
		return nil, cobra.ShellCompDirectiveDefault
	}
}

// runInstanceStart handles the instance start command
func runInstanceStart(cmd *cobra.Command, args []string) error {
	if err := ensureBaseImage(cmd.Context()); err != nil {
		return err
	}

	// Parse arguments: instance name and apptainer flags
	instanceName, apptainerFlags := ParseInstanceArgs("start")

	if instanceName == "" {
		return errors.New("instance name is required")
	}

	// Ensure only one positional argument (the instance name)
	if len(args) > 1 {
		return errors.New("only one positional argument (instance name) is allowed; use -o/--overlay for overlays")
	}

	// Use overlays from -o flag only
	overlayFinal := instanceFlags.Overlays

	// Resolve user-specified base image if provided
	var baseImageResolved string
	if instanceFlags.BaseImage != "" {
		if utils.FileExists(instanceFlags.BaseImage) {
			baseImageResolved = instanceFlags.BaseImage
		} else {
			resolvedBase, err := container.ResolveOverlayPaths([]string{instanceFlags.BaseImage})
			if err == nil && len(resolvedBase) > 0 {
				baseImageResolved = resolvedBase[0]
			} else {
				baseImageResolved = instanceFlags.BaseImage
			}
		}
	}

	resolvedOverlays, err := container.ResolveOverlayPaths(overlayFinal)
	if err != nil {
		return err
	}

	options := instance.Options{
		Name:           instanceName,
		Overlays:       resolvedOverlays,
		WritableImg:    instanceFlags.WritableImg,
		EnvSettings:    instanceFlags.EnvSettings,
		BindPaths:      instanceFlags.BindPaths,
		ApptainerFlags: apptainerFlags,
		Fakeroot:       instanceFlags.Fakeroot,
		BaseImage:      baseImageResolved,
		ApptainerBin:   config.Global.ApptainerBin,
	}

	if err := instance.Start(cmd.Context(), options); err != nil {
		if errors.Is(err, context.Canceled) || errors.Is(cmd.Context().Err(), context.Canceled) {
			return nil
		}
		return err
	}

	utils.PrintSuccess("Instance %s started successfully", utils.StyleInfo(instanceName))
	return nil
}

// runInstanceExec handles the instance exec command
func runInstanceExec(cmd *cobra.Command, args []string) error {
	// Parse args: unknown flags before positional args, then instance name, then command
	var apptainerFlags []string
	var instanceName string
	var command []string

	// Extract unknown flags and positional args
	for i := 0; i < len(args); i++ {
		arg := args[i]

		// If it's a flag (starts with -), add to apptainer flags
		if strings.HasPrefix(arg, "-") {
			apptainerFlags = append(apptainerFlags, arg)
			continue
		}

		// First non-flag arg is the instance name
		if instanceName == "" {
			instanceName = arg
			// Everything after instance name is the command
			command = args[i+1:]
			break
		}
	}

	if instanceName == "" {
		return errors.New("instance name is required")
	}

	if len(command) == 0 {
		command = []string{"bash"}
	}

	options := instance.ExecOptions{
		Name:           instanceName,
		Command:        command,
		EnvSettings:    instanceExecEnvSettings, // From --env flag (handled by Cobra)
		ApptainerFlags: apptainerFlags,
		ApptainerBin:   config.Global.ApptainerBin,
		PrintEnv:       len(command) <= 1, // Print env if command length <= 1
	}

	if err := instance.Exec(cmd.Context(), options); err != nil {
		if errors.Is(err, context.Canceled) || errors.Is(cmd.Context().Err(), context.Canceled) {
			return nil
		}
		return err
	}

	return nil
}
