package instance

import (
	"context"
	"fmt"
	"os"
	"sort"

	"log/slog"

	"github.com/Justype/condatainer/internal/apptainer"
	"github.com/Justype/condatainer/internal/config"
	"github.com/Justype/condatainer/internal/logging"
	"github.com/Justype/condatainer/internal/utils"
)

// ExecOptions holds configuration for executing commands in an instance
type ExecOptions struct {
	Name           string   // Instance name
	Command        []string // Command to execute
	EnvSettings    []string // Additional environment variables to set (format: "KEY=VALUE")
	ApptainerFlags []string // Additional apptainer flags to pass through
	ApptainerBin   string   // Path to apptainer binary
	PrintEnv       bool     // Whether to print environment variables (when command length <= 1)
}

// Exec executes a command in a running instance
func Exec(ctx context.Context, opts ExecOptions) error {
	if opts.ApptainerBin == "" {
		opts.ApptainerBin = config.Global.ApptainerBin
	}

	if err := apptainer.SetBin(opts.ApptainerBin); err != nil {
		return err
	}

	// Load instance state to get environment variables
	state, err := loadState(opts.Name)
	if err != nil {
		slog.Default().Debug("instance: could not load state", "err", err)
		logging.FromContext(ctx).Info("instance state not found, environment variables will not be set", "kind", "note")
	}

	// Print overlay environment notes for interactive shells with simple commands
	// (similar to exec behavior: only print when command length <= 1)
	if state != nil && utils.IsInteractiveShell() && opts.PrintEnv {
		log := logging.FromContext(ctx)
		if len(state.Notes) > 0 {
			log.Info("overlay envs:")
			sortedNotes := make([]string, 0, len(state.Notes))
			for key := range state.Notes {
				sortedNotes = append(sortedNotes, key)
			}
			sort.Strings(sortedNotes)
			for _, key := range sortedNotes {
				log.Info(fmt.Sprintf("  %s: %s", key, state.Notes[key]))
			}
		}
	}

	// Build environment variables list
	envList := opts.EnvSettings
	if state != nil {
		envList = append(state.Env, opts.EnvSettings...)
		slog.Default().Debug("instance: applied environment variables from state", "count", len(state.Env))
	}

	// Execute in instance using apptainer package
	execOpts := &apptainer.InstanceExecOptions{
		Env:        envList,
		Additional: opts.ApptainerFlags,
		Stdin:      os.Stdin,
	}

	return apptainer.InstanceExec(ctx, opts.Name, opts.Command, execOpts)
}
