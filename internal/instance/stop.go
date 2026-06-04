package instance

import (
	"context"
	"strings"

	"log/slog"

	"github.com/Justype/condatainer/internal/apptainer"
	"github.com/Justype/condatainer/internal/config"
)

// StopOptions holds configuration for stopping an instance
type StopOptions struct {
	ApptainerBin   string   // Path to apptainer binary
	ApptainerFlags []string // All arguments to pass to apptainer instance stop
	All            bool     // Whether --all flag is set
	InstanceNames  []string // Instance names/patterns to stop
}

// Stop stops a running instance
func Stop(ctx context.Context, opts StopOptions) error {
	if opts.ApptainerBin == "" {
		opts.ApptainerBin = config.Global.ApptainerBin
	}

	if err := apptainer.SetBin(opts.ApptainerBin); err != nil {
		return err
	}

	// Pass all flags directly to apptainer (includes instance name and all flags)
	if err := apptainer.InstanceStop(ctx, "", opts.ApptainerFlags); err != nil {
		return err
	}

	// Clean up state files based on the options
	if opts.All {
		// Delete all state files when --all flag is used
		if err := deleteAllStates(); err != nil {
			slog.Default().Debug("instance: failed to delete all states", "err", err)
		}
		return nil
	}

	// Otherwise, clean up state files for specific instances/patterns
	for _, name := range opts.InstanceNames {
		// Check if this is a wildcard pattern
		if strings.ContainsAny(name, "*?[]") {
			// Delete states matching the pattern
			if err := deleteStatesMatching(name); err != nil {
				slog.Default().Debug("instance: failed to delete states matching pattern", "pattern", name, "err", err)
			}
		} else {
			// Delete specific instance state
			if err := deleteState(name); err != nil {
				slog.Default().Debug("instance: failed to delete state", "name", name, "err", err)
			}
		}
	}

	return nil
}
