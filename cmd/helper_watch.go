package cmd

import (
	"context"
	"errors"
	"fmt"
	"os"
	"syscall"
	"time"

	ui "github.com/Justype/condatainer/cmd/internal/ui"
	"github.com/spf13/cobra"
	"golang.org/x/sys/unix"
)

var (
	helperWatchID       string
	helperWatchTTY      string
	helperWatchShellPID int
)

// helperWatchCmd is what `condatainer helper` leaves running after it submits:
// it prints to the launching terminal when the service is ready, then ends. It
// also ends when the launching shell does.
var helperWatchCmd = &cobra.Command{
	Use:    "_helper_watch",
	Hidden: true,
	Short:  "Watch a helper in the background (started by 'condatainer helper')",
	Args:   cobra.NoArgs,
	RunE: func(cmd *cobra.Command, args []string) error {
		cmd.SilenceUsage = true
		if helperWatchID == "" {
			return fmt.Errorf("--id is required")
		}
		if err := attachOutput(helperWatchTTY); err != nil {
			return err
		}
		ctx, cancel := context.WithCancel(cmd.Context())
		defer cancel()
		go endWithShell(ctx, cancel, helperWatchShellPID)
		return ui.MonitorHelper(ctx, helperWatchID, true)
	},
}

func init() {
	rootCmd.AddCommand(helperWatchCmd)
	helperWatchCmd.Flags().StringVar(&helperWatchID, "id", "", "Helper run ID")
	helperWatchCmd.Flags().StringVar(&helperWatchTTY, "tty", "", "Terminal device to print to")
	helperWatchCmd.Flags().IntVar(&helperWatchShellPID, "shell-pid", 0, "End when this process exits")
}

// endWithShell cancels the watcher once the launching shell's process is gone.
func endWithShell(ctx context.Context, cancel context.CancelFunc, pid int) {
	if pid <= 1 {
		return
	}
	for {
		select {
		case <-ctx.Done():
			return
		case <-time.After(2 * time.Second):
			if err := syscall.Kill(pid, 0); errors.Is(err, syscall.ESRCH) {
				cancel()
				return
			}
		}
	}
}

// attachOutput points stdout and stderr at tty, or at /dev/null when it cannot
// be opened, so the watcher's output reaches the launching terminal and nothing
// else.
func attachOutput(tty string) error {
	target := os.DevNull
	if tty != "" {
		target = tty
	}
	f, err := os.OpenFile(target, os.O_WRONLY|unix.O_NOCTTY, 0)
	if err != nil {
		if f, err = os.OpenFile(os.DevNull, os.O_WRONLY, 0); err != nil {
			return err
		}
	}
	defer f.Close()
	for _, fd := range []int{1, 2} {
		if err := unix.Dup3(int(f.Fd()), fd, 0); err != nil {
			return err
		}
	}
	return nil
}
