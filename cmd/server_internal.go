package cmd

import (
	"context"
	"fmt"
	"net"
	"os"
	"os/signal"
	"strconv"
	"syscall"
	"time"

	"github.com/Justype/condatainer/internal/helper"
	"github.com/Justype/condatainer/internal/image/producer"
	"github.com/Justype/condatainer/internal/server"
	"github.com/spf13/cobra"
)

func init() {
	rootCmd.AddCommand(serverReadyCmd)
	rootCmd.AddCommand(serverMessageCmd)
	rootCmd.AddCommand(serverDoneCmd)
	rootCmd.AddCommand(pickPortCmd)
	rootCmd.AddCommand(helperHoldCmd)
	rootCmd.AddCommand(helperStateDirCmd)
	rootCmd.AddCommand(serverDaemonCmd)
}

// ── _server_ready ──────────────────────────────────────────────────────────

var (
	serverReadyID          string
	serverReadyPort        int
	serverReadyNode        string
	serverReadyURLPath     string
	serverReadyExternalURL string
	serverReadyJobID       string
	serverReadyWait        int
)

var serverReadyCmd = &cobra.Command{
	Use:    "_server_ready",
	Hidden: true,
	Short:  "Signal that a helper service is ready (called from compute node)",
	RunE: func(cmd *cobra.Command, args []string) error {
		cmd.SilenceUsage = true
		id := serverReadyID
		if id == "" {
			id = os.Getenv("CNT_HELPER_ID")
		}
		if id == "" {
			return fmt.Errorf("--id or $CNT_HELPER_ID is required")
		}

		node := serverReadyNode
		if node == "" {
			node, _ = os.Hostname()
		}

		walltimeSec, _ := strconv.ParseInt(os.Getenv("CNT_HELPER_WALLTIME_SECS"), 10, 64)
		jobID := serverReadyJobID
		if jobID == "" {
			jobID = os.Getenv("CNT_HELPER_JOB_ID")
		}

		if serverReadyPort > 0 && serverReadyExternalURL == "" && serverReadyWait > 0 {
			timeout := time.Duration(serverReadyWait) * time.Second
			if !waitForPort(cmd.Context(), serverReadyPort, timeout) {
				fmt.Fprintf(os.Stderr, "warning: port %d is not accepting connections after %s; reporting ready anyway\n",
					serverReadyPort, timeout)
			}
		}

		ev := helper.StateEvent{
			Type:        "ready",
			Port:        serverReadyPort,
			Node:        node,
			WalltimeSec: walltimeSec,
			JobID:       jobID,
			URLPath:     serverReadyURLPath,
			ExternalURL: serverReadyExternalURL,
		}
		if err := helper.AppendEvent(id, ev); err != nil {
			return fmt.Errorf("writing ready event: %w", err)
		}
		return nil
	},
}

func init() {
	serverReadyCmd.Flags().StringVar(&serverReadyID, "id", "", "Helper run ID (defaults to $CNT_HELPER_ID)")
	serverReadyCmd.Flags().IntVar(&serverReadyPort, "port", 0, "Service port (0 = no tunnel)")
	serverReadyCmd.Flags().StringVar(&serverReadyNode, "node", "", "Compute node hostname (defaults to $HOSTNAME)")
	serverReadyCmd.Flags().StringVar(&serverReadyURLPath, "url-path", "", "URL path/query appended to proxy URL")
	serverReadyCmd.Flags().StringVar(&serverReadyExternalURL, "external-url", "", "External URL shown as-is")
	serverReadyCmd.Flags().StringVar(&serverReadyJobID, "job-id", "", "Scheduler job ID")
	serverReadyCmd.Flags().IntVar(&serverReadyWait, "wait", 60, "Seconds to wait for the port to accept connections before reporting ready (0 = don't wait)")
}

// waitForPort polls 127.0.0.1:port until it accepts a connection, timeout
// passes, or ctx ends, and reports whether it connected.
func waitForPort(ctx context.Context, port int, timeout time.Duration) bool {
	addr := net.JoinHostPort("127.0.0.1", strconv.Itoa(port))
	deadline := time.Now().Add(timeout)
	for {
		if conn, err := net.DialTimeout("tcp", addr, time.Second); err == nil {
			conn.Close() //nolint:errcheck
			return true
		}
		if time.Now().After(deadline) {
			return false
		}
		select {
		case <-ctx.Done():
			return false
		case <-time.After(250 * time.Millisecond):
		}
	}
}

// ── _server_message ────────────────────────────────────────────────────────

var (
	serverMessageID    string
	serverMessageLevel string
)

var serverMessageCmd = &cobra.Command{
	Use:    "_server_message <text>",
	Hidden: true,
	Short:  "Append a message to the helper messages file (called from compute node)",
	Args:   cobra.MinimumNArgs(1),
	RunE: func(cmd *cobra.Command, args []string) error {
		cmd.SilenceUsage = true
		id := serverMessageID
		if id == "" {
			id = os.Getenv("CNT_HELPER_ID")
		}
		if id == "" {
			return fmt.Errorf("--id or $CNT_HELPER_ID is required")
		}
		level := serverMessageLevel
		if level == "" {
			level = "info"
		}
		return helper.AppendEvent(id, helper.StateEvent{Type: "msg", Level: level, Text: args[0]})
	},
}

func init() {
	serverMessageCmd.Flags().StringVar(&serverMessageID, "id", "", "Helper run ID (defaults to $CNT_HELPER_ID)")
	serverMessageCmd.Flags().StringVar(&serverMessageLevel, "level", "info", "Message level: info|warn|error")
}

// ── _server_done ───────────────────────────────────────────────────────────

var (
	serverDoneID       string
	serverDoneExitCode int
)

var serverDoneCmd = &cobra.Command{
	Use:    "_server_done",
	Hidden: true,
	Short:  "Signal that a helper job has finished (called from the wrapper script)",
	RunE: func(cmd *cobra.Command, args []string) error {
		cmd.SilenceUsage = true
		id := serverDoneID
		if id == "" {
			id = os.Getenv("CNT_HELPER_ID")
		}
		if id == "" {
			return fmt.Errorf("--id or $CNT_HELPER_ID is required")
		}
		return helper.AppendEvent(id, helper.StateEvent{
			Type:     "done",
			ExitCode: &serverDoneExitCode,
		})
	},
}

func init() {
	serverDoneCmd.Flags().StringVar(&serverDoneID, "id", "", "Helper run ID (defaults to $CNT_HELPER_ID)")
	serverDoneCmd.Flags().IntVar(&serverDoneExitCode, "exit-code", 0, "Process exit code")
}

// ── _pick_port ─────────────────────────────────────────────────────────────

var pickPortCmd = &cobra.Command{
	Use:    "_pick_port",
	Hidden: true,
	Short:  "Print a free TCP port on this host",
	RunE: func(cmd *cobra.Command, args []string) error {
		cmd.SilenceUsage = true
		port, err := pickFreePort()
		if err != nil {
			return err
		}
		fmt.Println(port)
		return nil
	},
}

func pickFreePort() (int, error) {
	ln, err := net.Listen("tcp", "127.0.0.1:0")
	if err != nil {
		return 0, fmt.Errorf("finding free port: %w", err)
	}
	port := ln.Addr().(*net.TCPAddr).Port
	ln.Close()
	return port, nil
}

// ── _helper_hold ───────────────────────────────────────────────────────────

var helperHoldCmd = &cobra.Command{
	Use:    "_helper_hold <lock-file> <wrapper-pid>",
	Hidden: true,
	Short:  "Hold a headless helper's lock file until its wrapper exits (called from the wrapper)",
	Args:   cobra.ExactArgs(2),
	RunE: func(cmd *cobra.Command, args []string) error {
		cmd.SilenceUsage = true
		pid, err := strconv.Atoi(args[1])
		if err != nil || pid <= 0 {
			return fmt.Errorf("invalid wrapper pid %q", args[1])
		}
		// The wrapper's group TERM must not release the lock before it has
		// recorded done; the holder ends only when the wrapper does.
		signal.Ignore(syscall.SIGTERM, syscall.SIGINT, syscall.SIGHUP)
		lock, err := helper.HoldLock(args[0], producer.Info{
			Runner:    "local",
			Node:      producer.ShortHostname(),
			PID:       pid,
			CreatedAt: time.Now().Format(time.RFC3339),
		})
		if err != nil {
			return err
		}
		defer lock.Close()
		for os.Getppid() == pid {
			time.Sleep(time.Second)
		}
		return nil
	},
}

// ── _helper_state_dir ──────────────────────────────────────────────────────

var helperStateDirCmd = &cobra.Command{
	Use:    "_helper_state_dir <id>",
	Hidden: true,
	Short:  "Print the NFS state directory for a helper run ID",
	Args:   cobra.ExactArgs(1),
	RunE: func(cmd *cobra.Command, args []string) error {
		cmd.SilenceUsage = true
		dir := helper.StateDir(args[0])
		if dir == "" {
			return fmt.Errorf("cannot determine state directory")
		}
		fmt.Println(dir)
		return nil
	},
}

// ── _server_daemon ─────────────────────────────────────────────────────────

var (
	serverDaemonPort     int
	serverDaemonReportFd int
	serverDaemonWatchPID int
)

var serverDaemonCmd = &cobra.Command{
	Use:    "_server_daemon",
	Hidden: true,
	Short:  "Run the condatainer server daemon (invoked by 'server start')",
	RunE: func(cmd *cobra.Command, args []string) error {
		cmd.SilenceUsage = true

		port := serverDaemonPort
		if port == 0 {
			var err error
			port, err = pickFreePort()
			if err != nil {
				return err
			}
		}

		var reportPipe *os.File
		if serverDaemonReportFd > 0 {
			reportPipe = os.NewFile(uintptr(serverDaemonReportFd), "report")
		}

		return server.RunDaemon(port, reportPipe, serverDaemonWatchPID)
	},
}

func init() {
	serverDaemonCmd.Flags().IntVar(&serverDaemonPort, "port", 0, "Port to listen on (0 = use saved or free port)")
	serverDaemonCmd.Flags().IntVar(&serverDaemonReportFd, "report-fd", 0, "File descriptor to signal readiness to parent")
	serverDaemonCmd.Flags().IntVar(&serverDaemonWatchPID, "watch-pid", 0, "Exit when this PID dies (SSH session tracking)")
}
