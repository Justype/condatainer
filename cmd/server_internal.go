package cmd

import (
	"fmt"
	"net"
	"os"
	"strconv"

	"github.com/Justype/condatainer/internal/helper"
	"github.com/Justype/condatainer/internal/server"
	"github.com/spf13/cobra"
)

func init() {
	rootCmd.AddCommand(serverReadyCmd)
	rootCmd.AddCommand(serverMessageCmd)
	rootCmd.AddCommand(serverDoneCmd)
	rootCmd.AddCommand(pickPortCmd)
	rootCmd.AddCommand(helperStateDirCmd)
	rootCmd.AddCommand(serverDaemonCmd)
}

// ── _server_ready ──────────────────────────────────────────────────────────

var (
	serverReadyID          string
	serverReadyPort        int
	serverReadyLabel       string
	serverReadyNode        string
	serverReadyURLPath     string
	serverReadyExternalURL string
	serverReadyJobID       string
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

		ev := helper.StateEvent{
			Type:        "ready",
			Port:        serverReadyPort,
			Node:        node,
			Label:       serverReadyLabel,
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
	serverReadyCmd.Flags().StringVar(&serverReadyLabel, "label", "", "Human-readable service label")
	serverReadyCmd.Flags().StringVar(&serverReadyNode, "node", "", "Compute node hostname (defaults to $HOSTNAME)")
	serverReadyCmd.Flags().StringVar(&serverReadyURLPath, "url-path", "", "URL path/query appended to proxy URL")
	serverReadyCmd.Flags().StringVar(&serverReadyExternalURL, "external-url", "", "External URL shown as-is")
	serverReadyCmd.Flags().StringVar(&serverReadyJobID, "job-id", "", "Scheduler job ID")
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
