package cmd

import (
	"bytes"
	"fmt"
	"io"
	"os"
	"os/exec"
	"strconv"
	"strings"
	"syscall"
	"time"

	"github.com/Justype/condatainer/internal/config"
	"github.com/Justype/condatainer/internal/server"
	"github.com/Justype/condatainer/internal/utils"
	"github.com/spf13/cobra"
)

var serverCmd = &cobra.Command{
	Use:   "server",
	Short: "Manage the condatainer dashboard server",
	Long: `Manage the condatainer dashboard server.

The server runs on the login node and provides:
  - Web dashboard at http://localhost:{port}/
  - Reverse proxy to helper services on compute nodes via SSH tunnels
  - Live SSE log streaming

The port is persisted after the first start. Set up SSH once:
  LocalForward {port} localhost:{port}`,
}

var (
	serverStartPort   int
	serverStartDaemon bool
)

var serverStartCmd = &cobra.Command{
	Use:   "start",
	Short: "Start the dashboard server in the background",
	Long: `Start the dashboard server in the background.

By default the server exits when your SSH session ends.`,
	RunE: func(cmd *cobra.Command, args []string) error {
		cmd.SilenceUsage = true

		port := serverStartPort
		if port != 0 {
			// Explicit --port: persist it for future starts.
			if err := config.SaveServerPort(port); err != nil {
				utils.PrintDebug("saving server port: %v", err)
			}
		} else {
			port = config.ReadSavedServerPort()
			if port == 0 {
				var err error
				port, err = pickFreePort()
				if err != nil {
					return fmt.Errorf("finding free port: %w", err)
				}
			}
		}

		pidFile := config.GetServerPidFilePath()
		if server.IsAlive(pidFile) {
			if ss, err := server.ReadState(pidFile); err == nil {
				utils.PrintMessage("Server already running at http://localhost:%d/", ss.Port)
				utils.PrintMessage("Log: %s", config.GetServerLogPath())
				return nil
			}
		}

		return startServerDaemon(port, serverStartDaemon)
	},
}

var serverStopNodeFlag string

var serverStopCmd = &cobra.Command{
	Use:   "stop [node]",
	Short: "Stop the dashboard server",
	Long: `Stop the dashboard server.

Without arguments, stops the server on the current login node.
Pass a node name (or partial match) to stop a server on another login node.

Example:
  condatainer server stop              # stop local server
  condatainer server stop login2       # stop server on the login node matching "login2"
  condatainer server stop --node login2`,
	Args: cobra.MaximumNArgs(1),
	RunE: func(cmd *cobra.Command, args []string) error {
		cmd.SilenceUsage = true

		nodePattern := serverStopNodeFlag
		if nodePattern == "" && len(args) > 0 {
			nodePattern = args[0]
		}

		currentHost, _ := os.Hostname()

		// No pattern or pattern matches local host → stop local server.
		if nodePattern == "" || matchHost(nodePattern, currentHost) {
			pidFile := config.GetServerPidFilePath()
			ss, err := server.ReadState(pidFile)
			if err != nil {
				utils.PrintMessage("Server is not running.")
				return nil
			}
			p, err := os.FindProcess(ss.PID)
			if err != nil || p.Signal(syscall.Signal(0)) != nil {
				utils.PrintMessage("Server is not running (stale PID file).")
				os.Remove(pidFile)
				return nil
			}
			if err := p.Signal(syscall.SIGTERM); err != nil {
				return fmt.Errorf("sending SIGTERM to PID %d: %w", ss.PID, err)
			}
			utils.PrintSuccess("Server (PID %d) stopped.", ss.PID)
			return nil
		}

		// Pattern given → scan all PID files and stop on matching remote host.
		for _, pidFile := range config.ListServerPidFiles() {
			ss, err := server.ReadState(pidFile)
			if err != nil {
				continue
			}
			if !matchHost(nodePattern, ss.Host) {
				continue
			}
			if err := stopRemoteServer(ss); err != nil {
				return fmt.Errorf("stopping server on %s: %w", ss.Host, err)
			}
			utils.PrintSuccess("Server on %s (PID %d) stopped.", ss.Host, ss.PID)
			return nil
		}

		utils.PrintMessage("No server found matching %q.", nodePattern)
		return nil
	},
}

// matchHost reports whether pattern partially matches hostname.
// Exact match, or hostname starts with pattern+"." (subdomain prefix), or contains.
func matchHost(pattern, hostname string) bool {
	if hostname == pattern {
		return true
	}
	if strings.HasPrefix(hostname, pattern+".") {
		return true
	}
	return strings.Contains(hostname, pattern)
}

// stopRemoteServer sends SIGTERM to a server on a remote login node via SSH.
func stopRemoteServer(ss *config.ServerState) error {
	return exec.Command("ssh", ss.Host, "kill", "-TERM", strconv.Itoa(ss.PID)).Run()
}

var serverStatusCmd = &cobra.Command{
	Use:   "status",
	Short: "Show dashboard server status across all login nodes",
	RunE: func(cmd *cobra.Command, args []string) error {
		cmd.SilenceUsage = true

		currentHost, _ := os.Hostname()
		pidFiles := config.ListServerPidFiles()

		if len(pidFiles) == 0 {
			utils.PrintMessage("Server is not running.")
			return nil
		}

		found := false
		for _, pidFile := range pidFiles {
			ss, err := server.ReadState(pidFile)
			if err != nil {
				continue
			}
			found = true
			if ss.Host == currentHost {
				if server.IsAlive(pidFile) {
					utils.PrintSuccess("Running: %s  port:%d  PID:%d", ss.Host, ss.Port, ss.PID)
				} else {
					utils.PrintMessage("Not running: %s (stale PID file)", ss.Host)
				}
			} else {
				utils.PrintSuccess("Running: %s  port:%d  PID:%d", ss.Host, ss.Port, ss.PID)
			}
		}
		if !found {
			utils.PrintMessage("Server is not running.")
		}
		return nil
	},
}

var serverRestartPort int

var serverRestartCmd = &cobra.Command{
	Use:   "restart",
	Short: "Restart the dashboard server (stop + start)",
	RunE: func(cmd *cobra.Command, args []string) error {
		cmd.SilenceUsage = true

		pidFile := config.GetServerPidFilePath()
		ss, err := server.ReadState(pidFile)

		// Determine target port: flag > saved state > saved config > free port.
		port := serverRestartPort
		if port != 0 {
			if err := config.SaveServerPort(port); err != nil {
				utils.PrintDebug("saving server port: %v", err)
			}
		} else if ss != nil {
			port = ss.Port
		} else {
			port = config.ReadSavedServerPort()
		}

		if err != nil || ss == nil {
			utils.PrintMessage("Server is not running, starting fresh.")
			return startServerDaemon(port, serverStartDaemon)
		}

		p, err := os.FindProcess(ss.PID)
		if err != nil || p.Signal(syscall.Signal(0)) != nil {
			utils.PrintMessage("Server process not found, starting fresh.")
			return startServerDaemon(port, serverStartDaemon)
		}
		if err := p.Signal(syscall.SIGTERM); err != nil {
			return fmt.Errorf("sending SIGTERM to PID %d: %w", ss.PID, err)
		}
		utils.PrintMessage("Stopping server (PID %d)…", ss.PID)

		// Wait for the process to exit (up to 10s).
		for range 50 {
			time.Sleep(200 * time.Millisecond)
			if p.Signal(syscall.Signal(0)) != nil {
				break
			}
		}

		utils.PrintMessage("Starting server on port %d…", port)
		return startServerDaemon(port, serverStartDaemon)
	},
}

func init() {
	rootCmd.AddCommand(serverCmd)
	serverCmd.AddCommand(serverStartCmd, serverStopCmd, serverStatusCmd, serverRestartCmd)

	serverStartCmd.Flags().IntVarP(&serverStartPort, "port", "p", 0, "Port to listen on (saved to config for reuse)")
	serverStartCmd.Flags().BoolVar(&serverStartDaemon, "daemon", false, "Survive SSH logout (useful with a permanent SSH LocalForward)")
	serverStopCmd.Flags().StringVarP(&serverStopNodeFlag, "node", "n", "", "Login node to stop (partial match, e.g. login2)")
	serverRestartCmd.Flags().IntVarP(&serverRestartPort, "port", "p", 0, "Restart on a different port (saved to config for reuse)")
	serverRestartCmd.Flags().BoolVar(&serverStartDaemon, "daemon", false, "Survive SSH logout")
}

// startServerDaemon forks the server using the pipe-readiness protocol.
// The child always runs in a new session (Setsid) so sshd exits cleanly on logout.
// In non-daemon mode the child watches the parent shell PID and exits when it dies.
// stdout/stderr are redirected to the server log file.
func startServerDaemon(port int, daemon bool) error {
	exe, err := os.Executable()
	if err != nil {
		return fmt.Errorf("finding executable: %w", err)
	}

	// Open log file before forking so errors are visible to the user.
	logPath := config.GetServerLogPath()
	var logFile *os.File
	if logPath != "" {
		if utils.EnsureWritableDir(config.GetUserStateDir()) {
			logFile, _ = os.OpenFile(logPath, os.O_CREATE|os.O_APPEND|os.O_WRONLY, 0664)
		}
	}

	pr, pw, err := os.Pipe()
	if err != nil {
		if logFile != nil {
			logFile.Close()
		}
		return fmt.Errorf("creating pipe: %w", err)
	}

	args := []string{"_server_daemon", "--report-fd", "3", "--port", strconv.Itoa(port)}
	// Non-daemon: watch the SSH session leader so the server exits when the user
	// logs out, regardless of how many process layers deep this call originates
	// (e.g. auto-started from "condatainer helper run" via EnsureServer).
	if !daemon {
		if sid := sshSessionLeaderPID(); sid > 1 {
			args = append(args, "--watch-pid", strconv.Itoa(sid))
		}
	}
	daemonCmd := exec.Command(exe, args...)
	// Always create a new session so sshd exits cleanly on logout.
	daemonCmd.SysProcAttr = &syscall.SysProcAttr{Setsid: true}
	daemonCmd.ExtraFiles = []*os.File{pw}
	daemonCmd.Stdout = logFile
	daemonCmd.Stderr = logFile

	if err := daemonCmd.Start(); err != nil {
		pw.Close()
		pr.Close()
		if logFile != nil {
			logFile.Close()
		}
		return fmt.Errorf("starting server: %w", err)
	}
	pw.Close()
	if logFile != nil {
		logFile.Close()
	}
	daemonCmd.Process.Release() //nolint:errcheck

	msg, _ := io.ReadAll(pr)
	pr.Close()
	if len(msg) > 0 {
		return fmt.Errorf("%s", strings.TrimSpace(string(msg)))
	}

	utils.PrintSuccess("Sever Dashboard started: http://localhost:%d/", port)
	if logPath != "" {
		utils.PrintMessage("Log: %s", logPath)
	}
	return nil
}

// sshSessionLeaderPID returns the PID of the current session leader by reading
// /proc/self/stat. This is the login shell regardless of process nesting depth,
// so the server correctly exits when the user logs out even when auto-started
// from a helper command. Falls back to os.Getppid() if /proc is unavailable.
func sshSessionLeaderPID() int {
	data, err := os.ReadFile("/proc/self/stat")
	if err != nil {
		return os.Getppid()
	}
	// Format: pid (comm) state ppid pgrp session ...
	// Skip past comm which may contain spaces by finding the last ')'.
	idx := bytes.LastIndexByte(data, ')')
	if idx < 0 {
		return os.Getppid()
	}
	fields := bytes.Fields(data[idx+1:])
	// After ')': state(0) ppid(1) pgrp(2) session(3)
	if len(fields) < 4 {
		return os.Getppid()
	}
	sid, err := strconv.Atoi(string(fields[3]))
	if err != nil || sid <= 1 {
		return os.Getppid()
	}
	return sid
}
