package cmd

import (
	"fmt"
	"os"
	"os/exec"
	"strconv"
	"syscall"
	"time"

	"github.com/Justype/condatainer/internal/config"
	"github.com/Justype/condatainer/internal/proxy"
	"github.com/Justype/condatainer/internal/scheduler"
	"github.com/Justype/condatainer/internal/utils"
	"github.com/spf13/cobra"
)

var proxyCmd = &cobra.Command{
	Use:   "proxy",
	Short: "Manage the shared SOCKS5 proxy tunnel for compute nodes",
	Long: `Start and manage an SSH SOCKS5 proxy tunnel so that compute node jobs
can reach the internet through the login node.

Run "condatainer proxy start" once on the login node before submitting jobs.
All condatainer commands running inside jobs will automatically use the tunnel.`,
}

// proxy start flags
var (
	proxyStartHost string
	proxyStartPort int
)

var proxyStartCmd = &cobra.Command{
	Use:   "start",
	Short: "Start the SOCKS5 proxy tunnel (run on login node only)",
	RunE: func(cmd *cobra.Command, args []string) error {
		if scheduler.IsInsideJob() {
			return fmt.Errorf("proxy must be started on the login node, not inside a scheduler job")
		}
		if config.IsInsideContainer() {
			return fmt.Errorf("proxy must be started on the login node, not inside a container")
		}
		if scheduler.ActiveScheduler() == nil {
			return fmt.Errorf("no scheduler detected — proxy is only needed on HPC clusters")
		}

		// Check if already running
		if host, port, pid, err := proxy.ReadPidFile(); err == nil {
			if proxy.ProxyAlive(host, port) {
				utils.PrintMessage("Proxy already running on %s:%d (PID %d)", host, port, pid)
				return nil
			}
			// Stale PID file — clean it up
			proxy.RemovePidFile()
		}

		host := proxyStartHost
		if host == "" {
			host, _ = os.Hostname()
		}
		if host == "" {
			return fmt.Errorf("could not determine SSH host; use --host to specify")
		}

		port := proxyStartPort
		if port == 0 {
			var err error
			port, err = proxy.FreePort()
			if err != nil {
				return fmt.Errorf("failed to find free port: %w", err)
			}
		}

		// Double-fork: re-invoke ourselves as "_proxy_daemon" with Setsid so the
		// child is reparented to PID 1, surviving terminal close and KillUserProcesses.
		daemonCmd := exec.Command(os.Args[0], "_proxy_daemon", host, strconv.Itoa(port))
		daemonCmd.SysProcAttr = &syscall.SysProcAttr{Setsid: true}
		daemonCmd.Stdout = nil
		daemonCmd.Stderr = nil
		daemonCmd.Stdin = nil

		if err := daemonCmd.Start(); err != nil {
			return fmt.Errorf("failed to start proxy daemon: %w", err)
		}
		daemonCmd.Process.Release() //nolint:errcheck

		// Wait briefly for PID file to appear so we can confirm and print the port.
		pidFileReady := false
		for i := 0; i < 20; i++ {
			if _, _, _, err := proxy.ReadPidFile(); err == nil {
				pidFileReady = true
				break
			}
			time.Sleep(100 * time.Millisecond)
		}

		if pidFileReady {
			utils.PrintSuccess("Proxy started on %s:%d", host, port)
			utils.PrintMessage("Compute node jobs will use socks5://%s:%d", host, port)
		} else {
			utils.PrintWarning("Proxy daemon started but PID file not yet written — SSH may have failed")
			utils.PrintMessage("Check: condatainer proxy status")
		}

		// Try to enable linger so the proxy survives logout.
		// Safe when KillUserProcesses=no OR Linger=yes.
		if survives, _, _ := proxy.WillSurviveLogout(); !survives {
			user := os.Getenv("USER")
			lingerCmd := exec.Command("loginctl", "enable-linger", user)
			lingerCmd.Stderr = nil // suppress polkit "Authorization not available" noise
			if err := lingerCmd.Run(); err == nil {
				utils.PrintMessage("Enabled linger for %s — proxy will survive logout.", user)
			} else {
				utils.PrintWarning("Proxy may be killed on full logout (KillUserProcesses=yes, Linger=no, enable-linger failed).")
				utils.PrintMessage("Proxy will survive SSH disconnects but not complete logout.")
				utils.PrintMessage("Ask your sysadmin to run: loginctl enable-linger %s", user)
			}
		}
		return nil
	},
}

var proxyStopCmd = &cobra.Command{
	Use:   "stop",
	Short: "Stop the SOCKS5 proxy tunnel (works from any login node)",
	RunE: func(cmd *cobra.Command, args []string) error {
		host, _, pid, err := proxy.ReadPidFile()
		if err != nil {
			utils.PrintMessage("Proxy is not running (no PID file found)")
			return nil
		}

		currentHost, _ := os.Hostname()
		if host == currentHost {
			// Same node — kill directly
			proc, err := os.FindProcess(pid)
			if err != nil || proc.Signal(syscall.Signal(0)) == syscall.ESRCH {
				utils.PrintMessage("Proxy already stopped (stale PID file)")
				proxy.RemovePidFile()
				return nil
			}
			if err := proc.Signal(syscall.SIGTERM); err != nil {
				return fmt.Errorf("failed to stop proxy (PID %d): %w", pid, err)
			}
		} else {
			// Different login node — SSH kill
			out, err := exec.Command("ssh", host,
				"-o", "BatchMode=yes",
				"-o", "ConnectTimeout=5",
				"kill", strconv.Itoa(pid),
			).CombinedOutput()
			if err != nil {
				// If process not found, clean up stale file
				proxy.RemovePidFile()
				utils.PrintMessage("Proxy already stopped or unreachable: %s", string(out))
				return nil
			}
		}

		utils.PrintSuccess("Proxy stopped")
		return nil
	},
}

var proxyShowCmd = &cobra.Command{
	Use:   "show",
	Short: "Print the proxy URL for use by other programs",
	Run: func(cmd *cobra.Command, args []string) {
		host, port, _, err := proxy.ReadPidFile()
		if err != nil || !proxy.ProxyAlive(host, port) {
			return
		}
		fmt.Printf("socks5://%s:%d\n", host, port)
	},
}

var proxyStatusCmd = &cobra.Command{
	Use:   "status",
	Short: "Show the status of the SOCKS5 proxy tunnel",
	Run: func(cmd *cobra.Command, args []string) {
		host, port, pid, err := proxy.ReadPidFile()
		if err != nil {
			utils.PrintMessage("Proxy not running")
			return
		}

		if proxy.ProxyAlive(host, port) {
			utils.PrintSuccess("Proxy running on %s:%d (PID %d)", host, port, pid)
			utils.PrintMessage("Jobs will use: socks5://%s:%d", host, port)
		} else {
			utils.PrintWarning("Proxy dead — stale PID file at %s", proxy.PidFilePath())
			proxy.RemovePidFile()
			utils.PrintMessage("Run: condatainer proxy start")
		}
	},
}

// _proxy_daemon is an internal hidden subcommand invoked by "proxy start" after
// double-forking. It runs the SSH tunnel and blocks until the tunnel exits.
var proxyDaemonCmd = &cobra.Command{
	Use:    "_proxy_daemon",
	Hidden: true,
	Args:   cobra.ExactArgs(2),
	RunE: func(cmd *cobra.Command, args []string) error {
		host := args[0]
		port, err := strconv.Atoi(args[1])
		if err != nil {
			return fmt.Errorf("invalid port: %w", err)
		}
		return proxy.RunDaemon(host, port)
	},
}

func init() {
	proxyStartCmd.Flags().StringVar(&proxyStartHost, "host", "", "SSH host to tunnel through (default: current hostname)")
	proxyStartCmd.Flags().IntVar(&proxyStartPort, "port", 0, "Local port to listen on (default: OS-assigned free port)")

	proxyCmd.AddCommand(proxyStartCmd)
	proxyCmd.AddCommand(proxyStopCmd)
	proxyCmd.AddCommand(proxyStatusCmd)
	proxyCmd.AddCommand(proxyShowCmd)

	rootCmd.AddCommand(proxyCmd)
	rootCmd.AddCommand(proxyDaemonCmd)
}
