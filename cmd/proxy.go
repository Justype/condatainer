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
	Short: "Manage SOCKS5 proxy tunnels for compute nodes",
	Long: `Start and manage SSH SOCKS5 proxy tunnels so that compute node jobs
can reach the internet through the login node.

Shared mode (login node): run "condatainer proxy start" once on the login node.
All condatainer commands running inside jobs automatically use the tunnel.

Per-job mode (compute node): run "condatainer proxy start --via <login-node>" inside
a job, or set proxy_perjob=true in config for automatic per-job proxy startup.`,
}

// proxy start flags
var (
	proxyStartHost string
	proxyStartVia  string
	proxyStartPort int
)

var proxyStartCmd = &cobra.Command{
	Use:          "start",
	Short:        "Start a SOCKS5 proxy tunnel",
	SilenceUsage: true,
	RunE: func(cmd *cobra.Command, args []string) error {
		if config.IsInsideContainer() {
			return fmt.Errorf("proxy must be started outside a container")
		}
		if scheduler.ActiveScheduler() == nil {
			return fmt.Errorf("no scheduler detected — proxy is only needed on HPC clusters")
		}

		currentHost, _ := os.Hostname()
		isInsideJob := scheduler.ActiveScheduler().IsInsideJob()

		if !isInsideJob {
			// ── Login node ──────────────────────────────────────────────────
			targetHost := proxyStartHost
			if targetHost == "" {
				targetHost = currentHost
			}

			// Delegate to another login node if --host points elsewhere
			if targetHost != currentHost {
				utils.PrintMessage("Delegating proxy start to %s…", targetHost)
				return proxy.StartOnNode(targetHost, proxyStartVia, proxyStartPort)
			}

			// Start shared proxy on this node
			// Check if already running
			if ps, err := proxy.ReadPidFile(); err == nil {
				if proxy.ProxyAlive(ps.Host, ps.Port) {
					utils.PrintMessage("Proxy already running on %s:%d (PID %d)", ps.Host, ps.Port, ps.PID)
					return nil
				}
				proxy.RemovePidFile()
			}

			sshDest := proxyStartVia
			if sshDest == "" {
				sshDest = currentHost
			}

			port := proxyStartPort
			if port == 0 {
				var err error
				port, err = proxy.FreePort()
				if err != nil {
					return fmt.Errorf("failed to find free port: %w", err)
				}
			}

			return startProxyDaemon(sshDest, port, false)

		} else {
			// ── Compute node (inside job) ────────────────────────────────────
			// Check if already running (local or shared)
			if url, ok := proxy.FindActiveProxy(false); ok {
				utils.PrintMessage("Proxy already active: %s", url)
				return nil
			}

			viaHost, err := proxy.ResolveViaHost(proxyStartVia)
			if err != nil {
				return err
			}

			port := proxyStartPort
			if port == 0 {
				port, err = proxy.FreePort()
				if err != nil {
					return fmt.Errorf("failed to find free port: %w", err)
				}
			}

			return startProxyDaemon(viaHost, port, true)
		}
	},
}

// startProxyDaemon double-forks the proxy daemon and waits for the PID file to appear.
// localOnly=false: shared mode (0.0.0.0, NFS PID file).
// localOnly=true:  per-job mode (127.0.0.1, node-local PID file).
func startProxyDaemon(sshDest string, port int, localOnly bool) error {
	exe, err := os.Executable()
	if err != nil {
		return fmt.Errorf("failed to determine executable path: %w", err)
	}

	args := []string{"_proxy_daemon"}
	if localOnly {
		args = append(args, "--local")
	}
	args = append(args, sshDest, strconv.Itoa(port))

	daemonCmd := exec.Command(exe, args...)
	daemonCmd.SysProcAttr = &syscall.SysProcAttr{Setsid: true}

	if err := daemonCmd.Start(); err != nil {
		return fmt.Errorf("failed to start proxy daemon: %w", err)
	}
	daemonCmd.Process.Release() //nolint:errcheck

	readFn := proxy.ReadPidFile
	if localOnly {
		readFn = proxy.ReadLocalPidFile
	}

	pidFileReady := false
	for i := 0; i < 20; i++ {
		if _, err := readFn(); err == nil {
			pidFileReady = true
			break
		}
		time.Sleep(100 * time.Millisecond)
	}

	if !pidFileReady {
		utils.PrintWarning("Proxy daemon started but PID file not yet written — SSH may have failed")
		utils.PrintMessage("Check: condatainer proxy status")
		return nil
	}

	if localOnly {
		utils.PrintSuccess("Per-job proxy started on 127.0.0.1:%d (via %s)", port, sshDest)
		utils.PrintMessage("This job will use socks5://127.0.0.1:%d", port)
	} else {
		currentHost, _ := os.Hostname()
		utils.PrintSuccess("Proxy started on %s:%d", currentHost, port)
		utils.PrintMessage("Compute node jobs will use socks5://%s:%d", currentHost, port)

		// Try to enable linger so the proxy survives logout.
		if survives, _, _ := proxy.WillSurviveLogout(); !survives {
			user := os.Getenv("USER")
			lingerCmd := exec.Command("loginctl", "enable-linger", user)
			lingerCmd.Stderr = nil
			if err := lingerCmd.Run(); err == nil {
				utils.PrintMessage("Enabled linger for %s — proxy will survive logout.", user)
			} else {
				utils.PrintWarning("Proxy may be killed on full logout (KillUserProcesses=yes, Linger=no, enable-linger failed).")
				utils.PrintMessage("Proxy will survive SSH disconnects but not complete logout.")
				utils.PrintMessage("Ask your sysadmin to run: loginctl enable-linger %s", user)
			}
		}
	}
	return nil
}

var proxyStopCmd = &cobra.Command{
	Use:          "stop",
	Short:        "Stop the shared SOCKS5 proxy tunnel (works from any login node)",
	SilenceUsage: true,
	RunE: func(cmd *cobra.Command, args []string) error {
		ps, err := proxy.ReadPidFile()
		if err != nil {
			utils.PrintMessage("Proxy is not running (no PID file found)")
			return nil
		}

		currentHost, _ := os.Hostname()
		if ps.Host == currentHost {
			// Same node — kill directly
			proc, err := os.FindProcess(ps.PID)
			if err != nil || proc.Signal(syscall.Signal(0)) == syscall.ESRCH {
				utils.PrintMessage("Proxy already stopped (stale PID file)")
				proxy.RemovePidFile()
				return nil
			}
			if err := proc.Signal(syscall.SIGTERM); err != nil {
				return fmt.Errorf("failed to stop proxy (PID %d): %w", ps.PID, err)
			}
		} else {
			// Different login node — SSH kill
			out, err := exec.Command("ssh", ps.Host,
				"-o", "BatchMode=yes",
				"-o", "ConnectTimeout=5",
				"kill", strconv.Itoa(ps.PID),
			).CombinedOutput()
			if err != nil {
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
		if url, ok := proxy.FindActiveProxy(false); ok {
			fmt.Println(url)
		}
	},
}

var proxyStatusCmd = &cobra.Command{
	Use:   "status",
	Short: "Show the status of the SOCKS5 proxy tunnel",
	Run: func(cmd *cobra.Command, args []string) {
		// Per-job local proxy
		if ps, err := proxy.ReadLocalPidFile(); err == nil {
			if proxy.ProxyAlive(ps.Host, ps.Port) {
				if ps.Via != ps.Host {
					utils.PrintSuccess("Per-job proxy running on %s:%d (PID %d, via %s)", ps.Host, ps.Port, ps.PID, ps.Via)
				} else {
					utils.PrintSuccess("Per-job proxy running on %s:%d (PID %d)", ps.Host, ps.Port, ps.PID)
				}
				utils.PrintMessage("This job uses: socks5://%s:%d", ps.Host, ps.Port)
				return
			}
			utils.PrintWarning("Per-job proxy dead — stale PID file at %s", proxy.LocalPidFilePath())
			os.Remove(proxy.LocalPidFilePath()) //nolint:errcheck
		}

		// Shared NFS proxy
		ps, err := proxy.ReadPidFile()
		if err != nil {
			utils.PrintMessage("Proxy not running")
			return
		}
		if proxy.ProxyAlive(ps.Host, ps.Port) {
			if ps.Via != ps.Host {
				utils.PrintSuccess("Shared proxy running on %s:%d (PID %d, via %s)", ps.Host, ps.Port, ps.PID, ps.Via)
			} else {
				utils.PrintSuccess("Shared proxy running on %s:%d (PID %d)", ps.Host, ps.Port, ps.PID)
			}
			utils.PrintMessage("Jobs will use: socks5://%s:%d", ps.Host, ps.Port)
		} else {
			utils.PrintWarning("Proxy dead — stale PID file at %s", proxy.PidFilePath())
			proxy.RemovePidFile()
			utils.PrintMessage("Run: condatainer proxy start")
		}
	},
}

// _proxy_daemon is an internal hidden subcommand invoked by "proxy start" after
// double-forking. It runs the SSH tunnel and blocks until the tunnel exits.
// Without --local: shared mode (0.0.0.0, NFS PID file).
// With --local:    per-job mode (127.0.0.1, node-local PID file).
var proxyDaemonLocalMode bool

var proxyDaemonCmd = &cobra.Command{
	Use:    "_proxy_daemon",
	Hidden: true,
	Args:   cobra.ExactArgs(2),
	RunE: func(cmd *cobra.Command, args []string) error {
		sshDest := args[0]
		port, err := strconv.Atoi(args[1])
		if err != nil {
			return fmt.Errorf("invalid port: %w", err)
		}
		return proxy.RunDaemon(sshDest, port, proxyDaemonLocalMode)
	},
}

func init() {
	proxyStartCmd.Flags().StringVar(&proxyStartHost, "host", "", "Node to run the daemon on — delegates via SSH if different from current node (login node only)")
	proxyStartCmd.Flags().StringVar(&proxyStartVia, "via", "", "SSH server to tunnel through (default: current node for shared; auto-detected via CNT_PROXY_VIA for per-job)")
	proxyStartCmd.Flags().IntVar(&proxyStartPort, "port", 0, "Local port to listen on (default: OS-assigned free port)")

	proxyDaemonCmd.Flags().BoolVar(&proxyDaemonLocalMode, "local", false, "Per-job mode: bind 127.0.0.1 and write node-local PID file")

	proxyCmd.AddCommand(proxyStartCmd)
	proxyCmd.AddCommand(proxyStopCmd)
	proxyCmd.AddCommand(proxyStatusCmd)
	proxyCmd.AddCommand(proxyShowCmd)

	rootCmd.AddCommand(proxyCmd)
	rootCmd.AddCommand(proxyDaemonCmd)
}
