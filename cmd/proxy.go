package cmd

import (
	"fmt"
	"io"
	"os"
	"os/exec"
	"strconv"
	"strings"
	"syscall"

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
			if url, ok := proxy.FindActiveProxy(); ok {
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

// startProxyDaemon forks the proxy daemon and waits for it to signal readiness.
// The daemon reports back via a pipe: EOF = ready, non-empty message = error.
// localOnly=false: shared mode (0.0.0.0, NFS PID file).
// localOnly=true:  per-job mode (127.0.0.1, node-local PID file).
func startProxyDaemon(sshDest string, port int, localOnly bool) error {
	exe, err := os.Executable()
	if err != nil {
		return fmt.Errorf("failed to determine executable path: %w", err)
	}

	// Report pipe: daemon closes write end (EOF) on success or writes error message.
	pr, pw, err := os.Pipe()
	if err != nil {
		return fmt.Errorf("report pipe: %w", err)
	}

	args := []string{"_proxy_daemon", "--report-fd", "3"}
	if localOnly {
		args = append(args, "--local")
	}
	if utils.DebugMode {
		args = append(args, "--debug")
	}
	args = append(args, sshDest, strconv.Itoa(port))

	daemonCmd := exec.Command(exe, args...)
	daemonCmd.SysProcAttr = &syscall.SysProcAttr{Setsid: true}
	daemonCmd.ExtraFiles = []*os.File{pw} // ExtraFiles[0] → fd 3 in daemon
	if utils.DebugMode {
		daemonCmd.Stderr = os.Stderr
	}

	if err := daemonCmd.Start(); err != nil {
		pw.Close()
		pr.Close()
		return fmt.Errorf("failed to start proxy daemon: %w", err)
	}
	pw.Close()                  // parent closes its copy of the write end
	daemonCmd.Process.Release() //nolint:errcheck

	// Block until daemon signals ready (EOF) or reports an error.
	msg, _ := io.ReadAll(pr)
	pr.Close()
	if len(msg) > 0 {
		return fmt.Errorf("%s", strings.TrimSpace(string(msg)))
	}

	if localOnly {
		utils.PrintSuccess("Per-job proxy started on 127.0.0.1:%d (via %s)", port, sshDest)
		utils.PrintMessage("This job will use http://127.0.0.1:%d", port)
		return nil
	}

	currentHost, _ := os.Hostname()
	utils.PrintSuccess("Proxy started on %s:%d", currentHost, port)
	utils.PrintMessage("Compute node jobs will use http://%s:%d", currentHost, port)

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

var proxyShowExport bool

var proxyShowCmd = &cobra.Command{
	Use:   "show",
	Short: "Print the proxy URL for use by other programs",
	Long: `Print the active proxy URL.

Without flags: prints the HTTP CONNECT URL (http://host:PORT).

With --export: prints shell export statements for all proxy env vars.
Suitable for: eval $(condatainer proxy show --export)

If no proxy is active, prints nothing.`,
	Run: func(cmd *cobra.Command, args []string) {
		httpURL, ok := proxy.FindActiveProxy()
		if !ok {
			return
		}
		if !proxyShowExport {
			fmt.Println(httpURL)
			return
		}
		for _, kv := range proxy.ProxyEnvList(httpURL) {
			fmt.Printf("export %s\n", kv)
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
				utils.PrintMessage("This job uses: http://%s:%d", ps.Host, ps.Port)
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
			utils.PrintMessage("Jobs will use: http://%s:%d", ps.Host, ps.Port)
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
var proxyDaemonReportFd int

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
		return proxy.RunDaemon(sshDest, port, proxyDaemonLocalMode, proxyDaemonReportFd)
	},
}

func init() {
	proxyStartCmd.Flags().StringVar(&proxyStartHost, "host", "", "Node to run the daemon on — delegates via SSH if different from current node (login node only)")
	proxyStartCmd.Flags().StringVar(&proxyStartVia, "via", "", "SSH server to tunnel through (default: current node for shared mode; required for per-job mode inside a job)")
	proxyStartCmd.Flags().IntVar(&proxyStartPort, "port", 0, "Local port to listen on (default: OS-assigned free port)")

	proxyShowCmd.Flags().BoolVar(&proxyShowExport, "export", false, "Print shell export statements for all proxy env vars (eval-friendly)")
	proxyDaemonCmd.Flags().BoolVar(&proxyDaemonLocalMode, "local", false, "Per-job mode: bind 127.0.0.1 and write node-local PID file")
	proxyDaemonCmd.Flags().IntVar(&proxyDaemonReportFd, "report-fd", 0, "File descriptor to report startup result to parent")

	proxyCmd.AddCommand(proxyStartCmd)
	proxyCmd.AddCommand(proxyStopCmd)
	proxyCmd.AddCommand(proxyStatusCmd)
	proxyCmd.AddCommand(proxyShowCmd)

	rootCmd.AddCommand(proxyCmd)
	rootCmd.AddCommand(proxyDaemonCmd)
}
