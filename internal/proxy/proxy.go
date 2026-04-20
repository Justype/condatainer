package proxy

import (
	"encoding/json"
	"errors"
	"fmt"
	"net"
	"os"
	"os/exec"
	"path/filepath"
	"strconv"
	"strings"
	"sync"
	"syscall"
	"time"

	"github.com/Justype/condatainer/internal/config"
	"github.com/Justype/condatainer/internal/scheduler"
	"github.com/Justype/condatainer/internal/utils"
)

// PidFilePath returns the NFS-shared proxy PID file path (shared mode).
func PidFilePath() string {
	return filepath.Join(config.GetUserDataDir(), "proxy.pid")
}

// LocalPidFilePath returns the node-local proxy PID file path (per-job mode).
// Uses utils.GetTmpDir() which respects SLURM_TMPDIR, PBS_TMPDIR, etc.
func LocalPidFilePath() string {
	return filepath.Join(utils.GetTmpDir(), "proxy.pid")
}

// ProxyState holds the state written to a proxy PID file.
type ProxyState struct {
	Host string `json:"host"`
	Via  string `json:"via"`
	Port int    `json:"port"`
	PID  int    `json:"pid"`
}

// WritePidFileAt writes proxy state to an arbitrary PID file path as JSON.
func WritePidFileAt(path string, ps ProxyState) error {
	if err := os.MkdirAll(filepath.Dir(path), 0775); err != nil {
		return err
	}
	data, err := json.Marshal(ps)
	if err != nil {
		return err
	}
	return os.WriteFile(path, data, 0600)
}

// ReadPidFileAt reads and parses a JSON proxy PID file.
func ReadPidFileAt(path string) (*ProxyState, error) {
	data, err := os.ReadFile(path)
	if err != nil {
		return nil, err
	}
	var ps ProxyState
	if err := json.Unmarshal(data, &ps); err != nil {
		return nil, fmt.Errorf("malformed proxy PID file: %w", err)
	}
	return &ps, nil
}

// ReadPidFile reads the shared NFS proxy PID file.
func ReadPidFile() (*ProxyState, error) {
	return ReadPidFileAt(PidFilePath())
}

// ReadLocalPidFile reads the node-local per-job proxy PID file.
func ReadLocalPidFile() (*ProxyState, error) {
	return ReadPidFileAt(LocalPidFilePath())
}

// RemovePidFile removes the shared NFS proxy PID file.
func RemovePidFile() {
	os.Remove(PidFilePath()) //nolint:errcheck
}

// ProxyAlive checks if the SOCKS5 proxy is reachable at host:port.
func ProxyAlive(host string, port int) bool {
	conn, err := net.DialTimeout("tcp",
		fmt.Sprintf("%s:%d", host, port),
		500*time.Millisecond)
	if err != nil {
		return false
	}
	conn.Close()
	return true
}

// FreePort asks the OS for an available TCP port on all interfaces and returns it.
func FreePort() (int, error) {
	ln, err := net.Listen("tcp", "0.0.0.0:0")
	if err != nil {
		return 0, fmt.Errorf("failed to find free port: %w", err)
	}
	port := ln.Addr().(*net.TCPAddr).Port
	ln.Close()
	return port, nil
}

// ResolveViaHost resolves the SSH server to tunnel through for a per-job proxy.
// Priority: --via flag > CNT_PROXY_VIA env (baked at submit time by condatainer).
func ResolveViaHost(flagVia string) (string, error) {
	if flagVia != "" {
		return flagVia, nil
	}
	if h := os.Getenv("CNT_PROXY_VIA"); h != "" {
		return h, nil
	}
	return "", errors.New("cannot determine login node: use --via, or submit jobs via condatainer so CNT_PROXY_VIA is set automatically")
}

// StartOnNode SSHes to host and runs "condatainer proxy start" there (delegation).
// Passes --via and --port if set. Waits up to 10s for the NFS PID file to appear.
func StartOnNode(host, via string, port int) error {
	args := []string{
		host,
		"-o", "BatchMode=yes",
		"-o", "ConnectTimeout=10",
		"condatainer", "proxy", "start",
	}
	if via != "" {
		args = append(args, "--via", via)
	}
	if port > 0 {
		args = append(args, "--port", strconv.Itoa(port))
	}
	out, err := exec.Command("ssh", args...).CombinedOutput()
	if err != nil {
		return fmt.Errorf("failed to start proxy on %s: %w\n%s", host, err, strings.TrimSpace(string(out)))
	}
	// Wait up to 10s for NFS PID file to appear (NFS flush lag)
	for range 20 {
		if _, err := ReadPidFile(); err == nil {
			return nil
		}
		time.Sleep(500 * time.Millisecond)
	}
	return fmt.Errorf("proxy started on %s but PID file not yet visible — NFS lag?", host)
}

// FindActiveProxy returns the proxy URL to use for the current job.
// Checks the local per-job PID file first, then the shared NFS PID file.
// If autostart is true and no active proxy is found, attempts to start a
// per-job local proxy using ResolveViaHost.
func FindActiveProxy(autostart bool) (string, bool) {
	// 1. Local per-job proxy (127.0.0.1:PORT — this compute node only)
	if ps, err := ReadLocalPidFile(); err == nil && ProxyAlive(ps.Host, ps.Port) {
		return fmt.Sprintf("socks5h://%s:%d", ps.Host, ps.Port), true
	}
	// 2. Shared NFS proxy (loginNode:PORT — all compute nodes)
	if ps, err := ReadPidFile(); err == nil {
		if ProxyAlive(ps.Host, ps.Port) {
			return fmt.Sprintf("socks5h://%s:%d", ps.Host, ps.Port), true
		}
		utils.PrintWarning("Shared proxy on %s:%d is unreachable", ps.Host, ps.Port)
	}
	// 3. Auto-start per-job proxy if enabled and login node is known
	if autostart {
		if via, err := ResolveViaHost(""); err == nil {
			if url, ok := autoStartLocalProxy(via); ok {
				return url, true
			}
			utils.PrintWarning("Failed to auto-start per-job proxy via %s — SSH may have failed", via)
		}
	}
	return "", false
}

var (
	jobProxyOnce sync.Once
	jobProxyURL  string
)

// GetJobProxy returns the active proxy URL for injection into container envs and
// HTTP clients inside a scheduler job. Returns "", false on login nodes — they
// have direct internet access and need no proxy.
// Result is cached after the first call.
func GetJobProxy() (string, bool) {
	jobProxyOnce.Do(func() {
		if !scheduler.IsInsideJob() {
			return
		}
		if u, ok := FindActiveProxy(config.Global.ProxyPerJob); ok {
			jobProxyURL = u
		}
	})
	if jobProxyURL == "" {
		return "", false
	}
	return jobProxyURL, true
}

// autoStartLocalProxy double-forks a per-job daemon (127.0.0.1, localOnly=true)
// tunneling through via, waits for the local PID file, and returns the proxy URL.
func autoStartLocalProxy(via string) (string, bool) {
	port, err := FreePort()
	if err != nil {
		return "", false
	}

	exe, err := os.Executable()
	if err != nil {
		return "", false
	}

	daemonCmd := exec.Command(exe, "_proxy_daemon", "--local", via, strconv.Itoa(port))
	daemonCmd.SysProcAttr = &syscall.SysProcAttr{Setsid: true}

	if err := daemonCmd.Start(); err != nil {
		return "", false
	}
	daemonCmd.Process.Release() //nolint:errcheck

	// Wait up to 2s for local PID file
	for range 20 {
		if ps, err := ReadLocalPidFile(); err == nil && ProxyAlive(ps.Host, ps.Port) {
			return fmt.Sprintf("socks5h://%s:%d", ps.Host, ps.Port), true
		}
		time.Sleep(100 * time.Millisecond)
	}
	return "", false
}

// WillSurviveLogout returns true if the proxy daemon is expected to survive user logout.
// Safe when KillUserProcesses=no OR Linger=yes.
// Returns (survives, killUserProcesses, linger).
func WillSurviveLogout() (survives, killUserProcesses, linger bool) {
	killUserProcesses = readKillUserProcesses()
	if !killUserProcesses {
		return true, false, false
	}
	linger = readLinger()
	return linger, killUserProcesses, linger
}

// readKillUserProcesses returns true if systemd-logind is configured to kill user
// processes on logout. Tries loginctl first, then falls back to config files.
// Defaults to true (conservative) when the value cannot be determined.
func readKillUserProcesses() bool {
	// Try loginctl show (works on newer systemd)
	if out, err := exec.Command("loginctl", "show",
		"--property=KillUserProcesses", "--value").Output(); err == nil {
		v := strings.TrimSpace(string(out))
		if v == "yes" {
			return true
		}
		if v == "no" {
			return false
		}
		// empty — property not exposed by this systemd version, fall through
	}

	// Fall back: read config files directly
	for _, path := range []string{
		"/etc/systemd/logind.conf",
		"/etc/systemd/logind.conf.d/",
		"/usr/lib/systemd/logind.conf.d/",
	} {
		if v, found := grepLogindConfig(path, "KillUserProcesses"); found {
			return v == "yes"
		}
	}

	return true // conservative default: assume yes
}

// readLinger returns true if loginctl linger is enabled for the current user.
func readLinger() bool {
	user := os.Getenv("USER")
	if user == "" {
		return false
	}
	out, err := exec.Command("loginctl", "show-user", user,
		"--property=Linger", "--value").Output()
	if err != nil {
		return false
	}
	return strings.TrimSpace(string(out)) == "yes"
}

// grepLogindConfig reads KillUserProcesses from a logind config file or directory.
func grepLogindConfig(path, key string) (value string, found bool) {
	entries, err := os.ReadDir(path)
	if err != nil {
		// Try as a plain file
		return grepLogindFile(path, key)
	}
	// Directory: scan all .conf files
	for _, e := range entries {
		if e.IsDir() || !strings.HasSuffix(e.Name(), ".conf") {
			continue
		}
		if v, ok := grepLogindFile(filepath.Join(path, e.Name()), key); ok {
			return v, true
		}
	}
	return "", false
}

func grepLogindFile(path, key string) (value string, found bool) {
	data, err := os.ReadFile(path)
	if err != nil {
		return "", false
	}
	for line := range strings.SplitSeq(string(data), "\n") {
		line = strings.TrimSpace(line)
		if strings.HasPrefix(line, "#") {
			continue
		}
		if v, ok := strings.CutPrefix(line, key+"="); ok {
			return v, true
		}
	}
	return "", false
}
