package proxy

import (
	"fmt"
	"net"
	"os"
	"os/exec"
	"path/filepath"
	"strconv"
	"strings"
	"syscall"
	"time"

	"github.com/Justype/condatainer/internal/config"
)

// PidFilePath returns the path to the proxy PID file in the user data directory.
func PidFilePath() string {
	return filepath.Join(config.GetUserDataDir(), "proxy.pid")
}

// WritePidFile writes host, port, and pid to the proxy PID file (NFS-visible).
func WritePidFile(host string, port, pid int) error {
	path := PidFilePath()
	if err := os.MkdirAll(filepath.Dir(path), 0775); err != nil {
		return err
	}
	content := fmt.Sprintf("%s\n%d\n%d\n", host, port, pid)
	return os.WriteFile(path, []byte(content), 0600)
}

// ReadPidFile reads the proxy PID file and returns host, port, pid.
// Returns an error if the file does not exist or cannot be parsed.
func ReadPidFile() (host string, port, pid int, err error) {
	data, err := os.ReadFile(PidFilePath())
	if err != nil {
		return "", 0, 0, err
	}
	parts := strings.SplitN(strings.TrimSpace(string(data)), "\n", 3)
	if len(parts) != 3 {
		return "", 0, 0, fmt.Errorf("malformed proxy PID file")
	}
	host = strings.TrimSpace(parts[0])
	port, err = strconv.Atoi(strings.TrimSpace(parts[1]))
	if err != nil {
		return "", 0, 0, fmt.Errorf("invalid port in proxy PID file: %w", err)
	}
	pid, err = strconv.Atoi(strings.TrimSpace(parts[2]))
	if err != nil {
		return "", 0, 0, fmt.Errorf("invalid PID in proxy PID file: %w", err)
	}
	return host, port, pid, nil
}

// RemovePidFile removes the proxy PID file.
func RemovePidFile() {
	os.Remove(PidFilePath()) //nolint:errcheck
}

// ProxyAlive checks if the SOCKS5 proxy is reachable at host:port.
func ProxyAlive(host string, port int) bool {
	conn, err := net.DialTimeout("tcp",
		fmt.Sprintf("%s:%d", host, port),
		2*time.Second)
	if err != nil {
		return false
	}
	conn.Close()
	return true
}

// ProcessAlive checks if a process with the given PID is still running on the current node.
func ProcessAlive(pid int) bool {
	proc, err := os.FindProcess(pid)
	if err != nil {
		return false
	}
	err = proc.Signal(syscall.Signal(0))
	return err == nil || err == syscall.EPERM
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

// RestartProxy SSHes to host and runs "condatainer proxy start" to bring the tunnel
// back up on the same port. Waits up to 3 seconds for the proxy to become reachable.
// Returns true if the proxy is alive after the restart attempt.
func RestartProxy(host string, port int) bool {
	exec.Command("ssh", host, //nolint:errcheck
		"-o", "BatchMode=yes",
		"-o", "ConnectTimeout=5",
		"condatainer", "proxy", "start",
		"--host", host,
		"--port", strconv.Itoa(port),
	).Run()

	for range 6 {
		time.Sleep(500 * time.Millisecond)
		if ProxyAlive(host, port) {
			return true
		}
	}
	return false
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
	for _, line := range strings.Split(string(data), "\n") {
		line = strings.TrimSpace(line)
		if strings.HasPrefix(line, "#") {
			continue
		}
		if strings.HasPrefix(line, key+"=") {
			return strings.TrimPrefix(line, key+"="), true
		}
	}
	return "", false
}
