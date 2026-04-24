package proxy

import (
	"context"
	"fmt"
	"os"
	"os/exec"
	"os/signal"
	"strconv"
	"strings"
	"syscall"
)

// RunDaemon starts a dual-protocol (SOCKS5 + HTTP CONNECT) proxy on the visible port.
// It creates an internal SSH SOCKS5 tunnel on a loopback-only port and runs a Go
// dual-protocol proxy in front of it on the published port.
//
// localOnly=false: shared mode — visible proxy binds 0.0.0.0:port, NFS PID file.
// localOnly=true:  per-job mode — visible proxy binds 127.0.0.1:port, node-local PID file.
//
// The daemon blocks until SIGTERM/SIGINT, then cleans up the PID file.
func RunDaemon(sshDest string, port int, localOnly bool) error {
	// Allocate an internal loopback port for the SSH SOCKS5 tunnel.
	// This port is never published — it's only used by the Go dual-proxy internally.
	hiddenPort, err := FreePort()
	if err != nil {
		return fmt.Errorf("failed to allocate internal port: %w", err)
	}

	// Start SSH SOCKS5 tunnel on loopback only.
	ssh := exec.Command("ssh",
		"-N",
		"-D", fmt.Sprintf("127.0.0.1:%d", hiddenPort),
		sshDest,
		"-o", "BatchMode=yes",
		"-o", "StrictHostKeyChecking=no",
		"-o", "ServerAliveInterval=60",
		"-o", "ServerAliveCountMax=3",
	)
	if err := ssh.Start(); err != nil {
		return fmt.Errorf("failed to start SSH tunnel: %w", err)
	}
	sshDone := make(chan error, 1)
	go func() { sshDone <- ssh.Wait() }()

	bind := "0.0.0.0"
	if localOnly {
		bind = "127.0.0.1"
	}

	var pidHost, pidPath string
	if localOnly {
		pidHost = "127.0.0.1"
		pidPath = LocalPidFilePath()
	} else {
		pidHost, _ = os.Hostname()
		pidPath = PidFilePath()
	}

	if err := WritePidFileAt(pidPath, ProxyState{
		Host: pidHost,
		Via:  sshDest,
		Port: port,
		PID:  os.Getpid(),
	}); err != nil {
		ssh.Process.Kill() //nolint:errcheck
		return fmt.Errorf("failed to write PID file: %w", err)
	}
	defer os.Remove(pidPath)

	ctx, cancel := context.WithCancel(context.Background())
	sig := make(chan os.Signal, 1)
	signal.Notify(sig, syscall.SIGTERM, syscall.SIGINT)
	go func() {
		select {
		case <-sig:
		case <-sshDone: // SSH tunnel exited unexpectedly
		}
		cancel()
		ssh.Process.Kill() //nolint:errcheck
	}()

	socks5Addr := "127.0.0.1:" + strconv.Itoa(hiddenPort)
	return RunProxy(ctx, fmt.Sprintf("%s:%d", bind, port), socks5Addr)
}

// TestSSH verifies SSH connectivity to sshDest using the system ssh binary.
// Uses BatchMode so it never prompts — mirrors exactly what the daemon tunnel does.
func TestSSH(sshDest string) error {
	out, err := exec.Command("ssh",
		"-o", "BatchMode=yes",
		"-o", "ConnectTimeout=10",
		"-o", "StrictHostKeyChecking=no",
		sshDest, "true",
	).CombinedOutput()
	if err != nil {
		msg := strings.TrimSpace(string(out))
		if msg == "" {
			return fmt.Errorf("SSH connection to %s failed", sshDest)
		}
		return fmt.Errorf("SSH connection to %s failed: %s", sshDest, msg)
	}
	return nil
}
