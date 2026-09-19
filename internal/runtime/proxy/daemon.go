package proxy

import (
	"bytes"
	"context"
	"fmt"
	"io"
	"net"
	"os"
	"os/exec"
	"os/signal"
	"strconv"
	"strings"
	"sync"
	"syscall"
	"time"

	"log/slog"
)

// RunDaemon starts a dual-protocol (SOCKS5 + HTTP CONNECT) proxy on the visible port.
// It creates an SSH tunnel and runs a Go dual-protocol proxy in front of it.
//
// localOnly=false: shared mode — visible proxy binds 0.0.0.0:port, NFS PID file.
// localOnly=true:  per-job mode — visible proxy binds 127.0.0.1:port, node-local PID file.
//
// reportFd, if > 0, is a file descriptor for a pipe back to the parent process:
// the daemon closes it (EOF) once the tunnel and PID file are ready, or writes an
// error message before closing if startup fails. The parent blocks on this pipe
// instead of polling the PID file.
//
// The daemon blocks until SIGTERM/SIGINT, then cleans up the PID file and socket.
func RunDaemon(sshDest string, port int, localOnly bool, reportFd int) error {
	var sockPath, pidPath, pidHost string
	if localOnly {
		sockPath = LocalSockPath()
		pidPath = LocalPidFilePath()
		pidHost = "127.0.0.1"
	} else {
		sockPath = SockPath()
		pidPath = PidFilePath()
		pidHost, _ = os.Hostname()
	}

	// report writes msg to the parent pipe and closes it. Empty msg = success (EOF).
	var reportPipe *os.File
	if reportFd > 0 {
		reportPipe = os.NewFile(uintptr(reportFd), "report")
	}
	report := func(msg string) {
		if reportPipe == nil {
			return
		}
		if msg != "" {
			fmt.Fprint(reportPipe, msg) //nolint:errcheck
		}
		reportPipe.Close()
		reportPipe = nil
	}

	dial, tunnelStop, tunnelDone, method, err := EstablishTunnel(sshDest, sockPath)
	if err != nil {
		report(fmt.Sprintf("no working tunnel to %s: %v", sshDest, err))
		return fmt.Errorf("no working tunnel to %s: %w", sshDest, err)
	}
	defer tunnelStop()
	slog.Default().Debug("proxy tunnel established", "method", method)

	if err := WritePidFileAt(pidPath, ProxyState{
		Host: pidHost,
		Via:  sshDest,
		Port: port,
		PID:  os.Getpid(),
	}); err != nil {
		report(fmt.Sprintf("failed to write PID file: %v", err))
		return fmt.Errorf("failed to write PID file: %w", err)
	}
	defer os.Remove(pidPath)

	// Signal parent: tunnel up and PID file written.
	report("")

	ctx, cancel := context.WithCancel(context.Background())
	sig := make(chan os.Signal, 1)
	signal.Notify(sig, syscall.SIGTERM, syscall.SIGINT)
	go func() {
		select {
		case <-sig:
		case <-tunnelDone: // tunnel exited unexpectedly
		}
		cancel()
		tunnelStop()
	}()

	bind := "0.0.0.0"
	if localOnly {
		bind = "127.0.0.1"
	}
	return RunProxy(ctx, fmt.Sprintf("%s:%d", bind, port), dial)
}

// tunnelStartTimeout is how long a system ssh -D gets to open its listener.
const tunnelStartTimeout = 10 * time.Second

// startSocksSSH runs `ssh -D listen <sshArgs>` and waits for the listener at
// (network, addr). It gives up as soon as ssh exits, naming what ssh printed,
// rather than waiting out the timeout. cleanup runs on stop and on failure.
func startSocksSSH(listen, network, addr string, sshArgs []string, cleanup func()) (DialFunc, func(), <-chan struct{}, error) {
	var stderr bytes.Buffer
	cmd := exec.Command("ssh", append([]string{"-D", listen}, sshArgs...)...)
	// Pdeathsig: if the daemon itself dies (crash, OOM-kill) before stop()
	// runs, the kernel kills this ssh process directly instead of leaving
	// it holding the tunnel open. ssh never escalates privilege, so
	// unlike freeze.MountedRun's namespaced mount this needs no sentinel.
	cmd.SysProcAttr = &syscall.SysProcAttr{Pdeathsig: syscall.SIGKILL}
	cmd.Stderr = &stderr
	if err := cmd.Start(); err != nil {
		return nil, nil, nil, err
	}
	exited := make(chan struct{})
	go func() { cmd.Wait(); close(exited) }() //nolint:errcheck

	var once sync.Once
	stop := func() {
		once.Do(func() {
			cmd.Process.Kill() //nolint:errcheck
			cleanup()
		})
	}

	deadline := time.Now().Add(tunnelStartTimeout)
	for {
		if c, err := net.DialTimeout(network, addr, 200*time.Millisecond); err == nil {
			c.Close()
			return dialViaSocks5(network, addr), stop, exited, nil
		}
		if time.Now().After(deadline) {
			stop()
			return nil, nil, nil, fmt.Errorf("listener did not open within %s", tunnelStartTimeout)
		}
		select {
		case <-exited:
			stop()
			return nil, nil, nil, fmt.Errorf("ssh exited: %s", strings.TrimSpace(stderr.String()))
		case <-time.After(100 * time.Millisecond):
		}
	}
}

// EstablishTunnel attempts three methods to create a tunnel through sshDest,
// returning a DialFunc, a cleanup stop(), a done channel (closed when the tunnel
// exits unexpectedly), and a label for logging. sockPath is where the
// `ssh -D` Unix socket is created, so concurrent tunnels need distinct paths.
//
// Priority:
//  1. Go SSH library with non-interactive auth (hostbased or publickey)
//  2. System ssh -D Unix socket (current approach, needs SSH keys)
//  3. System ssh -D TCP port (127.0.0.1, fallback for older OpenSSH)
func EstablishTunnel(sshDest, sockPath string) (DialFunc, func(), <-chan struct{}, string, error) {
	var errs []string

	// Option 1: Go SSH with non-interactive auth
	slog.Default().Debug("proxy tunnel: trying go-ssh", "dest", sshDest)
	if dial, stop, done, err := DialGoSSH(sshDest); err == nil {
		return dial, stop, done, "go-ssh", nil
	} else {
		slog.Default().Debug("proxy tunnel: go-ssh failed", "err", err)
		errs = append(errs, "go-ssh: "+err.Error())
	}

	sshArgs := []string{
		"-N",
		"-o", "BatchMode=yes",
		"-o", "StrictHostKeyChecking=no",
		"-o", "ServerAliveInterval=60",
		"-o", "ServerAliveCountMax=3",
		sshDest,
	}

	// Option 2: system ssh -D Unix socket
	slog.Default().Debug("proxy tunnel: trying ssh/unix-sock", "dest", sshDest, "sock", sockPath)
	os.Remove(sockPath)                                                                                         //nolint:errcheck
	dial, stop, done, err := startSocksSSH(sockPath, "unix", sockPath, sshArgs, func() { os.Remove(sockPath) }) //nolint:errcheck
	if err == nil {
		return dial, stop, done, "ssh/unix-sock", nil
	}
	slog.Default().Debug("proxy tunnel: ssh/unix-sock failed", "err", err)
	errs = append(errs, "ssh/unix-sock: "+err.Error())

	// Option 3: system ssh -D TCP port (127.0.0.1)
	slog.Default().Debug("proxy tunnel: trying ssh/tcp", "dest", sshDest)
	innerPort, err := FreePort()
	if err != nil {
		errs = append(errs, "ssh/tcp free port: "+err.Error())
		return nil, nil, nil, "", fmt.Errorf("all tunnel methods failed: %s", strings.Join(errs, "; "))
	}
	innerAddr := fmt.Sprintf("127.0.0.1:%d", innerPort)
	dial, stop, done, err = startSocksSSH(innerAddr, "tcp", innerAddr, sshArgs, func() {})
	if err != nil {
		slog.Default().Debug("proxy tunnel: ssh/tcp failed", "err", err)
		errs = append(errs, "ssh/tcp: "+err.Error())
		return nil, nil, nil, "", fmt.Errorf("all tunnel methods failed: %s", strings.Join(errs, "; "))
	}
	return dial, stop, done, "ssh/tcp", nil
}

// dialViaSocks5 returns a DialFunc that reaches targets through a SOCKS5 proxy.
// sockNet is "unix" for a Unix socket or "tcp" for a TCP address.
func dialViaSocks5(sockNet, sockAddr string) DialFunc {
	return func(ctx context.Context, network, target string) (net.Conn, error) {
		host, portStr, err := net.SplitHostPort(target)
		if err != nil {
			return nil, err
		}
		p, err := strconv.ParseUint(portStr, 10, 16)
		if err != nil {
			return nil, fmt.Errorf("invalid port %q: %w", portStr, err)
		}

		conn, err := net.Dial(sockNet, sockAddr)
		if err != nil {
			return nil, err
		}

		// SOCKS5 no-auth greeting
		if _, err := conn.Write([]byte{0x05, 0x01, 0x00}); err != nil {
			conn.Close()
			return nil, err
		}
		resp := make([]byte, 2)
		if _, err := io.ReadFull(conn, resp); err != nil || resp[0] != 0x05 || resp[1] != 0x00 {
			conn.Close()
			return nil, fmt.Errorf("socks5 handshake failed")
		}

		// CONNECT request with domain name
		req := []byte{0x05, 0x01, 0x00, 0x03, byte(len(host))}
		req = append(req, []byte(host)...)
		port := uint16(p)
		req = append(req, byte(port>>8), byte(port))
		if _, err := conn.Write(req); err != nil {
			conn.Close()
			return nil, err
		}

		// Response: 4-byte header + variable bound address
		hdr := make([]byte, 4)
		if _, err := io.ReadFull(conn, hdr); err != nil || hdr[1] != 0x00 {
			conn.Close()
			return nil, fmt.Errorf("socks5 CONNECT failed: status %d", hdr[1])
		}
		switch hdr[3] {
		case 0x01:
			io.ReadFull(conn, make([]byte, 6)) //nolint:errcheck
		case 0x03:
			l := make([]byte, 1)
			io.ReadFull(conn, l)                         //nolint:errcheck
			io.ReadFull(conn, make([]byte, int(l[0])+2)) //nolint:errcheck
		case 0x04:
			io.ReadFull(conn, make([]byte, 18)) //nolint:errcheck
		}
		return conn, nil
	}
}
