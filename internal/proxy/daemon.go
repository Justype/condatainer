package proxy

import (
	"os/exec"
	"strconv"
)

// RunDaemon starts the SSH SOCKS5 tunnel, writes the PID file, and blocks until
// the tunnel exits (either via SIGTERM from "proxy stop" or a network failure).
// Called by the daemonized child process spawned by "condatainer proxy start".
// Removes the PID file on exit so status checks see a clean state.
func RunDaemon(host string, port int) error {
	ssh := exec.Command("ssh",
		"-N",
		"-D", "0.0.0.0:"+strconv.Itoa(port),
		host,
		"-o", "BatchMode=yes",
		"-o", "StrictHostKeyChecking=no",
		"-o", "ServerAliveInterval=60",
		"-o", "ServerAliveCountMax=3",
	)

	if err := ssh.Start(); err != nil {
		return err
	}

	if err := WritePidFile(host, port, ssh.Process.Pid); err != nil {
		ssh.Process.Kill() //nolint:errcheck
		return err
	}

	// Block until SSH exits — triggered by SIGTERM (proxy stop) or network death.
	ssh.Wait() //nolint:errcheck
	RemovePidFile()
	return nil
}
