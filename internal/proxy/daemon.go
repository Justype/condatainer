package proxy

import (
	"os"
	"os/exec"
	"strconv"
)

// RunDaemon starts the SSH SOCKS5 tunnel and blocks until it exits.
// localOnly=false: shared mode (0.0.0.0, NFS PID file, host=hostname).
// localOnly=true:  per-job mode (127.0.0.1, local PID file, host=127.0.0.1).
func RunDaemon(sshDest string, port int, localOnly bool) error {
	bind := "0.0.0.0"
	if localOnly {
		bind = "127.0.0.1"
	}

	ssh := exec.Command("ssh",
		"-N",
		"-D", bind+":"+strconv.Itoa(port),
		sshDest,
		"-o", "BatchMode=yes",
		"-o", "StrictHostKeyChecking=no",
		"-o", "ServerAliveInterval=60",
		"-o", "ServerAliveCountMax=3",
	)

	if err := ssh.Start(); err != nil {
		return err
	}

	var pidHost, pidPath string
	if localOnly {
		pidHost = "127.0.0.1"
		pidPath = LocalPidFilePath()
	} else {
		pidHost, _ = os.Hostname() // listener hostname, not the SSH destination
		pidPath = PidFilePath()
	}

	if err := WritePidFileAt(pidPath, ProxyState{Host: pidHost, Via: sshDest, Port: port, PID: ssh.Process.Pid}); err != nil {
		ssh.Process.Kill() //nolint:errcheck
		return err
	}

	// Block until SSH exits — triggered by SIGTERM (proxy stop) or network death.
	ssh.Wait()         //nolint:errcheck
	os.Remove(pidPath) //nolint:errcheck
	return nil
}
