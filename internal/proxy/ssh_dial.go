package proxy

import (
	"context"
	"fmt"
	"io"
	"net"
	"os"
	"os/exec"
	"os/user"
	"path/filepath"
	"sync"
	"time"

	"golang.org/x/crypto/ssh"
	"golang.org/x/crypto/ssh/agent"

	"github.com/Justype/condatainer/internal/utils"
)

// rsaSHA2Signer wraps an ssh.AlgorithmSigner to advertise and sign with
// rsa-sha2-256. Go crypto/ssh v0.50.0 offers "ssh-rsa" (SHA-1) but signs with
// SHA-256, causing a server-side mismatch. This wrapper fixes it.
type rsaSHA2Signer struct{ ssh.AlgorithmSigner }

func (r *rsaSHA2Signer) PublicKey() ssh.PublicKey {
	return &rsaSHA2PublicKey{r.AlgorithmSigner.PublicKey()}
}

func (r *rsaSHA2Signer) Sign(rand io.Reader, data []byte) (*ssh.Signature, error) {
	return r.AlgorithmSigner.SignWithAlgorithm(rand, data, ssh.KeyAlgoRSASHA256)
}

type rsaSHA2PublicKey struct{ ssh.PublicKey }

func (k *rsaSHA2PublicKey) Type() string { return ssh.KeyAlgoRSASHA256 }

// findHostKey reads the first available host public key from /etc/ssh/.
// Certificates are preferred over plain keys (many clusters require cert-based
// hostbased auth). Within each group: ed25519 > ecdsa > rsa.
func findHostKey() (ssh.PublicKey, error) {
	for _, path := range []string{
		"/etc/ssh/ssh_host_ed25519_key-cert.pub",
		"/etc/ssh/ssh_host_ecdsa_key-cert.pub",
		"/etc/ssh/ssh_host_rsa_key-cert.pub",
		"/etc/ssh/ssh_host_ed25519_key.pub",
		"/etc/ssh/ssh_host_ecdsa_key.pub",
		"/etc/ssh/ssh_host_rsa_key.pub",
	} {
		data, err := os.ReadFile(path)
		if err != nil {
			continue
		}
		key, _, _, _, err := ssh.ParseAuthorizedKey(data)
		if err != nil {
			continue
		}
		return key, nil
	}
	return nil, fmt.Errorf("no readable host public key in /etc/ssh/")
}

// findSSHKeysign returns the path to the setuid ssh-keysign binary.
func findSSHKeysign() (string, error) {
	for _, p := range []string{
		"/usr/libexec/openssh/ssh-keysign",
		"/usr/lib/openssh/ssh-keysign",
		"/usr/libexec/ssh-keysign",
		"/usr/local/libexec/ssh-keysign",
	} {
		if _, err := os.Stat(p); err == nil {
			return p, nil
		}
	}
	if p, err := exec.LookPath("ssh-keysign"); err == nil {
		return p, nil
	}
	return "", fmt.Errorf("ssh-keysign not found")
}

// buildAuthMethods returns all non-interactive ssh.AuthMethod values available
// on this machine (hostbased, agent publickey, file publickey). Never prompts.
func buildAuthMethods() []ssh.AuthMethod {
	var methods []ssh.AuthMethod

	// Hostbased: requires ssh-keysign + a local host key in /etc/ssh/
	if keysign, err := findSSHKeysign(); err == nil {
		if hostKey, err := findHostKey(); err == nil {
			methods = append(methods, ssh.Hostbased(hostKey, keysign))
			utils.PrintDebug("proxy go-ssh auth: hostbased (%s via %s)", hostKey.Type(), keysign)
		} else {
			utils.PrintDebug("proxy go-ssh auth: hostbased skip — %v", err)
		}
	} else {
		utils.PrintDebug("proxy go-ssh auth: hostbased skip — %v", err)
	}

	// SSH agent: handles RSA sha-256 correctly on its own
	if sock := os.Getenv("SSH_AUTH_SOCK"); sock != "" {
		if conn, err := net.Dial("unix", sock); err == nil {
			methods = append(methods, ssh.PublicKeysCallback(agent.NewClient(conn).Signers))
			utils.PrintDebug("proxy go-ssh auth: agent (%s)", sock)
		} else {
			utils.PrintDebug("proxy go-ssh auth: agent skip — dial %s: %v", sock, err)
		}
	} else {
		utils.PrintDebug("proxy go-ssh auth: agent skip — SSH_AUTH_SOCK not set")
	}

	// Key files: ed25519/ecdsa (no algorithm negotiation issues)
	var signers []ssh.Signer
	home := sshHomeDir()
	for _, name := range []string{"id_ed25519", "id_ecdsa", "id_ecdsa_sk", "id_ed25519_sk"} {
		path := filepath.Join(home, ".ssh", name)
		data, err := os.ReadFile(path)
		if err != nil {
			continue
		}
		s, err := ssh.ParsePrivateKey(data)
		if err != nil {
			utils.PrintDebug("proxy go-ssh auth: %s skip — encrypted or unreadable", name)
			continue
		}
		signers = append(signers, s)
		utils.PrintDebug("proxy go-ssh auth: key file %s", path)
	}
	// RSA with sha2-256 wrapper to fix Go crypto/ssh algorithm mismatch
	if data, err := os.ReadFile(filepath.Join(home, ".ssh", "id_rsa")); err == nil {
		if s, err := ssh.ParsePrivateKey(data); err == nil {
			if alg, ok := s.(ssh.AlgorithmSigner); ok {
				signers = append(signers, &rsaSHA2Signer{alg})
				utils.PrintDebug("proxy go-ssh auth: key file id_rsa (rsa-sha2-256 wrapper)")
			}
		} else {
			utils.PrintDebug("proxy go-ssh auth: id_rsa skip — encrypted or unreadable")
		}
	}
	if len(signers) > 0 {
		methods = append(methods, ssh.PublicKeys(signers...))
	}

	return methods
}

// DialGoSSH dials sshDest using the Go SSH library with all available
// non-interactive auth methods. Returns a DialFunc (wrapping client.Dial),
// a stop closer, and a done channel closed when the connection drops.
func DialGoSSH(sshDest string) (DialFunc, func(), <-chan struct{}, error) {
	auths := buildAuthMethods()
	if len(auths) == 0 {
		utils.PrintDebug("proxy go-ssh: no non-interactive auth methods available — skipping")
		return nil, nil, nil, fmt.Errorf("no non-interactive SSH auth methods available")
	}
	utils.PrintDebug("proxy go-ssh: dialing %s with %d auth method(s)", sshDest, len(auths))

	addr := sshDest
	if _, _, err := net.SplitHostPort(sshDest); err != nil {
		addr = sshDest + ":22"
	}

	username := sshUsername()
	cfg := &ssh.ClientConfig{
		User:            username,
		Auth:            auths,
		HostKeyCallback: ssh.InsecureIgnoreHostKey(), //nolint:gosec
		HostKeyAlgorithms: []string{
			ssh.KeyAlgoED25519,
			ssh.KeyAlgoRSASHA256,
			ssh.KeyAlgoRSASHA512,
		},
		Timeout: 15 * time.Second,
	}

	client, err := ssh.Dial("tcp", addr, cfg)
	if err != nil {
		utils.PrintDebug("proxy go-ssh: dial %s failed: %v", sshDest, err)
		return nil, nil, nil, fmt.Errorf("go-ssh dial %s: %w", sshDest, err)
	}
	utils.PrintDebug("proxy go-ssh: connected to %s as %s", sshDest, username)

	done := make(chan struct{})
	go func() {
		client.Wait() //nolint:errcheck
		close(done)
	}()

	// Send keepalives so firewalls don't silently drop idle SSH connections.
	go func() {
		t := time.NewTicker(60 * time.Second)
		defer t.Stop()
		for {
			select {
			case <-t.C:
				client.SendRequest("keepalive@openssh.com", true, nil) //nolint:errcheck
			case <-done:
				return
			}
		}
	}()

	dial := func(ctx context.Context, network, tgt string) (net.Conn, error) {
		type result struct {
			conn net.Conn
			err  error
		}
		ch := make(chan result, 1)
		go func() {
			c, e := client.Dial(network, tgt)
			ch <- result{c, e}
		}()
		select {
		case r := <-ch:
			return r.conn, r.err
		case <-ctx.Done():
			return nil, ctx.Err()
		}
	}

	var once sync.Once
	stop := func() {
		once.Do(func() { client.Close() }) //nolint:errcheck
	}
	return dial, stop, done, nil
}

func sshUsername() string {
	if u, err := user.Current(); err == nil && u.Username != "" {
		return u.Username
	}
	return os.Getenv("USER")
}

func sshHomeDir() string {
	if u, err := user.Current(); err == nil {
		return u.HomeDir
	}
	return os.Getenv("HOME")
}
