package proxy

import (
	"os"
	"path/filepath"
	"strings"
	"testing"
	"time"
)

// A failed login is reported as soon as ssh exits, with what ssh printed,
// instead of after the listener timeout.
func TestStartSocksSSHFailsFastWhenSSHExits(t *testing.T) {
	bin := t.TempDir()
	script := "#!/bin/sh\necho 'Permission denied (publickey)' >&2\nexit 255\n"
	if err := os.WriteFile(filepath.Join(bin, "ssh"), []byte(script), 0o755); err != nil {
		t.Fatal(err)
	}
	t.Setenv("PATH", bin+":"+os.Getenv("PATH"))

	sock := filepath.Join(t.TempDir(), "s.sock")
	start := time.Now()
	_, _, _, err := startSocksSSH(sock, "unix", sock, []string{"-N", "node"}, func() {})
	if err == nil || !strings.Contains(err.Error(), "Permission denied") {
		t.Fatalf("err = %v, want ssh's message", err)
	}
	if elapsed := time.Since(start); elapsed > tunnelStartTimeout/2 {
		t.Fatalf("took %s, should not wait out the timeout", elapsed)
	}
}
