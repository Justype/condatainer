package freeze

import (
	"context"
	"os"
	"os/exec"
	"path/filepath"
	"testing"
	"time"

	execpkg "github.com/Justype/condatainer/internal/runtime/exec"
)

// A ctx cancelled while work is still running must not leave the FUSE
// process behind: MountedRun puts the whole unshare/bash/FUSE tree in one
// process group (Setpgid) specifically so Cancel can kill all of it, not
// just the top-level bash Go tracks.
func TestMountedRunCancelKillsFuseProcess(t *testing.T) {
	squashfuse, err := findSquashfuse()
	if err != nil {
		t.Skip("no squashfuse available")
	}
	if _, err := exec.LookPath("unshare"); err != nil {
		t.Skip("unshare not available")
	}

	dir := t.TempDir()
	src := filepath.Join(dir, "src")
	if err := os.MkdirAll(src, 0o755); err != nil {
		t.Fatal(err)
	}
	if err := os.WriteFile(filepath.Join(src, "f.txt"), []byte("hello"), 0o644); err != nil {
		t.Fatal(err)
	}
	sqf := filepath.Join(dir, "t.sqf")
	mksquashfsBin, err := exec.LookPath("mksquashfs")
	if err != nil {
		t.Skip("mksquashfs not available")
	}
	if out, err := exec.Command(mksquashfsBin, src, sqf, "-no-progress").CombinedOutput(); err != nil {
		t.Fatalf("mksquashfs: %v\n%s", err, out)
	}

	mnt := filepath.Join(dir, "mnt")
	if err := os.MkdirAll(mnt, 0o755); err != nil {
		t.Fatal(err)
	}

	ctx, cancel := context.WithCancel(t.Context())
	done := make(chan error, 1)
	go func() {
		done <- MountedRun(ctx, squashfuse, []string{sqf}, mnt, "sleep 30", execpkg.IO{})
	}()

	// Give the mount time to appear before cancelling mid-"work".
	time.Sleep(1 * time.Second)
	cancel()

	select {
	case <-done:
	case <-time.After(10 * time.Second):
		t.Fatal("MountedRun did not return after cancellation")
	}

	// The kernel needs a moment to reap the killed process group.
	time.Sleep(1 * time.Second)

	out, _ := exec.Command("pgrep", "-f", sqf).Output()
	if len(out) != 0 {
		t.Errorf("FUSE process for %s survived cancellation:\n%s", sqf, out)
	}
}
