package freeze

import (
	"context"
	"fmt"
	"os"
	"os/exec"
	"path/filepath"
	"testing"
	"time"

	execpkg "github.com/Justype/condatainer/internal/runtime/exec"
)

// testHarnessPath is the compiled testdata/mountharness binary, built once
// in TestMain and used for the rest of this file: MountedRun (called by
// nearly every test in this package, directly or via Pack/Unfreeze) re-execs
// os.Executable() into the hidden _mount_sentinel command, and the `go test`
// binary can't play that role since it has no cmd/cobra dispatch at all.
var testHarnessPath string

func TestMain(m *testing.M) {
	if _, err := exec.LookPath("unshare"); err == nil {
		if bin, err := buildMountHarness(); err != nil {
			fmt.Fprintln(os.Stderr, "building test harness:", err)
		} else {
			testHarnessPath = bin
			SetExecutablePathForTest(bin)
		}
	}
	os.Exit(m.Run())
}

func buildMountHarness() (string, error) {
	dir, err := os.MkdirTemp("", "mountharness")
	if err != nil {
		return "", err
	}
	bin := filepath.Join(dir, "mountharness")
	out, err := exec.Command("go", "build", "-o", bin, "./testdata/mountharness").CombinedOutput()
	if err != nil {
		return "", fmt.Errorf("%w\n%s", err, out)
	}
	return bin, nil
}

func buildTestSqf(t *testing.T) (sqf, mnt string) {
	t.Helper()
	dir := t.TempDir()
	src := filepath.Join(dir, "src")
	if err := os.MkdirAll(src, 0o755); err != nil {
		t.Fatal(err)
	}
	if err := os.WriteFile(filepath.Join(src, "f.txt"), []byte("hello"), 0o644); err != nil {
		t.Fatal(err)
	}
	sqf = filepath.Join(dir, "t.sqf")
	mksquashfsBin, err := exec.LookPath("mksquashfs")
	if err != nil {
		t.Skip("mksquashfs not available")
	}
	if out, err := exec.Command(mksquashfsBin, src, sqf, "-no-progress").CombinedOutput(); err != nil {
		t.Fatalf("mksquashfs: %v\n%s", err, out)
	}
	mnt = filepath.Join(dir, "mnt")
	if err := os.MkdirAll(mnt, 0o755); err != nil {
		t.Fatal(err)
	}
	return sqf, mnt
}

// A ctx cancelled while work is still running must not leave the FUSE
// process behind: the sentinel joins bash/FUSE into its own process group
// specifically so Cancel can kill all of it, not just the sentinel itself.
func TestMountedRunCancelKillsFuseProcess(t *testing.T) {
	squashfuse, err := findSquashfuse()
	if err != nil {
		t.Skip("no squashfuse available")
	}
	if testHarnessPath == "" {
		t.Skip("unshare not available (test harness was not built)")
	}

	sqf, mnt := buildTestSqf(t)

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

// If the process that called MountedRun dies abruptly (SIGKILL, standing in
// for a crash or an OOM-kill) with no chance to run cmd.Cancel at all, the
// FUSE process must still not be left behind. This is exactly the gap
// Pdeathsig can't close directly on the namespaced process (entering the
// user namespace clears it) -- RunSentinel exists to close it one process
// out instead. See sentinel.go's doc comment.
func TestMountedRunSurvivesAbruptDeath(t *testing.T) {
	if _, err := findSquashfuse(); err != nil {
		t.Skip("no squashfuse available")
	}
	if testHarnessPath == "" {
		t.Skip("unshare not available (test harness was not built)")
	}

	sqf, mnt := buildTestSqf(t)

	// "drive" makes the harness call MountedRun for real, as its own OS
	// process -- standing in for condatainer -- so it can be killed from
	// outside without taking this test binary down with it.
	driver := exec.Command(testHarnessPath, "drive", "squashfuse", sqf, mnt, "sleep 30")
	if err := driver.Start(); err != nil {
		t.Fatalf("starting driver: %v", err)
	}
	driverPid := driver.Process.Pid

	// Wait for the FUSE process to actually come up (the mount lives in a
	// private namespace, invisible to us, so watch for the process instead).
	// Scoped to this test's own sqf path, never a bare tool name.
	deadline := time.Now().Add(10 * time.Second)
	for time.Now().Before(deadline) {
		out, _ := exec.Command("pgrep", "-f", "squashfuse.*"+sqf).Output()
		if len(out) != 0 {
			break
		}
		time.Sleep(100 * time.Millisecond)
	}
	out, _ := exec.Command("pgrep", "-f", "squashfuse.*"+sqf).Output()
	if len(out) == 0 {
		driver.Process.Kill() //nolint:errcheck
		t.Fatal("FUSE process never appeared")
	}

	// Kill only the exact driver pid -- nothing else.
	if err := driver.Process.Kill(); err != nil {
		t.Fatalf("killing driver: %v", err)
	}
	driver.Wait() //nolint:errcheck

	// Give the kernel a moment to deliver Pdeathsig to the sentinel and let
	// it run the group kill.
	time.Sleep(2 * time.Second)

	survivor, _ := exec.Command("pgrep", "-f", "squashfuse.*"+sqf).Output()
	if len(survivor) != 0 {
		exec.Command("pkill", "-KILL", "-f", "squashfuse.*"+sqf).Run() //nolint:errcheck
		t.Errorf("FUSE process for %s survived the driver's abrupt death (driver pid %d):\n%s", sqf, driverPid, survivor)
	}
}
