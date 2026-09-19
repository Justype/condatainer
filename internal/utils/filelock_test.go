package utils

import (
	"bufio"
	"errors"
	"os"
	"os/exec"
	"path/filepath"
	"testing"

	"golang.org/x/sys/unix"
)

func lockTarget(t *testing.T) string {
	t.Helper()
	path := filepath.Join(t.TempDir(), "image")
	if err := os.WriteFile(path, make([]byte, 4096), 0o644); err != nil {
		t.Fatal(err)
	}
	return path
}

// Shared holders coexist; an exclusive holder excludes everyone, including
// another acquisition in the same process; Close releases.
func TestAcquireFileLockConflicts(t *testing.T) {
	path := lockTarget(t)

	r1, err := AcquireFileLock(path, false)
	if err != nil {
		t.Fatal(err)
	}
	r2, err := AcquireFileLock(path, false)
	if err != nil {
		t.Fatalf("second shared lock: %v", err)
	}
	if _, err := AcquireFileLock(path, true); !errors.Is(err, ErrLockConflict) {
		t.Fatalf("exclusive over shared: err = %v, want ErrLockConflict", err)
	}
	r1.Close()
	r2.Close()

	w, err := AcquireFileLock(path, true)
	if err != nil {
		t.Fatalf("exclusive after release: %v", err)
	}
	for _, write := range []bool{false, true} {
		if _, err := AcquireFileLock(path, write); !errors.Is(err, ErrLockConflict) {
			t.Errorf("write=%v over exclusive: err = %v, want ErrLockConflict", write, err)
		}
	}
	w.Close()
	if _, err := AcquireFileLock(path, true); err != nil {
		t.Errorf("after Close: %v", err)
	}
}

// A plain fcntl lock in another process — what apptainer holds on a mounted
// ext3 image — is seen.
func TestAcquireFileLockSeesAnotherProcessesPosixLock(t *testing.T) {
	if os.Getenv("CNT_TEST_HOLD_LOCK") != "" {
		f, err := os.OpenFile(os.Getenv("CNT_TEST_HOLD_LOCK"), os.O_RDWR, 0)
		if err != nil {
			os.Exit(2)
		}
		lock := unix.Flock_t{Type: unix.F_WRLCK}
		if unix.FcntlFlock(f.Fd(), unix.F_SETLK, &lock) != nil {
			os.Exit(3)
		}
		os.Stdout.WriteString("held\n")
		bufio.NewReader(os.Stdin).ReadString('\n') // parent closes stdin to release
		os.Exit(0)
	}

	path := lockTarget(t)
	child := exec.Command(os.Args[0], "-test.run=^TestAcquireFileLockSeesAnotherProcessesPosixLock$")
	child.Env = append(os.Environ(), "CNT_TEST_HOLD_LOCK="+path)
	stdin, _ := child.StdinPipe()
	stdout, _ := child.StdoutPipe()
	if err := child.Start(); err != nil {
		t.Fatal(err)
	}
	defer child.Wait()
	defer stdin.Close()
	if line, _ := bufio.NewReader(stdout).ReadString('\n'); line != "held\n" {
		t.Fatalf("child did not take its lock: %q", line)
	}

	for _, write := range []bool{false, true} {
		if _, err := AcquireFileLock(path, write); !errors.Is(err, ErrLockConflict) {
			t.Errorf("write=%v: err = %v, want ErrLockConflict", write, err)
		}
	}
}

// A file that cannot be opened is reported as such, not as a lock conflict.
func TestAcquireFileLockOpenFailureIsNotAConflict(t *testing.T) {
	_, err := AcquireFileLock(filepath.Join(t.TempDir(), "missing"), true)
	if !errors.Is(err, os.ErrNotExist) || errors.Is(err, ErrLockConflict) {
		t.Errorf("err = %v, want not-exist and not a conflict", err)
	}
}
