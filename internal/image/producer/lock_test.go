package producer

import (
	"errors"
	"os"
	"path/filepath"
	"testing"
)

// asJob makes the process look like the given scheduler job for one test.
func asJob(t *testing.T, id string) {
	t.Helper()
	previous := currentJobID
	currentJobID = func() string { return id }
	t.Cleanup(func() { currentJobID = previous })
}

func TestAcquireLocalSerializesTarget(t *testing.T) {
	target := filepath.Join(t.TempDir(), "demo.sqf")
	first, err := AcquireLocal(target)
	if err != nil {
		t.Fatal(err)
	}
	defer first.Release() //nolint:errcheck

	_, err = AcquireLocal(target)
	var producing *ProducingError
	if !errors.As(err, &producing) {
		t.Fatalf("second AcquireLocal error = %v, want *ProducingError", err)
	}
	if producing.Info.PID != os.Getpid() {
		t.Fatalf("ProducingError names PID %d, want %d", producing.Info.PID, os.Getpid())
	}
	if err := first.Release(); err != nil {
		t.Fatal(err)
	}
	second, err := AcquireLocal(target)
	if err != nil {
		t.Fatalf("AcquireLocal after release: %v", err)
	}
	second.Release() //nolint:errcheck
}

func TestAcquireLocalClearsStaleLockAndPartial(t *testing.T) {
	target := filepath.Join(t.TempDir(), "demo.sqf")
	stale := Info{}
	if err := Acquire(Path(target), stale); err != nil {
		t.Fatal(err)
	}
	partial := PreparedPath(target, stale)
	if err := os.WriteFile(partial, []byte("partial"), 0o600); err != nil {
		t.Fatal(err)
	}

	guard, err := AcquireLocal(target)
	if err != nil {
		t.Fatal(err)
	}
	defer guard.Release() //nolint:errcheck
	if _, err := os.Stat(partial); !os.IsNotExist(err) {
		t.Fatalf("stale partial remains: %v", err)
	}
}

func TestTagDistinguishesLocalAndSchedulerOwners(t *testing.T) {
	local := Tag(Info{Runner: "local", Node: "node-a", PID: 42})
	job := Tag(Info{Runner: "slurm", JobID: "42"})
	if local == job || local != "local-node-a-42" || job != "slurm-42" {
		t.Fatalf("tags = %q and %q", local, job)
	}
}

// A submitted job finds a lock its own submitter created and must adopt it, not
// report itself as another producer. A job ID names one job, so a lock carrying
// ours is ours.
func TestAcquireLocalAdoptsItsOwnJobLock(t *testing.T) {
	asJob(t, "12345")
	target := filepath.Join(t.TempDir(), "demo.sqf")
	submitted := Info{Runner: "slurm", JobID: "12345", Node: "login01", CreatedAt: "2026-01-01T00:00:00Z"}
	if err := Acquire(Path(target), submitted); err != nil {
		t.Fatal(err)
	}

	guard, err := AcquireLocal(target)
	if err != nil {
		t.Fatalf("a job did not adopt its own lock: %v", err)
	}
	if guard.Info().JobID != "12345" {
		t.Errorf("adopted info = %#v, want the submitter's record", guard.Info())
	}
	// Releasing hands the lock back, which is what lets the next restore see the
	// artifact as no longer in flight.
	if err := guard.Release(); err != nil {
		t.Fatal(err)
	}
	if _, err := os.Stat(Path(target)); !os.IsNotExist(err) {
		t.Errorf("the adopted lock survived release: %v", err)
	}
}

// Another job's live lock is still a conflict: only a matching job ID is ours.
func TestAcquireLocalRefusesAnotherJobsLock(t *testing.T) {
	asJob(t, "12345")
	target := filepath.Join(t.TempDir(), "demo.sqf")
	other := Info{Runner: "slurm", JobID: "99999", Node: "login01", CreatedAt: "2026-01-01T00:00:00Z"}
	if err := Acquire(Path(target), other); err != nil {
		t.Fatal(err)
	}
	if _, err := AcquireLocal(target); err == nil {
		t.Fatal("another job's lock was adopted")
	}
}
