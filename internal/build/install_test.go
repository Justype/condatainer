package build

import (
	"os"
	"path/filepath"
	"strings"
	"testing"
)

// The prepared path must sit beside the target, or the final rename crosses a
// filesystem boundary and stops being atomic.
func TestPreparedPathIsBesideTarget(t *testing.T) {
	target := "/images/samtools--1.21.sqf"
	got := preparedPathFor(target, BuildLockInfo{Runner: "local", Node: "cn001", PID: 4242})

	if filepath.Dir(got) != filepath.Dir(target) {
		t.Errorf("prepared %q is not in the target's directory %q", got, filepath.Dir(target))
	}
	if got == target {
		t.Error("prepared path equals the target; a build would write over the installed image")
	}
	if !strings.HasSuffix(got, preparedSuffix) {
		t.Errorf("prepared %q does not end in %q", got, preparedSuffix)
	}
}

// The owner has to be recoverable from the filename, since that is how an orphan
// left by a killed process is attributed and cleaned up.
func TestPreparedPathNamesItsOwner(t *testing.T) {
	target := "/images/x--1.sqf"

	local := preparedPathFor(target, BuildLockInfo{Runner: "local", Node: "cn001", PID: 4242})
	if !strings.Contains(local, "cn001") || !strings.Contains(local, "4242") {
		t.Errorf("local prepared path %q does not name node and pid", local)
	}

	job := preparedPathFor(target, BuildLockInfo{Runner: "slurm", JobID: "98765", Node: "cn002", PID: 7})
	if !strings.Contains(job, "slurm") || !strings.Contains(job, "98765") {
		t.Errorf("scheduler prepared path %q does not name runner and job", job)
	}
	if local == job {
		t.Error("two different owners derived the same prepared path")
	}
}

// Recomputable from the lock alone: stale-lock cleanup has no record of the path
// the dead owner used, so it must derive the identical one.
func TestPreparedPathIsDeterministic(t *testing.T) {
	info := BuildLockInfo{Runner: "slurm", JobID: "42", Node: "cn003", PID: 9}
	first := preparedPathFor("/images/y--2.sqf", info)
	second := preparedPathFor("/images/y--2.sqf", info)
	if first != second {
		t.Errorf("derived two different paths for one owner: %q vs %q", first, second)
	}
}

// A job ID from a scheduler is not guaranteed to be filename-safe.
func TestPreparedPathSanitizesOwner(t *testing.T) {
	got := preparedPathFor("/images/z--1.sqf", BuildLockInfo{Runner: "lsf", JobID: "job/1 [2]"})
	base := filepath.Base(got)
	if strings.ContainsAny(base, "/ []") {
		t.Errorf("prepared basename %q contains characters that are unsafe in a filename", base)
	}
	if filepath.Dir(got) != "/images" {
		t.Errorf("a slash in the job ID escaped the target directory: %q", got)
	}
}

// The installed image must survive a failed install, which is the whole reason
// the old target is not removed first.
func TestAtomicInstallReplacesWithoutRemovingFirst(t *testing.T) {
	dir := t.TempDir()
	target := filepath.Join(dir, "img.sqf")
	prepared := filepath.Join(dir, "img.sqf.local-cn-1.part")

	if err := os.WriteFile(target, []byte("old"), 0o644); err != nil {
		t.Fatal(err)
	}
	if err := os.WriteFile(prepared, []byte("new"), 0o644); err != nil {
		t.Fatal(err)
	}

	if err := atomicInstall(prepared, target); err != nil {
		t.Fatalf("atomicInstall: %v", err)
	}

	data, err := os.ReadFile(target)
	if err != nil {
		t.Fatal(err)
	}
	if string(data) != "new" {
		t.Errorf("target holds %q, want the prepared content", data)
	}
	if _, err := os.Stat(prepared); !os.IsNotExist(err) {
		t.Error("prepared path still exists after install")
	}
}

// A failed install must leave the previously installed image untouched.
func TestAtomicInstallFailureKeepsInstalledImage(t *testing.T) {
	dir := t.TempDir()
	target := filepath.Join(dir, "img.sqf")
	if err := os.WriteFile(target, []byte("installed"), 0o644); err != nil {
		t.Fatal(err)
	}

	// A prepared path that was never produced: rename fails.
	missing := filepath.Join(dir, "img.sqf.local-cn-2.part")
	if err := atomicInstall(missing, target); err == nil {
		t.Fatal("installing a missing prepared path succeeded")
	}

	data, err := os.ReadFile(target)
	if err != nil {
		t.Fatalf("installed image was destroyed by a failed install: %v", err)
	}
	if string(data) != "installed" {
		t.Errorf("installed image now holds %q", data)
	}
}

// Cleanup on failure removes this build's own output and nothing else. The
// previous behaviour deleted the target itself outside update mode, which with
// prepared paths would delete a working image on an unrelated failure.
func TestCleanupRemovesOnlyPreparedOutput(t *testing.T) {
	dir := t.TempDir()
	target := filepath.Join(dir, "img.sqf")
	prepared := filepath.Join(dir, "img.sqf.local-cn-3.part")
	for _, p := range []string{target, prepared} {
		if err := os.WriteFile(p, []byte("x"), 0o644); err != nil {
			t.Fatal(err)
		}
	}

	b := &BuildObject{tgt: Target{Path: target, Prepared: prepared}}
	if err := b.Cleanup(true); err != nil {
		t.Fatalf("Cleanup: %v", err)
	}

	if _, err := os.Stat(prepared); !os.IsNotExist(err) {
		t.Error("failed build left its partial output behind")
	}
	if _, err := os.Stat(target); err != nil {
		t.Errorf("failed build removed the installed image: %v", err)
	}
}

// A killed process leaves both a lock and a partial output; clearing the lock
// without the output would strand the bytes forever.
func TestStaleLockCleanupRemovesItsOrphan(t *testing.T) {
	dir := t.TempDir()
	target := filepath.Join(dir, "img.sqf")
	info := BuildLockInfo{Runner: "local", Node: "cn009", PID: 31337}
	orphan := preparedPathFor(target, info)
	if err := os.WriteFile(orphan, []byte("half an image"), 0o644); err != nil {
		t.Fatal(err)
	}

	b := &BuildObject{tgt: targetFor(target)}
	b.removeOrphanedOutput(t.Context(), info)

	if _, err := os.Stat(orphan); !os.IsNotExist(err) {
		t.Error("orphaned output survived stale-lock cleanup")
	}
}

// Every constructor has to carry update through to the object, or the backend
// skips the rebuild its caller asked for. FromExternalSource did not take one at
// all, which made `build -f script.sh --update` a silent no-op.
func TestExternalSourceCarriesUpdate(t *testing.T) {
	dir := t.TempDir()
	script := filepath.Join(dir, "demo.sh")
	if err := os.WriteFile(script, []byte("#!/usr/bin/env bash\n#DESC:demo\ntrue\n"), 0o755); err != nil {
		t.Fatal(err)
	}

	for _, want := range []bool{true, false} {
		obj, err := FromExternalSource(t.Context(), filepath.Join(dir, "demo"), script, false, dir, want)
		if err != nil {
			t.Fatalf("FromExternalSource(update=%v): %v", want, err)
		}
		if obj.update != want {
			t.Errorf("update = %v, want %v", obj.update, want)
		}
		if obj.Update() != want {
			t.Errorf("Update() = %v, want %v", obj.Update(), want)
		}
	}
}
