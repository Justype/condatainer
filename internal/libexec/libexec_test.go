package libexec

import (
	"os"
	"path/filepath"
	"testing"

	"github.com/Justype/condatainer/internal/config"
	"github.com/Justype/condatainer/internal/utils"
)

// withScratchTier points only the (uncached) scratch tier at a fresh temp
// directory and forces config's lazily-cached search paths to recompute.
// CNT_ROOT/CNT_EXTRA_ROOT are left empty throughout this file rather than
// varied per test: internal/config memoizes GetRootDir/GetExtraRootDir behind
// a package-private sync.Once this package cannot reset, so a test that
// changed them would silently keep reading whichever value a prior test in
// this binary saw first. internal/config's own tests already cover
// tier-ordering; this package only needs one uncached tier to exercise
// Dir/BinDir's own marker-checking logic.
func withScratchTier(t *testing.T) string {
	t.Helper()
	scratch := filepath.Join(t.TempDir(), "condatainer")
	t.Setenv("SCRATCH", filepath.Dir(scratch))
	t.Setenv("XDG_DATA_HOME", "")
	t.Setenv("CNT_EXTRA_ROOT", "")
	t.Setenv("CNT_ROOT", "")
	config.InitDataPaths()
	return scratch
}

// provisionedStub drops a fake bin/apptainer and a lock sentinel under
// tier/libexec, just enough to satisfy Dir's marker check and let
// AcquireUse/Update's own locking be exercised without a real bootstrap.
func provisionedStub(t *testing.T, tier string) {
	t.Helper()
	bin := filepath.Join(tier, "libexec", "bin")
	if err := utils.MkdirAllShared(bin); err != nil {
		t.Fatalf("failed to create stub bin dir: %v", err)
	}
	if err := os.WriteFile(filepath.Join(bin, "apptainer"), []byte("#!/bin/sh\n"), 0755); err != nil {
		t.Fatalf("failed to write stub apptainer: %v", err)
	}
	if err := os.WriteFile(filepath.Join(tier, "libexec", lockFileName), nil, 0644); err != nil {
		t.Fatalf("failed to write stub lock sentinel: %v", err)
	}
}

func TestDirFindsAProvisionedTier(t *testing.T) {
	scratch := withScratchTier(t)
	provisionedStub(t, scratch)

	dir, ok := Dir()
	if !ok {
		t.Fatal("Dir() = not found, want the scratch tier")
	}
	if want := filepath.Join(scratch, "libexec"); dir != want {
		t.Errorf("Dir() = %s, want %s", dir, want)
	}
}

func TestDirIgnoresUnprovisionedTiers(t *testing.T) {
	withScratchTier(t)
	if dir, ok := Dir(); ok {
		t.Fatalf("Dir() = %s, want not found (no tier provisioned)", dir)
	}
}

func TestDirIgnoresAnEmptyLibexecDir(t *testing.T) {
	scratch := withScratchTier(t)
	// libexec/ exists but was never actually provisioned (no bin/apptainer) —
	// e.g. a half-finished bootstrap or a stale .new left by a crash.
	if err := utils.MkdirAllShared(filepath.Join(scratch, "libexec")); err != nil {
		t.Fatalf("failed to create empty libexec dir: %v", err)
	}
	if dir, ok := Dir(); ok {
		t.Fatalf("Dir() = %s, want not found (empty libexec dir)", dir)
	}
}

func TestApptainerPathAndMksquashfsPath(t *testing.T) {
	scratch := withScratchTier(t)
	provisionedStub(t, scratch)

	path, ok := ApptainerPath()
	if !ok {
		t.Fatal("ApptainerPath() = not found")
	}
	if want := filepath.Join(scratch, "libexec", "bin", "apptainer"); path != want {
		t.Errorf("ApptainerPath() = %s, want %s", path, want)
	}

	// mksquashfs was never written by the stub; MksquashfsPath only names the
	// path within an already-provisioned bin/, it does not check this one file.
	msPath, ok := MksquashfsPath()
	if !ok {
		t.Fatal("MksquashfsPath() = not found")
	}
	if want := filepath.Join(scratch, "libexec", "bin", "mksquashfs"); msPath != want {
		t.Errorf("MksquashfsPath() = %s, want %s", msPath, want)
	}
}

func TestAcquireUseWithNothingProvisioned(t *testing.T) {
	withScratchTier(t)
	lock, err := AcquireUse()
	if lock != nil || err != nil {
		t.Errorf("AcquireUse() = %v, %v; want nil, nil (nothing to protect)", lock, err)
	}
}

func TestAcquireUseSucceedsWhenFree(t *testing.T) {
	scratch := withScratchTier(t)
	provisionedStub(t, scratch)

	lock, err := AcquireUse()
	if err != nil {
		t.Fatalf("AcquireUse() failed: %v", err)
	}
	if lock == nil {
		t.Fatal("AcquireUse() = nil lock, want a held one")
	}
	defer lock.Close()

	// Another concurrent reader is fine — shared locks do not conflict with
	// each other, only with Update's exclusive one.
	second, err := AcquireUse()
	if err != nil {
		t.Fatalf("a second concurrent reader should not conflict: %v", err)
	}
	second.Close()
}

func TestAcquireUseFailsWhileUpdateHoldsTheLock(t *testing.T) {
	scratch := withScratchTier(t)
	provisionedStub(t, scratch)

	writer, err := utils.AcquireFlock(filepath.Join(scratch, "libexec", lockFileName), true)
	if err != nil {
		t.Fatalf("failed to simulate Update's exclusive lock: %v", err)
	}
	defer writer.Close()

	if lock, err := AcquireUse(); err == nil {
		if lock != nil {
			lock.Close()
		}
		t.Error("AcquireUse() succeeded while the exclusive lock was held, want a refusal")
	}
}
