package libexec

import (
	"context"
	"os"
	"path/filepath"
	"runtime"
	"testing"

	"github.com/Justype/condatainer/internal/utils"
)

func TestParseVersion(t *testing.T) {
	cases := map[string]string{
		"apptainer version 1.5.3":        "1.5.3",
		"1.4.0":                          "1.4.0",
		"singularity-ce version 4.1.1\n": "4.1.1",
		"":                               "",
		"no version here":                "",
	}
	for input, want := range cases {
		if got := parseVersion(input); got != want {
			t.Errorf("parseVersion(%q) = %q, want %q", input, got, want)
		}
	}
}

func TestMeetsFloor(t *testing.T) {
	cases := map[string]bool{
		"1.5.3": true,
		"1.4.0": true,
		"1.3.9": false,
		"1.0.0": false,
		"2.0.0": true,
		"0.9.0": false,
		"1":     false, // no minor component to compare
	}
	for version, want := range cases {
		if got := meetsFloor(version, apptainerZstdFloorMajor, apptainerZstdFloorMinor); got != want {
			t.Errorf("meetsFloor(%q, %d, %d) = %v, want %v", version, apptainerZstdFloorMajor, apptainerZstdFloorMinor, got, want)
		}
	}
}

// mksquashfs and unsquashfs only recognize their own single-dash "-version",
// not the GNU-style "--version" every other provisioned tool takes.
func TestVersionFlag(t *testing.T) {
	cases := map[string]string{
		"mksquashfs": "-version",
		"unsquashfs": "-version",
		"apptainer":  "--version",
		"squashfuse": "--version",
		"micromamba": "--version",
	}
	for name, want := range cases {
		if got := versionFlag(name); got != want {
			t.Errorf("versionFlag(%q) = %q, want %q", name, got, want)
		}
	}
}

func TestMicromambaAssetName(t *testing.T) {
	if runtime.GOOS != "linux" {
		t.Skip("the self-provisioned toolchain is linux-only")
	}
	asset, err := micromambaAssetName()
	switch runtime.GOARCH {
	case "amd64":
		if err != nil || asset != "micromamba-linux-64" {
			t.Errorf("micromambaAssetName() = %q, %v; want micromamba-linux-64, nil", asset, err)
		}
	case "arm64":
		if err != nil || asset != "micromamba-linux-aarch64" {
			t.Errorf("micromambaAssetName() = %q, %v; want micromamba-linux-aarch64, nil", asset, err)
		}
	default:
		if err == nil {
			t.Errorf("micromambaAssetName() = %q, nil; want an error for unsupported arch %s", asset, runtime.GOARCH)
		}
	}
}

func TestUpdateRefusesWhileInUse(t *testing.T) {
	scratch := withScratchTier(t)
	provisionedStub(t, scratch)

	reader, err := utils.AcquireFlock(filepath.Join(scratch, "libexec", lockFileName), false)
	if err != nil {
		t.Fatalf("failed to simulate an active reader: %v", err)
	}
	defer reader.Close()

	// Must fail before ever reaching the network (no bootstrap download, no
	// micromamba invocation) — the lock check is the very first thing Update
	// does once it finds a live generation.
	if err := Update(context.Background()); err == nil {
		t.Fatal("Update() succeeded while a reader held the lock, want a refusal")
	}

	// The live generation must be untouched — Update bailed out before
	// creating a staging directory or touching anything.
	if _, ok := Dir(); !ok {
		t.Error("Update() left no provisioned toolchain after refusing")
	}
	if info, err := os.Lstat(filepath.Join(scratch, stagingName)); err == nil && info.Mode()&os.ModeSymlink == 0 {
		t.Error("Update() created a staging directory despite refusing up front")
	}
}

func TestVersionsRefusesWithNothingProvisioned(t *testing.T) {
	withScratchTier(t)
	if _, err := Versions(context.Background()); err != ErrNotProvisioned {
		t.Errorf("Versions() error = %v, want ErrNotProvisioned", err)
	}
}

func TestVersionsParsesEachBinary(t *testing.T) {
	scratch := withScratchTier(t)
	provisionedStub(t, scratch)
	bin := filepath.Join(scratch, "libexec", "bin")

	// apptainer already exists (provisionedStub's marker); give it a real
	// version-shaped banner instead of the empty stub, and add the other
	// three, one printing nothing --version-shaped to prove that case
	// reports "unknown" rather than failing the rest.
	writeFakeVersionBin(t, bin, "apptainer", "apptainer version 1.5.2\n")
	writeFakeVersionBin(t, bin, "mksquashfs", "mksquashfs version 4.6.1 (2023-08-19)\n")
	writeFakeVersionBin(t, bin, "squashfuse", "squashfuse\n") // no version in output
	writeFakeVersionBin(t, bin, "micromamba", "2.0.5\n")

	versions, err := Versions(context.Background())
	if err != nil {
		t.Fatalf("Versions(): %v", err)
	}
	want := map[string]string{
		"apptainer":  "1.5.2",
		"mksquashfs": "4.6.1",
		"squashfuse": "unknown",
		"micromamba": "2.0.5",
	}
	got := map[string]string{}
	for _, v := range versions {
		got[v.Name] = v.Version
	}
	for name, wantVersion := range want {
		if got[name] != wantVersion {
			t.Errorf("%s version = %q, want %q", name, got[name], wantVersion)
		}
	}
}

// writeStubApptainer drops a fake bin/apptainer directly under dir, enough to
// satisfy isProvisioned/Dir's own marker check.
func writeStubApptainer(t *testing.T, dir string) {
	t.Helper()
	bin := filepath.Join(dir, "bin")
	if err := utils.MkdirAllShared(bin); err != nil {
		t.Fatalf("failed to create stub bin dir: %v", err)
	}
	if err := os.WriteFile(filepath.Join(bin, "apptainer"), []byte("#!/bin/sh\n"), 0755); err != nil {
		t.Fatalf("failed to write stub apptainer: %v", err)
	}
}

// provisionedVerifiableDir writes fake apptainer/mksquashfs/unsquashfs/
// squashfuse binaries under dir/bin that pass verifyToolchain outright: each
// prints its own name for its version flag, and the two floor checks clear
// their thresholds.
func provisionedVerifiableDir(t *testing.T, dir string) {
	t.Helper()
	bin := filepath.Join(dir, "bin")
	if err := utils.MkdirAllShared(bin); err != nil {
		t.Fatalf("failed to create stub bin dir: %v", err)
	}
	writeFakeVersionBin(t, bin, "apptainer", "apptainer version 1.5.0\n")
	writeFakeVersionBin(t, bin, "mksquashfs", "mksquashfs version 4.5.0 (2023-08-19)\n")
	writeFakeVersionBin(t, bin, "unsquashfs", "unsquashfs version 4.5.0\n")
	writeFakeVersionBin(t, bin, "squashfuse", "squashfuse version 0.5.0\n")
}

func TestEnsureCompatSymlinkCreatesAndRepairs(t *testing.T) {
	dir := t.TempDir()
	target := filepath.Join(dir, "libexec")
	staging := filepath.Join(dir, stagingName)
	writeStubApptainer(t, target)

	if err := ensureCompatSymlink(target, staging); err != nil {
		t.Fatalf("ensureCompatSymlink() (create): %v", err)
	}
	if link, err := os.Readlink(staging); err != nil || link != target {
		t.Fatalf("staging -> %q, %v; want %q", link, err, target)
	}

	// Idempotent when already correct.
	if err := ensureCompatSymlink(target, staging); err != nil {
		t.Fatalf("ensureCompatSymlink() (idempotent): %v", err)
	}

	// Repairs a symlink pointing somewhere stale.
	os.Remove(staging)
	if err := os.Symlink(filepath.Join(dir, "somewhere-else"), staging); err != nil {
		t.Fatalf("failed to seed a stale symlink: %v", err)
	}
	if err := ensureCompatSymlink(target, staging); err != nil {
		t.Fatalf("ensureCompatSymlink() (repair): %v", err)
	}
	if link, err := os.Readlink(staging); err != nil || link != target {
		t.Fatalf("staging -> %q, %v; want %q after repair", link, err, target)
	}
}

func TestEnsureFusermount3InBinLinksFromSbin(t *testing.T) {
	prefix := t.TempDir()
	sbin := filepath.Join(prefix, "sbin")
	if err := utils.MkdirAllShared(sbin); err != nil {
		t.Fatalf("failed to create sbin: %v", err)
	}
	if err := os.WriteFile(filepath.Join(sbin, "fusermount3"), []byte("#!/bin/sh\n"), 0755); err != nil {
		t.Fatalf("failed to write stub fusermount3: %v", err)
	}
	if err := utils.MkdirAllShared(filepath.Join(prefix, "bin")); err != nil {
		t.Fatalf("failed to create bin: %v", err)
	}

	if err := ensureFusermount3InBin(prefix); err != nil {
		t.Fatalf("ensureFusermount3InBin(): %v", err)
	}

	link, err := os.Readlink(filepath.Join(prefix, "bin", "fusermount3"))
	if err != nil {
		t.Fatalf("bin/fusermount3 not created as a symlink: %v", err)
	}
	if want := filepath.Join("..", "sbin", "fusermount3"); link != want {
		t.Errorf("bin/fusermount3 -> %q, want %q", link, want)
	}
}

func TestEnsureFusermount3InBinNoOpWhenAlreadyInBin(t *testing.T) {
	prefix := t.TempDir()
	bin := filepath.Join(prefix, "bin")
	if err := utils.MkdirAllShared(bin); err != nil {
		t.Fatalf("failed to create bin: %v", err)
	}
	if err := os.WriteFile(filepath.Join(bin, "fusermount3"), []byte("#!/bin/sh\n"), 0755); err != nil {
		t.Fatalf("failed to write stub fusermount3: %v", err)
	}

	if err := ensureFusermount3InBin(prefix); err != nil {
		t.Fatalf("ensureFusermount3InBin(): %v", err)
	}
	info, err := os.Lstat(filepath.Join(bin, "fusermount3"))
	if err != nil {
		t.Fatalf("bin/fusermount3 disappeared: %v", err)
	}
	if info.Mode()&os.ModeSymlink != 0 {
		t.Error("ensureFusermount3InBin() replaced an existing bin/fusermount3 with a symlink")
	}
}

func TestEnsureFusermount3InBinNoOpWhenNotInSbinEither(t *testing.T) {
	prefix := t.TempDir()
	if err := utils.MkdirAllShared(filepath.Join(prefix, "bin")); err != nil {
		t.Fatalf("failed to create bin: %v", err)
	}
	if err := ensureFusermount3InBin(prefix); err != nil {
		t.Fatalf("ensureFusermount3InBin(): %v", err)
	}
	if _, err := os.Lstat(filepath.Join(prefix, "bin", "fusermount3")); !os.IsNotExist(err) {
		t.Error("ensureFusermount3InBin() created bin/fusermount3 out of nothing")
	}
}

func TestRecoverStaleWithNothingUsable(t *testing.T) {
	tier := t.TempDir()
	target := filepath.Join(tier, "libexec")

	recoverStale(context.Background(), target)

	if isProvisioned(target) {
		t.Fatal("recoverStale() produced a provisioned target from nothing")
	}
	if _, err := os.Lstat(filepath.Join(tier, stagingName)); !os.IsNotExist(err) {
		t.Error("recoverStale() created a compat symlink with nothing to point at")
	}
}

func TestRecoverStaleRestoresOldWhenStagingIsMissing(t *testing.T) {
	tier := t.TempDir()
	target := filepath.Join(tier, "libexec")
	old := filepath.Join(tier, staleName)
	writeStubApptainer(t, old)

	recoverStale(context.Background(), target)

	if !isProvisioned(target) {
		t.Fatal("recoverStale() did not restore the outgoing generation")
	}
	if _, err := os.Stat(old); !os.IsNotExist(err) {
		t.Error("recoverStale() left .libexec.old behind after restoring it")
	}
	if link, err := os.Readlink(filepath.Join(tier, stagingName)); err != nil || link != target {
		t.Errorf("compat symlink = %q, %v; want -> %q", link, err, target)
	}
}

func TestRecoverStaleFinishesAVerifiedActivation(t *testing.T) {
	tier := t.TempDir()
	target := filepath.Join(tier, "libexec")
	staging := filepath.Join(tier, stagingName)
	old := filepath.Join(tier, staleName)
	provisionedVerifiableDir(t, staging)
	writeStubApptainer(t, old)

	recoverStale(context.Background(), target)

	if !isProvisioned(target) {
		t.Fatal("recoverStale() did not finish activating the verified staging build")
	}
	if info, err := os.Lstat(staging); err != nil || info.Mode()&os.ModeSymlink == 0 {
		t.Errorf("recoverStale() left the staging directory unactivated instead of the compat symlink: %v", err)
	}
	if _, err := os.Stat(old); !os.IsNotExist(err) {
		t.Error("recoverStale() left the outgoing generation behind after activating the new one")
	}
}

func TestRecoverStaleDiscardsUnverifiedStagingInFavorOfOld(t *testing.T) {
	tier := t.TempDir()
	target := filepath.Join(tier, "libexec")
	staging := filepath.Join(tier, stagingName)
	old := filepath.Join(tier, staleName)
	if err := utils.MkdirAllShared(staging); err != nil { // real dir, but no binaries in it
		t.Fatalf("failed to create incomplete staging dir: %v", err)
	}
	writeStubApptainer(t, old)

	recoverStale(context.Background(), target)

	if !isProvisioned(target) {
		t.Fatal("recoverStale() did not fall back to the outgoing generation")
	}
	if info, err := os.Lstat(staging); err != nil || info.Mode()&os.ModeSymlink == 0 {
		t.Errorf("recoverStale() kept the unverified staging directory around: %v", err)
	}
}

func TestRecoverStaleCleansUpDebrisAroundALiveGeneration(t *testing.T) {
	tier := t.TempDir()
	target := filepath.Join(tier, "libexec")
	staging := filepath.Join(tier, stagingName)
	old := filepath.Join(tier, staleName)
	writeStubApptainer(t, target)
	if err := utils.MkdirAllShared(staging); err != nil {
		t.Fatalf("failed to create stray staging dir: %v", err)
	}
	writeStubApptainer(t, old)

	recoverStale(context.Background(), target)

	if _, err := os.Stat(old); !os.IsNotExist(err) {
		t.Error("recoverStale() left stray .libexec.old behind a live generation")
	}
	info, err := os.Lstat(staging)
	if err != nil || info.Mode()&os.ModeSymlink == 0 {
		t.Errorf("recoverStale() did not replace stray staging debris with the compat symlink: %v", err)
	}
}

func writeFakeVersionBin(t *testing.T, bin, name, output string) {
	t.Helper()
	script := "#!/bin/sh\ncat <<'EOF'\n" + output + "EOF\n"
	if err := os.WriteFile(filepath.Join(bin, name), []byte(script), 0755); err != nil {
		t.Fatalf("failed to write stub %s: %v", name, err)
	}
}
