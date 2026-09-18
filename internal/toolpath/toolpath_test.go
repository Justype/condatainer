package toolpath

import (
	"errors"
	"os"
	"path/filepath"
	"runtime"
	"testing"

	"github.com/Justype/condatainer/internal/config"
	"github.com/Justype/condatainer/internal/utils"
)

// writeFakeExecutable writes an empty, executable file at dir/name.
func writeFakeExecutable(t *testing.T, dir, name string) string {
	t.Helper()
	p := filepath.Join(dir, name)
	if err := os.WriteFile(p, []byte("#!/bin/sh\n"), 0o755); err != nil {
		t.Fatalf("write fake executable: %v", err)
	}
	return p
}

// withProvisionedLibexec points the scratch tier at a fresh temp directory
// and drops a stub bin/apptainer + name there, satisfying libexec.Dir's
// marker check. Mirrors internal/libexec's own withScratchTier/
// provisionedStub, duplicated here rather than imported: libexec's test
// helpers are unexported.
func withProvisionedLibexec(t *testing.T, name string) string {
	t.Helper()
	scratch := filepath.Join(t.TempDir(), "condatainer")
	t.Setenv("SCRATCH", filepath.Dir(scratch))
	t.Setenv("XDG_DATA_HOME", "")
	t.Setenv("CNT_EXTRA_ROOT", "")
	t.Setenv("CNT_ROOT", "")
	config.InitDataPaths()

	bin := filepath.Join(scratch, "libexec", "bin")
	if err := utils.MkdirAllShared(bin); err != nil {
		t.Fatalf("failed to create stub bin dir: %v", err)
	}
	if err := os.WriteFile(filepath.Join(bin, "micromamba"), []byte("#!/bin/sh\n"), 0o755); err != nil {
		t.Fatalf("failed to write stub micromamba: %v", err)
	}
	if name != "" {
		writeFakeExecutable(t, bin, name)
	}
	return bin
}

func TestResolveFindsOnPath(t *testing.T) {
	if runtime.GOOS != "linux" {
		t.Skip("PATH/exec-bit semantics are POSIX-specific")
	}
	dir := t.TempDir()
	writeFakeExecutable(t, dir, "fake-tool")
	t.Setenv("PATH", dir)

	p, err := Resolve("fake-tool")
	if err != nil {
		t.Fatalf("Resolve: %v", err)
	}
	if filepath.Dir(p) != dir {
		t.Fatalf("resolved %s, want a binary under %s", p, dir)
	}
}

func TestResolveFallsBackToFHSDirs(t *testing.T) {
	if runtime.GOOS != "linux" {
		t.Skip("PATH/exec-bit semantics are POSIX-specific")
	}
	t.Setenv("PATH", t.TempDir()) // empty: nothing findable on PATH

	fhsDir := t.TempDir()
	want := writeFakeExecutable(t, fhsDir, "fake-tool")

	orig := fhsFallbackDirs
	fhsFallbackDirs = []string{fhsDir}
	t.Cleanup(func() { fhsFallbackDirs = orig })

	p, err := Resolve("fake-tool")
	if err != nil {
		t.Fatalf("Resolve: %v", err)
	}
	if p != want {
		t.Fatalf("Resolve() = %s, want %s", p, want)
	}
}

func TestResolveMissingReturnsErrToolMissing(t *testing.T) {
	if runtime.GOOS != "linux" {
		t.Skip("PATH/exec-bit semantics are POSIX-specific")
	}
	t.Setenv("PATH", t.TempDir())

	orig := fhsFallbackDirs
	fhsFallbackDirs = []string{t.TempDir()}
	t.Cleanup(func() { fhsFallbackDirs = orig })

	_, err := Resolve("does-not-exist-anywhere")
	if !errors.Is(err, ErrToolMissing) {
		t.Fatalf("Resolve() error = %v, want ErrToolMissing", err)
	}
}

// The provisioned toolchain wins even when the same name is also on PATH —
// it is required, not merely tried first among equals.
func TestResolvePrefersLibexecOverPath(t *testing.T) {
	if runtime.GOOS != "linux" {
		t.Skip("PATH/exec-bit semantics are POSIX-specific")
	}
	pathDir := t.TempDir()
	writeFakeExecutable(t, pathDir, "mksquashfs")
	t.Setenv("PATH", pathDir)

	libexecBin := withProvisionedLibexec(t, "mksquashfs")

	p, err := Resolve("mksquashfs")
	if err != nil {
		t.Fatalf("Resolve: %v", err)
	}
	if want := filepath.Join(libexecBin, "mksquashfs"); p != want {
		t.Errorf("Resolve(mksquashfs) = %s, want the provisioned %s", p, want)
	}
}

// A name libexec never provisions (debugfs, e2fsck, ...) isn't found under
// its bin/, and falls straight through to PATH unaffected.
func TestResolveFallsThroughForANameLibexecNeverProvisions(t *testing.T) {
	if runtime.GOOS != "linux" {
		t.Skip("PATH/exec-bit semantics are POSIX-specific")
	}
	pathDir := t.TempDir()
	want := writeFakeExecutable(t, pathDir, "debugfs")
	t.Setenv("PATH", pathDir)

	withProvisionedLibexec(t, "") // provisioned, but never wrote a debugfs

	p, err := Resolve("debugfs")
	if err != nil {
		t.Fatalf("Resolve: %v", err)
	}
	if p != want {
		t.Fatalf("Resolve(debugfs) = %s, want %s", p, want)
	}
}
