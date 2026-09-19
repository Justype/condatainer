package toolpath

import (
	"errors"
	"os"
	"path/filepath"
	"runtime"
	"strings"
	"sync"
	"testing"
	"time"

	"github.com/Justype/condatainer/internal/config"
	"github.com/Justype/condatainer/internal/utils"
)

// TestMain keeps every test off the real per-user cache.
func TestMain(m *testing.M) {
	cachePath = func() string { return "" }
	os.Exit(m.Run())
}

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
// and drops a stub bin/micromamba + name there, satisfying libexec.Dir's
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

// writeScript writes an executable shell script printing output at dir/name.
func writeScript(t *testing.T, dir, name, body string) string {
	t.Helper()
	p := filepath.Join(dir, name)
	if err := os.WriteFile(p, []byte("#!/bin/sh\n"+body+"\n"), 0o755); err != nil {
		t.Fatal(err)
	}
	return p
}

// noHostApptainer keeps the host's own apptainer out of a test's search.
func noHostApptainer(t *testing.T) {
	t.Helper()
	prev := config.Global.Build.SystemApptainer
	config.Global.Build.SystemApptainer = filepath.Join(t.TempDir(), "no-apptainer")
	t.Cleanup(func() { config.Global.Build.SystemApptainer = prev })
}

// A host mksquashfs below the version floor is skipped, and the error says why.
func TestResolveSkipsAHostToolBelowItsFloor(t *testing.T) {
	if runtime.GOOS != "linux" {
		t.Skip("PATH/exec-bit semantics are POSIX-specific")
	}
	withScratchWithoutLibexec(t)
	noHostApptainer(t)
	orig := fhsFallbackDirs
	fhsFallbackDirs = []string{t.TempDir()}
	t.Cleanup(func() { fhsFallbackDirs = orig })

	old := t.TempDir()
	writeScript(t, old, "mksquashfs", `echo "mksquashfs version 4.3 (2014/05/12)"`)
	t.Setenv("PATH", old)
	_, err := Resolve("mksquashfs")
	if !errors.Is(err, ErrToolMissing) || !strings.Contains(err.Error(), "4.3") {
		t.Fatalf("Resolve(mksquashfs) error = %v, want ErrToolMissing naming version 4.3", err)
	}

	current := t.TempDir()
	want := writeScript(t, current, "mksquashfs", `echo "mksquashfs version 4.6.1 (2023/03/25)"`)
	t.Setenv("PATH", old+string(os.PathListSeparator)+current)
	got, err := Resolve("mksquashfs")
	if err != nil || got != want {
		t.Errorf("Resolve(mksquashfs) = %q, %v; want %q (the old one skipped)", got, err, want)
	}
}

// With nothing on PATH, a tool bundled with the host apptainer is found in
// the directory its `buildcfg` reports.
func TestResolveFindsAToolBundledWithApptainer(t *testing.T) {
	if runtime.GOOS != "linux" {
		t.Skip("PATH/exec-bit semantics are POSIX-specific")
	}
	withScratchWithoutLibexec(t)
	root := t.TempDir()
	bundleBin := filepath.Join(root, "libexec", "apptainer", "bin")
	if err := os.MkdirAll(bundleBin, 0o755); err != nil {
		t.Fatal(err)
	}
	want := writeScript(t, bundleBin, "squashfuse_ll", "true")
	apptainerDir := t.TempDir()
	// The binary's base name selects <libexecdir>/<name>/bin, so it is called apptainer.
	fake := writeScript(t, apptainerDir, "apptainer", `echo "LIBEXECDIR=`+filepath.Join(root, "libexec")+`"`)

	prev := config.Global.Build.SystemApptainer
	config.Global.Build.SystemApptainer = fake
	t.Cleanup(func() { config.Global.Build.SystemApptainer = prev })
	orig := fhsFallbackDirs
	fhsFallbackDirs = []string{t.TempDir()}
	t.Cleanup(func() { fhsFallbackDirs = orig })
	t.Setenv("PATH", t.TempDir())

	got, err := Resolve("squashfuse_ll")
	if err != nil || got != want {
		t.Errorf("Resolve(squashfuse_ll) = %q, %v; want the bundled %q", got, err, want)
	}
}

// withScratchWithoutLibexec points the search at an empty scratch tier.
func withScratchWithoutLibexec(t *testing.T) {
	t.Helper()
	scratch := filepath.Join(t.TempDir(), "condatainer")
	t.Setenv("SCRATCH", filepath.Dir(scratch))
	t.Setenv("XDG_DATA_HOME", "")
	t.Setenv("CNT_EXTRA_ROOT", "")
	t.Setenv("CNT_ROOT", "")
	config.InitDataPaths()
}

// The host apptainer is only asked where its bundled tools are when PATH did
// not already have the tool.
func TestResolveDoesNotRunApptainerWhenPathHasTheTool(t *testing.T) {
	if runtime.GOOS != "linux" {
		t.Skip("PATH/exec-bit semantics are POSIX-specific")
	}
	withScratchWithoutLibexec(t)
	marker := filepath.Join(t.TempDir(), "asked")
	fake := writeScript(t, t.TempDir(), "apptainer", "touch "+marker)
	prev := config.Global.Build.SystemApptainer
	config.Global.Build.SystemApptainer = fake
	t.Cleanup(func() { config.Global.Build.SystemApptainer = prev })

	dir := t.TempDir()
	want := writeFakeExecutable(t, dir, "fuse2fs")
	t.Setenv("PATH", dir)

	got, err := Resolve("fuse2fs")
	if err != nil || got != want {
		t.Fatalf("Resolve(fuse2fs) = %q, %v; want %q", got, err, want)
	}
	if _, err := os.Stat(marker); err == nil {
		t.Error("Resolve ran the host apptainer although PATH had the tool")
	}
}

// resetMemoryCaches makes the next Resolve behave like a new process.
func resetMemoryCaches() {
	for _, m := range []*sync.Map{&floors, &bundled} {
		m.Range(func(k, _ any) bool { m.Delete(k); return true })
	}
}

// A new process reuses what an earlier one worked out, until the binary changes.
func TestResolveReusesTheOnDiskCache(t *testing.T) {
	if runtime.GOOS != "linux" {
		t.Skip("PATH/exec-bit semantics are POSIX-specific")
	}
	withScratchWithoutLibexec(t)
	cache := filepath.Join(t.TempDir(), "cache", cacheFileName)
	cachePath = func() string { return cache }
	t.Cleanup(func() { cachePath = func() string { return "" }; resetMemoryCaches() })

	root := t.TempDir()
	bundleBin := filepath.Join(root, "libexec", "apptainer", "bin")
	if err := os.MkdirAll(bundleBin, 0o755); err != nil {
		t.Fatal(err)
	}
	askedBuildcfg := filepath.Join(root, "buildcfg-runs")
	askedVersion := filepath.Join(root, "version-runs")
	writeScript(t, bundleBin, "mksquashfs", `echo x >> `+askedVersion+`; echo "mksquashfs version 4.7.5 (2026/03/01)"`)
	fake := writeScript(t, t.TempDir(), "apptainer", `echo x >> `+askedBuildcfg+`; echo "LIBEXECDIR=`+filepath.Join(root, "libexec")+`"`)

	prev := config.Global.Build.SystemApptainer
	config.Global.Build.SystemApptainer = fake
	t.Cleanup(func() { config.Global.Build.SystemApptainer = prev })
	orig := fhsFallbackDirs
	fhsFallbackDirs = []string{t.TempDir()}
	t.Cleanup(func() { fhsFallbackDirs = orig })
	t.Setenv("PATH", t.TempDir())

	runs := func(file string) int {
		data, _ := os.ReadFile(file)
		return len(data) / 2
	}
	for i := 0; i < 3; i++ {
		resetMemoryCaches() // a new process each time
		if _, err := Resolve("mksquashfs"); err != nil {
			t.Fatalf("Resolve #%d: %v", i, err)
		}
	}
	if got := runs(askedBuildcfg); got != 1 {
		t.Errorf("apptainer buildcfg ran %d times across 3 processes, want 1", got)
	}
	if got := runs(askedVersion); got != 1 {
		t.Errorf("mksquashfs -version ran %d times across 3 processes, want 1", got)
	}

	// A changed binary is not trusted: it is asked again.
	future := time.Now().Add(time.Hour)
	if err := os.Chtimes(fake, future, future); err != nil {
		t.Fatal(err)
	}
	resetMemoryCaches()
	if _, err := Resolve("mksquashfs"); err != nil {
		t.Fatal(err)
	}
	if got := runs(askedBuildcfg); got != 2 {
		t.Errorf("apptainer buildcfg ran %d times after the binary changed, want 2", got)
	}
}

// Inside a container there is no on-disk cache to read or write.
func TestNoOnDiskCacheInsideAContainer(t *testing.T) {
	t.Setenv("XDG_CACHE_HOME", t.TempDir())
	config.InitDataPaths()
	prev := insideContainer
	t.Cleanup(func() { insideContainer = prev })

	insideContainer = func() bool { return true }
	if got := defaultCachePath(); got != "" {
		t.Errorf("defaultCachePath() inside a container = %q, want none", got)
	}
	insideContainer = func() bool { return false }
	if got := defaultCachePath(); got == "" || filepath.Base(got) != cacheFileName {
		t.Errorf("defaultCachePath() on a host = %q, want the personal cache file", got)
	}
}

// Remember runs compute once per unchanged binary, and never stores an error.
func TestRememberComputesOnceAndNotErrors(t *testing.T) {
	cache := filepath.Join(t.TempDir(), cacheFileName)
	cachePath = func() string { return cache }
	t.Cleanup(func() { cachePath = func() string { return "" } })
	bin := writeFakeExecutable(t, t.TempDir(), "tool")

	runs := 0
	compute := func() (string, error) { runs++; return "1.5.3", nil }
	for i := 0; i < 2; i++ {
		if got, err := Remember("version", bin, compute); err != nil || got != "1.5.3" {
			t.Fatalf("Remember = %q, %v", got, err)
		}
	}
	if runs != 1 {
		t.Errorf("compute ran %d times, want 1", runs)
	}

	failed := errors.New("boom")
	other := writeFakeExecutable(t, t.TempDir(), "tool")
	for i := 0; i < 2; i++ {
		if _, err := Remember("version", other, func() (string, error) { runs++; return "", failed }); !errors.Is(err, failed) {
			t.Fatalf("err = %v, want %v", err, failed)
		}
	}
	if runs != 3 {
		t.Errorf("a failure was stored: compute ran %d times in total, want 3", runs)
	}
}
