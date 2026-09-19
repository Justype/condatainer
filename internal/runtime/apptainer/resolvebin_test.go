package apptainer

import (
	"bytes"
	"context"
	"errors"
	"os"
	"path/filepath"
	"strings"
	"testing"

	"github.com/Justype/condatainer/internal/config"
)

// withLibexecTier points the libexec search path at one scratch directory,
// so a test decides for itself whether the self-provisioned toolchain looks
// provisioned.
func withLibexecTier(t *testing.T) string {
	t.Helper()
	dir := filepath.Join(t.TempDir(), "libexec")
	prev := config.GlobalDataPaths
	config.GlobalDataPaths.LibexecDirs = []string{dir}
	t.Cleanup(func() { config.GlobalDataPaths = prev })
	return dir
}

// writeFakeBin writes an executable shell script at dir/bin/name that prints
// output verbatim to stdout, standing in for a real apptainer/singularity
// binary so Version's "--version" subprocess has something to read.
func writeFakeBin(t *testing.T, dir, name, output string) string {
	t.Helper()
	bin := filepath.Join(dir, "bin")
	if err := os.MkdirAll(bin, 0o755); err != nil {
		t.Fatal(err)
	}
	// micromamba is what marks a libexec tier provisioned.
	if err := os.WriteFile(filepath.Join(bin, "micromamba"), []byte("#!/bin/sh\n"), 0o755); err != nil {
		t.Fatal(err)
	}
	path := filepath.Join(bin, name)
	script := "#!/bin/sh\necho '" + output + "'\n"
	if err := os.WriteFile(path, []byte(script), 0o755); err != nil {
		t.Fatal(err)
	}
	return path
}

// resetApptainerState clears the remembered binary between tests and pretends
// to be on a host, so one test's resolution does not leak into the next.
func resetApptainerState(t *testing.T) {
	t.Helper()
	prevLast, prevInside := last, insideContainer
	last, insideContainer = Bin{}, func() bool { return false }
	t.Cleanup(func() { last, insideContainer = prevLast, prevInside })
}

// systemApptainer points the configured system/module binary at path.
func systemApptainer(t *testing.T, path string) {
	t.Helper()
	prev := config.Global.Build.SystemApptainer
	config.Global.Build.SystemApptainer = path
	t.Cleanup(func() { config.Global.Build.SystemApptainer = prev })
}

// With no installed libexec apptainer and no usable system one, Normal refuses
// and names both ways out.
func TestNormalRefusesWithNothingUsable(t *testing.T) {
	withLibexecTier(t)
	resetApptainerState(t)
	systemApptainer(t, filepath.Join(t.TempDir(), "no-apptainer"))

	_, err := Normal()
	if err == nil {
		t.Fatal("Normal succeeded with no apptainer anywhere")
	}
	for _, want := range []string{"update --libexec apptainer", "module"} {
		if !strings.Contains(err.Error(), want) {
			t.Errorf("err = %v, want it to mention %q", err, want)
		}
	}
}

// Inside a container the same refusal points at nested running, not at a fix
// that is refused there.
func TestNormalInsideContainerPointsAtNestedRun(t *testing.T) {
	withLibexecTier(t)
	resetApptainerState(t)
	insideContainer = func() bool { return true }
	systemApptainer(t, "")

	_, err := Normal()
	if err == nil || !strings.Contains(err.Error(), "nested_run") || strings.Contains(err.Error(), "update --libexec") {
		t.Fatalf("Normal error = %v, want it to mention nested_run and not update --libexec", err)
	}
}

// Without libexec's apptainer, a system one that clears the floor is used.
func TestNormalFallsBackToSystem(t *testing.T) {
	withLibexecTier(t)
	resetApptainerState(t)
	sysBin := writeFakeBin(t, t.TempDir(), "apptainer", "apptainer version 1.5.3")
	systemApptainer(t, sysBin)

	bin, err := Normal()
	if err != nil {
		t.Fatalf("Normal: %v", err)
	}
	if bin.Path != sysBin || bin.Libexec {
		t.Errorf("Normal = %+v, want the system binary %q", bin, sysBin)
	}
}

// A system apptainer below the zstd floor is refused, with the reason.
func TestNormalRefusesOldSystemApptainer(t *testing.T) {
	withLibexecTier(t)
	resetApptainerState(t)
	systemApptainer(t, writeFakeBin(t, t.TempDir(), "apptainer", "apptainer version 1.3.9"))

	if _, err := Normal(); err == nil || !strings.Contains(err.Error(), "zstd") {
		t.Fatalf("Normal error = %v, want it to name zstd", err)
	}
}

// An apptainer installed in libexec wins over the system binary.
func TestNormalUsesLibexec(t *testing.T) {
	dir := withLibexecTier(t)
	resetApptainerState(t)
	binPath := writeFakeBin(t, dir, "apptainer", "apptainer version 1.5.3")
	systemApptainer(t, "/should/not/be/used")

	bin, err := Normal()
	if err != nil {
		t.Fatalf("Normal: %v", err)
	}
	if bin.Path != binPath || !bin.Libexec {
		t.Errorf("Normal = %+v, want the libexec binary %q", bin, binPath)
	}
}

// Fakeroot and ForBuild use the system/module binary whatever the toolchain
// state: only its setuid starter can escalate, which libexec's non-setuid
// apptainer cannot.
func TestFakerootAndForBuildUseSystemBinary(t *testing.T) {
	dir := withLibexecTier(t)
	resetApptainerState(t)
	writeFakeBin(t, dir, "apptainer", "apptainer version 1.5.3")
	sysBin := writeFakeBin(t, t.TempDir(), "apptainer", "apptainer version 1.5.3")
	systemApptainer(t, sysBin)

	for name, resolve := range map[string]func() (Bin, error){"Fakeroot": Fakeroot, "ForBuild": ForBuild} {
		bin, err := resolve()
		if err != nil {
			t.Fatalf("%s: %v", name, err)
		}
		if bin.Path != sysBin || bin.Libexec {
			t.Errorf("%s = %+v, want the system binary %q", name, bin, sysBin)
		}
	}
}

// A system apptainer below the zstd floor cannot mount what condatainer packs,
// and singularity cannot at any version, so a fakeroot exec refuses rather than
// failing later inside Apptainer's own mount step.
func TestFakerootRefusesOldApptainerAndSingularity(t *testing.T) {
	resetApptainerState(t)
	for _, c := range []struct{ name, output, want string }{
		{"apptainer", "apptainer version 1.3.9", "zstd"},
		{"singularity", "singularity-ce version 4.1.1", "singularity"},
	} {
		systemApptainer(t, writeFakeBin(t, t.TempDir(), c.name, c.output))
		if _, err := Fakeroot(); err == nil || !strings.Contains(err.Error(), c.want) {
			t.Errorf("Fakeroot with %s: error = %v, want it to name %q", c.name, err, c.want)
		}
	}
}

// A definition build has no version check: its output is a sandbox.
func TestForBuildAcceptsAnySystemBinary(t *testing.T) {
	resetApptainerState(t)
	systemApptainer(t, writeFakeBin(t, t.TempDir(), "singularity", "singularity-ce version 3.0.0"))

	if _, err := ForBuild(); err != nil {
		t.Fatalf("ForBuild: %v", err)
	}
}

// A container has no setuid starter, so both refuse at once.
func TestFakerootAndForBuildRefuseInsideContainer(t *testing.T) {
	resetApptainerState(t)
	insideContainer = func() bool { return true }
	systemApptainer(t, writeFakeBin(t, t.TempDir(), "apptainer", "apptainer version 1.5.3"))

	for name, resolve := range map[string]func() (Bin, error){"Fakeroot": Fakeroot, "ForBuild": ForBuild} {
		if _, err := resolve(); !errors.Is(err, ErrNeedsHost) {
			t.Errorf("%s error = %v, want ErrNeedsHost", name, err)
		}
	}
}

// Apptainer shells out to unsquashfs/mksquashfs for some of its own
// operations (extracting a .sqf to build a sandbox, or as a mount fallback),
// on its own PATH — not anything a launched container's environment
// controls. When the binary is the self-provisioned libexec copy, its own
// bin/ (where those tools are provisioned alongside it) must be on that PATH
// too, or Apptainer cannot find them itself.
func TestRunApptainerPrependsLibexecBinToPath(t *testing.T) {
	dir := withLibexecTier(t)
	resetApptainerState(t)
	binPath := writeFakeBin(t, dir, "apptainer", "unused")
	if err := os.WriteFile(binPath, []byte("#!/bin/sh\necho \"$PATH\"\n"), 0o755); err != nil {
		t.Fatal(err)
	}

	bin, err := Normal()
	if err != nil {
		t.Fatalf("Normal: %v", err)
	}

	var out bytes.Buffer
	if err := runApptainerWithOutput(context.Background(), bin, "exec", "", false, nil, &out, &out, nil); err != nil {
		t.Fatalf("runApptainerWithOutput: %v", err)
	}

	libexecBin := filepath.Dir(binPath)
	if !strings.Contains(out.String(), libexecBin) {
		t.Errorf("subprocess PATH = %q, want it to contain %q", out.String(), libexecBin)
	}
}

// A binary's version is read from the per-user cache before the binary is run,
// so a new process does not run it again while it is unchanged.
func TestVersionUsesThePerUserCache(t *testing.T) {
	if insideContainer() {
		t.Skip("the cache is off inside a container")
	}
	counter := filepath.Join(t.TempDir(), "runs")
	path := filepath.Join(t.TempDir(), "apptainer")
	script := "#!/bin/sh\necho x >> '" + counter + "'\necho 'apptainer version 1.5.3'\n"
	if err := os.WriteFile(path, []byte(script), 0o755); err != nil {
		t.Fatal(err)
	}
	bin := Bin{Path: path}

	for i := 0; i < 2; i++ {
		versions.Delete(path) // a new process has no in-memory answer
		if v, err := bin.Version(); err != nil || v != "1.5.3" {
			t.Fatalf("Version = %q, %v", v, err)
		}
	}
	data, _ := os.ReadFile(counter)
	if runs := strings.Count(string(data), "x"); runs != 1 {
		t.Errorf("the binary ran %d times, want 1", runs)
	}
}
