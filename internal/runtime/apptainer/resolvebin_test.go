package apptainer

import (
	"bytes"
	"context"
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
// binary so GetVersion's "--version" subprocess has something to read.
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

// resetApptainerState clears the package's cached binary/version between
// tests, so one test's SetBin does not leak into the next.
func resetApptainerState(t *testing.T) {
	t.Helper()
	prevCmd, prevVersion := apptainerCmd, cachedVersion
	apptainerCmd, cachedVersion = "", ""
	t.Cleanup(func() { apptainerCmd, cachedVersion = prevCmd, prevVersion })
}

// systemApptainer points the configured system/module binary at path.
func systemApptainer(t *testing.T, path string) {
	t.Helper()
	prev := config.Global.Build.SystemApptainer
	config.Global.Build.SystemApptainer = path
	t.Cleanup(func() { config.Global.Build.SystemApptainer = prev })
}

// With no installed libexec apptainer and no usable system one, a
// non-fakeroot resolution refuses and names both ways out.
func TestResolveBinNonFakerootRefusesWithNothingUsable(t *testing.T) {
	withLibexecTier(t)
	resetApptainerState(t)
	systemApptainer(t, filepath.Join(t.TempDir(), "no-apptainer"))

	err := ResolveBin(false)
	if err == nil {
		t.Fatal("ResolveBin(false) succeeded with no apptainer anywhere")
	}
	for _, want := range []string{"update --libexec apptainer", "module"} {
		if !strings.Contains(err.Error(), want) {
			t.Errorf("err = %v, want it to mention %q", err, want)
		}
	}
}

// Without libexec's apptainer, a system one that clears the floor is used.
func TestResolveBinNonFakerootFallsBackToSystem(t *testing.T) {
	withLibexecTier(t)
	resetApptainerState(t)
	sysBin := writeFakeBin(t, t.TempDir(), "apptainer", "apptainer version 1.5.3")
	systemApptainer(t, sysBin)

	if err := ResolveBin(false); err != nil {
		t.Fatalf("ResolveBin(false): %v", err)
	}
	if apptainerCmd != sysBin {
		t.Errorf("configured binary = %q, want the system one %q", apptainerCmd, sysBin)
	}
}

// A system apptainer below the zstd floor is refused for a non-fakeroot exec
// too, with the reason.
func TestResolveBinNonFakerootRefusesOldSystemApptainer(t *testing.T) {
	withLibexecTier(t)
	resetApptainerState(t)
	systemApptainer(t, writeFakeBin(t, t.TempDir(), "apptainer", "apptainer version 1.3.9"))

	err := ResolveBin(false)
	if err == nil || !strings.Contains(err.Error(), "zstd") {
		t.Fatalf("ResolveBin(false) error = %v, want it to name zstd", err)
	}
}

// An apptainer installed in libexec wins over the system binary.
func TestResolveBinNonFakerootUsesLibexec(t *testing.T) {
	dir := withLibexecTier(t)
	resetApptainerState(t)
	binPath := writeFakeBin(t, dir, "apptainer", "apptainer version 1.5.3")

	prevBin := config.Global.Build.SystemApptainer
	config.Global.Build.SystemApptainer = "/should/not/be/used"
	t.Cleanup(func() { config.Global.Build.SystemApptainer = prevBin })

	if err := ResolveBin(false); err != nil {
		t.Fatalf("ResolveBin(false): %v", err)
	}
	if apptainerCmd != binPath {
		t.Errorf("configured binary = %q, want the libexec one %q", apptainerCmd, binPath)
	}
}

// Fakeroot always resolves to the system/module binary, whatever the
// toolchain state — it can mount its own overlays via a setuid starter,
// which libexec's non-setuid apptainer cannot.
func TestResolveBinFakerootUsesSystemBinary(t *testing.T) {
	resetApptainerState(t)
	sysDir := t.TempDir()
	sysBin := writeFakeBin(t, sysDir, "apptainer", "apptainer version 1.5.3")

	prevBin := config.Global.Build.SystemApptainer
	config.Global.Build.SystemApptainer = sysBin
	t.Cleanup(func() { config.Global.Build.SystemApptainer = prevBin })

	if err := ResolveBin(true); err != nil {
		t.Fatalf("ResolveBin(true): %v", err)
	}
	if apptainerCmd != sysBin {
		t.Errorf("configured binary = %q, want the system one %q", apptainerCmd, sysBin)
	}
}

// A system apptainer below the zstd floor cannot mount what condatainer
// packs, so a fakeroot exec refuses rather than failing later inside
// Apptainer's own mount step with a message naming neither cause.
func TestResolveBinFakerootRefusesOldApptainer(t *testing.T) {
	resetApptainerState(t)
	sysDir := t.TempDir()
	sysBin := writeFakeBin(t, sysDir, "apptainer", "apptainer version 1.3.9")

	prevBin := config.Global.Build.SystemApptainer
	config.Global.Build.SystemApptainer = sysBin
	t.Cleanup(func() { config.Global.Build.SystemApptainer = prevBin })

	err := ResolveBin(true)
	if err == nil {
		t.Fatal("ResolveBin(true) accepted an apptainer below the zstd floor")
	}
	if !strings.Contains(err.Error(), "zstd") {
		t.Errorf("err = %v, want it to name zstd as the reason", err)
	}
}

// Singularity is refused outright for fakeroot, regardless of version: it
// cannot mount a zstd-compressed SquashFS the way Apptainer >= 1.4 can.
func TestResolveBinFakerootRefusesSingularity(t *testing.T) {
	resetApptainerState(t)
	sysDir := t.TempDir()
	sysBin := writeFakeBin(t, sysDir, "singularity", "singularity-ce version 4.1.1")

	prevBin := config.Global.Build.SystemApptainer
	config.Global.Build.SystemApptainer = sysBin
	t.Cleanup(func() { config.Global.Build.SystemApptainer = prevBin })

	err := ResolveBin(true)
	if err == nil {
		t.Fatal("ResolveBin(true) accepted singularity")
	}
	if !strings.Contains(err.Error(), "singularity") {
		t.Errorf("err = %v, want it to name singularity as the reason", err)
	}
}

// Apptainer shells out to unsquashfs/mksquashfs for some of its own
// operations (extracting a .sqf to build a sandbox, or as a mount fallback),
// on its own PATH — not anything a launched container's environment
// controls. When the resolved binary is the self-provisioned libexec copy,
// its own bin/ (where those tools are provisioned alongside it) must be on
// that PATH too, or Apptainer cannot find them itself.
func TestRunApptainerPrependsLibexecBinToPath(t *testing.T) {
	dir := withLibexecTier(t)
	resetApptainerState(t)
	binPath := writeFakeBin(t, dir, "apptainer", "unused")
	if err := os.WriteFile(binPath, []byte("#!/bin/sh\necho \"$PATH\"\n"), 0o755); err != nil {
		t.Fatal(err)
	}

	if err := ResolveBin(false); err != nil {
		t.Fatalf("ResolveBin(false): %v", err)
	}

	var out bytes.Buffer
	if err := runApptainerWithOutput(context.Background(), "exec", "", false, nil, &out, &out, nil); err != nil {
		t.Fatalf("runApptainerWithOutput: %v", err)
	}

	libexecBin := filepath.Dir(binPath)
	if !strings.Contains(out.String(), libexecBin) {
		t.Errorf("subprocess PATH = %q, want it to contain %q", out.String(), libexecBin)
	}
}
