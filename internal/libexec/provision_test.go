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
	if _, err := os.Stat(filepath.Join(scratch, "libexec.new")); !os.IsNotExist(err) {
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

func writeFakeVersionBin(t *testing.T, bin, name, output string) {
	t.Helper()
	script := "#!/bin/sh\ncat <<'EOF'\n" + output + "EOF\n"
	if err := os.WriteFile(filepath.Join(bin, name), []byte(script), 0755); err != nil {
		t.Fatalf("failed to write stub %s: %v", name, err)
	}
}
