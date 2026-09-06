package container

import (
	"os"
	"os/exec"
	"path/filepath"
	"testing"

	"github.com/Justype/condatainer/internal/artifact/meta"
)

// packEnvSnapshot builds a real env-typed .sqf at path, with one conda-meta
// package file so PairedPackages has something to find.
func packEnvSnapshot(t *testing.T, path string) {
	t.Helper()
	root := t.TempDir()
	dir := filepath.Join(root, meta.DirName)
	if err := meta.StageRuntime(dir, envRuntime()); err != nil {
		t.Fatalf("StageRuntime: %v", err)
	}
	metaFile := filepath.Join(root, "cnt_env", "conda-meta", "numpy-1.24.0-py311h1234567_0.json")
	if err := os.MkdirAll(filepath.Dir(metaFile), 0o755); err != nil {
		t.Fatal(err)
	}
	if err := os.WriteFile(metaFile, []byte("{}"), 0o644); err != nil {
		t.Fatal(err)
	}
	cmd := exec.Command("mksquashfs", root, path, "-no-progress", "-noappend", "-quiet", "-no-xattrs")
	if output, err := cmd.CombinedOutput(); err != nil {
		t.Fatalf("mksquashfs: %v\n%s", err, output)
	}
}

// PairedPackages unions a paired snapshot's packages with the .img's own,
// which is what makes #IMG_PACKAGES: checks pass against a thin .img sitting
// on a snapshot instead of falsely reporting everything as not installed.
func TestPairedPackagesUnionsSnapshotAndImg(t *testing.T) {
	requireSquashfsTools(t)
	t.Setenv("USER", "alice")
	dir := t.TempDir()
	snapshot := filepath.Join(dir, "env-alice.sqf")
	packEnvSnapshot(t, snapshot)

	// A placeholder .img: PairedPackages's own conda.ListCondaPackages call
	// against it will fail to find debugfs content (no real ext3 filesystem
	// here), which is fine — this test exercises the union with a snapshot
	// found via LookupSnapshot, already covered end to end.
	img := filepath.Join(dir, "env-alice.img")
	if err := os.WriteFile(img, []byte{}, 0o644); err != nil {
		t.Fatal(err)
	}

	pkgs, err := PairedPackages(img)
	if err != nil {
		t.Fatalf("PairedPackages: %v", err)
	}
	if pkgs["numpy"] != "1.24.0" {
		t.Fatalf("pkgs = %v, want numpy=1.24.0 from the snapshot", pkgs)
	}
}

// A bare .sqf (no paired .img at all) is checked directly — the third
// overlay state, where ResolveEnvOverlayInDir resolves EnvImg to the
// snapshot itself.
func TestPairedPackagesBareSqf(t *testing.T) {
	requireSquashfsTools(t)
	dir := t.TempDir()
	snapshot := filepath.Join(dir, "env.sqf")
	packEnvSnapshot(t, snapshot)

	pkgs, err := PairedPackages(snapshot)
	if err != nil {
		t.Fatalf("PairedPackages: %v", err)
	}
	if pkgs["numpy"] != "1.24.0" {
		t.Fatalf("pkgs = %v, want numpy=1.24.0 from the bare snapshot", pkgs)
	}
}
