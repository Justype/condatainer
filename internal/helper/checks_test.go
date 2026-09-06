package helper

import (
	"context"
	"os"
	"os/exec"
	"path/filepath"
	"testing"

	"github.com/Justype/condatainer/catalog"
	"github.com/Justype/condatainer/internal/artifact/meta"
)

func requireSquashfsTools(t *testing.T) {
	t.Helper()
	for _, bin := range []string{"mksquashfs", "unsquashfs"} {
		if _, err := exec.LookPath(bin); err != nil {
			t.Skipf("%s not available", bin)
		}
	}
}

// packEnvSnapshot builds a real env-typed .sqf at path, with one conda-meta
// package file so ListCondaPackagesSqf has something to find.
func packEnvSnapshot(t *testing.T, path string) {
	t.Helper()
	root := t.TempDir()
	dir := filepath.Join(root, meta.DirName)
	rt := meta.Runtime{
		SchemaVersion: meta.SchemaVersion,
		Name:          meta.EnvName,
		Type:          catalog.TypeEnv,
		Platform:      meta.NativePlatform(),
		Prefix:        meta.EnvPrefix,
	}
	if err := meta.StageRuntime(dir, rt); err != nil {
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

// CheckEnv reports the pair's combined size and names the paired snapshot,
// rather than the thin .img's own small size alone.
func TestCheckEnvReportsPairedSnapshot(t *testing.T) {
	requireSquashfsTools(t)
	t.Setenv("USER", "alice")
	dir := t.TempDir()
	snapshot := filepath.Join(dir, "env-alice.sqf")
	packEnvSnapshot(t, snapshot)
	snapInfo, err := os.Stat(snapshot)
	if err != nil {
		t.Fatal(err)
	}

	img := filepath.Join(dir, "env-alice.img")
	if err := os.WriteFile(img, make([]byte, 4096), 0o644); err != nil {
		t.Fatal(err)
	}

	st, err := CheckEnv(context.Background(), img)
	if err != nil {
		t.Fatalf("CheckEnv: %v", err)
	}
	if st.Snapshot != snapshot {
		t.Fatalf("Snapshot = %q, want %q", st.Snapshot, snapshot)
	}
	wantSizeMB := (int64(4096)/(1024*1024) + snapInfo.Size()/(1024*1024))
	if st.SizeMB != wantSizeMB {
		t.Fatalf("SizeMB = %d, want %d (img + snapshot)", st.SizeMB, wantSizeMB)
	}
}

// A .img with no paired snapshot reports its own size only, and no Snapshot.
func TestCheckEnvNoSnapshot(t *testing.T) {
	t.Setenv("USER", "alice")
	dir := t.TempDir()
	img := filepath.Join(dir, "env-alice.img")
	if err := os.WriteFile(img, make([]byte, 4096), 0o644); err != nil {
		t.Fatal(err)
	}

	st, err := CheckEnv(context.Background(), img)
	if err != nil {
		t.Fatalf("CheckEnv: %v", err)
	}
	if st.Snapshot != "" {
		t.Fatalf("Snapshot = %q, want none", st.Snapshot)
	}
}

// The third overlay state: no .img on disk at all, but a paired snapshot
// sitting where one would go. CheckEnv reports Exists=false (nothing to edit)
// alongside a non-empty Snapshot (not the ordinary "start from scratch" case
// either) — the dashboard tells the two apart by that combination.
func TestCheckEnvThirdStateImgMissingSnapshotPresent(t *testing.T) {
	requireSquashfsTools(t)
	t.Setenv("USER", "alice")
	dir := t.TempDir()
	snapshot := filepath.Join(dir, "env-alice.sqf")
	packEnvSnapshot(t, snapshot)
	img := filepath.Join(dir, "env-alice.img")

	st, err := CheckEnv(context.Background(), img)
	if err != nil {
		t.Fatalf("CheckEnv: %v", err)
	}
	if st.Exists {
		t.Fatal("Exists = true, want false: the .img is not on disk")
	}
	if st.Snapshot != snapshot {
		t.Fatalf("Snapshot = %q, want %q", st.Snapshot, snapshot)
	}
}

// FindEnvSnapshot finds a shared snapshot with no .img yet.
func TestFindEnvSnapshotSharedLine(t *testing.T) {
	requireSquashfsTools(t)
	t.Setenv("USER", "alice")
	dir := t.TempDir()
	shared := filepath.Join(dir, "env.sqf")
	packEnvSnapshot(t, shared)

	if got := FindEnvSnapshot(dir); got != shared {
		t.Fatalf("FindEnvSnapshot = %q, want %q", got, shared)
	}
}

// A personal snapshot line is checked before the shared one.
func TestFindEnvSnapshotPersonalLine(t *testing.T) {
	requireSquashfsTools(t)
	t.Setenv("USER", "alice")
	dir := t.TempDir()
	personal := filepath.Join(dir, "env-alice.sqf")
	packEnvSnapshot(t, personal)

	if got := FindEnvSnapshot(dir); got != personal {
		t.Fatalf("FindEnvSnapshot = %q, want %q", got, personal)
	}
}

// No snapshot in cwd means nothing found.
func TestFindEnvSnapshotNone(t *testing.T) {
	dir := t.TempDir()
	if got := FindEnvSnapshot(dir); got != "" {
		t.Fatalf("FindEnvSnapshot = %q, want \"\"", got)
	}
}
