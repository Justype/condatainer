package config

import (
	"path/filepath"
	"sync"
	"testing"
)

// withAllTiers points each of the four data tiers at its own directory and
// resets the sync.Once caches that GetRootDir/GetExtraRootDir memoize.
func withAllTiers(t *testing.T) (scratch, user, extraRoot, root string) {
	t.Helper()
	base := t.TempDir()
	scratch = filepath.Join(base, "scratch")
	user = filepath.Join(base, "xdg")
	extraRoot = filepath.Join(base, "lab")
	root = filepath.Join(base, "install")

	t.Setenv("SCRATCH", scratch)
	t.Setenv("XDG_DATA_HOME", user)
	t.Setenv("CNT_EXTRA_ROOT", extraRoot)
	t.Setenv("CNT_ROOT", root)

	reset := func() {
		extraRootOnce, extraRootCache = sync.Once{}, ""
		rootDirOnce, rootDirCache = sync.Once{}, ""
	}
	reset()
	t.Cleanup(reset)

	// GetScratchDataDir/GetUserDataDir append "condatainer" to their env roots.
	return filepath.Join(scratch, "condatainer"), filepath.Join(user, "condatainer"), extraRoot, root
}

// Reads resolve nearest-first so a personal build shadows a shared one for the
// person who made it; writes go furthest-first so one copy serves the group.
// The two orders are deliberately opposite — pinning both here because nothing
// else in the tree states the contract executably.
func TestSearchAndWriteOrdersAreOpposite(t *testing.T) {
	scratch, user, extraRoot, root := withAllTiers(t)

	wantRead := []string{
		filepath.Join(scratch, "images"),
		filepath.Join(user, "images"),
		filepath.Join(extraRoot, "images"),
		filepath.Join(root, "images"),
	}
	got := searchPaths("images")
	if len(got) != len(wantRead) {
		t.Fatalf("searchPaths = %v, want %v", got, wantRead)
	}
	for i := range wantRead {
		if got[i] != wantRead[i] {
			t.Errorf("read order [%d] = %s, want %s", i, got[i], wantRead[i])
		}
	}

	wantWrite := []string{
		filepath.Join(extraRoot, "images"),
		filepath.Join(root, "images"),
		filepath.Join(scratch, "images"),
		filepath.Join(user, "images"),
	}
	dirs := imageWriteDirs()
	if len(dirs) != len(wantWrite) {
		t.Fatalf("imageWriteDirs = %v, want %v", dirs, wantWrite)
	}
	for i := range wantWrite {
		if dirs[i].Path != wantWrite[i] {
			t.Errorf("write order [%d] = %s, want %s", i, dirs[i].Path, wantWrite[i])
		}
	}

	// The tiers swap halves; order *within* a tier does not reverse. Assert the
	// swap directly so a future "just reverse the slice" refactor fails here:
	// reads must start personal and end shared, writes the other way round.
	personal := map[string]bool{
		filepath.Join(scratch, "images"): true,
		filepath.Join(user, "images"):    true,
	}
	if !personal[got[0]] || !personal[got[1]] || personal[got[2]] || personal[got[3]] {
		t.Errorf("read order does not put the personal tier first: %v", got)
	}
	if personal[dirs[0].Path] || personal[dirs[1].Path] || !personal[dirs[2].Path] || !personal[dirs[3].Path] {
		t.Errorf("write order does not put the shared tier first: %v", dirs)
	}
}

// Helper scripts use the same tier order as images.
func TestHelperScriptsShareTheImageOrder(t *testing.T) {
	scratch, user, extraRoot, root := withAllTiers(t)

	want := []string{
		filepath.Join(scratch, "helper-scripts"),
		filepath.Join(user, "helper-scripts"),
		filepath.Join(extraRoot, "helper-scripts"),
		filepath.Join(root, "helper-scripts"),
	}
	got := helperScriptSearchPaths()
	if len(got) != len(want) {
		t.Fatalf("helperScriptSearchPaths = %v, want %v", got, want)
	}
	for i := range want {
		if got[i] != want[i] {
			t.Errorf("[%d] = %s, want %s", i, got[i], want[i])
		}
	}
}

// The self-provisioned toolchain uses the same tier order as images: reads
// nearest-first, writes furthest-first, so one copy serves the whole group.
func TestLibexecSharesTheImageOrder(t *testing.T) {
	scratch, user, extraRoot, root := withAllTiers(t)

	wantRead := []string{
		filepath.Join(scratch, "libexec"),
		filepath.Join(user, "libexec"),
		filepath.Join(extraRoot, "libexec"),
		filepath.Join(root, "libexec"),
	}
	got := libexecSearchPaths()
	if len(got) != len(wantRead) {
		t.Fatalf("libexecSearchPaths = %v, want %v", got, wantRead)
	}
	for i := range wantRead {
		if got[i] != wantRead[i] {
			t.Errorf("read order [%d] = %s, want %s", i, got[i], wantRead[i])
		}
	}

	wantWrite := []string{
		filepath.Join(extraRoot, "libexec"),
		filepath.Join(root, "libexec"),
		filepath.Join(scratch, "libexec"),
		filepath.Join(user, "libexec"),
	}
	dirs := libexecWriteDirs()
	if len(dirs) != len(wantWrite) {
		t.Fatalf("libexecWriteDirs = %v, want %v", dirs, wantWrite)
	}
	for i := range wantWrite {
		if dirs[i].Path != wantWrite[i] {
			t.Errorf("write order [%d] = %s, want %s", i, dirs[i].Path, wantWrite[i])
		}
	}
}

// A directory serving two tiers (CNT_ROOT == $SCRATCH/condatainer) appears once.
// It is read at the nearer position, and still classified as the shared tier,
// which is what decides whether a write may land there.
func TestSharedAndPersonalTierCollapseToOnePath(t *testing.T) {
	base := t.TempDir()
	t.Setenv("SCRATCH", base)
	t.Setenv("XDG_DATA_HOME", filepath.Join(base, "xdg"))
	t.Setenv("CNT_EXTRA_ROOT", "")
	t.Setenv("CNT_ROOT", filepath.Join(base, "condatainer"))

	extraRootOnce, extraRootCache = sync.Once{}, ""
	rootDirOnce, rootDirCache = sync.Once{}, ""
	t.Cleanup(func() {
		extraRootOnce, extraRootCache = sync.Once{}, ""
		rootDirOnce, rootDirCache = sync.Once{}, ""
	})

	shared := filepath.Join(base, "condatainer", "images")
	got := searchPaths("images")
	if got[0] != shared {
		t.Errorf("read order [0] = %s, want the collapsed path %s", got[0], shared)
	}
	for i, p := range got[1:] {
		if p == shared {
			t.Errorf("collapsed path repeated at [%d]", i+1)
		}
	}
	if layer := ClassifyDataDir(shared); layer != LayerAppRoot {
		t.Errorf("ClassifyDataDir = %s, want %s (writes depend on the shared label)", layer, LayerAppRoot)
	}
}
