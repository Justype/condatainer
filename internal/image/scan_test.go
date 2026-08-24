package image

import (
	"os"
	"path/filepath"
	"testing"

	"github.com/Justype/condatainer/internal/config"
)

// writeImages creates dir and touches one file per name.
func writeImages(t *testing.T, dir string, names ...string) string {
	t.Helper()
	if err := os.MkdirAll(dir, 0o755); err != nil {
		t.Fatal(err)
	}
	for _, name := range names {
		if err := os.WriteFile(filepath.Join(dir, name), nil, 0o644); err != nil {
			t.Fatal(err)
		}
	}
	return dir
}

func withBase(t *testing.T, base string) {
	t.Helper()
	prev := config.Global.Base
	config.Global.Base = base
	t.Cleanup(func() { config.Global.Base = prev })
}

// A filename is the name with -- for /, and only .sqf/.img count. A .sif is the
// container root, not an overlay, so it must not appear.
func TestScanNormalizesNamesAndSkipsNonOverlays(t *testing.T) {
	dir := writeImages(t, t.TempDir(),
		"samtools--1.23.1.sqf", "grch38--genome--gencode.sqf", "env.img",
		"ubuntu24--base.sif", "notes.txt")

	scan, err := ScanOverlays(ScanOptions{Dirs: []string{dir}})
	if err != nil {
		t.Fatal(err)
	}

	want := []string{"samtools/1.23.1", "grch38/genome/gencode", "env"}
	for _, name := range want {
		if len(scan[name]) != 1 {
			t.Errorf("%q: got %v, want one copy", name, scan[name])
		}
	}
	if len(scan) != len(want) {
		t.Errorf("scan has %d names, want %d: %v", len(scan), len(want), scan)
	}
}

// Search order is priority order: the first directory holding a name owns it,
// and the shadowed copies stay reachable behind it.
func TestScanKeepsEveryCopyInSearchOrder(t *testing.T) {
	root := t.TempDir()
	high := writeImages(t, filepath.Join(root, "high"), "samtools--1.23.1.sqf")
	low := writeImages(t, filepath.Join(root, "low"), "samtools--1.23.1.sqf")

	scan, err := ScanOverlays(ScanOptions{Dirs: []string{high, low}})
	if err != nil {
		t.Fatal(err)
	}

	copies := scan["samtools/1.23.1"]
	if len(copies) != 2 {
		t.Fatalf("got %d copies, want 2: %v", len(copies), copies)
	}
	if copies[0] != filepath.Join(high, "samtools--1.23.1.sqf") {
		t.Errorf("highest-priority copy is %q, want the one in %q", copies[0], high)
	}
	if got := FirstPaths(scan)["samtools/1.23.1"]; got != copies[0] {
		t.Errorf("FirstPaths chose %q, want %q", got, copies[0])
	}
}

// A distro overlay answers to its bare name, but only when asked: a #DEP: or a
// removal names an image exactly, and must not hit the alias.
func TestScanAliasesAreOptional(t *testing.T) {
	withBase(t, "ubuntu24")
	dir := writeImages(t, t.TempDir(), "ubuntu24--build-essential.sqf")

	plain, err := ScanOverlays(ScanOptions{Dirs: []string{dir}})
	if err != nil {
		t.Fatal(err)
	}
	if _, ok := plain["build-essential"]; ok {
		t.Error("bare name resolved without Aliases")
	}

	aliased, err := ScanOverlays(ScanOptions{Dirs: []string{dir}, Aliases: true})
	if err != nil {
		t.Fatal(err)
	}
	if got, want := FirstPaths(aliased)["build-essential"], filepath.Join(dir, "ubuntu24--build-essential.sqf"); got != want {
		t.Errorf("alias resolved to %q, want %q", got, want)
	}
	if _, ok := aliased["ubuntu24/build-essential"]; !ok {
		t.Error("alias replaced the full name instead of adding to it")
	}
}

// An image named exactly like the alias wins it: the alias is a fallback, not an
// override of a real overlay.
func TestScanAliasNeverShadowsRealName(t *testing.T) {
	withBase(t, "ubuntu24")
	dir := writeImages(t, t.TempDir(), "ubuntu24--build-essential.sqf", "build-essential.sqf")

	scan, err := ScanOverlays(ScanOptions{Dirs: []string{dir}, Aliases: true})
	if err != nil {
		t.Fatal(err)
	}
	if got, want := FirstPaths(scan)["build-essential"], filepath.Join(dir, "build-essential.sqf"); got != want {
		t.Errorf("bare name resolved to %q, want the real image %q", got, want)
	}
}

// A directory that is missing is not a problem; one that cannot be read is
// reported, and the overlays found elsewhere still come back.
func TestScanReportsUnreadableDirButKeepsScanning(t *testing.T) {
	root := t.TempDir()
	good := writeImages(t, filepath.Join(root, "good"), "samtools--1.23.1.sqf")
	bad := writeImages(t, filepath.Join(root, "bad"))
	if err := os.Chmod(bad, 0o000); err != nil {
		t.Skip("cannot make a directory unreadable here")
	}
	t.Cleanup(func() { os.Chmod(bad, 0o755) }) //nolint:errcheck
	if _, err := os.ReadDir(bad); err == nil {
		t.Skip("directory is still readable (running as root?)")
	}

	scan, err := ScanOverlays(ScanOptions{Dirs: []string{bad, good, filepath.Join(root, "absent")}})
	if err == nil {
		t.Error("unreadable directory was not reported")
	}
	if len(scan["samtools/1.23.1"]) != 1 {
		t.Errorf("overlays after the unreadable directory were dropped: %v", scan)
	}
}

func TestNamesIsTheKeySet(t *testing.T) {
	names := Names(map[string][]string{"a/1": {"/x"}, "b/2": {"/y", "/z"}})
	if len(names) != 2 || !names["a/1"] || !names["b/2"] {
		t.Errorf("Names = %v, want the two keys", names)
	}
}

// A store entry is addressed by identity and reached through the store
// commands. The overlay scan is the by-name view, so it stays flat: descending
// store/ would make one name answer to several identities here.
func TestScanDoesNotDescendTheStore(t *testing.T) {
	root := t.TempDir()
	writeImages(t, root, "samtools--1.23.1.sqf")
	writeImages(t, filepath.Join(root, "store"), "samtools--1.23.1@abcdef012345.sqf")

	scan, err := ScanOverlays(ScanOptions{Dirs: []string{root}})
	if err != nil {
		t.Fatal(err)
	}
	if len(scan) != 1 {
		t.Fatalf("scan = %v, want only the flat image", scan)
	}
	copies := scan["samtools/1.23.1"]
	if len(copies) != 1 || filepath.Dir(copies[0]) != root {
		t.Fatalf("copies = %v, want only the flat one in %s", copies, root)
	}
}
