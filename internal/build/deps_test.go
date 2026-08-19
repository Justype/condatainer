package build

import (
	"os"
	"path/filepath"
	"slices"
	"testing"
)

// withInstalled pins the installed-overlay set for one test, so nothing depends
// on what the machine running it happens to have.
func withInstalled(t *testing.T, names ...string) {
	t.Helper()
	set := make(map[string]bool, len(names))
	for _, name := range names {
		set[name] = true
	}
	previous := cachedInstalledOverlays
	cachedInstalledOverlays = set
	t.Cleanup(func() { cachedInstalledOverlays = previous })
}

// An overlay path is satisfied by the file existing. The installed map is keyed
// by name, so a path compared against it would look permanently missing and
// send the caller off to build an artifact called "./overlays/genome.sqf".
func TestGetMissingDependenciesResolvesOverlayPathsByExistence(t *testing.T) {
	withInstalled(t, "samtools/1.21")

	dir := t.TempDir()
	present := filepath.Join(dir, "genome.sqf")
	if err := os.WriteFile(present, []byte("payload"), 0o644); err != nil {
		t.Fatalf("write fixture: %v", err)
	}
	absent := filepath.Join(dir, "index.sqf")

	b := &BuildObject{spec: Spec{Dependencies: []string{
		present,
		absent,
		"samtools/1.21",
		"cutadapt/5.0",
	}}}

	missing, err := b.GetMissingDependencies()
	if err != nil {
		t.Fatalf("GetMissingDependencies: %v", err)
	}

	if slices.Contains(missing, present) {
		t.Errorf("an existing overlay path was reported missing: %v", missing)
	}
	if !slices.Contains(missing, absent) {
		t.Errorf("an absent overlay path was not reported missing: %v", missing)
	}
	if slices.Contains(missing, "samtools/1.21") {
		t.Errorf("an installed name was reported missing: %v", missing)
	}
	if !slices.Contains(missing, "cutadapt/5.0") {
		t.Errorf("an uninstalled name was not reported missing: %v", missing)
	}
}

// A writable .img dependency follows the same rule: it is a path, not a name.
func TestGetMissingDependenciesAcceptsAnExistingImg(t *testing.T) {
	withInstalled(t)

	env := filepath.Join(t.TempDir(), "env.img")
	if err := os.WriteFile(env, []byte("payload"), 0o644); err != nil {
		t.Fatalf("write fixture: %v", err)
	}

	b := &BuildObject{spec: Spec{Dependencies: []string{env}}}
	missing, err := b.GetMissingDependencies()
	if err != nil {
		t.Fatalf("GetMissingDependencies: %v", err)
	}
	if len(missing) != 0 {
		t.Errorf("an existing .img was reported missing: %v", missing)
	}
}
