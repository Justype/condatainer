package project

import (
	"os/exec"
	"path/filepath"
	"testing"

	"github.com/Justype/condatainer/catalog"
	"github.com/Justype/condatainer/internal/artifact/meta"
	"github.com/Justype/condatainer/internal/helperhistory"
	"github.com/Justype/condatainer/internal/project/lock"
)

func requireSquashfsTools(t *testing.T) {
	t.Helper()
	for _, bin := range []string{"mksquashfs"} {
		if _, err := exec.LookPath(bin); err != nil {
			t.Skipf("%s not available", bin)
		}
	}
}

// packEnvSqf builds a .sqf at path carrying an env-typed runtime.json, the
// way `overlay freeze` stages one — enough for container.LookupSnapshot to
// recognize it.
func packEnvSqf(t *testing.T, path string) {
	t.Helper()
	staging := t.TempDir()
	dir := filepath.Join(staging, meta.DirName)
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
	cmd := exec.Command("mksquashfs", staging, path, "-no-progress", "-noappend", "-quiet", "-no-xattrs")
	if output, err := cmd.CombinedOutput(); err != nil {
		t.Fatalf("mksquashfs: %v\n%s", err, output)
	}
}

// A frozen, shared env.sqf at a location folds into the index as an
// ordinary KindPath candidate — no different from a helper's own recorded
// overlay at that path.
func TestFrozenEnvKeyResolvesTheSharedSnapshot(t *testing.T) {
	requireSquashfsTools(t)
	root := projectRoot(t)
	packEnvSqf(t, filepath.Join(root, "env.sqf"))

	key, ok := frozenEnvKey(root, ".")
	if !ok {
		t.Fatal("frozenEnvKey did not find the shared env.sqf")
	}
	if want := lock.PathPrefix + "env.sqf"; key != want {
		t.Fatalf("key = %q, want %q", key, want)
	}
}

// A personal env-$USER.sqf never appears in the project-wide index — it has
// no audience beyond whoever created it.
func TestFrozenEnvKeyExcludesThePersonalVariant(t *testing.T) {
	requireSquashfsTools(t)
	t.Setenv("USER", "alice")
	root := projectRoot(t)
	packEnvSqf(t, filepath.Join(root, "env-alice.sqf"))

	if _, ok := frozenEnvKey(root, "."); ok {
		t.Fatal("frozenEnvKey must not surface a personal env-$USER.sqf")
	}
}

// No env overlay at all, and nothing to report.
func TestFrozenEnvKeyNoneFound(t *testing.T) {
	root := projectRoot(t)
	if _, ok := frozenEnvKey(root, "."); ok {
		t.Fatal("frozenEnvKey found something where nothing exists")
	}
}

// usageIndex folds a location's frozen env.sqf into the aggregate right
// alongside its recorded overlay, once some helper has a combination
// recorded there — env.sqf's mere presence contributes nothing on its own,
// since usageIndex only ever walks the locations recorded combinations
// name.
func TestUsageIndexIncludesFrozenEnvOnceAHelperHasRunThere(t *testing.T) {
	requireSquashfsTools(t)
	root := projectRoot(t)
	packEnvSqf(t, filepath.Join(root, "env.sqf"))
	if err := helperhistory.RecordUsed(root, "rstudio-server", ".", []string{"build-essential"}); err != nil {
		t.Fatal(err)
	}

	unpinned, err := UnpinnedHelperOverlays(root, lock.New())
	if err != nil {
		t.Fatal(err)
	}
	want := map[string]bool{"build-essential": true, lock.PathPrefix + "env.sqf": true}
	if len(unpinned) != len(want) {
		t.Fatalf("unpinned = %v, want %v", unpinned, want)
	}
	for _, key := range unpinned {
		if !want[key] {
			t.Fatalf("unexpected key %q in %v", key, unpinned)
		}
	}
}
