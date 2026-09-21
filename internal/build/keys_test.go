package build

import (
	"context"
	"os"
	"os/exec"
	"path/filepath"
	"testing"

	"github.com/Justype/condatainer/catalog"
	"github.com/Justype/condatainer/internal/artifact/key"
	"github.com/Justype/condatainer/internal/artifact/meta"
	"github.com/Justype/condatainer/internal/config"
)

// packDependency builds a real .sqf a dependency edge can be read out of. A
// locked rebuild supplies every dependency as a path, so a path is what the
// interesting cases hand to dependencyKeys.
func packDependency(t *testing.T, filename, name, recipe string) (string, meta.Manifest) {
	t.Helper()
	if _, err := exec.LookPath("mksquashfs"); err != nil {
		t.Skip("mksquashfs not available")
	}

	manifest := meta.Manifest{
		SchemaVersion: meta.SchemaVersion,
		Name:          name,
		Type:          catalog.TypeApp,
		BuildType:     "script",
		Platform:      meta.NativePlatform(),
		Source:        meta.Source{Files: []string{meta.RecipeFileName}},
	}
	derived, err := key.Generate(manifest, key.Sources{meta.RecipeFileName: []byte(recipe)})
	if err != nil {
		t.Fatal(err)
	}
	manifest.Keys = derived.Keys()

	root := t.TempDir()
	dir := filepath.Join(root, meta.DirName)
	if err := meta.StageRuntime(dir, meta.Runtime{
		SchemaVersion: meta.SchemaVersion,
		Name:          name,
		Type:          catalog.TypeApp,
		Platform:      meta.NativePlatform(),
		Prefix:        "/cnt/" + name,
	}); err != nil {
		t.Fatal(err)
	}
	if err := meta.StageManifest(dir, manifest); err != nil {
		t.Fatal(err)
	}
	if err := meta.StageBytes(dir, meta.RecipeFileName, []byte(recipe)); err != nil {
		t.Fatal(err)
	}

	out := filepath.Join(t.TempDir(), filename)
	cmd := exec.Command("mksquashfs", root, out, "-no-progress", "-noappend", "-quiet", "-no-xattrs")
	if output, err := cmd.CombinedOutput(); err != nil {
		t.Fatalf("mksquashfs: %v\n%s", err, output)
	}
	return out, manifest
}

// A path is recognized before it is parsed, never after. catalog.ParseDep
// *succeeds* on an overlay path and normalizes it into a name, turning
// .../hello--1.0.sqf into ".../hello/1.0.sqf" — which names nothing, so no
// manifest is read and the edge carries no keys at all. A locked rebuild
// supplies every dependency as a path, so that is its whole edge set.
func TestDependencyKeysReadsAPathDependencysKeys(t *testing.T) {
	path, manifest := packDependency(t, "zlib--1.3.sqf", "zlib/1.3", "#TYPE:app\necho zlib\n")

	b := &BuildObject{spec: Spec{Dependencies: []string{path}}}
	deps := b.dependencyKeys(context.Background())

	if len(deps) != 1 {
		t.Fatalf("got %d edges, want 1", len(deps))
	}
	if deps[0].Name != manifest.Name {
		t.Errorf("edge name = %q, want the name the image records, %q", deps[0].Name, manifest.Name)
	}
	if deps[0].Identity != manifest.Keys.Identity {
		t.Errorf("edge identity = %+v, want %+v", deps[0].Identity, manifest.Keys.Identity)
	}
	if deps[0].Equiv != manifest.Keys.Equiv {
		t.Errorf("edge equivalence = %+v, want %+v", deps[0].Equiv, manifest.Keys.Equiv)
	}
	if got := b.depImagePaths[manifest.Name]; got != path {
		t.Errorf("capsule would be composed from %q, want %q", got, path)
	}
}

// The path must survive as written even when nothing can be read from it. Under
// the parse-first order it came back mangled, which is what made the failure
// silent: an edge named for a path that does not exist, carrying no keys.
func TestDependencyKeysKeepsAnUnreadablePathIntact(t *testing.T) {
	path := filepath.Join(t.TempDir(), "hello--1.0.sqf")

	b := &BuildObject{spec: Spec{Dependencies: []string{path}}}
	deps := b.dependencyKeys(context.Background())

	if len(deps) != 1 {
		t.Fatalf("got %d edges, want 1", len(deps))
	}
	if deps[0].Name != path {
		t.Errorf("edge name = %q, want the path verbatim, %q", deps[0].Name, path)
	}
}

// A name dependency is still normalized: a constraint is a build-time request,
// never part of the edge an artifact records.
func TestDependencyKeysNormalizesANameDependency(t *testing.T) {
	b := &BuildObject{spec: Spec{Dependencies: []string{"samtools/1.21>=1.19"}}}
	deps := b.dependencyKeys(context.Background())

	if len(deps) != 1 {
		t.Fatalf("got %d edges, want 1", len(deps))
	}
	if deps[0].Name != "samtools/1.21" {
		t.Errorf("edge name = %q, want the resolved name/version", deps[0].Name)
	}
}

// A constraint picks the installed version the build mounts, so the edge names
// that version and carries its keys — not the preferred version, which may not be
// installed at all.
func TestDependencyKeysRecordsTheVersionAConstraintResolvesTo(t *testing.T) {
	path, manifest := packDependency(t, "samtools--1.20.sqf", "samtools/1.20", "#TYPE:app\necho samtools\n")

	root := t.TempDir()
	images := filepath.Join(root, "condatainer", "images")
	if err := os.MkdirAll(images, 0o755); err != nil {
		t.Fatal(err)
	}
	if err := os.Rename(path, filepath.Join(images, "samtools--1.20.sqf")); err != nil {
		t.Fatal(err)
	}
	t.Setenv("HOME", t.TempDir())
	t.Setenv("SCRATCH", root)
	t.Setenv("XDG_DATA_HOME", "")
	t.Setenv("CNT_EXTRA_ROOT", "")
	t.Setenv("CNT_ROOT", "")
	config.InitDataPaths()
	InvalidateInstalledOverlays()
	t.Cleanup(func() { config.InitDataPaths(); InvalidateInstalledOverlays() })

	b := &BuildObject{spec: Spec{Dependencies: []string{"samtools/1.24>=1.10"}}}
	deps := b.dependencyKeys(context.Background())

	if len(deps) != 1 {
		t.Fatalf("got %d edges, want 1", len(deps))
	}
	if deps[0].Name != "samtools/1.20" || deps[0].Identity != manifest.Keys.Identity {
		t.Errorf("edge = %s %+v, want the installed samtools/1.20 and its identity", deps[0].Name, deps[0].Identity)
	}
}

func TestDependencyKeysReturnsNothingForNoDependencies(t *testing.T) {
	b := &BuildObject{}
	if deps := b.dependencyKeys(context.Background()); deps != nil {
		t.Errorf("got %v, want nil", deps)
	}
}
