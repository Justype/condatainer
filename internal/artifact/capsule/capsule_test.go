package capsule

import (
	"errors"
	"os"
	"os/exec"
	"path/filepath"
	"slices"
	"strings"
	"testing"

	"github.com/Justype/condatainer/catalog"
	"github.com/Justype/condatainer/internal/artifact/key"
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

// depImage packs an image carrying the metadata a dependency would, plus any
// capsule entries of its own.
func depImage(t *testing.T, name string, files map[string]string, inherited map[string]string, complete bool) (string, meta.KeyRef) {
	t.Helper()
	root := t.TempDir()
	cnt := filepath.Join(root, meta.DirName)
	if err := os.MkdirAll(cnt, 0o755); err != nil {
		t.Fatal(err)
	}

	manifest := meta.Manifest{
		SchemaVersion:      meta.SchemaVersion,
		Name:               name,
		Type:               catalog.TypeData,
		BuildType:          "script",
		Platform:           meta.NativePlatform(),
		Source:             meta.Source{Files: []string{meta.RecipeFileName}},
		ProvenanceComplete: &complete,
	}
	recipe, ok := files[meta.RecipeFileName]
	if !ok {
		t.Fatal("dependency fixture has no recipe")
	}
	derived, err := key.Generate(manifest, key.Sources{meta.RecipeFileName: []byte(recipe)})
	if err != nil {
		t.Fatal(err)
	}
	manifest.Keys = derived.Keys()
	if err := meta.StageManifest(cnt, manifest); err != nil {
		t.Fatal(err)
	}
	if err := meta.StageRuntime(cnt, meta.Runtime{
		SchemaVersion: meta.SchemaVersion,
		Name:          name,
		Type:          catalog.TypeData,
		Platform:      meta.NativePlatform(),
		Prefix:        "/cnt/" + name,
	}); err != nil {
		t.Fatal(err)
	}
	for file, body := range files {
		if err := os.WriteFile(filepath.Join(cnt, file), []byte(body), 0o644); err != nil {
			t.Fatal(err)
		}
	}
	for path, body := range inherited {
		full := filepath.Join(cnt, DirName, path)
		if err := os.MkdirAll(filepath.Dir(full), 0o755); err != nil {
			t.Fatal(err)
		}
		if err := os.WriteFile(full, []byte(body), 0o644); err != nil {
			t.Fatal(err)
		}
	}

	// A payload, so the image is not metadata alone.
	payload := filepath.Join(root, "cnt", name)
	if err := os.MkdirAll(payload, 0o755); err != nil {
		t.Fatal(err)
	}
	if err := os.WriteFile(filepath.Join(payload, "data"), []byte("x"), 0o644); err != nil {
		t.Fatal(err)
	}

	out := filepath.Join(t.TempDir(), "dep.sqf")
	cmd := exec.Command("mksquashfs", root, out, "-no-progress", "-noappend", "-quiet", "-no-xattrs")
	if output, err := cmd.CombinedOutput(); err != nil {
		t.Fatalf("mksquashfs: %v\n%s", err, output)
	}
	return out, manifest.Keys.Identity
}

func TestEntryName(t *testing.T) {
	tests := []struct {
		name, identity, want string
	}{
		{"grch38/gtf-gencode/49", "sha256:41ab1c2d3e4f5a6b7c8d", "grch38--gtf-gencode--49@41ab1c2d3e4f"},
		{"star/2.7.11b", "sha256:0c19a7f34b02aaaa", "star--2.7.11b@0c19a7f34b02"},
		// A name with slashes must never become directory levels.
		{"a/b/c/d", "sha256:abcdef012345678", "a--b--c--d@abcdef012345"},
		// Short digests are taken as they are; truncation is addressing only.
		{"x/1", "sha256:abc", "x--1@abc"},
	}
	for _, tt := range tests {
		if got := EntryName(tt.name, tt.identity); got != tt.want {
			t.Errorf("EntryName(%q, %q) = %q, want %q", tt.name, tt.identity, got, tt.want)
		}
		if strings.Contains(EntryName(tt.name, tt.identity), "/") {
			t.Errorf("%q produced a nested path", tt.name)
		}
	}
}

// The capsule is a dependency's records, plus that dependency's own capsule
// copied across unchanged. Nothing is re-derived and no payload comes with it.
func TestComposeUnionsRecordsAndInheritedEntries(t *testing.T) {
	requireSquashfsTools(t)

	dep, identity := depImage(t, "grch38/gtf-gencode/49",
		map[string]string{
			meta.RecipeFileName: "#DESC:gtf\necho build\n",
		},
		map[string]string{
			"grch38--genome--gencode@aaaabbbbcccc/" + meta.FileName:       "{}\n",
			"grch38--genome--gencode@aaaabbbbcccc/" + meta.RecipeFileName: "echo genome\n",
		}, true)

	metaDir := t.TempDir()
	complete, err := Compose(metaDir, []Dep{{
		Name:      "grch38/gtf-gencode/49",
		Identity:  identity,
		ImagePath: dep,
	}})
	if err != nil {
		t.Fatalf("Compose: %v", err)
	}
	if !complete {
		t.Error("a fully recorded dependency produced an incomplete capsule")
	}

	entries, err := Entries(filepath.Join(metaDir, DirName))
	if err != nil {
		t.Fatalf("Entries: %v", err)
	}
	var dirs []string
	for _, e := range entries {
		dirs = append(dirs, e.Dir)
	}
	want := []string{"grch38--genome--gencode@aaaabbbbcccc", EntryName("grch38/gtf-gencode/49", identity.Digest())}
	if !slices.Equal(dirs, want) {
		t.Fatalf("entries = %v, want %v", dirs, want)
	}

	// The direct dependency's manifest and recipe came across.
	direct := entries[1]
	for _, file := range []string{meta.FileName, meta.RecipeFileName} {
		if !hasFile(direct.Files, file) {
			t.Errorf("%s is missing %s (has %v)", direct.Dir, file, direct.Files)
		}
	}
	// ...and its runtime did not. A capsule entry exists to rebuild an artifact,
	// never to mount one.
	if hasFile(direct.Files, meta.RuntimeFileName) {
		t.Errorf("%s copied runtime.json", direct.Dir)
	}
	// The name round-trips out of the directory name.
	if direct.Name != "grch38/gtf-gencode/49" {
		t.Errorf("name = %q", direct.Name)
	}

	// No payload rides along.
	if _, err := os.Stat(filepath.Join(metaDir, DirName, direct.Dir, "data")); err == nil {
		t.Error("a dependency payload reached the capsule")
	}
}

// A diamond stores the shared dependency once: deduplication is by directory
// name, which is (name, identity).
func TestComposeDeduplicatesADiamond(t *testing.T) {
	requireSquashfsTools(t)

	shared := map[string]string{
		"grch38--genome--gencode@aaaabbbbcccc/" + meta.FileName: "{}\n",
	}
	left, leftIdentity := depImage(t, "grch38/gtf/49", map[string]string{meta.RecipeFileName: "echo left\n"}, shared, true)
	right, rightIdentity := depImage(t, "grch38/vcf/49", map[string]string{meta.RecipeFileName: "echo right\n"}, shared, true)

	metaDir := t.TempDir()
	if _, err := Compose(metaDir, []Dep{
		{Name: "grch38/gtf/49", Identity: leftIdentity, ImagePath: left},
		{Name: "grch38/vcf/49", Identity: rightIdentity, ImagePath: right},
	}); err != nil {
		t.Fatalf("Compose: %v", err)
	}

	entries, err := Entries(filepath.Join(metaDir, DirName))
	if err != nil {
		t.Fatal(err)
	}
	if len(entries) != 3 {
		var dirs []string
		for _, e := range entries {
			dirs = append(dirs, e.Dir)
		}
		t.Errorf("entries = %v, want the shared dependency stored once", dirs)
	}
}

// An unrecorded dependency has nothing to copy and nothing to name it by, so it
// contributes no entry — and makes the closure incomplete, which is what the
// manifest says out loud.
func TestComposeWithAnUnrecordedDependency(t *testing.T) {
	metaDir := t.TempDir()
	complete, err := Compose(metaDir, []Dep{{Name: "samtools/1.23.1"}})
	if err != nil {
		t.Fatalf("Compose: %v", err)
	}
	if complete {
		t.Error("an unrecorded dependency produced a complete capsule")
	}
	if entries, err := Entries(filepath.Join(metaDir, DirName)); err != nil || len(entries) != 0 {
		t.Errorf("entries = %v, err = %v", entries, err)
	}
}

// Incompleteness is inherited one level up: if a dependency's own manifest says
// its closure is incomplete, so is everything built on it.
func TestComposeInheritsIncompleteness(t *testing.T) {
	requireSquashfsTools(t)

	out, identity := depImage(t, "grch38/gtf/49",
		map[string]string{meta.RecipeFileName: "echo gtf\n"}, nil, false)

	complete, err := Compose(t.TempDir(), []Dep{{
		Name: "grch38/gtf/49", Identity: identity, ImagePath: out,
	}})
	if err != nil {
		t.Fatalf("Compose: %v", err)
	}
	if complete {
		t.Error("a dependency with an incomplete closure did not propagate")
	}
}

func TestValidate(t *testing.T) {
	write := func(t *testing.T, dir string, files map[string]string) {
		t.Helper()
		for path, body := range files {
			full := filepath.Join(dir, path)
			if err := os.MkdirAll(filepath.Dir(full), 0o755); err != nil {
				t.Fatal(err)
			}
			if err := os.WriteFile(full, []byte(body), 0o644); err != nil {
				t.Fatal(err)
			}
		}
	}
	writeEntry := func(t *testing.T, dir, name string) string {
		t.Helper()
		recipe := []byte("echo build\n")
		manifest := meta.Manifest{
			SchemaVersion: meta.SchemaVersion,
			Name:          name, Type: catalog.TypeApp, BuildType: "script",
			Platform: meta.NativePlatform(),
			Source:   meta.Source{Files: []string{meta.RecipeFileName}},
		}
		derived, err := key.Generate(manifest, key.Sources{meta.RecipeFileName: recipe})
		if err != nil {
			t.Fatal(err)
		}
		manifest.Keys = derived.Keys()
		entry := filepath.Join(dir, EntryName(name, manifest.Keys.Identity.Digest()))
		if err := meta.StageManifest(entry, manifest); err != nil {
			t.Fatal(err)
		}
		if err := meta.StageBytes(entry, meta.RecipeFileName, recipe); err != nil {
			t.Fatal(err)
		}
		return manifest.Keys.Identity.Digest()
	}

	t.Run("a complete capsule passes", func(t *testing.T) {
		dir := t.TempDir()
		writeEntry(t, dir, "star/2.7.11b")
		if err := Validate(dir, "grch38/star/2.7.11b/gencode49", "sha256:41ab1c2d3e4f"); err != nil {
			t.Errorf("Validate: %v", err)
		}
	})

	t.Run("a self-reference is rejected", func(t *testing.T) {
		dir := t.TempDir()
		identity := writeEntry(t, dir, "grch38/star")
		err := Validate(dir, "grch38/star", identity)
		if !errors.Is(err, ErrInvalid) {
			t.Errorf("err = %v, want ErrInvalid", err)
		}
	})

	t.Run("a truncated entry is rejected", func(t *testing.T) {
		dir := t.TempDir()
		write(t, dir, map[string]string{"star--2.7.11b@0c19a7f34b02/" + meta.RecipeFileName: "echo\n"})
		err := Validate(dir, "x/1", "sha256:aaaa")
		if !errors.Is(err, ErrInvalid) {
			t.Errorf("err = %v, want ErrInvalid", err)
		}
	})

	t.Run("a malformed entry name is rejected", func(t *testing.T) {
		dir := t.TempDir()
		write(t, dir, map[string]string{"no-identity-here/" + meta.FileName: "{}\n"})
		if err := Validate(dir, "x/1", "sha256:aaaa"); !errors.Is(err, ErrInvalid) {
			t.Errorf("err = %v, want ErrInvalid", err)
		}
	})

	t.Run("an absent capsule is not an error", func(t *testing.T) {
		if err := Validate(filepath.Join(t.TempDir(), "absent"), "x/1", "sha256:aaaa"); err != nil {
			t.Errorf("Validate: %v", err)
		}
	})
}
