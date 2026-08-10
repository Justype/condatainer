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
func depImage(t *testing.T, name string, files map[string]string, inherited map[string]string) string {
	t.Helper()
	root := t.TempDir()
	cnt := filepath.Join(root, meta.DirName)
	if err := os.MkdirAll(cnt, 0o755); err != nil {
		t.Fatal(err)
	}

	complete := true
	if err := meta.StageManifest(cnt, meta.Manifest{
		SchemaVersion:      meta.SchemaVersion,
		Name:               name,
		Type:               catalog.TypeData,
		BuildType:          "script",
		Platform:           meta.NativePlatform(),
		ProvenanceComplete: &complete,
	}); err != nil {
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
	return out
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

	dep := depImage(t, "grch38/gtf-gencode/49",
		map[string]string{
			meta.IdentityFileName: "cnt-identity-v1\ntype=data\n",
			meta.EquivFileName:    "cnt-equiv-v1\ntype=data\n",
			meta.RecipeFileName:   "#DESC:gtf\necho build\n",
		},
		map[string]string{
			"grch38--genome--gencode@aaaabbbbcccc/" + meta.FileName:       "{}\n",
			"grch38--genome--gencode@aaaabbbbcccc/" + meta.RecipeFileName: "echo genome\n",
		})

	metaDir := t.TempDir()
	complete, err := Compose(metaDir, []Dep{{
		Name:      "grch38/gtf-gencode/49",
		Identity:  "sha256:41ab1c2d3e4f5a6b",
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
	want := []string{"grch38--genome--gencode@aaaabbbbcccc", "grch38--gtf-gencode--49@41ab1c2d3e4f"}
	if !slices.Equal(dirs, want) {
		t.Fatalf("entries = %v, want %v", dirs, want)
	}

	// The direct dependency's own records came across...
	direct := entries[1]
	for _, file := range []string{meta.FileName, meta.IdentityFileName, meta.EquivFileName, meta.RecipeFileName} {
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
	left := depImage(t, "grch38/gtf/49", map[string]string{meta.RecipeFileName: "echo left\n"}, shared)
	right := depImage(t, "grch38/vcf/49", map[string]string{meta.RecipeFileName: "echo right\n"}, shared)

	metaDir := t.TempDir()
	if _, err := Compose(metaDir, []Dep{
		{Name: "grch38/gtf/49", Identity: "sha256:1111222233334444", ImagePath: left},
		{Name: "grch38/vcf/49", Identity: "sha256:5555666677778888", ImagePath: right},
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

	root := t.TempDir()
	cnt := filepath.Join(root, meta.DirName)
	if err := os.MkdirAll(cnt, 0o755); err != nil {
		t.Fatal(err)
	}
	incomplete := false
	if err := meta.StageManifest(cnt, meta.Manifest{
		SchemaVersion:      meta.SchemaVersion,
		Name:               "grch38/gtf/49",
		Type:               catalog.TypeData,
		BuildType:          "script",
		Platform:           meta.NativePlatform(),
		ProvenanceComplete: &incomplete,
	}); err != nil {
		t.Fatal(err)
	}
	out := filepath.Join(t.TempDir(), "dep.sqf")
	if output, err := exec.Command("mksquashfs", root, out, "-no-progress", "-noappend", "-quiet", "-no-xattrs").CombinedOutput(); err != nil {
		t.Fatalf("mksquashfs: %v\n%s", err, output)
	}

	complete, err := Compose(t.TempDir(), []Dep{{
		Name: "grch38/gtf/49", Identity: "sha256:1111222233334444", ImagePath: out,
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

	t.Run("a complete capsule passes", func(t *testing.T) {
		dir := t.TempDir()
		write(t, dir, map[string]string{"star--2.7.11b@0c19a7f34b02/" + meta.FileName: "{}\n"})
		if err := Validate(dir, "grch38/star/2.7.11b/gencode49", "sha256:41ab1c2d3e4f"); err != nil {
			t.Errorf("Validate: %v", err)
		}
	})

	t.Run("a self-reference is rejected", func(t *testing.T) {
		dir := t.TempDir()
		write(t, dir, map[string]string{"grch38--star@41ab1c2d3e4f/" + meta.FileName: "{}\n"})
		err := Validate(dir, "grch38/star", "sha256:41ab1c2d3e4f5a6b")
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
