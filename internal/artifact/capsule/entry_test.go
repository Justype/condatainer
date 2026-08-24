package capsule

import (
	"errors"
	"os"
	"path/filepath"
	"strings"
	"testing"

	"github.com/Justype/condatainer/catalog"
	"github.com/Justype/condatainer/internal/artifact/key"
	"github.com/Justype/condatainer/internal/artifact/meta"
)

// stageEntry lays down one entry directory the way both a capsule and a lock
// hold it, and returns its correct directory name.
func stageEntry(t *testing.T, parent, name, recipe string) string {
	t.Helper()
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

	dir := EntryName(name, manifest.Keys.Identity.Digest())
	if err := meta.StageManifest(filepath.Join(parent, dir), manifest); err != nil {
		t.Fatal(err)
	}
	if err := meta.StageBytes(filepath.Join(parent, dir), meta.RecipeFileName, []byte(recipe)); err != nil {
		t.Fatal(err)
	}
	return dir
}

func TestReadRecordRegeneratesBothKeys(t *testing.T) {
	root := t.TempDir()
	dir := stageEntry(t, root, "samtools/1.23.1", "#TYPE:app\necho hi\n")

	record, err := ReadRecord(filepath.Join(root, dir), PlainReader)
	if err != nil {
		t.Fatalf("ReadRecord: %v", err)
	}
	if record.Manifest.Name != "samtools/1.23.1" {
		t.Fatalf("name is %q", record.Manifest.Name)
	}
	if record.Derived.Identity.Ref.Empty() || record.Derived.Equiv.Ref.Empty() {
		t.Fatal("both keys should be regenerated from the sources")
	}
	if err := CheckDirName(dir, record); err != nil {
		t.Fatalf("CheckDirName on the name ReadRecord's own keys produce: %v", err)
	}
}

func TestReadRecordRefusesEditedSources(t *testing.T) {
	root := t.TempDir()
	dir := stageEntry(t, root, "samtools/1.23.1", "#TYPE:app\necho hi\n")

	recipe := filepath.Join(root, dir, meta.RecipeFileName)
	if err := os.WriteFile(recipe, []byte("#TYPE:app\necho tampered\n"), 0o644); err != nil {
		t.Fatal(err)
	}
	_, err := ReadRecord(filepath.Join(root, dir), PlainReader)
	if err == nil || !strings.Contains(err.Error(), "regenerate") {
		t.Fatalf("edited sources should fail regeneration, got %v", err)
	}
}

// The reader is the caller's policy, and it is the only thing that differs
// between an image's entry and a checkout's. A bounded reader must therefore be
// able to refuse bytes without the shared code having to know why.
func TestReadRecordSurfacesTheCallersReaderError(t *testing.T) {
	root := t.TempDir()
	dir := stageEntry(t, root, "samtools/1.23.1", "#TYPE:app\necho hi\n")

	tooBig := errors.New("over the limit")
	bounded := func(path string) ([]byte, error) {
		if filepath.Base(path) == meta.RecipeFileName {
			return nil, tooBig
		}
		return os.ReadFile(path)
	}
	_, err := ReadRecord(filepath.Join(root, dir), bounded)
	if !errors.Is(err, tooBig) {
		t.Fatalf("want the reader's own error, got %v", err)
	}
}

func TestFileNamesIsTheManifestPlusItsSources(t *testing.T) {
	m := meta.Manifest{Source: meta.Source{Files: []string{"a", "b"}}}
	got := FileNames(m)
	want := []string{meta.FileName, "a", "b"}
	if len(got) != len(want) {
		t.Fatalf("got %v, want %v", got, want)
	}
	for i := range want {
		if got[i] != want[i] {
			t.Fatalf("got %v, want %v", got, want)
		}
	}
}
