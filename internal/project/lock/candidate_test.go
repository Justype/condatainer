package lock

import (
	"encoding/json"
	"errors"
	"os"
	"os/exec"
	"path/filepath"
	"testing"

	"github.com/Justype/condatainer/catalog"
	"github.com/Justype/condatainer/internal/artifact/capsule"
	"github.com/Justype/condatainer/internal/artifact/key"
	"github.com/Justype/condatainer/internal/artifact/meta"
)

// packArtifact builds a real .sqf under whatever filename the caller asks for.
// The rest of this package stages an extracted /.cnt instead, which is enough
// for the vendoring half; candidateFromPath is the half that reads a file, and
// what it reads has to be a real archive.
func packArtifact(t *testing.T, filename, name, recipe string) (string, meta.Manifest) {
	t.Helper()
	for _, bin := range []string{"mksquashfs", "unsquashfs"} {
		if _, err := exec.LookPath(bin); err != nil {
			t.Skipf("%s not available", bin)
		}
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

// A path selection exists so a project can call the file whatever it likes. A
// directory scan applies the flat-name rule — the filename must encode the
// artifact name — which is right where the filename is the address and wrong
// here, and it would refuse overlays/combined.sqf for not being
// testdata--combined--1.0.sqf.
func TestCandidateFromPathIgnoresWhatTheFileIsCalled(t *testing.T) {
	path, manifest := packArtifact(t, "combined.sqf", "testdata/combined/1.0", "#TYPE:app\necho one\n")

	got, err := candidateFromPath(path)
	if err != nil {
		t.Fatalf("candidateFromPath: %v", err)
	}
	if got.Name != manifest.Name {
		t.Errorf("name = %q, want %q", got.Name, manifest.Name)
	}
	if got.Identity != manifest.Keys.Identity {
		t.Errorf("identity = %+v, want %+v", got.Identity, manifest.Keys.Identity)
	}
	if got.Equiv != manifest.Keys.Equiv {
		t.Errorf("equivalence = %+v, want %+v", got.Equiv, manifest.Keys.Equiv)
	}
}

// This and project.LookupAt read the same file by the same rule. If they
// disagree a lock can be written and then never resolve, so the two are pinned
// to the same answer here.
func TestCandidateFromPathAgreesWithTheResolverThatWillReadIt(t *testing.T) {
	path, manifest := packArtifact(t, "overlays.sqf", "testdata/combined/1.0", "#TYPE:app\necho one\n")

	got, err := candidateFromPath(path)
	if err != nil {
		t.Fatalf("candidateFromPath: %v", err)
	}
	// The lock records what this returns; the resolver later reads the same
	// file and must produce keys equal to the recorded ones.
	if got.Identity != manifest.Keys.Identity || got.Equiv != manifest.Keys.Equiv {
		t.Fatalf("selection would record keys the resolver cannot reproduce: %+v", got)
	}
	if got.Path != path {
		t.Errorf("path = %q, want the absolute selected path %q", got.Path, path)
	}
}

// A relative path is stored and resolved from the project root, so the
// candidate has to carry the absolute one.
func TestCandidateFromPathReturnsAnAbsolutePath(t *testing.T) {
	path, _ := packArtifact(t, "combined.sqf", "testdata/combined/1.0", "#TYPE:app\necho one\n")
	dir := filepath.Dir(path)
	t.Chdir(dir)

	got, err := candidateFromPath("combined.sqf")
	if err != nil {
		t.Fatalf("candidateFromPath: %v", err)
	}
	if !filepath.IsAbs(got.Path) {
		t.Errorf("path = %q, want an absolute path", got.Path)
	}
}

// A writable .img has no identity to pin and a .sif is a container root, so
// neither can be selected.
func TestCandidateFromPathRefusesWhatCannotBeSelected(t *testing.T) {
	for _, name := range []string{"env.img", "base.sif", "notes.txt"} {
		path := filepath.Join(t.TempDir(), name)
		if err := os.WriteFile(path, []byte("x"), 0o664); err != nil {
			t.Fatal(err)
		}
		if _, err := candidateFromPath(path); !errors.Is(err, ErrInvalid) {
			t.Errorf("%s: want ErrInvalid, got %v", name, err)
		}
	}
}

// An artifact with no verifiable keys is nothing to pin. That is a distinct
// answer from "the path is wrong": the file is there and says nothing.
func TestCandidateFromPathRefusesAnUnkeyedArtifact(t *testing.T) {
	for _, bin := range []string{"mksquashfs", "unsquashfs"} {
		if _, err := exec.LookPath(bin); err != nil {
			t.Skipf("%s not available", bin)
		}
	}
	root := t.TempDir()
	dir := filepath.Join(root, meta.DirName)
	if err := meta.StageRuntime(dir, meta.Runtime{
		SchemaVersion: meta.SchemaVersion,
		Name:          "testdata/combined/1.0",
		Type:          catalog.TypeApp,
		Platform:      meta.NativePlatform(),
		Prefix:        "/cnt/testdata/combined/1.0",
	}); err != nil {
		t.Fatal(err)
	}
	out := filepath.Join(t.TempDir(), "combined.sqf")
	cmd := exec.Command("mksquashfs", root, out, "-no-progress", "-noappend", "-quiet", "-no-xattrs")
	if output, err := cmd.CombinedOutput(); err != nil {
		t.Fatalf("mksquashfs: %v\n%s", err, output)
	}

	if _, err := candidateFromPath(out); !errors.Is(err, ErrNoCandidate) {
		t.Errorf("want ErrNoCandidate, got %v", err)
	}
}

func TestCandidateFromPathRefusesASymlink(t *testing.T) {
	path, _ := packArtifact(t, "combined.sqf", "testdata/combined/1.0", "#TYPE:app\necho one\n")
	link := filepath.Join(t.TempDir(), "link.sqf")
	if err := os.Symlink(path, link); err != nil {
		t.Fatal(err)
	}

	if _, err := candidateFromPath(link); !errors.Is(err, ErrInvalid) {
		t.Errorf("want ErrInvalid, got %v", err)
	}
}

// An entry directory is one format with two producers: capsule.Compose writes it
// into an image, StageEntry writes it into a checkout, and both are read by
// capsule.ReadRecord. They are pinned to each other here because a drift means a
// lock that verifies against an image that does not, or the reverse — and each
// half's own tests would stay green through it.
func TestAnEntryStagedIntoALockReadsBackAsAnImagesOwn(t *testing.T) {
	root := t.TempDir()
	manifest := meta.Manifest{
		SchemaVersion: meta.SchemaVersion,
		Name:          "zlib/1.3",
		Type:          catalog.TypeApp,
		BuildType:     "script",
		Platform:      meta.NativePlatform(),
		Source:        meta.Source{Files: []string{meta.RecipeFileName}},
	}
	recipe := []byte("#TYPE:app\necho zlib\n")
	derived, err := key.Generate(manifest, key.Sources{meta.RecipeFileName: recipe})
	if err != nil {
		t.Fatal(err)
	}
	manifest.Keys = derived.Keys()
	body, err := json.MarshalIndent(manifest, "", "  ")
	if err != nil {
		t.Fatal(err)
	}

	entryName := capsule.EntryName(manifest.Name, manifest.Keys.Identity.Digest())
	relative, err := StageEntry(root, entryName, map[string][]byte{
		meta.FileName: body, meta.RecipeFileName: recipe,
	})
	if err != nil {
		t.Fatalf("StageEntry: %v", err)
	}

	// Read it the way an image's own provenance is read: no bound, no lock
	// knowledge, just the shared reader.
	record, err := capsule.ReadRecord(filepath.Join(Dir(root), relative), capsule.PlainReader)
	if err != nil {
		t.Fatalf("a staged entry must read back as a capsule entry: %v", err)
	}
	if err := capsule.CheckDirName(filepath.Base(relative), record); err != nil {
		t.Errorf("StageEntry named the directory something ReadRecord disagrees with: %v", err)
	}
	if record.Derived.Identity.Ref != manifest.Keys.Identity {
		t.Errorf("regenerated identity %+v, want %+v", record.Derived.Identity.Ref, manifest.Keys.Identity)
	}
}

// The other direction: what vendoring copies out of an image's capsule must be
// exactly what an entry holds, so nothing is dropped on the way into the lock
// and nothing extra rides along.
func TestVendoringCarriesExactlyTheEntryFileSet(t *testing.T) {
	manifest := meta.Manifest{
		SchemaVersion: meta.SchemaVersion,
		Name:          "zlib/1.3",
		BuildType:     "script",
		Source:        meta.Source{Files: []string{meta.RecipeFileName}},
	}
	got := capsule.FileNames(manifest)
	if len(got) != 2 || got[0] != meta.FileName || got[1] != meta.RecipeFileName {
		t.Fatalf("entry file set = %v, want [%s %s]", got, meta.FileName, meta.RecipeFileName)
	}
	for _, name := range got {
		if name == meta.RuntimeFileName {
			t.Fatal("runtime.json is not part of any hashed preimage and must never be vendored")
		}
	}
}
