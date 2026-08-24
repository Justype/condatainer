package project

import (
	"os"
	"os/exec"
	"path/filepath"
	"testing"

	"github.com/Justype/condatainer/catalog"
	"github.com/Justype/condatainer/internal/artifact/key"
	"github.com/Justype/condatainer/internal/artifact/meta"
)

// packArtifact builds a real .sqf carrying a full /.cnt, under whatever filename
// the caller asks for. Every other test in this package injects a fake in place
// of LookupAt, which is exactly where a bug in it hides: the seam it crosses —
// a prefixed digest into a bare-hex KeyRef — fails by never matching, so a fake
// that always matches proves nothing.
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
	if err := os.MkdirAll(filepath.Dir(out), 0o775); err != nil {
		t.Fatal(err)
	}
	cmd := exec.Command("mksquashfs", root, out, "-no-progress", "-noappend", "-quiet", "-no-xattrs")
	if output, err := cmd.CombinedOutput(); err != nil {
		t.Fatalf("mksquashfs: %v\n%s", err, output)
	}
	return out, manifest
}

func TestLookupAtAcceptsTheExactIdentity(t *testing.T) {
	path, manifest := packArtifact(t, "combined.sqf", "testdata/combined/1.0", "#TYPE:app\necho one\n")

	got, ok := LookupAt(path, manifest.Name, manifest.Keys, MatchIdentity)
	if !ok {
		t.Fatal("the artifact at the declared path should satisfy its own identity")
	}
	if got.Path != path {
		t.Errorf("candidate path = %q, want %q", got.Path, path)
	}
	if got.Identity != manifest.Keys.Identity {
		t.Errorf("candidate identity = %+v, want %+v", got.Identity, manifest.Keys.Identity)
	}
	if got.Equiv != manifest.Keys.Equiv {
		t.Errorf("candidate equivalence = %+v, want %+v", got.Equiv, manifest.Keys.Equiv)
	}
}

// The filename is an address here, never a name. A project path artifact takes
// its name from the script's #TARGET:, deliberately unrelated to what the file
// is called, so decoding the filename would reject every artifact that used one.
func TestLookupAtIgnoresWhatTheFileIsCalled(t *testing.T) {
	path, manifest := packArtifact(t, "overlays.sqf", "testdata/combined/1.0", "#TYPE:app\necho one\n")

	if _, ok := LookupAt(path, manifest.Name, manifest.Keys, MatchIdentity); !ok {
		t.Fatal("a filename that does not encode the artifact name should still resolve")
	}
}

// A different build at the declared path is not the locked one. Under
// MatchIdentity that is a miss, whatever else about it agrees.
func TestLookupAtRefusesADifferentBuildUnderMatchIdentity(t *testing.T) {
	_, locked := packArtifact(t, "locked.sqf", "testdata/combined/1.0", "#TYPE:app\necho one\n")
	path, other := packArtifact(t, "combined.sqf", "testdata/combined/1.0", "#TYPE:app\necho two\n")
	if locked.Keys.Identity == other.Keys.Identity {
		t.Fatal("fixture is wrong: the two builds should differ")
	}

	if _, ok := LookupAt(path, locked.Name, locked.Keys, MatchIdentity); ok {
		t.Error("a different build should not satisfy MatchIdentity")
	}
}

// The name still has to agree: the path says which file, the manifest says what
// it is, and an artifact that is something else is not a substitute for it.
func TestLookupAtRefusesAnArtifactRecordingAnotherName(t *testing.T) {
	path, manifest := packArtifact(t, "combined.sqf", "testdata/other/1.0", "#TYPE:app\necho one\n")

	if _, ok := LookupAt(path, "testdata/combined/1.0", manifest.Keys, MatchIdentity); ok {
		t.Error("an artifact recording another name should not resolve")
	}
}

func TestLookupAtReportsAnAbsentPath(t *testing.T) {
	missing := filepath.Join(t.TempDir(), "combined.sqf")
	if _, ok := LookupAt(missing, "testdata/combined/1.0", meta.Keys{}, MatchIdentity); ok {
		t.Error("an absent path should not resolve")
	}
}

// A symlink is not a regular file. The declared path is where restore writes the
// artifact, so following a link there would read something the project does not
// own.
func TestLookupAtRefusesASymlink(t *testing.T) {
	path, manifest := packArtifact(t, "combined.sqf", "testdata/combined/1.0", "#TYPE:app\necho one\n")
	link := filepath.Join(t.TempDir(), "link.sqf")
	if err := os.Symlink(path, link); err != nil {
		t.Fatal(err)
	}

	if _, ok := LookupAt(link, manifest.Name, manifest.Keys, MatchIdentity); ok {
		t.Error("a symlink should not resolve")
	}
}
