package lock

import (
	"encoding/json"
	"os"
	"path/filepath"
	"strings"
	"testing"

	"github.com/Justype/condatainer/catalog"
	"github.com/Justype/condatainer/internal/artifact/capsule"
	"github.com/Justype/condatainer/internal/artifact/key"
	"github.com/Justype/condatainer/internal/artifact/meta"
)

// recipeArtifact builds a manifest whose keys really are what its recipe
// derives, so verification exercises regeneration rather than a stub.
//
// Only data may declare dependencies — an app is self-contained — and the key
// schemes enforce that at generation, so a fixture with edges is data.
func recipeArtifact(t *testing.T, name, recipe string, deps ...meta.Dependency) (meta.Manifest, map[string][]byte) {
	t.Helper()
	typ := catalog.TypeApp
	if len(deps) > 0 {
		typ = catalog.TypeData
	}
	manifest := meta.Manifest{
		SchemaVersion: meta.SchemaVersion,
		Name:          name,
		Type:          typ,
		BuildType:     "script",
		Platform:      meta.Platform{OS: "linux", Arch: "amd64"},
		Source:        meta.Source{Files: []string{meta.RecipeFileName}},
		Dependencies:  deps,
	}
	sources := key.Sources{meta.RecipeFileName: []byte(recipe)}
	derived, err := key.Generate(manifest, sources)
	if err != nil {
		t.Fatalf("generate keys for %s: %v", name, err)
	}
	manifest.Keys = derived.Keys()
	return manifest, map[string][]byte{meta.RecipeFileName: []byte(recipe)}
}

// vendor writes an artifact into cnt-lock/artifacts/ and returns its relative path.
func vendor(t *testing.T, root string, manifest meta.Manifest, files map[string][]byte) string {
	t.Helper()
	data, err := json.MarshalIndent(manifest, "", "  ")
	if err != nil {
		t.Fatal(err)
	}
	all := map[string][]byte{meta.FileName: data}
	for name, body := range files {
		all[name] = body
	}
	relative, err := StageArtifact(root, capsule.EntryName(manifest.Name, manifest.Keys.Identity.Digest()), all)
	if err != nil {
		t.Fatal(err)
	}
	return relative
}

func edge(manifest meta.Manifest) meta.Dependency {
	return meta.Dependency{Name: manifest.Name, Type: manifest.Type,
		Identity: manifest.Keys.Identity, Equiv: manifest.Keys.Equiv, Role: meta.RoleApp}
}

func projectRoot(t *testing.T) string {
	t.Helper()
	root := t.TempDir()
	if err := os.MkdirAll(filepath.Join(root, DirName), 0o775); err != nil {
		t.Fatal(err)
	}
	return root
}

func problemText(problems []Problem) string {
	var b strings.Builder
	for _, problem := range problems {
		b.WriteString(problem.String())
		b.WriteString("\n")
	}
	return b.String()
}

func TestVerifyAcceptsACompleteClosure(t *testing.T) {
	root := projectRoot(t)
	dep, depFiles := recipeArtifact(t, "zlib/1.3", "echo zlib\n")
	depPath := vendor(t, root, dep, depFiles)
	app, appFiles := recipeArtifact(t, "star/2.7.11b", "echo star\n", edge(dep))
	appPath := vendor(t, root, app, appFiles)

	l := New()
	l.Selections["star/2.7.11b"] = Selection{Artifact: appPath}

	verified, problems := Verify(root, l)
	if len(problems) != 0 {
		t.Fatalf("problems:\n%s", problemText(problems))
	}
	if !verified.Reachable[appPath] || !verified.Reachable[depPath] {
		t.Fatalf("closure = %#v, want both entries reachable", verified.Reachable)
	}
	if verified.Entries[appPath].Identity != app.Keys.Identity {
		t.Errorf("identity was not regenerated to the recorded value")
	}
}

// The lock is a rebuild specification, so a missing source is fatal even when
// the manifest and directory name look right.
func TestVerifyRejectsAMissingRecipe(t *testing.T) {
	root := projectRoot(t)
	app, appFiles := recipeArtifact(t, "star/2.7.11b", "echo star\n")
	appPath := vendor(t, root, app, appFiles)
	if err := os.Remove(filepath.Join(Dir(root), appPath, meta.RecipeFileName)); err != nil {
		t.Fatal(err)
	}

	l := New()
	l.Selections["star/2.7.11b"] = Selection{Artifact: appPath}
	if _, problems := Verify(root, l); len(problems) == 0 {
		t.Fatal("a vendored artifact with no recipe verified")
	}
}

// Removing any transitive dependency's directory breaks the closure, even
// though the selected artifact itself is intact.
func TestVerifyRejectsAnUnvendoredDependency(t *testing.T) {
	root := projectRoot(t)
	dep, depFiles := recipeArtifact(t, "zlib/1.3", "echo zlib\n")
	depPath := vendor(t, root, dep, depFiles)
	app, appFiles := recipeArtifact(t, "star/2.7.11b", "echo star\n", edge(dep))
	appPath := vendor(t, root, app, appFiles)
	if err := os.RemoveAll(filepath.Join(Dir(root), depPath)); err != nil {
		t.Fatal(err)
	}

	l := New()
	l.Selections["star/2.7.11b"] = Selection{Artifact: appPath}
	_, problems := Verify(root, l)
	if !strings.Contains(problemText(problems), "is not vendored") {
		t.Fatalf("problems:\n%s", problemText(problems))
	}
}

// An edge with no recorded keys cannot be followed, so a lock containing one is
// not a rebuild specification.
func TestVerifyRejectsAnUnrecordedEdge(t *testing.T) {
	root := projectRoot(t)
	app, appFiles := recipeArtifact(t, "star/2.7.11b", "echo star\n",
		meta.Dependency{Name: "mystery/1.0", Type: catalog.TypeApp, Role: meta.RoleApp, Records: "unrecorded"})
	appPath := vendor(t, root, app, appFiles)

	l := New()
	l.Selections["star/2.7.11b"] = Selection{Artifact: appPath}
	_, problems := Verify(root, l)
	if !strings.Contains(problemText(problems), "no complete keys") {
		t.Fatalf("problems:\n%s", problemText(problems))
	}
}

// The directory name is addressing, not identity: renaming it must not make a
// different artifact answer to the name.
func TestVerifyRejectsARenamedDirectory(t *testing.T) {
	root := projectRoot(t)
	app, appFiles := recipeArtifact(t, "star/2.7.11b", "echo star\n")
	appPath := vendor(t, root, app, appFiles)
	moved := ArtifactPath("star--9.9.9@aaaaaaaaaaaa")
	if err := os.Rename(filepath.Join(Dir(root), appPath), filepath.Join(Dir(root), moved)); err != nil {
		t.Fatal(err)
	}

	l := New()
	l.Selections["star/2.7.11b"] = Selection{Artifact: moved}
	_, problems := Verify(root, l)
	if !strings.Contains(problemText(problems), "directory name does not match") {
		t.Fatalf("problems:\n%s", problemText(problems))
	}
}

// Editing a vendored recipe changes what it derives, which is the whole point
// of regenerating rather than trusting the manifest.
func TestVerifyRejectsAnEditedRecipe(t *testing.T) {
	root := projectRoot(t)
	app, appFiles := recipeArtifact(t, "star/2.7.11b", "echo star\n")
	appPath := vendor(t, root, app, appFiles)
	if err := os.WriteFile(filepath.Join(Dir(root), appPath, meta.RecipeFileName), []byte("echo tampered\n"), 0o664); err != nil {
		t.Fatal(err)
	}

	l := New()
	l.Selections["star/2.7.11b"] = Selection{Artifact: appPath}
	_, problems := Verify(root, l)
	if !strings.Contains(problemText(problems), "do not regenerate") {
		t.Fatalf("problems:\n%s", problemText(problems))
	}
}

// A constrained selection key cannot come from a scan any more, but a
// hand-edited lock could still carry one, so verification refuses it rather
// than resolving a range it has no way to disambiguate.
func TestVerifyRefusesAConstrainedSelectionKey(t *testing.T) {
	root := projectRoot(t)
	app, appFiles := recipeArtifact(t, "star/2.7.11b", "echo star\n")
	appPath := vendor(t, root, app, appFiles)

	l := New()
	l.Selections["star/2.7.11b>=2.7.0"] = Selection{Artifact: appPath}
	_, problems := Verify(root, l)
	if len(problems) == 0 {
		t.Fatal("a constrained selection key was accepted")
	}
	if !strings.Contains(problemText(problems), "build recipe") {
		t.Errorf("problem does not say where a range belongs:\n%s", problemText(problems))
	}
}

func TestVerifyReportsUnreachableEntries(t *testing.T) {
	root := projectRoot(t)
	app, appFiles := recipeArtifact(t, "star/2.7.11b", "echo star\n")
	appPath := vendor(t, root, app, appFiles)
	orphan, orphanFiles := recipeArtifact(t, "cutadapt/5.0", "echo cutadapt\n")
	vendor(t, root, orphan, orphanFiles)

	l := New()
	l.Selections["star/2.7.11b"] = Selection{Artifact: appPath}
	_, problems := Verify(root, l)
	if !strings.Contains(problemText(problems), "not reachable") {
		t.Fatalf("problems:\n%s", problemText(problems))
	}
}

// A diamond stores one directory per (name, identity) and both parents reach it.
func TestVerifyHandlesADiamond(t *testing.T) {
	root := projectRoot(t)
	shared, sharedFiles := recipeArtifact(t, "zlib/1.3", "echo zlib\n")
	sharedPath := vendor(t, root, shared, sharedFiles)
	left, leftFiles := recipeArtifact(t, "star/2.7.11b", "echo star\n", edge(shared))
	leftPath := vendor(t, root, left, leftFiles)
	right, rightFiles := recipeArtifact(t, "cutadapt/5.0", "echo cutadapt\n", edge(shared))
	rightPath := vendor(t, root, right, rightFiles)

	l := New()
	l.Selections["star/2.7.11b"] = Selection{Artifact: leftPath}
	l.Selections["cutadapt/5.0"] = Selection{Artifact: rightPath}

	verified, problems := Verify(root, l)
	if len(problems) != 0 {
		t.Fatalf("problems:\n%s", problemText(problems))
	}
	if len(verified.Reachable) != 3 || !verified.Reachable[sharedPath] {
		t.Fatalf("closure = %#v", verified.Reachable)
	}
}

func TestVerifyRejectsAnOriginForAnAbsentArtifact(t *testing.T) {
	root := projectRoot(t)
	app, appFiles := recipeArtifact(t, "star/2.7.11b", "echo star\n")
	appPath := vendor(t, root, app, appFiles)

	l := New()
	l.Selections["star/2.7.11b"] = Selection{Artifact: appPath}
	if err := l.AddOrigin("artifacts/ghost--1.0@aaaaaaaaaaaa", Origin{Repository: "ghcr.io/x/y", ManifestDigest: digestA}); err != nil {
		t.Fatal(err)
	}
	_, problems := Verify(root, l)
	if !strings.Contains(problemText(problems), "origin refers to an artifact that is not vendored") {
		t.Fatalf("problems:\n%s", problemText(problems))
	}
}

// Verification is checkout-local: no images, no catalog, no config, no network.
func TestVerifyNeedsNoPayload(t *testing.T) {
	root := projectRoot(t)
	app, appFiles := recipeArtifact(t, "star/2.7.11b", "echo star\n")
	appPath := vendor(t, root, app, appFiles)

	l := New()
	l.Selections["star/2.7.11b"] = Selection{Artifact: appPath}
	t.Setenv("CNT_ROOT", filepath.Join(root, "nonexistent"))
	t.Setenv("SCRATCH", filepath.Join(root, "nonexistent"))
	if _, problems := Verify(root, l); len(problems) != 0 {
		t.Fatalf("verification touched something outside the checkout:\n%s", problemText(problems))
	}
}
