package lock

import (
	"encoding/json"
	"errors"
	"os"
	"path/filepath"
	"strings"
	"testing"

	"github.com/Justype/condatainer/catalog"
	"github.com/Justype/condatainer/internal/artifact/capsule"
	"github.com/Justype/condatainer/internal/artifact/meta"
	"github.com/Justype/condatainer/internal/store"
)

// fakeImage is an extracted /.cnt: what an image would yield, without needing
// squashfs tooling in a unit test.
type fakeImage struct {
	dir      string
	manifest meta.Manifest
}

// stageImage writes a /.cnt tree for one artifact plus the capsule entries it
// would carry, and returns an extractor that serves it.
func stageImage(t *testing.T, manifest meta.Manifest, files map[string][]byte, closure ...fakeImage) fakeImage {
	t.Helper()
	root := t.TempDir()
	metaDir := filepath.Join(root, meta.DirName)
	if err := os.MkdirAll(metaDir, 0o775); err != nil {
		t.Fatal(err)
	}
	data, err := json.MarshalIndent(manifest, "", "  ")
	if err != nil {
		t.Fatal(err)
	}
	writeAll(t, metaDir, map[string][]byte{meta.FileName: data})
	writeAll(t, metaDir, files)
	// A runtime document is present in a real image and must never be vendored.
	writeAll(t, metaDir, map[string][]byte{meta.RuntimeFileName: []byte(`{"schema_version":1}`)})

	if len(closure) > 0 {
		capsuleDir := filepath.Join(metaDir, capsule.DirName)
		for _, dep := range closure {
			entry := filepath.Join(capsuleDir, capsule.EntryName(dep.manifest.Name, dep.manifest.Keys.Identity.Digest()))
			if err := os.MkdirAll(entry, 0o775); err != nil {
				t.Fatal(err)
			}
			source := filepath.Join(dep.dir, meta.DirName)
			names, err := os.ReadDir(source)
			if err != nil {
				t.Fatal(err)
			}
			for _, name := range names {
				if name.IsDir() || name.Name() == meta.RuntimeFileName {
					continue
				}
				body, err := os.ReadFile(filepath.Join(source, name.Name()))
				if err != nil {
					t.Fatal(err)
				}
				writeAll(t, entry, map[string][]byte{name.Name(): body})
			}
		}
	}
	return fakeImage{dir: root, manifest: manifest}
}

func writeAll(t *testing.T, dir string, files map[string][]byte) {
	t.Helper()
	for name, body := range files {
		if err := os.WriteFile(filepath.Join(dir, name), body, 0o664); err != nil {
			t.Fatal(err)
		}
	}
}

// storeCandidate is what candidate resolution would have produced for a staged
// image, so vendoring can be exercised without squashfs tooling.
func storeCandidate(manifest meta.Manifest) store.Candidate {
	return store.Candidate{Name: manifest.Name, Path: "/fake/" + manifest.Name + ".sqf",
		Layout: store.LayoutFlat, Identity: manifest.Keys.Identity, Equiv: manifest.Keys.Equiv}
}

// vendorOnly drives the capture half directly. Candidate resolution needs a real
// verifiable .sqf, which a unit test has no way to build.
func vendorOnly(t *testing.T, root string, image fakeImage) ([]string, error) {
	t.Helper()
	staging := t.TempDir()
	if err := os.CopyFS(staging, os.DirFS(image.dir)); err != nil {
		t.Fatal(err)
	}
	return vendorClosure(root, filepath.Join(staging, meta.DirName), storeCandidate(image.manifest))
}

func TestVendorClosureCapturesTheArtifactAndItsCapsule(t *testing.T) {
	root := projectRoot(t)
	dep, depFiles := recipeArtifact(t, "zlib/1.3", "echo zlib\n")
	depImage := stageImage(t, dep, depFiles)
	data, dataFiles := recipeArtifact(t, "index/1.0", "echo index\n", edge(dep))
	dataImage := stageImage(t, data, dataFiles, depImage)

	vendored, err := vendorOnly(t, root, dataImage)
	if err != nil {
		t.Fatal(err)
	}
	if len(vendored) != 2 {
		t.Fatalf("vendored = %v, want the artifact and its dependency", vendored)
	}

	l := New()
	l.Pins["index/1.0"] = PinEntry{Artifact: EntryPath(capsule.EntryName(data.Name, data.Keys.Identity.Digest()))}
	if _, problems := Verify(root, l); len(problems) != 0 {
		t.Fatalf("a freshly vendored closure does not verify:\n%s", problemText(problems))
	}
}

// runtime.json is not part of any hashed preimage, so a lock carrying it would
// be carrying something it cannot verify.
func TestVendorClosureNeverCapturesRuntimeJSON(t *testing.T) {
	root := projectRoot(t)
	app, appFiles := recipeArtifact(t, "star/2.7.11b", "echo star\n")
	image := stageImage(t, app, appFiles)

	vendored, err := vendorOnly(t, root, image)
	if err != nil {
		t.Fatal(err)
	}
	entries, err := os.ReadDir(filepath.Join(Dir(root), vendored[0]))
	if err != nil {
		t.Fatal(err)
	}
	for _, entry := range entries {
		if entry.Name() == meta.RuntimeFileName {
			t.Fatalf("runtime.json was vendored into %s", vendored[0])
		}
	}
}

// An artifact built over an unrecorded dependency cannot be rebuilt from the
// checkout, whatever its own sources say.
func TestVendorClosureRefusesIncompleteProvenance(t *testing.T) {
	root := projectRoot(t)
	incomplete := false
	data, dataFiles := recipeArtifact(t, "index/1.0", "echo index\n",
		meta.Dependency{Name: "mystery/1.0", Type: catalog.TypeApp, Role: meta.RoleApp, Records: "unrecorded"})
	data.ProvenanceComplete = &incomplete
	image := stageImage(t, data, dataFiles)

	if _, err := vendorOnly(t, root, image); !errors.Is(err, ErrInvalid) {
		t.Fatalf("error = %v, want ErrInvalid", err)
	}
}

func TestVendorClosureRefusesAKeylessArtifact(t *testing.T) {
	root := projectRoot(t)
	app, appFiles := recipeArtifact(t, "star/2.7.11b", "echo star\n")
	app.Keys = meta.Keys{}
	image := stageImage(t, app, appFiles)

	if _, err := vendorOnly(t, root, image); !errors.Is(err, ErrInvalid) {
		t.Fatalf("error = %v, want ErrInvalid", err)
	}
}

// A diamond stores one directory per (name, identity), so re-selecting a second
// parent confirms the shared entry rather than duplicating it.
func TestVendorClosureDeduplicatesADiamond(t *testing.T) {
	root := projectRoot(t)
	shared, sharedFiles := recipeArtifact(t, "zlib/1.3", "echo zlib\n")
	sharedImage := stageImage(t, shared, sharedFiles)
	left, leftFiles := recipeArtifact(t, "index/1.0", "echo index\n", edge(shared))
	right, rightFiles := recipeArtifact(t, "other/1.0", "echo other\n", edge(shared))

	if _, err := vendorOnly(t, root, stageImage(t, left, leftFiles, sharedImage)); err != nil {
		t.Fatal(err)
	}
	if _, err := vendorOnly(t, root, stageImage(t, right, rightFiles, sharedImage)); err != nil {
		t.Fatal(err)
	}

	entries, err := os.ReadDir(ProvenancePath(root))
	if err != nil {
		t.Fatal(err)
	}
	if len(entries) != 3 {
		names := make([]string, 0, len(entries))
		for _, entry := range entries {
			names = append(names, entry.Name())
		}
		t.Fatalf("artifacts = %v, want three distinct entries", names)
	}
}

// Apply verifies before publishing, so a lock is never left pointing at
// something that does not validate.
func TestApplyRestoresTheLockWhenVerificationFails(t *testing.T) {
	root := projectRoot(t)
	app, appFiles := recipeArtifact(t, "star/2.7.11b", "echo star\n")
	appPath := vendor(t, root, app, appFiles)

	l := New()
	l.Pins["star/2.7.11b"] = PinEntry{Artifact: appPath}
	if err := Publish(root, l); err != nil {
		t.Fatal(err)
	}

	err := Apply(root, l, &Pinned{Request: "ghost/1.0", Artifact: EntryPath("ghost--1.0@aaaaaaaaaaaa")})
	if !errors.Is(err, ErrInvalid) {
		t.Fatalf("error = %v, want ErrInvalid", err)
	}
	if _, present := l.Pins["ghost/1.0"]; present {
		t.Error("a failed Apply left its pin behind")
	}
	loaded, err := Load(root)
	if err != nil {
		t.Fatal(err)
	}
	if len(loaded.Pins) != 1 || loaded.Pins["star/2.7.11b"].Artifact != appPath {
		t.Fatalf("a failed Apply changed the published lock: %#v", loaded.Pins)
	}
}

func TestApplyPublishesAValidPin(t *testing.T) {
	root := projectRoot(t)
	app, appFiles := recipeArtifact(t, "star/2.7.11b", "echo star\n")
	appPath := vendor(t, root, app, appFiles)

	l := New()
	if err := Apply(root, l, &Pinned{Request: "star/2.7.11b", Artifact: appPath}); err != nil {
		t.Fatal(err)
	}
	loaded, err := Load(root)
	if err != nil {
		t.Fatal(err)
	}
	if loaded.Pins["star/2.7.11b"].Artifact != appPath {
		t.Fatalf("loaded = %#v", loaded.Pins)
	}
}

// Reconcile drops what nothing requests and reports what has no pin. It
// never invents one: choosing an artifact is an explicit act.
func TestReconcileDropsStaleAndReportsWhatNeedsPinning(t *testing.T) {
	root := projectRoot(t)
	app, appFiles := recipeArtifact(t, "star/2.7.11b", "echo star\n")
	appPath := vendor(t, root, app, appFiles)

	l := New()
	l.Pins["star/2.7.11b"] = PinEntry{Artifact: appPath}
	l.Pins["gone/1.0"] = PinEntry{Artifact: appPath}

	write(t, root, "run.sh", "#DEP: star/2.7.11b\n#DEP: cutadapt/5.0\n#DEP: env.img\nrun\n")
	result := scan(t, root)

	needPin := Reconcile(root, l, result)
	if _, stale := l.Pins["gone/1.0"]; stale {
		t.Error("a pin nothing requests survived")
	}
	if len(needPin) != 1 || needPin[0].Key != "cutadapt/5.0" {
		keys := make([]string, 0, len(needPin))
		for _, request := range needPin {
			keys = append(keys, request.Key)
		}
		t.Fatalf("needPin = %v, want only cutadapt/5.0", keys)
	}
}

// A pin whose artifact stopped verifying is dropped and re-reported, so a
// lock never keeps pointing at something broken.
func TestReconcileDropsAnInvalidPin(t *testing.T) {
	root := projectRoot(t)
	app, appFiles := recipeArtifact(t, "star/2.7.11b", "echo star\n")
	appPath := vendor(t, root, app, appFiles)
	if err := os.Remove(filepath.Join(Dir(root), appPath, meta.RecipeFileName)); err != nil {
		t.Fatal(err)
	}

	l := New()
	l.Pins["star/2.7.11b"] = PinEntry{Artifact: appPath}
	write(t, root, "run.sh", "#DEP: star/2.7.11b\nrun\n")

	needPin := Reconcile(root, l, scan(t, root))
	if _, kept := l.Pins["star/2.7.11b"]; kept {
		t.Error("an invalid pin was kept")
	}
	if len(needPin) != 1 {
		t.Fatalf("needPin = %#v", needPin)
	}
}

// A pin names an identity, so a file is refused whatever its extension.
// Nothing downstream could use one: restore and push both find a selected
// artifact through the image roots and the store, so a file elsewhere would
// record an identity only a rebuild could satisfy.
func TestPinRefusesAFileAsTheIdentity(t *testing.T) {
	root := projectRoot(t)
	for _, target := range []string{
		"/shared/overlays/star.sqf", "./star.sqf", "/tmp/env.img", "/tmp/base.sif", "base.sif",
	} {
		_, err := Pin(root, "star/2.7.11b", target, PinOptions{})
		if !errors.Is(err, ErrInvalid) || !strings.Contains(err.Error(), "names an identity") {
			t.Errorf("Pin(%s) error = %v", target, err)
		}
	}
}

func TestPinRefusesAnEmptyRequest(t *testing.T) {
	if _, err := Pin(projectRoot(t), "  ", "abcdef", PinOptions{}); !errors.Is(err, ErrInvalid) {
		t.Fatalf("error = %v, want ErrInvalid", err)
	}
}

// Selecting an unpinnable request is refused with the marker as the remedy,
// rather than failing later inside Validate with a path complaint.
func TestPinRefusesUnpinnableRequests(t *testing.T) {
	root := t.TempDir()
	for _, request := range []string{
		PathPrefix + "/shared/lab/genome.sqf",
		PathPrefix + "../outside/tool.sqf",
		PathPrefix + "env.img",
	} {
		_, err := Pin(root, request, "sha256:"+strings.Repeat("a", 64), PinOptions{})
		if !errors.Is(err, ErrInvalid) {
			t.Errorf("Pin(%q) = %v, want ErrInvalid", request, err)
			continue
		}
		// The refusal has to say what closes the gap, not only that it refused.
		if !strings.Contains(err.Error(), "instead") && !strings.Contains(err.Error(), "pin that") {
			t.Errorf("Pin(%q) error does not say what to do instead: %v", request, err)
		}
	}
}

// MatchPin tries the literal key first, so a pin already recorded there keeps
// meaning exactly that.
func TestMatchPinPrefersTheLiteralKeyWhenBothArePinned(t *testing.T) {
	l := New()
	l.Pins[PathPrefix+"xxx.sqf"] = PinEntry{Artifact: "provenance/root--1@abc"}
	l.Pins[PathPrefix+"steps1/xxx.sqf"] = PinEntry{Artifact: "provenance/nested--1@def"}
	request := Request{Key: PathPrefix + "xxx.sqf", Kind: KindPath, Path: "xxx.sqf", Scripts: []string{"steps1/run.sh"}}

	key, entry, ok := MatchPin(l, request)
	if !ok || key != PathPrefix+"xxx.sqf" || entry.Artifact != "provenance/root--1@abc" {
		t.Errorf("MatchPin = (%q, %#v, %v), want the literal key preferred", key, entry, ok)
	}
}

// A pin recorded only under a declaring script's own directory still answers
// the literal, root-relative declaration — the fallback MatchPin and
// pinFirstCandidate share.
func TestMatchPinFallsBackToADeclaringScriptsDirectory(t *testing.T) {
	l := New()
	l.Pins[PathPrefix+"steps1/xxx.sqf"] = PinEntry{Artifact: "provenance/tool--1@abc"}
	request := Request{Key: PathPrefix + "xxx.sqf", Kind: KindPath, Path: "xxx.sqf", Scripts: []string{"steps1/run.sh"}}

	key, entry, ok := MatchPin(l, request)
	if !ok || key != PathPrefix+"steps1/xxx.sqf" || entry.Artifact != "provenance/tool--1@abc" {
		t.Errorf("MatchPin = (%q, %#v, %v), want the fallback candidate", key, entry, ok)
	}
}

// Nothing pinned under any candidate is an ordinary miss, not a crash.
func TestMatchPinFalseWhenNothingAnswers(t *testing.T) {
	l := New()
	request := Request{Key: PathPrefix + "xxx.sqf", Kind: KindPath, Path: "xxx.sqf", Scripts: []string{"steps1/run.sh"}}

	if _, _, ok := MatchPin(l, request); ok {
		t.Fatal("MatchPin found a pin that was never set")
	}
}

// A KindName request has no ambiguity to offer: MatchPin degenerates to an
// ordinary single-key lookup.
func TestMatchPinHasOneCandidateForAName(t *testing.T) {
	l := New()
	l.Pins["star/2.7.11b"] = PinEntry{Artifact: "provenance/star--2.7.11b@abc"}
	request := Request{Key: "star/2.7.11b", Kind: KindName}

	key, _, ok := MatchPin(l, request)
	if !ok || key != "star/2.7.11b" {
		t.Errorf("MatchPin = (%q, ok=%v), want the name's own key", key, ok)
	}
}

// A loose name request (no pin under its own key) matches an existing pin
// whose name agrees and whose version the request's Dep admits — the newest
// one, when more than one qualifies — so a second script's looser reference
// to an already-pinned floating dependency reuses that pin instead of
// reporting unpinned.
func TestMatchPinLooseNameMatchesExistingSatisfyingPin(t *testing.T) {
	l := New()
	l.Pins["openjdk/17.0.15"] = PinEntry{Artifact: "provenance/openjdk--17.0.15@abc"}
	l.Pins["openjdk/17.0.18"] = PinEntry{Artifact: "provenance/openjdk--17.0.18@def"}
	l.Pins["openjdk/16.0.2"] = PinEntry{Artifact: "provenance/openjdk--16.0.2@ghi"}
	dep, err := catalog.ParseDep("openjdk/17")
	if err != nil {
		t.Fatal(err)
	}
	request := Request{Key: "openjdk/17", Kind: KindName, Dep: dep}

	key, entry, ok := MatchPin(l, request)
	if !ok || key != "openjdk/17.0.18" || entry.Artifact != "provenance/openjdk--17.0.18@def" {
		t.Errorf("MatchPin = (%q, %#v, %v), want the newest satisfying pin openjdk/17.0.18", key, entry, ok)
	}
}

// An exact name request with nothing pinned under its own key gains nothing
// from a same-named pin at a different, non-admitted version.
func TestMatchPinExactNameDoesNotMatchADifferentVersion(t *testing.T) {
	l := New()
	l.Pins["star/2.7.10"] = PinEntry{Artifact: "provenance/star--2.7.10@abc"}
	dep, err := catalog.ParseDep("star/2.7.11b")
	if err != nil {
		t.Fatal(err)
	}
	request := Request{Key: "star/2.7.11b", Kind: KindName, Dep: dep}

	if _, _, ok := MatchPin(l, request); ok {
		t.Error("MatchPin should not match an exact request against an unrelated version")
	}
}

// pinFirstCandidate resolves a script-relative file when the root-relative
// one does not exist, and records the resolved key, not the one as scanned.
func TestPinFirstCandidateFallsBackToADeclaringScriptsDirectory(t *testing.T) {
	root := t.TempDir()
	built, manifest := packArtifact(t, "xxx.sqf", "testdata/tool", "recipe")
	dest := filepath.Join(root, "steps1", "xxx.sqf")
	if err := os.MkdirAll(filepath.Dir(dest), 0o775); err != nil {
		t.Fatal(err)
	}
	data, err := os.ReadFile(built)
	if err != nil {
		t.Fatal(err)
	}
	if err := os.WriteFile(dest, data, 0o664); err != nil {
		t.Fatal(err)
	}

	request := Request{Key: PathPrefix + "xxx.sqf", Kind: KindPath, Path: "xxx.sqf", Scripts: []string{"steps1/run.sh"}}
	pinned, err := pinFirstCandidate(root, request, PinOptions{})
	if err != nil {
		t.Fatal(err)
	}
	if pinned.Request != PathPrefix+"steps1/xxx.sqf" {
		t.Errorf("Request = %q, want the resolved script-relative key", pinned.Request)
	}
	if pinned.Name != manifest.Name {
		t.Errorf("Name = %q, want %q", pinned.Name, manifest.Name)
	}
}

// The root-relative candidate wins when a real file answers there too.
func TestPinFirstCandidatePrefersTheRootRelativeFile(t *testing.T) {
	root := t.TempDir()
	rootBuilt, rootManifest := packArtifact(t, "xxx.sqf", "testdata/root-tool", "root-recipe")
	nestedBuilt, _ := packArtifact(t, "xxx.sqf", "testdata/nested-tool", "nested-recipe")

	for _, pair := range []struct{ built, dest string }{
		{rootBuilt, filepath.Join(root, "xxx.sqf")},
		{nestedBuilt, filepath.Join(root, "steps1", "xxx.sqf")},
	} {
		if err := os.MkdirAll(filepath.Dir(pair.dest), 0o775); err != nil {
			t.Fatal(err)
		}
		data, err := os.ReadFile(pair.built)
		if err != nil {
			t.Fatal(err)
		}
		if err := os.WriteFile(pair.dest, data, 0o664); err != nil {
			t.Fatal(err)
		}
	}

	request := Request{Key: PathPrefix + "xxx.sqf", Kind: KindPath, Path: "xxx.sqf", Scripts: []string{"steps1/run.sh"}}
	pinned, err := pinFirstCandidate(root, request, PinOptions{})
	if err != nil {
		t.Fatal(err)
	}
	if pinned.Request != PathPrefix+"xxx.sqf" {
		t.Errorf("Request = %q, want the literal root-relative key preferred", pinned.Request)
	}
	if pinned.Name != rootManifest.Name {
		t.Errorf("Name = %q, want the root file's %q", pinned.Name, rootManifest.Name)
	}
}
