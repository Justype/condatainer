package lock

import (
	"os"
	"path/filepath"
	"strings"
	"testing"

	"github.com/Justype/condatainer/internal/artifact/meta"
)

func TestPublishWritesTheLockAtomically(t *testing.T) {
	root := projectRoot(t)
	app, appFiles := recipeArtifact(t, "star/2.7.11b", "echo star\n")
	appPath := vendor(t, root, app, appFiles)

	l := New()
	l.Pins["star/2.7.11b"] = PinEntry{Artifact: appPath}
	if err := Publish(root, l); err != nil {
		t.Fatal(err)
	}

	loaded, err := Load(root)
	if err != nil {
		t.Fatal(err)
	}
	if loaded.Pins["star/2.7.11b"].Artifact != appPath {
		t.Fatalf("loaded = %#v", loaded.Pins)
	}
	leftovers, err := os.ReadDir(Dir(root))
	if err != nil {
		t.Fatal(err)
	}
	for _, entry := range leftovers {
		if entry.Name() != ignoreFile && strings.HasPrefix(entry.Name(), stagingPrefix) {
			t.Errorf("staging file survived publication: %s", entry.Name())
		}
	}
}

// Pruning follows dependency edges, so a closure entry is kept even though no
// pin names it directly.
func TestPruneKeepsTheClosureAndRemovesOrphans(t *testing.T) {
	root := projectRoot(t)
	dep, depFiles := recipeArtifact(t, "zlib/1.3", "echo zlib\n")
	depPath := vendor(t, root, dep, depFiles)
	data, dataFiles := recipeArtifact(t, "index/1.0", "echo index\n", edge(dep))
	dataPath := vendor(t, root, data, dataFiles)
	orphan, orphanFiles := recipeArtifact(t, "cutadapt/5.0", "echo cutadapt\n")
	orphanPath := vendor(t, root, orphan, orphanFiles)

	l := New()
	l.Pins["index/1.0"] = PinEntry{Artifact: dataPath}
	if err := Publish(root, l); err != nil {
		t.Fatal(err)
	}

	if _, err := os.Stat(filepath.Join(Dir(root), depPath)); err != nil {
		t.Errorf("pruning deleted a closure entry: %v", err)
	}
	if _, err := os.Stat(filepath.Join(Dir(root), dataPath)); err != nil {
		t.Errorf("pruning deleted a selected entry: %v", err)
	}
	if _, err := os.Stat(filepath.Join(Dir(root), orphanPath)); !os.IsNotExist(err) {
		t.Errorf("an unreachable entry survived pruning: %v", err)
	}
}

// An entry that cannot be read cannot be shown to be unreachable, so deleting
// it would turn a corrupt file into data loss.
func TestPruneLeavesUnreadableEntriesAlone(t *testing.T) {
	root := projectRoot(t)
	app, appFiles := recipeArtifact(t, "star/2.7.11b", "echo star\n")
	appPath := vendor(t, root, app, appFiles)
	broken := EntryPath("broken--1.0@aaaaaaaaaaaa")
	if err := os.MkdirAll(filepath.Join(Dir(root), broken), 0o775); err != nil {
		t.Fatal(err)
	}
	if err := os.WriteFile(filepath.Join(Dir(root), broken, meta.FileName), []byte("{not json"), 0o664); err != nil {
		t.Fatal(err)
	}

	l := New()
	l.Pins["star/2.7.11b"] = PinEntry{Artifact: appPath}
	if err := Publish(root, l); err != nil {
		t.Fatal(err)
	}
	if _, err := os.Stat(filepath.Join(Dir(root), broken)); err != nil {
		t.Errorf("an unreadable entry was deleted: %v", err)
	}
}

// A lock that fails to marshal must not reach the filesystem, and must not
// trigger a prune against a document that was never published.
func TestPublishRefusesAnInvalidLockAndChangesNothing(t *testing.T) {
	root := projectRoot(t)
	app, appFiles := recipeArtifact(t, "star/2.7.11b", "echo star\n")
	appPath := vendor(t, root, app, appFiles)

	good := New()
	good.Pins["star/2.7.11b"] = PinEntry{Artifact: appPath}
	if err := Publish(root, good); err != nil {
		t.Fatal(err)
	}

	bad := New()
	bad.Pins["star/2.7.11b"] = PinEntry{Artifact: "/etc/passwd"}
	if err := Publish(root, bad); err == nil {
		t.Fatal("Publish accepted an absolute artifact path")
	}

	loaded, err := Load(root)
	if err != nil {
		t.Fatal(err)
	}
	if loaded.Pins["star/2.7.11b"].Artifact != appPath {
		t.Fatalf("a failed publish replaced the authoritative lock: %#v", loaded.Pins)
	}
	if _, err := os.Stat(filepath.Join(Dir(root), appPath)); err != nil {
		t.Errorf("a failed publish pruned artifacts: %v", err)
	}
}

// One (name, identity) has one set of records, so re-staging is a no-op rather
// than a rewrite.
func TestStageEntryIsIdempotent(t *testing.T) {
	root := projectRoot(t)
	app, appFiles := recipeArtifact(t, "star/2.7.11b", "echo star\n")
	first := vendor(t, root, app, appFiles)
	second := vendor(t, root, app, appFiles)
	if first != second {
		t.Fatalf("second stage produced %q, want %q", second, first)
	}
}

func TestStageEntryRejectsUnsafeNames(t *testing.T) {
	root := projectRoot(t)
	for _, name := range []string{"../escape", "nested/entry", ".", ".."} {
		if _, err := StageEntry(root, name, map[string][]byte{meta.FileName: []byte("{}")}); err == nil {
			t.Errorf("StageEntry accepted %q", name)
		}
	}
	if _, err := StageEntry(root, "star--2.7@aaaaaaaaaaaa", map[string][]byte{"../escape": []byte("x")}); err == nil {
		t.Error("StageEntry accepted a traversing source name")
	}
}

// A staging directory is a write in progress, never an artifact: it must not
// show up as a validation problem, and a sweep clears it.
func TestStagingDirectoriesAreInvisibleAndSweepable(t *testing.T) {
	root := projectRoot(t)
	app, appFiles := recipeArtifact(t, "star/2.7.11b", "echo star\n")
	appPath := vendor(t, root, app, appFiles)
	abandoned := filepath.Join(ProvenancePath(root), stagingName("staging"))
	if err := os.MkdirAll(abandoned, 0o775); err != nil {
		t.Fatal(err)
	}

	l := New()
	l.Pins["star/2.7.11b"] = PinEntry{Artifact: appPath}
	if _, problems := Verify(root, l); len(problems) != 0 {
		t.Fatalf("an abandoned staging directory was read as an artifact:\n%s", problemText(problems))
	}
	if err := SweepStaging(root); err != nil {
		t.Fatal(err)
	}
	if _, err := os.Stat(abandoned); !os.IsNotExist(err) {
		t.Errorf("sweep left the staging directory: %v", err)
	}
	if _, err := os.Stat(filepath.Join(Dir(root), appPath)); err != nil {
		t.Errorf("sweep removed a real artifact: %v", err)
	}
}

// A remote addresses a vendored artifact, so it must go when that artifact does
// — and it must go from the *published bytes*, not just from the in-memory lock.
// Verify calls a remote naming an absent artifact a problem, so leaving one
// behind makes the project fail its own validate after a re-pin.
func TestPublishDropsRemotesForUnreachableArtifacts(t *testing.T) {
	root := projectRoot(t)
	l := New()
	app, appFiles := recipeArtifact(t, "star/2.7.11b", "echo star\n")
	kept := vendor(t, root, app, appFiles)
	l.Pins["star/2.7.11b"] = PinEntry{Artifact: kept}

	digest := "sha256:" + strings.Repeat("a", 64)
	if err := l.AddRemote(kept, Remote{Repository: "ghcr.io/lab/p", ManifestDigest: digest}); err != nil {
		t.Fatal(err)
	}
	orphan := EntryPath("cutadapt--5.0@" + strings.Repeat("b", 12))
	if err := l.AddRemote(orphan, Remote{Repository: "ghcr.io/lab/p", ManifestDigest: digest}); err != nil {
		t.Fatal(err)
	}

	if err := Publish(root, l); err != nil {
		t.Fatalf("Publish: %v", err)
	}
	data, err := os.ReadFile(FilePath(root))
	if err != nil {
		t.Fatal(err)
	}
	published, err := Unmarshal(data)
	if err != nil {
		t.Fatalf("the published lock does not parse: %v", err)
	}
	if _, still := published.Remotes[orphan]; still {
		t.Errorf("the published lock still records a remote for %s:\n%s", orphan, data)
	}
	if len(published.Remotes[kept]) != 1 {
		t.Errorf("the reachable artifact lost its remote:\n%s", data)
	}
	if _, still := l.Remotes[orphan]; still {
		t.Error("the in-memory lock still records the orphaned remote")
	}
}

// Everything is reachable, so nothing is touched — and a lock with no remotes at
// all must not gain an empty map that would change its bytes.
func TestPublishKeepsRemotesItStillReaches(t *testing.T) {
	root := projectRoot(t)
	l := New()
	app, appFiles := recipeArtifact(t, "star/2.7.11b", "echo star\n")
	kept := vendor(t, root, app, appFiles)
	l.Pins["star/2.7.11b"] = PinEntry{Artifact: kept}
	if err := Publish(root, l); err != nil {
		t.Fatalf("Publish: %v", err)
	}
	if l.Remotes != nil {
		t.Errorf("Publish invented a remotes map: %#v", l.Remotes)
	}
}

func TestEnsureIgnoreWritesOnceAndKeepsEdits(t *testing.T) {
	root := projectRoot(t)
	path := filepath.Join(Dir(root), ignoreFile)
	if err := EnsureIgnore(root); err != nil {
		t.Fatal(err)
	}
	if data, _ := os.ReadFile(path); string(data) != "/.*/\n" {
		t.Fatalf("got %q", data)
	}
	if err := os.WriteFile(path, []byte("edited\n"), 0o644); err != nil {
		t.Fatal(err)
	}
	if err := EnsureIgnore(root); err != nil {
		t.Fatal(err)
	}
	if data, _ := os.ReadFile(path); string(data) != "edited\n" {
		t.Fatalf("an existing file was overwritten: %q", data)
	}
}
