package store

import (
	"errors"
	"os"
	"path/filepath"
	"strings"
	"testing"

	"github.com/Justype/condatainer/internal/artifact/compare"
	"github.com/Justype/condatainer/internal/artifact/meta"
)

func transactionArtifact(name, identity string) compare.Artifact {
	return compare.Artifact{Name: name, IdentityScheme: "identity-v1", Identity: "sha256:" + identity,
		EquivScheme: "equiv-v1", Equiv: "sha256:" + strings.Repeat("e", 64)}
}

// occupyFlat puts an unreadable file at a flat name, which is enough to make it
// a conflict: it is not the exact identity being installed either way.
func occupyFlat(t *testing.T, root, filename string) {
	t.Helper()
	if err := os.WriteFile(filepath.Join(root, filename), []byte("incumbent"), 0o664); err != nil {
		t.Fatal(err)
	}
}

func TestTransactionPublishesPreparedSiblingAndReleasesLock(t *testing.T) {
	root := t.TempDir()
	identity := strings.Repeat("a", 64)
	read := func(path string) (compare.Artifact, error) {
		if filepath.Ext(path) == ".part" {
			return transactionArtifact("star/2.7", identity), nil
		}
		return compare.Artifact{}, errors.New("not an artifact")
	}
	tx, err := begin("star/2.7", meta.KeyRef{Scheme: "identity-v1", SHA256: identity}, BeginOptions{ImagesDir: root, SearchDirs: []string{root}}, read)
	if err != nil {
		t.Fatal(err)
	}
	if filepath.Dir(tx.Prepared) != filepath.Dir(tx.TargetPath) || tx.TargetPath+".lock" == tx.Prepared {
		t.Fatalf("prepared %q is not a private sibling of %q", tx.Prepared, tx.TargetPath)
	}
	if err := os.WriteFile(tx.Prepared, []byte("complete"), 0o664); err != nil {
		t.Fatal(err)
	}
	candidate, err := tx.Commit()
	if err != nil {
		t.Fatal(err)
	}
	if candidate.Path != tx.TargetPath {
		t.Fatalf("published path = %q", candidate.Path)
	}
	if _, err := os.Stat(tx.TargetPath + ".lock"); !os.IsNotExist(err) {
		t.Fatalf("producer lock remains: %v", err)
	}
}

func TestBeginLengthensOccupiedPrefixAndAbortCleans(t *testing.T) {
	root := t.TempDir()
	storeDir := filepath.Join(root, DirName)
	if err := os.Mkdir(storeDir, 0o775); err != nil {
		t.Fatal(err)
	}
	// The flat name must be held by something else, or this identity would take
	// it and never reach the store.
	occupyFlat(t, root, "star--2.7.sqf")
	want := "aaaaaaaaaaaa" + strings.Repeat("b", 52)
	occupied := filepath.Join(storeDir, "star--2.7@aaaaaaaaaaaa.sqf")
	if err := os.WriteFile(occupied, []byte("other"), 0o664); err != nil {
		t.Fatal(err)
	}
	read := func(path string) (compare.Artifact, error) {
		if path == occupied {
			return transactionArtifact("star/2.7", strings.Repeat("a", 64)), nil
		}
		return compare.Artifact{}, errors.New("not an artifact")
	}
	tx, err := begin("star/2.7", meta.KeyRef{Scheme: "identity-v1", SHA256: want}, BeginOptions{ImagesDir: root, SearchDirs: []string{root}}, read)
	if err != nil {
		t.Fatal(err)
	}
	if !strings.Contains(filepath.Base(tx.TargetPath), "@aaaaaaaaaaaab.sqf") {
		t.Fatalf("collision did not lengthen target: %s", tx.TargetPath)
	}
	if err := os.WriteFile(tx.Prepared, []byte("partial"), 0o664); err != nil {
		t.Fatal(err)
	}
	tx.Abort()
	if _, err := os.Stat(tx.Prepared); !os.IsNotExist(err) {
		t.Fatalf("prepared remains: %v", err)
	}
	if _, err := os.Stat(tx.TargetPath + ".lock"); !os.IsNotExist(err) {
		t.Fatalf("lock remains: %v", err)
	}
}

func TestBeginRejectsSameDigestUnderDifferentScheme(t *testing.T) {
	root := t.TempDir()
	storeDir := filepath.Join(root, DirName)
	if err := os.Mkdir(storeDir, 0o775); err != nil {
		t.Fatal(err)
	}
	occupyFlat(t, root, "star--2.7.sqf")
	identity := strings.Repeat("c", 64)
	occupied := filepath.Join(storeDir, "star--2.7@cccccccccccc.sqf")
	if err := os.WriteFile(occupied, []byte("other"), 0o664); err != nil {
		t.Fatal(err)
	}
	read := func(path string) (compare.Artifact, error) {
		if path != occupied {
			return compare.Artifact{}, errors.New("not an artifact")
		}
		artifact := transactionArtifact("star/2.7", identity)
		artifact.IdentityScheme = "other-identity-v1"
		return artifact, nil
	}
	_, err := begin("star/2.7", meta.KeyRef{Scheme: "identity-v1", SHA256: identity}, BeginOptions{ImagesDir: root, SearchDirs: []string{root}}, read)
	if !errors.Is(err, ErrIdentityCollision) {
		t.Fatalf("Begin error = %v, want ErrIdentityCollision", err)
	}
}

func TestCommitMismatchCleansPreparedAndLock(t *testing.T) {
	root := t.TempDir()
	want := strings.Repeat("d", 64)
	read := func(path string) (compare.Artifact, error) {
		if filepath.Ext(path) == ".part" {
			return transactionArtifact("wrong/1", strings.Repeat("e", 64)), nil
		}
		return compare.Artifact{}, errors.New("not an artifact")
	}
	tx, err := begin("star/2.7", meta.KeyRef{Scheme: "identity-v1", SHA256: want}, BeginOptions{ImagesDir: root, SearchDirs: []string{root}}, read)
	if err != nil {
		t.Fatal(err)
	}
	if err := os.WriteFile(tx.Prepared, []byte("mismatch"), 0o664); err != nil {
		t.Fatal(err)
	}
	if _, err := tx.Commit(); err == nil {
		t.Fatal("Commit accepted mismatched artifact")
	}
	for _, path := range []string{tx.Prepared, tx.TargetPath, tx.TargetPath + ".lock"} {
		if _, err := os.Lstat(path); !os.IsNotExist(err) {
			t.Errorf("failed commit left %s: %v", path, err)
		}
	}
}

func TestProducerContentionThenExactAdoption(t *testing.T) {
	root := t.TempDir()
	identity := strings.Repeat("f", 64)
	read := func(path string) (compare.Artifact, error) {
		if _, err := os.Lstat(path); err == nil {
			return transactionArtifact("star/2.7", identity), nil
		}
		return compare.Artifact{}, errors.New("not an artifact")
	}
	opts := BeginOptions{ImagesDir: root, SearchDirs: []string{root}}
	key := meta.KeyRef{Scheme: "identity-v1", SHA256: identity}
	first, err := begin("star/2.7", key, opts, read)
	if err != nil {
		t.Fatal(err)
	}
	if _, err := begin("star/2.7", key, opts, read); err == nil {
		t.Fatal("second producer acquired the same target")
	}
	if err := os.WriteFile(first.Prepared, []byte("complete"), 0o664); err != nil {
		t.Fatal(err)
	}
	published, err := first.Commit()
	if err != nil {
		t.Fatal(err)
	}
	second, err := begin("star/2.7", key, opts, read)
	if err != nil {
		t.Fatal(err)
	}
	if second.Adopted == nil || second.Adopted.Path != published.Path || second.Prepared != "" {
		t.Fatalf("second transaction did not adopt: %#v", second)
	}
}

func TestCommitRejectsPreparedSymlink(t *testing.T) {
	root := t.TempDir()
	identity := strings.Repeat("1", 64)
	tx, err := begin("star/2.7", meta.KeyRef{Scheme: "identity-v1", SHA256: identity}, BeginOptions{ImagesDir: root, SearchDirs: []string{root}}, func(string) (compare.Artifact, error) {
		return compare.Artifact{}, errors.New("must not read symlink")
	})
	if err != nil {
		t.Fatal(err)
	}
	target := filepath.Join(root, "payload")
	if err := os.WriteFile(target, []byte("payload"), 0o664); err != nil {
		t.Fatal(err)
	}
	if err := os.Symlink(target, tx.Prepared); err != nil {
		t.Fatal(err)
	}
	if _, err := tx.Commit(); err == nil {
		t.Fatal("Commit accepted a prepared symlink")
	}
	if _, err := os.Stat(tx.TargetPath + ".lock"); !os.IsNotExist(err) {
		t.Fatalf("producer lock remains: %v", err)
	}
}

// A free bare name is taken, not routed around: a flat artifact answers to
// `exec -o` and `list`, so the store is reserved for genuine conflicts.
func TestBeginTakesTheFreeFlatNameAndCreatesNoStore(t *testing.T) {
	root := t.TempDir()
	identity := strings.Repeat("a", 64)
	read := func(path string) (compare.Artifact, error) {
		if filepath.Ext(path) == ".part" {
			return transactionArtifact("star/2.7", identity), nil
		}
		return compare.Artifact{}, errors.New("not an artifact")
	}
	tx, err := begin("star/2.7", meta.KeyRef{Scheme: "identity-v1", SHA256: identity},
		BeginOptions{ImagesDir: root, SearchDirs: []string{root}}, read)
	if err != nil {
		t.Fatal(err)
	}
	if tx.Layout != LayoutFlat {
		t.Errorf("layout = %q, want flat", tx.Layout)
	}
	if want := filepath.Join(root, "star--2.7.sqf"); tx.TargetPath != want {
		t.Errorf("target = %q, want %q", tx.TargetPath, want)
	}
	if err := os.WriteFile(tx.Prepared, []byte("complete"), 0o664); err != nil {
		t.Fatal(err)
	}
	candidate, err := tx.Commit()
	if err != nil {
		t.Fatal(err)
	}
	if candidate.Layout != LayoutFlat || candidate.Root != root {
		t.Errorf("published %q in %q, want flat in %q", candidate.Layout, candidate.Root, root)
	}
	if _, err := os.Stat(filepath.Join(root, DirName)); !os.IsNotExist(err) {
		t.Errorf("a free name created a store directory: %v", err)
	}
}

// The name held by a different identity is the only case that overflows, and
// the incumbent is never touched.
func TestBeginOverflowsToStoreOnConflictAndKeepsIncumbent(t *testing.T) {
	root := t.TempDir()
	flat := filepath.Join(root, "star--2.7.sqf")
	occupyFlat(t, root, "star--2.7.sqf")
	identity := strings.Repeat("b", 64)
	read := func(path string) (compare.Artifact, error) {
		if flat == path {
			return transactionArtifact("star/2.7", strings.Repeat("a", 64)), nil
		}
		if filepath.Ext(path) == ".part" {
			return transactionArtifact("star/2.7", identity), nil
		}
		return compare.Artifact{}, errors.New("not an artifact")
	}
	tx, err := begin("star/2.7", meta.KeyRef{Scheme: "identity-v1", SHA256: identity},
		BeginOptions{ImagesDir: root, SearchDirs: []string{root}}, read)
	if err != nil {
		t.Fatal(err)
	}
	if tx.Layout != LayoutStored {
		t.Fatalf("layout = %q, want store", tx.Layout)
	}
	if filepath.Dir(tx.TargetPath) != filepath.Join(root, DirName) {
		t.Errorf("target %q is not in the store", tx.TargetPath)
	}
	if err := os.WriteFile(tx.Prepared, []byte("complete"), 0o664); err != nil {
		t.Fatal(err)
	}
	if _, err := tx.Commit(); err != nil {
		t.Fatal(err)
	}
	kept, err := os.ReadFile(flat)
	if err != nil || string(kept) != "incumbent" {
		t.Errorf("incumbent flat artifact was modified: %q %v", kept, err)
	}
}

// "Free" spans every readable root. A differing copy in a nearer read-only root
// would shadow a new flat install, so it forces the store too.
func TestBeginJudgesFlatOccupancyAcrossAllRoots(t *testing.T) {
	nearer, dest := t.TempDir(), t.TempDir()
	occupyFlat(t, nearer, "star--2.7.sqf")
	identity := strings.Repeat("c", 64)
	read := func(path string) (compare.Artifact, error) {
		if filepath.Ext(path) == ".part" {
			return transactionArtifact("star/2.7", identity), nil
		}
		return compare.Artifact{}, errors.New("not an artifact")
	}
	tx, err := begin("star/2.7", meta.KeyRef{Scheme: "identity-v1", SHA256: identity},
		BeginOptions{ImagesDir: dest, SearchDirs: []string{nearer, dest}}, read)
	if err != nil {
		t.Fatal(err)
	}
	if tx.Layout != LayoutStored {
		t.Fatalf("layout = %q, want store: a nearer root holds the name", tx.Layout)
	}
	if _, err := os.Stat(filepath.Join(dest, "star--2.7.sqf")); !os.IsNotExist(err) {
		t.Errorf("a shadowed flat install was reserved anyway: %v", err)
	}
}

// An exact flat copy is adopted with no write at all.
func TestBeginAdoptsAnExactFlatCopyWithoutWriting(t *testing.T) {
	root := t.TempDir()
	flat := filepath.Join(root, "star--2.7.sqf")
	occupyFlat(t, root, "star--2.7.sqf")
	identity := strings.Repeat("d", 64)
	read := func(path string) (compare.Artifact, error) {
		if path == flat {
			return transactionArtifact("star/2.7", identity), nil
		}
		return compare.Artifact{}, errors.New("not an artifact")
	}
	tx, err := begin("star/2.7", meta.KeyRef{Scheme: "identity-v1", SHA256: identity},
		BeginOptions{ImagesDir: root, SearchDirs: []string{root}}, read)
	if err != nil {
		t.Fatal(err)
	}
	if tx.Adopted == nil || tx.Adopted.Path != flat || tx.Adopted.Layout != LayoutFlat {
		t.Fatalf("did not adopt the exact flat copy: %#v", tx.Adopted)
	}
	if tx.Prepared != "" {
		t.Errorf("adoption reserved a prepared file: %q", tx.Prepared)
	}
	if _, err := os.Stat(filepath.Join(root, DirName)); !os.IsNotExist(err) {
		t.Errorf("adoption created a store directory: %v", err)
	}
}
