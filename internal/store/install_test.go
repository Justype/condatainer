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
