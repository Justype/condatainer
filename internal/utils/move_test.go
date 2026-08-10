package utils

import (
	"context"
	"os"
	"path/filepath"
	"testing"
)

// Which branch MoveFile takes must not decide permissions. Only the rename branch
// is reachable here — a temp dir is one filesystem — and it is the branch that
// used to skip the fix, since a rename carries the source mode across.
func TestMoveFileSharesWithAGroupWritableParent(t *testing.T) {
	root := t.TempDir()
	srcDir := filepath.Join(root, "tmp")
	dstDir := filepath.Join(root, "images")
	for _, d := range []string{srcDir, dstDir} {
		if err := os.Mkdir(d, 0o755); err != nil {
			t.Fatal(err)
		}
	}
	// The shared-install case: a group-writable destination.
	if err := os.Chmod(dstDir, 0o775); err != nil {
		t.Fatal(err)
	}

	src := filepath.Join(srcDir, "overlay.img")
	if err := os.WriteFile(src, []byte("payload"), 0o600); err != nil {
		t.Fatal(err)
	}
	// WriteFile is umask-subject, so say what the test means.
	if err := os.Chmod(src, 0o600); err != nil {
		t.Fatal(err)
	}

	dst := filepath.Join(dstDir, "overlay.img")
	copied, err := MoveFile(context.Background(), src, dst, false)
	if err != nil {
		t.Fatalf("MoveFile: %v", err)
	}
	if copied {
		t.Error("a same-filesystem move copied instead of renaming")
	}

	fi, err := os.Stat(dst)
	if err != nil {
		t.Fatalf("the destination is missing: %v", err)
	}
	if mode := fi.Mode().Perm(); mode&0o060 != 0o060 {
		t.Errorf("mode = %#o, want group read-write inside a group-writable parent", mode)
	}
	if _, err := os.Stat(src); !os.IsNotExist(err) {
		t.Error("the source survived the move")
	}
	if data, err := os.ReadFile(dst); err != nil || string(data) != "payload" {
		t.Errorf("contents = %q, %v", data, err)
	}
}

// A personal install is umask-default and must stay that way.
func TestMoveFileLeavesAPrivateParentAlone(t *testing.T) {
	root := t.TempDir()
	dstDir := filepath.Join(root, "images")
	if err := os.Mkdir(dstDir, 0o755); err != nil {
		t.Fatal(err)
	}
	if err := os.Chmod(dstDir, 0o755); err != nil { // not group-writable
		t.Fatal(err)
	}

	src := filepath.Join(root, "overlay.img")
	if err := os.WriteFile(src, []byte("payload"), 0o600); err != nil {
		t.Fatal(err)
	}
	if err := os.Chmod(src, 0o600); err != nil {
		t.Fatal(err)
	}

	dst := filepath.Join(dstDir, "overlay.img")
	if _, err := MoveFile(context.Background(), src, dst, false); err != nil {
		t.Fatalf("MoveFile: %v", err)
	}
	fi, err := os.Stat(dst)
	if err != nil {
		t.Fatal(err)
	}
	if mode := fi.Mode().Perm(); mode != 0o600 {
		t.Errorf("mode = %#o, want 0600 preserved under a private parent", mode)
	}
}

// A missing destination directory is created rather than reported.
func TestMoveFileCreatesTheDestinationDirectory(t *testing.T) {
	root := t.TempDir()
	src := filepath.Join(root, "x.img")
	if err := os.WriteFile(src, []byte("x"), 0o644); err != nil {
		t.Fatal(err)
	}
	dst := filepath.Join(root, "a", "b", "x.img")
	if _, err := MoveFile(context.Background(), src, dst, false); err != nil {
		t.Fatalf("MoveFile: %v", err)
	}
	if _, err := os.Stat(dst); err != nil {
		t.Errorf("the destination is missing: %v", err)
	}
}
