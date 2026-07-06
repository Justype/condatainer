package server

import (
	"context"
	"os"
	"path/filepath"
	"syscall"
	"testing"
	"time"
)

func mkTree(t *testing.T, root string) {
	t.Helper()
	for _, d := range []string{"a/b", "c"} {
		if err := os.MkdirAll(filepath.Join(root, d), 0o755); err != nil {
			t.Fatal(err)
		}
	}
	for _, f := range []string{"a/f1", "a/b/f2", "c/f3", "f4"} {
		if err := os.WriteFile(filepath.Join(root, f), []byte("x"), 0o644); err != nil {
			t.Fatal(err)
		}
	}
	if err := os.Symlink(filepath.Join(root, "f4"), filepath.Join(root, "link")); err != nil {
		t.Fatal(err)
	}
}

func TestRemoveTree(t *testing.T) {
	root := t.TempDir()
	tree := filepath.Join(root, "tree")
	mkTree(t, tree)

	total, err := countTreeEntries(tree)
	if err != nil {
		t.Fatal(err)
	}
	// tree + a + a/b + c + 4 files + link = 9
	if total != 9 {
		t.Fatalf("countTreeEntries = %d, want 9", total)
	}

	n := 0
	if err := removeTree(context.Background(), tree, func() { n++ }); err != nil {
		t.Fatal(err)
	}
	if n != total {
		t.Fatalf("onEntry fired %d times, want %d", n, total)
	}
	if _, err := os.Lstat(tree); !os.IsNotExist(err) {
		t.Fatalf("tree still exists: %v", err)
	}
	// the symlink target outside... f4 was inside tree; check root only has nothing left
}

func TestRemoveTreeCancel(t *testing.T) {
	root := t.TempDir()
	tree := filepath.Join(root, "tree")
	mkTree(t, tree)

	ctx, cancel := context.WithCancel(context.Background())
	n := 0
	err := removeTree(ctx, tree, func() {
		n++
		if n == 2 {
			cancel() // cancel mid-delete
		}
	})
	if err != context.Canceled {
		t.Fatalf("err = %v, want context.Canceled", err)
	}
	if _, statErr := os.Lstat(tree); statErr != nil {
		t.Fatalf("tree root should still exist after cancel: %v", statErr)
	}
}

func TestRemoveTreeSymlinkNotFollowed(t *testing.T) {
	root := t.TempDir()
	outside := filepath.Join(root, "outside")
	if err := os.MkdirAll(outside, 0o755); err != nil {
		t.Fatal(err)
	}
	if err := os.WriteFile(filepath.Join(outside, "keep"), []byte("x"), 0o644); err != nil {
		t.Fatal(err)
	}
	tree := filepath.Join(root, "tree")
	if err := os.MkdirAll(tree, 0o755); err != nil {
		t.Fatal(err)
	}
	if err := os.Symlink(outside, filepath.Join(tree, "dirlink")); err != nil {
		t.Fatal(err)
	}
	if err := removeTree(context.Background(), tree, func() {}); err != nil {
		t.Fatal(err)
	}
	if _, err := os.Stat(filepath.Join(outside, "keep")); err != nil {
		t.Fatalf("symlink target was deleted: %v", err)
	}
}

func TestCopyTreeCancelMidFile(t *testing.T) {
	root := t.TempDir()
	src := filepath.Join(root, "src")
	if err := os.MkdirAll(src, 0o755); err != nil {
		t.Fatal(err)
	}
	// 8 MiB file so the chunked copy loops several times
	big := make([]byte, 8<<20)
	if err := os.WriteFile(filepath.Join(src, "big"), big, 0o644); err != nil {
		t.Fatal(err)
	}

	ctx, cancel := context.WithCancel(context.Background())
	cancel() // already cancelled: copy must abort without copying the file
	err := copyTree(ctx, src, filepath.Join(root, "dst"), false, func() {})
	if err != context.Canceled {
		t.Fatalf("err = %v, want context.Canceled", err)
	}

	// and a normal copy still works end-to-end
	n := 0
	if err := copyTree(context.Background(), src, filepath.Join(root, "dst2"), false, func() { n++ }); err != nil {
		t.Fatal(err)
	}
	got, err := os.ReadFile(filepath.Join(root, "dst2", "big"))
	if err != nil || len(got) != len(big) {
		t.Fatalf("copied file wrong: err=%v len=%d", err, len(got))
	}
	if n != 2 { // src dir + big
		t.Fatalf("onEntry fired %d times, want 2", n)
	}
}

// TestCopyTreeSkipsFifoAndKeepsExecBit checks that non-regular files are
// skipped (an opened FIFO would block forever) and that a copied
// executable keeps its exec bit.
func TestCopyTreeSkipsFifoAndKeepsExecBit(t *testing.T) {
	root := t.TempDir()
	src := filepath.Join(root, "src")
	if err := os.MkdirAll(src, 0o755); err != nil {
		t.Fatal(err)
	}
	if err := os.WriteFile(filepath.Join(src, "tool.sh"), []byte("#!/bin/sh\n"), 0o755); err != nil {
		t.Fatal(err)
	}
	if err := os.WriteFile(filepath.Join(src, "data"), []byte("x"), 0o644); err != nil {
		t.Fatal(err)
	}
	if err := syscall.Mkfifo(filepath.Join(src, "pipe"), 0o644); err != nil {
		t.Skipf("mkfifo not supported here: %v", err)
	}

	dst := filepath.Join(root, "dst")
	n := 0
	if err := copyTree(context.Background(), src, dst, false, func() { n++ }); err != nil {
		t.Fatal(err)
	}
	if n != 4 { // src dir + tool.sh + data + pipe (skipped but counted)
		t.Fatalf("onEntry fired %d times, want 4", n)
	}
	if _, err := os.Lstat(filepath.Join(dst, "pipe")); !os.IsNotExist(err) {
		t.Fatalf("fifo was copied: %v", err)
	}
	st, err := os.Stat(filepath.Join(dst, "tool.sh"))
	if err != nil {
		t.Fatal(err)
	}
	if st.Mode()&0o111 == 0 {
		t.Fatalf("exec bit lost: mode = %v", st.Mode())
	}
	st, err = os.Stat(filepath.Join(dst, "data"))
	if err != nil {
		t.Fatal(err)
	}
	if st.Mode()&0o111 != 0 {
		t.Fatalf("non-executable gained exec bit: mode = %v", st.Mode())
	}
}

// TestCopyTreePreservesPermsAndMoveTimes checks cp/mv semantics: exact
// permission bits are always preserved (including on a read-only dir),
// mtimes only with preserveTimes (the cross-filesystem move fallback).
func TestCopyTreePreservesPermsAndMoveTimes(t *testing.T) {
	root := t.TempDir()
	src := filepath.Join(root, "src")
	roDir := filepath.Join(src, "ro")
	if err := os.MkdirAll(roDir, 0o755); err != nil {
		t.Fatal(err)
	}
	if err := os.WriteFile(filepath.Join(roDir, "f"), []byte("x"), 0o600); err != nil {
		t.Fatal(err)
	}
	mtime := time.Date(2020, 1, 2, 3, 4, 5, 0, time.Local)
	if err := os.Chtimes(filepath.Join(roDir, "f"), mtime, mtime); err != nil {
		t.Fatal(err)
	}
	if err := os.Chmod(roDir, 0o555); err != nil {
		t.Fatal(err)
	}
	defer os.Chmod(roDir, 0o755) //nolint:errcheck // so t.TempDir cleanup can remove it

	// preserveTimes=false (copy): perms exact, mtime fresh
	dst := filepath.Join(root, "dst")
	if err := copyTree(context.Background(), src, dst, false, func() {}); err != nil {
		t.Fatal(err)
	}
	defer os.Chmod(filepath.Join(dst, "ro"), 0o755) //nolint:errcheck
	st, err := os.Stat(filepath.Join(dst, "ro"))
	if err != nil {
		t.Fatal(err)
	}
	if st.Mode().Perm() != 0o555 {
		t.Fatalf("dir perm = %v, want 0555", st.Mode().Perm())
	}
	st, err = os.Stat(filepath.Join(dst, "ro", "f"))
	if err != nil {
		t.Fatal(err)
	}
	if st.Mode().Perm() != 0o600 {
		t.Fatalf("file perm = %v, want 0600", st.Mode().Perm())
	}
	if st.ModTime().Equal(mtime) {
		t.Fatalf("plain copy should not preserve mtime")
	}

	// preserveTimes=true (move fallback): mtime preserved too
	dst2 := filepath.Join(root, "dst2")
	if err := copyTree(context.Background(), src, dst2, true, func() {}); err != nil {
		t.Fatal(err)
	}
	defer os.Chmod(filepath.Join(dst2, "ro"), 0o755) //nolint:errcheck
	st, err = os.Stat(filepath.Join(dst2, "ro", "f"))
	if err != nil {
		t.Fatal(err)
	}
	if !st.ModTime().Equal(mtime) {
		t.Fatalf("move should preserve mtime: got %v, want %v", st.ModTime(), mtime)
	}
}

// TestRemoveTreeToleratesVanishingEntries checks that an entry disappearing
// mid-delete doesn't fail the whole delete (os.RemoveAll parity).
func TestRemoveTreeToleratesVanishingEntries(t *testing.T) {
	root := t.TempDir()
	tree := filepath.Join(root, "tree")
	if err := os.MkdirAll(tree, 0o755); err != nil {
		t.Fatal(err)
	}
	for _, f := range []string{"a", "b"} {
		if err := os.WriteFile(filepath.Join(tree, f), []byte("x"), 0o644); err != nil {
			t.Fatal(err)
		}
	}
	first := true
	err := removeTree(context.Background(), tree, func() {
		if first { // after "a" is removed, yank "b" out from under the walk
			first = false
			os.Remove(filepath.Join(tree, "b"))
		}
	})
	if err != nil {
		t.Fatalf("removeTree failed on vanished entry: %v", err)
	}
	if _, err := os.Lstat(tree); !os.IsNotExist(err) {
		t.Fatalf("tree still exists: %v", err)
	}
}
