package lock

import (
	"errors"
	"os"
	"path/filepath"
	"testing"
)

// RootAbove finds a project from its own root, same as RootAt.
func TestRootAboveFindsTheRootItself(t *testing.T) {
	root := projectRoot(t)

	got, err := RootAbove(root)
	if err != nil {
		t.Fatal(err)
	}
	if got != root {
		t.Fatalf("RootAbove(root) = %q, want %q", got, root)
	}
}

// A subdirectory of a project is not a project for RootAt, but RootAbove
// walks up and finds the same root.
func TestRootAboveWalksUpFromASubdirectory(t *testing.T) {
	root := projectRoot(t)
	sub := filepath.Join(root, "steps1")
	if err := os.MkdirAll(sub, 0o775); err != nil {
		t.Fatal(err)
	}

	if _, err := RootAt(sub); !errors.Is(err, ErrNoProject) {
		t.Fatalf("RootAt(sub) = %v, want ErrNoProject", err)
	}
	got, err := RootAbove(sub)
	if err != nil {
		t.Fatal(err)
	}
	if got != root {
		t.Fatalf("RootAbove(sub) = %q, want %q", got, root)
	}
}

// No cnt-lock/ at or above dir is ErrNoProject, not a crash at the
// filesystem root.
func TestRootAboveErrorsWithNoProjectAnywhereAbove(t *testing.T) {
	dir := t.TempDir()

	if _, err := RootAbove(dir); !errors.Is(err, ErrNoProject) {
		t.Fatalf("RootAbove(dir) = %v, want ErrNoProject", err)
	}
}
