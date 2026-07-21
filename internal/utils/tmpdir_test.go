package utils

import (
	"os"
	"path/filepath"
	"testing"
)

func TestEnsureTmpSubdir_WorldWritableParentIsPrivate(t *testing.T) {
	withUmask(t, 0022)

	parent := t.TempDir()
	if err := os.Chmod(parent, 0777); err != nil { // o+w, like /tmp
		t.Fatalf("chmod parent: %v", err)
	}
	leaf := filepath.Join(parent, "cnt-user")
	if err := EnsureTmpSubdir(leaf); err != nil {
		t.Fatalf("EnsureTmpSubdir: %v", err)
	}
	if mode := statMode(t, leaf); mode != 0700 {
		t.Errorf("leaf under world-writable parent = %o, want 0700", mode)
	}
}

func TestEnsureTmpSubdir_SharedParentInheritsGroupWrite(t *testing.T) {
	withUmask(t, 0022)

	parent := t.TempDir()
	if err := os.Chmod(parent, 0775); err != nil { // group-writable, not world-writable
		t.Fatalf("chmod parent: %v", err)
	}
	leaf := filepath.Join(parent, "cnt-user")
	if err := EnsureTmpSubdir(leaf); err != nil {
		t.Fatalf("EnsureTmpSubdir: %v", err)
	}
	if mode := statMode(t, leaf); mode != 0775 {
		t.Errorf("leaf under shared parent = %o, want 0775", mode)
	}
}

func TestEnsureTmpSubdir_PrivateParentKeepsUmaskDefault(t *testing.T) {
	withUmask(t, 0022)

	parent := t.TempDir()
	if err := os.Chmod(parent, 0755); err != nil { // neither group- nor world-writable
		t.Fatalf("chmod parent: %v", err)
	}
	leaf := filepath.Join(parent, "cnt-user")
	if err := EnsureTmpSubdir(leaf); err != nil {
		t.Fatalf("EnsureTmpSubdir: %v", err)
	}
	if mode := statMode(t, leaf); mode != 0755 {
		t.Errorf("leaf under private parent = %o, want 0755", mode)
	}
}
