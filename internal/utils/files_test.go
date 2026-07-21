package utils

import (
	"os"
	"path/filepath"
	"syscall"
	"testing"
)

// withUmask sets the process umask for the duration of the test and restores
// the previous value afterward. Umask is process-global, so this test must not
// run in parallel with others that depend on it.
func withUmask(t *testing.T, mask int) {
	t.Helper()
	old := syscall.Umask(mask)
	t.Cleanup(func() { syscall.Umask(old) })
}

func TestShareWithParentGroup_GroupWritableParent(t *testing.T) {
	withUmask(t, 0022)

	parent := t.TempDir()
	if err := os.Chmod(parent, 0775); err != nil {
		t.Fatalf("chmod parent: %v", err)
	}

	// New file created under umask 022 gets 0644 (0664 &^ 022), not group-writable yet.
	filePath := filepath.Join(parent, "file")
	if err := os.WriteFile(filePath, []byte("x"), PermFile); err != nil {
		t.Fatalf("write file: %v", err)
	}
	if mode := statMode(t, filePath); mode != 0644 {
		t.Fatalf("precondition: expected file mode 0644 before sharing, got %o", mode)
	}

	ShareWithParentGroup(filePath)
	if mode := statMode(t, filePath); mode != 0664 {
		t.Errorf("file mode after ShareWithParentGroup = %o, want 0664", mode)
	}

	// New dir created under umask 022 gets 0755, not group-writable yet.
	dirPath := filepath.Join(parent, "dir")
	if err := os.Mkdir(dirPath, PermDir); err != nil {
		t.Fatalf("mkdir: %v", err)
	}
	if mode := statMode(t, dirPath); mode != 0755 {
		t.Fatalf("precondition: expected dir mode 0755 before sharing, got %o", mode)
	}

	ShareWithParentGroup(dirPath)
	if mode := statMode(t, dirPath); mode != 0775 {
		t.Errorf("dir mode after ShareWithParentGroup = %o, want 0775", mode)
	}
}

func TestShareWithParentGroup_NonGroupWritableParent(t *testing.T) {
	withUmask(t, 0022)

	parent := t.TempDir()
	if err := os.Chmod(parent, 0755); err != nil {
		t.Fatalf("chmod parent: %v", err)
	}

	filePath := filepath.Join(parent, "file")
	if err := os.WriteFile(filePath, []byte("x"), PermFile); err != nil {
		t.Fatalf("write file: %v", err)
	}
	if mode := statMode(t, filePath); mode != 0644 {
		t.Fatalf("precondition: expected file mode 0644 before sharing, got %o", mode)
	}

	ShareWithParentGroup(filePath)
	if mode := statMode(t, filePath); mode != 0644 {
		t.Errorf("file mode after ShareWithParentGroup = %o, want unchanged 0644", mode)
	}
}

func TestShareWithParentGroup_ExecutableGainsGroupExecOnlyWhenShared(t *testing.T) {
	withUmask(t, 0022)

	// Group-writable parent: executable gains g+w and g+x.
	sharedParent := t.TempDir()
	if err := os.Chmod(sharedParent, 0775); err != nil {
		t.Fatalf("chmod parent: %v", err)
	}
	execPath := filepath.Join(sharedParent, "exe")
	if err := os.WriteFile(execPath, []byte("x"), 0744); err != nil {
		t.Fatalf("write file: %v", err)
	}
	ShareWithParentGroup(execPath)
	if mode := statMode(t, execPath); mode != 0774 {
		t.Errorf("executable mode under shared parent = %o, want 0774", mode)
	}

	// Non-group-writable parent: executable is left untouched.
	privateParent := t.TempDir()
	if err := os.Chmod(privateParent, 0755); err != nil {
		t.Fatalf("chmod parent: %v", err)
	}
	privExecPath := filepath.Join(privateParent, "exe")
	if err := os.WriteFile(privExecPath, []byte("x"), 0744); err != nil {
		t.Fatalf("write file: %v", err)
	}
	ShareWithParentGroup(privExecPath)
	if mode := statMode(t, privExecPath); mode != 0744 {
		t.Errorf("executable mode under private parent = %o, want unchanged 0744", mode)
	}
}

func TestMakeExecutable_PersonalParentRespectsUmask(t *testing.T) {
	withUmask(t, 0022)

	// Private (non-group-writable) parent: MakeExecutable must add x mirroring the
	// read bits only, never forcing group/other-write like os.Chmod(_, PermExec) would.
	parent := t.TempDir()
	if err := os.Chmod(parent, 0755); err != nil {
		t.Fatalf("chmod parent: %v", err)
	}

	// File created under umask 022 is 0644; MakeExecutable -> 0755 (no group-write leak).
	filePath := filepath.Join(parent, "script")
	f, err := CreateFileWritable(filePath)
	if err != nil {
		t.Fatalf("create file: %v", err)
	}
	f.Close()
	if mode := statMode(t, filePath); mode != 0644 {
		t.Fatalf("precondition: expected 0644 before MakeExecutable, got %o", mode)
	}
	if err := MakeExecutable(filePath); err != nil {
		t.Fatalf("MakeExecutable: %v", err)
	}
	if mode := statMode(t, filePath); mode != 0755 {
		t.Errorf("file mode after MakeExecutable = %o, want 0755", mode)
	}
}

func TestMakeExecutable_SharedParentGainsGroupWriteExec(t *testing.T) {
	withUmask(t, 0022)

	// Group-writable (2775-style) parent: children must become group-writable+exec.
	parent := t.TempDir()
	if err := os.Chmod(parent, 0775); err != nil {
		t.Fatalf("chmod parent: %v", err)
	}

	filePath := filepath.Join(parent, "script")
	f, err := CreateFileWritable(filePath) // ShareWithParentGroup already made it 0664
	if err != nil {
		t.Fatalf("create file: %v", err)
	}
	f.Close()
	if err := MakeExecutable(filePath); err != nil {
		t.Fatalf("MakeExecutable: %v", err)
	}
	if mode := statMode(t, filePath); mode != 0775 {
		t.Errorf("file mode after MakeExecutable under shared parent = %o, want 0775", mode)
	}
}

func TestMakeExecutable_TightUmaskNoWorldExec(t *testing.T) {
	withUmask(t, 0027)

	// umask 027 strips other-read; MakeExecutable must not grant other-exec.
	parent := t.TempDir()
	if err := os.Chmod(parent, 0750); err != nil {
		t.Fatalf("chmod parent: %v", err)
	}
	filePath := filepath.Join(parent, "script")
	f, err := CreateFileWritable(filePath)
	if err != nil {
		t.Fatalf("create file: %v", err)
	}
	f.Close()
	if mode := statMode(t, filePath); mode != 0640 {
		t.Fatalf("precondition: expected 0640 before MakeExecutable, got %o", mode)
	}
	if err := MakeExecutable(filePath); err != nil {
		t.Fatalf("MakeExecutable: %v", err)
	}
	if mode := statMode(t, filePath); mode != 0750 {
		t.Errorf("file mode after MakeExecutable = %o, want 0750 (no world-exec)", mode)
	}
}

func TestMkdirAllShared_InheritsGroupWriteFromParent(t *testing.T) {
	withUmask(t, 0022)

	// Shared (2775-style) parent: created dir must be group-writable (0775).
	shared := t.TempDir()
	if err := os.Chmod(shared, 0775); err != nil {
		t.Fatalf("chmod parent: %v", err)
	}
	sharedChild := filepath.Join(shared, "a", "b") // also exercises parent creation
	if err := MkdirAllShared(sharedChild); err != nil {
		t.Fatalf("MkdirAllShared: %v", err)
	}
	if mode := statMode(t, sharedChild); mode != 0775 {
		t.Errorf("dir under shared parent = %o, want 0775", mode)
	}

	// Personal (non-group-writable) parent: created dir keeps umask default (0755).
	personal := t.TempDir()
	if err := os.Chmod(personal, 0755); err != nil {
		t.Fatalf("chmod parent: %v", err)
	}
	personalChild := filepath.Join(personal, "c")
	if err := MkdirAllShared(personalChild); err != nil {
		t.Fatalf("MkdirAllShared: %v", err)
	}
	if mode := statMode(t, personalChild); mode != 0755 {
		t.Errorf("dir under personal parent = %o, want unchanged 0755", mode)
	}
}

func statMode(t *testing.T, path string) os.FileMode {
	t.Helper()
	info, err := os.Stat(path)
	if err != nil {
		t.Fatalf("stat %s: %v", path, err)
	}
	return info.Mode().Perm()
}
