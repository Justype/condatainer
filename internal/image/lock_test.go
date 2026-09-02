package image

import (
	"errors"
	"io/fs"
	"os"
	"path/filepath"
	"testing"
)

// A protected image must be refused as protected, never as in use: the two send
// a reader to entirely different places.
func TestWriteLockTellsProtectedFromInUse(t *testing.T) {
	if os.Geteuid() == 0 {
		t.Skip("root ignores the write bit")
	}
	path := filepath.Join(t.TempDir(), "pinned.sqf")
	if err := os.WriteFile(path, []byte("payload"), 0o444); err != nil {
		t.Fatal(err)
	}

	_, err := AcquireLock(path, true)
	if err == nil {
		t.Fatal("locking a read-only image for writing should fail")
	}
	if !errors.Is(err, ErrProtected) {
		t.Errorf("want ErrProtected, got %v", err)
	}
	if errors.Is(err, ErrInUse) {
		t.Errorf("a protected image is not in use: %v", err)
	}

	// The same image still reads, which is the point of pinning it.
	lock, err := AcquireLock(path, false)
	if err != nil {
		t.Fatalf("a protected image must stay readable: %v", err)
	}
	lock.Close()
}

// Restoring the write bit releases the pin, so protection stays a reversible
// marker rather than a state a user cannot get out of.
func TestClearingProtectionAllowsTheWriteLock(t *testing.T) {
	if os.Geteuid() == 0 {
		t.Skip("root ignores the write bit")
	}
	path := filepath.Join(t.TempDir(), "pinned.sqf")
	if err := os.WriteFile(path, []byte("payload"), 0o444); err != nil {
		t.Fatal(err)
	}
	if err := os.Chmod(path, 0o644); err != nil {
		t.Fatal(err)
	}
	lock, err := AcquireLock(path, true)
	if err != nil {
		t.Fatalf("chmod +w should allow the write lock: %v", err)
	}
	lock.Close()
}

// A held exclusive lock is the one cause that is "in use".
func TestConflictingLockReportsInUse(t *testing.T) {
	path := filepath.Join(t.TempDir(), "busy.sqf")
	if err := os.WriteFile(path, []byte("payload"), 0o644); err != nil {
		t.Fatal(err)
	}
	held, err := AcquireLock(path, true)
	if err != nil {
		t.Fatal(err)
	}
	defer held.Close()

	if err := CheckAvailable(path, false); !errors.Is(err, ErrInUse) {
		t.Errorf("want ErrInUse while an exclusive lock is held, got %v", err)
	}
}

// A missing image is neither protected nor in use.
func TestMissingImageIsNeitherProtectedNorInUse(t *testing.T) {
	path := filepath.Join(t.TempDir(), "absent.sqf")
	err := CheckAvailable(path, true)
	if !errors.Is(err, fs.ErrNotExist) && err == nil {
		t.Fatal("a missing image must fail")
	}
	if errors.Is(err, ErrProtected) || errors.Is(err, ErrInUse) {
		t.Errorf("missing image misreported: %v", err)
	}
}

// A freeze holds a shared lock while it packs: it reads, so it needs no write
// bit, and shared already conflicts with the exclusive lock a writer takes.
func TestSharedLockHeldWhileReadingExcludesAWriter(t *testing.T) {
	path := filepath.Join(t.TempDir(), "env.img")
	if err := os.WriteFile(path, []byte("payload"), 0o444); err != nil {
		t.Fatal(err)
	}

	held, err := AcquireLock(path, false)
	if err != nil {
		t.Fatalf("a protected overlay must stay readable: %v", err)
	}
	defer held.Close()

	// Another reader is fine; a writer is not.
	if err := CheckAvailable(path, false); err != nil {
		t.Errorf("a second reader was refused: %v", err)
	}
	if err := os.Chmod(path, 0o644); err != nil {
		t.Fatal(err)
	}
	if err := CheckAvailable(path, true); !errors.Is(err, ErrInUse) {
		t.Errorf("a writer was allowed in during a read: %v", err)
	}
}
