package producer

import (
	"os"
	"path/filepath"
	"strings"
	"testing"
)

func TestAcquireLocalSerializesTarget(t *testing.T) {
	target := filepath.Join(t.TempDir(), "demo.sqf")
	first, err := AcquireLocal(target)
	if err != nil {
		t.Fatal(err)
	}
	defer first.Release() //nolint:errcheck

	if _, err := AcquireLocal(target); err == nil || !strings.Contains(err.Error(), "another build or pull") {
		t.Fatalf("second AcquireLocal error = %v", err)
	}
	if err := first.Release(); err != nil {
		t.Fatal(err)
	}
	second, err := AcquireLocal(target)
	if err != nil {
		t.Fatalf("AcquireLocal after release: %v", err)
	}
	second.Release() //nolint:errcheck
}

func TestAcquireLocalClearsStaleLockAndPartial(t *testing.T) {
	target := filepath.Join(t.TempDir(), "demo.sqf")
	stale := Info{}
	if err := Acquire(Path(target), stale); err != nil {
		t.Fatal(err)
	}
	partial := PreparedPath(target, stale)
	if err := os.WriteFile(partial, []byte("partial"), 0o600); err != nil {
		t.Fatal(err)
	}

	guard, err := AcquireLocal(target)
	if err != nil {
		t.Fatal(err)
	}
	defer guard.Release() //nolint:errcheck
	if _, err := os.Stat(partial); !os.IsNotExist(err) {
		t.Fatalf("stale partial remains: %v", err)
	}
}
