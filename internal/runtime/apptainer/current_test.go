package apptainer

import (
	"errors"
	"testing"
)

// Current reads back what the latest resolver returned; it must never resolve
// on its own, so with nothing resolved it refuses rather than searching PATH.
func TestCurrentRefusesWithNothingResolved(t *testing.T) {
	resetApptainerState(t)

	if _, _, err := Current(); !errors.Is(err, ErrNotResolved) {
		t.Errorf("err = %v, want %v", err, ErrNotResolved)
	}
}

// Current reports whatever binary was already resolved, without re-deciding
// which one that should be.
func TestCurrentReadsBackTheResolvedBinary(t *testing.T) {
	resetApptainerState(t)
	systemApptainer(t, writeFakeBin(t, t.TempDir(), "singularity", "singularity-ce version 4.1.1"))
	if _, err := ForBuild(); err != nil {
		t.Fatalf("ForBuild: %v", err)
	}

	implementation, version, err := Current()
	if err != nil {
		t.Fatalf("Current: %v", err)
	}
	if implementation != "singularity" {
		t.Errorf("implementation = %q, want singularity", implementation)
	}
	if version != "4.1.1" {
		t.Errorf("version = %q, want 4.1.1", version)
	}
}
