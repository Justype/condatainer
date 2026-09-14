package apptainer

import (
	"errors"
	"testing"
)

// Current reads back whatever SetBin/ResolveBin/EnsureApptainer already
// decided; it must never resolve on its own, so with nothing configured it
// refuses rather than falling back to PATH.
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
	dir := t.TempDir()
	binPath := writeFakeBin(t, dir, "singularity", "singularity-ce version 4.1.1")

	if err := SetBin(binPath); err != nil {
		t.Fatalf("SetBin: %v", err)
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
