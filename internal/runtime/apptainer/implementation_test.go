package apptainer

import "testing"

func TestImplementation(t *testing.T) {
	previous := apptainerCmd
	t.Cleanup(func() { apptainerCmd = previous })

	apptainerCmd = "/usr/bin/apptainer"
	if got := Implementation(); got != "apptainer" {
		t.Errorf("Apptainer implementation = %q", got)
	}

	apptainerCmd = "/opt/singularity/bin/singularity"
	if got := Implementation(); got != "singularity" {
		t.Errorf("Singularity implementation = %q", got)
	}
}
