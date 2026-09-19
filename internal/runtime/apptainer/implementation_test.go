package apptainer

import "testing"

func TestImplementation(t *testing.T) {
	if got := (Bin{Path: "/usr/bin/apptainer"}).Implementation(); got != "apptainer" {
		t.Errorf("Apptainer implementation = %q", got)
	}
	if got := (Bin{Path: "/opt/singularity/bin/singularity"}).Implementation(); got != "singularity" {
		t.Errorf("Singularity implementation = %q", got)
	}
}
