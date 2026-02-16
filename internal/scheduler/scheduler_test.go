package scheduler

import (
	"os"
	"testing"
	"time"
)

func TestSpecDefaultsInitialValues(t *testing.T) {
	defaults := GetSpecDefaults()
	if defaults.Ncpus != 2 {
		t.Errorf("Ncpus = %d; want 2", defaults.Ncpus)
	}
	if defaults.Nodes != 1 {
		t.Errorf("Nodes = %d; want 1", defaults.Nodes)
	}
	if defaults.Ntasks != 1 {
		t.Errorf("Ntasks = %d; want 1", defaults.Ntasks)
	}
	if defaults.MemMB != 8192 {
		t.Errorf("MemMB = %d; want 8192", defaults.MemMB)
	}
	if defaults.Time != 4*time.Hour {
		t.Errorf("Time = %v; want 4h", defaults.Time)
	}
}

func TestSpecDefaultsSetAndGet(t *testing.T) {
	// Save original and restore after test
	original := GetSpecDefaults()
	defer SetSpecDefaults(original)

	custom := SpecDefaults{
		Ncpus:  8,
		MemMB:  16384,
		Time:   4 * time.Hour,
		Nodes:  2,
		Ntasks: 4,
	}
	SetSpecDefaults(custom)

	got := GetSpecDefaults()
	if got.Ncpus != 8 {
		t.Errorf("Ncpus = %d; want 8", got.Ncpus)
	}
	if got.MemMB != 16384 {
		t.Errorf("MemMB = %d; want 16384", got.MemMB)
	}
	if got.Time != 4*time.Hour {
		t.Errorf("Time = %v; want 4h", got.Time)
	}
	if got.Nodes != 2 {
		t.Errorf("Nodes = %d; want 2", got.Nodes)
	}
	if got.Ntasks != 4 {
		t.Errorf("Ntasks = %d; want 4", got.Ntasks)
	}
}

func TestSpecDefaultsAffectsReadScriptSpecs(t *testing.T) {
	// Save original and restore after test
	original := GetSpecDefaults()
	defer SetSpecDefaults(original)

	SetSpecDefaults(SpecDefaults{
		Ncpus:  16,
		MemMB:  32768,
		Time:   8 * time.Hour,
		Nodes:  1,
		Ntasks: 1,
	})

	// Parse a script with no resource directives — should pick up custom defaults
	slurm := newTestSlurmScheduler()
	tmpDir := t.TempDir()
	scriptPath := tmpDir + "/test.sh"
	if err := writeTestScript(scriptPath, "#!/bin/bash\n#SBATCH --job-name=test\necho hello"); err != nil {
		t.Fatal(err)
	}

	specs, err := slurm.ReadScriptSpecs(scriptPath)
	if err != nil {
		t.Fatalf("ReadScriptSpecs failed: %v", err)
	}

	if specs.Ncpus != 16 {
		t.Errorf("Ncpus = %d; want 16 (from custom defaults)", specs.Ncpus)
	}
	if specs.MemMB != 32768 {
		t.Errorf("MemMB = %d; want 32768 (from custom defaults)", specs.MemMB)
	}
	if specs.Time != 8*time.Hour {
		t.Errorf("Time = %v; want 8h (from custom defaults)", specs.Time)
	}
}

func writeTestScript(path, content string) error {
	return os.WriteFile(path, []byte(content), 0644)
}
