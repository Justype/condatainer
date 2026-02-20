package scheduler

import (
	"os"
	"testing"
	"time"
)

func TestSpecDefaultsInitialValues(t *testing.T) {
	defaults := GetSpecDefaults()
	if defaults.CpusPerTask != 2 {
		t.Errorf("CpusPerTask = %d; want 2", defaults.CpusPerTask)
	}
	if defaults.Nodes != 1 {
		t.Errorf("Nodes = %d; want 1", defaults.Nodes)
	}
	if defaults.TasksPerNode != 1 {
		t.Errorf("TasksPerNode = %d; want 1", defaults.TasksPerNode)
	}
	if defaults.MemPerNodeMB != 8192 {
		t.Errorf("MemPerNodeMB = %d; want 8192", defaults.MemPerNodeMB)
	}
	if defaults.Time != 4*time.Hour {
		t.Errorf("Time = %v; want 4h", defaults.Time)
	}
}

func TestSpecDefaultsSetAndGet(t *testing.T) {
	// Save original and restore after test
	original := GetSpecDefaults()
	defer SetSpecDefaults(original)

	custom := ResourceSpec{
		CpusPerTask:  8,
		MemPerNodeMB: 16384,
		Time:         4 * time.Hour,
		Nodes:        2,
		TasksPerNode: 4,
	}
	SetSpecDefaults(custom)

	got := GetSpecDefaults()
	if got.CpusPerTask != 8 {
		t.Errorf("CpusPerTask = %d; want 8", got.CpusPerTask)
	}
	if got.MemPerNodeMB != 16384 {
		t.Errorf("MemPerNodeMB = %d; want 16384", got.MemPerNodeMB)
	}
	if got.Time != 4*time.Hour {
		t.Errorf("Time = %v; want 4h", got.Time)
	}
	if got.Nodes != 2 {
		t.Errorf("Nodes = %d; want 2", got.Nodes)
	}
	if got.TasksPerNode != 4 {
		t.Errorf("TasksPerNode = %d; want 4", got.TasksPerNode)
	}
}

func TestSpecDefaultsAffectsReadScriptSpecs(t *testing.T) {
	// Save original and restore after test
	original := GetSpecDefaults()
	defer SetSpecDefaults(original)

	SetSpecDefaults(ResourceSpec{
		CpusPerTask:  16,
		MemPerNodeMB: 32768,
		Time:         8 * time.Hour,
		Nodes:        1,
		TasksPerNode: 1,
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

	if specs.Spec == nil {
		t.Fatal("Spec is nil; expected non-nil ResourceSpec")
	}
	if specs.Spec.CpusPerTask != 16 {
		t.Errorf("CpusPerTask = %d; want 16 (from custom defaults)", specs.Spec.CpusPerTask)
	}
	if specs.Spec.MemPerNodeMB != 32768 {
		t.Errorf("MemPerNodeMB = %d; want 32768 (from custom defaults)", specs.Spec.MemPerNodeMB)
	}
	if specs.Spec.Time != 8*time.Hour {
		t.Errorf("Time = %v; want 8h (from custom defaults)", specs.Spec.Time)
	}
}

func writeTestScript(path, content string) error {
	return os.WriteFile(path, []byte(content), 0644)
}
