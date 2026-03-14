package scheduler

import (
	"testing"
	"time"
)

func TestSpecDefaultsInitialValues(t *testing.T) {
	defaults := GetSpecDefaults()
	if defaults.CpusPerTask != 1 {
		t.Errorf("CpusPerTask = %d; want 1", defaults.CpusPerTask)
	}
	if defaults.Nodes != 0 {
		t.Errorf("Nodes = %d; want 0", defaults.Nodes)
	}
	if defaults.TasksPerNode != 0 {
		t.Errorf("TasksPerNode = %d; want 0", defaults.TasksPerNode)
	}
	if defaults.MemPerNodeMB != 1024 {
		t.Errorf("MemPerNodeMB = %d; want 1024", defaults.MemPerNodeMB)
	}
	if defaults.Time != 2*time.Hour {
		t.Errorf("Time = %v; want 2h", defaults.Time)
	}
}
