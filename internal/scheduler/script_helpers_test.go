package scheduler

import (
	"strings"
	"testing"
	"time"
)

func TestFormatHMSTime(t *testing.T) {
	tests := []struct {
		input time.Duration
		want  string
	}{
		{2*time.Hour + 30*time.Minute, "02:30:00"},
		{time.Hour + 30*time.Second, "01:00:30"},
		{90 * time.Minute, "01:30:00"},
		{168 * time.Hour, "168:00:00"},
		{0, "00:00:00"},
	}

	for _, tt := range tests {
		t.Run(tt.want, func(t *testing.T) {
			got := formatHMSTime(tt.input)
			if got != tt.want {
				t.Errorf("formatHMSTime(%v) = %q; want %q", tt.input, got, tt.want)
			}
		})
	}
}

func TestResourceEnvVars(t *testing.T) {
	tests := []struct {
		name         string
		rs           *ResourceSpec
		wantVars     map[string]string
		mustNotExist []string
	}{
		{
			name: "OpenMP: 1 task, 8 CPUs, 32GB per node",
			rs: &ResourceSpec{
				Nodes:        1,
				Ntasks:       1,
				TasksPerNode: 1,
				CpusPerTask:  8,
				MemPerNodeMB: 32768,
			},
			wantVars: map[string]string{
				"NNODES":          "1",
				"NTASKS":          "1",
				"NCPUS":           "8",
				"OMP_NUM_THREADS": "8",
				"NTASKS_PER_NODE": "1",
				"MEM_PER_CPU":     "4096", // 32768 / 8
				"MEM_PER_CPU_MB":  "4096",
				"MEM_PER_CPU_GB":  "4",
				"MEM":             "32768",
				"MEM_MB":          "32768",
				"MEM_GB":          "32",
			},
		},
		{
			name: "MPI: 2 nodes, 16 tasks, 8 tasks/node, 1 CPU/task, 4GB per CPU",
			rs: &ResourceSpec{
				Nodes:        2,
				Ntasks:       16,
				TasksPerNode: 8,
				CpusPerTask:  1,
				MemPerCpuMB:  4096,
			},
			wantVars: map[string]string{
				"NNODES":          "2",
				"NTASKS":          "16",
				"NCPUS":           "1",
				"OMP_NUM_THREADS": "1",
				"NTASKS_PER_NODE": "8",
				"MEM_PER_CPU":     "4096",
				"MEM_PER_CPU_MB":  "4096",
				"MEM_PER_CPU_GB":  "4",
				"MEM":             "32768", // 4096 * 1 * 8
				"MEM_MB":          "32768",
				"MEM_GB":          "32",
			},
		},
		{
			name: "Hybrid: 2 nodes, 8 tasks, 4 tasks/node, 4 CPUs/task, 64GB per node",
			rs: &ResourceSpec{
				Nodes:        2,
				Ntasks:       8,
				TasksPerNode: 4,
				CpusPerTask:  4,
				MemPerNodeMB: 65536,
			},
			wantVars: map[string]string{
				"NNODES":          "2",
				"NTASKS":          "8",
				"NCPUS":           "4",
				"OMP_NUM_THREADS": "4",
				"NTASKS_PER_NODE": "4",
				"MEM_PER_CPU":     "4096", // 65536 / (4 * 4)
				"MEM_PER_CPU_MB":  "4096",
				"MEM_PER_CPU_GB":  "4",
				"MEM":             "65536",
				"MEM_MB":          "65536",
				"MEM_GB":          "64",
			},
		},
		{
			name: "Manually-constructed spec with Nodes but no TasksPerNode: no MEM variables",
			rs: &ResourceSpec{
				Nodes:        2, // Nodes specified
				Ntasks:       7,
				TasksPerNode: 0, // Free distribution (manual edge case; parsing would set this to 4)
				CpusPerTask:  1,
				MemPerCpuMB:  8192,
			},
			wantVars: map[string]string{
				"NNODES":          "2",
				"NTASKS":          "7",
				"NCPUS":           "1",
				"OMP_NUM_THREADS": "1",
				"MEM_PER_CPU":     "8192",
				"MEM_PER_CPU_MB":  "8192",
				"MEM_PER_CPU_GB":  "8",
				// MEM* NOT emitted: cannot guarantee per-node distribution without TasksPerNode
			},
			mustNotExist: []string{"NTASKS_PER_NODE", "MEM", "MEM_MB", "MEM_GB"},
		},
		{
			name: "Free-Distribution MPI without Nodes: 7 tasks only, no node info, 1 CPU/task, 8GB per CPU",
			rs: &ResourceSpec{
				Nodes:        0, // Not set - no node info
				Ntasks:       7,
				TasksPerNode: 0, // Free distribution
				CpusPerTask:  1,
				MemPerCpuMB:  8192,
			},
			wantVars: map[string]string{
				"NNODES":          "1", // Defaults to 1
				"NTASKS":          "7",
				"NCPUS":           "1",
				"OMP_NUM_THREADS": "1",
				"MEM_PER_CPU":     "8192",
				"MEM_PER_CPU_MB":  "8192",
				"MEM_PER_CPU_GB":  "8",
			},
			mustNotExist: []string{"NTASKS_PER_NODE", "MEM", "MEM_MB", "MEM_GB"},
		},
		{
			name: "OpenMP with MemPerCpuMB: Nodes=0, 1 task, 16 CPUs, 2GB per CPU",
			rs: &ResourceSpec{
				Nodes:        0, // Not set - single node implied
				Ntasks:       1,
				TasksPerNode: 1,
				CpusPerTask:  16,
				MemPerCpuMB:  2048,
			},
			wantVars: map[string]string{
				"NNODES":          "1", // Defaults to 1
				"NTASKS":          "1",
				"NCPUS":           "16",
				"OMP_NUM_THREADS": "16",
				"NTASKS_PER_NODE": "1",
				"MEM_PER_CPU":     "2048",
				"MEM_PER_CPU_MB":  "2048",
				"MEM_PER_CPU_GB":  "2",
				"MEM":             "32768", // 2048 * 16 * 1
				"MEM_MB":          "32768",
				"MEM_GB":          "32",
			},
		},
		{
			name: "No memory specified",
			rs: &ResourceSpec{
				Nodes:        2,
				Ntasks:       8,
				TasksPerNode: 4,
				CpusPerTask:  2,
			},
			wantVars: map[string]string{
				"NNODES":          "2",
				"NTASKS":          "8",
				"NCPUS":           "2",
				"OMP_NUM_THREADS": "2",
				"NTASKS_PER_NODE": "4",
			},
			mustNotExist: []string{"MEM_PER_CPU", "MEM_PER_CPU_MB", "MEM_PER_CPU_GB", "MEM", "MEM_MB", "MEM_GB"},
		},
		{
			name: "Minimal defaults: nil ResourceSpec",
			rs:   nil,
			wantVars: map[string]string{
				"NNODES":          "1",
				"NTASKS":          "1",
				"NCPUS":           "1",
				"OMP_NUM_THREADS": "1",
			},
			mustNotExist: []string{"NTASKS_PER_NODE", "MEM_PER_CPU", "MEM"},
		},
		{
			name: "Ntasks from topology: 2 nodes, 4 tasks/node (Ntasks not set)",
			rs: &ResourceSpec{
				Nodes:        2,
				TasksPerNode: 4,
				CpusPerTask:  2,
			},
			wantVars: map[string]string{
				"NNODES":          "2",
				"NTASKS":          "8", // 2 * 4
				"NCPUS":           "2",
				"OMP_NUM_THREADS": "2",
				"NTASKS_PER_NODE": "4",
			},
		},
	}

	for _, tt := range tests {
		t.Run(tt.name, func(t *testing.T) {
			envVars := ResourceEnvVars(tt.rs)

			// Convert slice to map for easier checking
			gotMap := make(map[string]string)
			for _, env := range envVars {
				parts := strings.SplitN(env, "=", 2)
				if len(parts) == 2 {
					gotMap[parts[0]] = parts[1]
				}
			}

			// Check expected variables exist with correct values
			for key, wantVal := range tt.wantVars {
				gotVal, exists := gotMap[key]
				if !exists {
					t.Errorf("Expected variable %s not found in output", key)
					continue
				}
				if gotVal != wantVal {
					t.Errorf("Variable %s = %q; want %q", key, gotVal, wantVal)
				}
			}

			// Check that unwanted variables don't exist
			for _, key := range tt.mustNotExist {
				if val, exists := gotMap[key]; exists {
					t.Errorf("Variable %s should not exist, but found with value %q", key, val)
				}
			}
		})
	}
}
