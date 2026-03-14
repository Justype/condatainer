package scheduler

import (
	"fmt"
	"os"
	"path/filepath"
	"regexp"
	"strings"
	"testing"
	"time"
)

func newTestPbsScheduler() *PbsScheduler {
	return &PbsScheduler{
		qsubBin:     "/usr/bin/qsub",
		directiveRe: regexp.MustCompile(`^\s*#PBS\s+(.+)$`),
		jobIDRe:     regexp.MustCompile(`^(\d+\..*|^\d+)$`),
	}
}

// Test that scheduler uses outputDir for log files
func TestPbsCreateScriptUsesOutputDir(t *testing.T) {
	tmpOutputDir := t.TempDir()

	pbs := newTestPbsScheduler()

	jobSpec := &JobSpec{
		Name:    "pbstest/job",
		Command: "echo 'hello'",
		Specs: &ScriptSpecs{
			RemainingFlags: []string{},
			Spec: &ResourceSpec{
				CpusPerTask: 1,
			},
		},
	}

	scriptPath, err := pbs.CreateScriptWithSpec(jobSpec, tmpOutputDir)
	if err != nil {
		t.Fatalf("CreateScriptWithSpec failed: %v", err)
	}

	// Ensure the generated script uses outputDir for logs
	logPath := filepath.Join(tmpOutputDir, "pbstest--job.log")
	content, err := os.ReadFile(scriptPath)
	if err != nil {
		t.Fatalf("Failed to read generated script: %v", err)
	}
	if !strings.Contains(string(content), fmt.Sprintf("-o %s", logPath)) {
		t.Errorf("Generated script does not contain expected output path %s\nScript:\n%s", logPath, string(content))
	}
}

func TestPbsEmailParsing(t *testing.T) {
	tests := []struct {
		name         string
		scriptLines  []string
		wantBegin    bool
		wantEnd      bool
		wantFail     bool
		wantMailUser string
	}{
		{
			name: "all notifications (abe)",
			scriptLines: []string{
				"#!/bin/bash",
				"#PBS -m abe",
				"#PBS -M test@example.com",
			},
			wantBegin:    true,
			wantEnd:      true,
			wantFail:     true,
			wantMailUser: "test@example.com",
		},
		{
			name: "begin and end only (be)",
			scriptLines: []string{
				"#!/bin/bash",
				"#PBS -m be",
				"#PBS -M user@domain.org",
			},
			wantBegin:    true,
			wantEnd:      true,
			wantFail:     false,
			wantMailUser: "user@domain.org",
		},
		{
			name: "abort/fail only (a)",
			scriptLines: []string{
				"#!/bin/bash",
				"#PBS -m a",
			},
			wantBegin: false,
			wantEnd:   false,
			wantFail:  true,
		},
		{
			name: "end only (e)",
			scriptLines: []string{
				"#!/bin/bash",
				"#PBS -m e",
			},
			wantBegin: false,
			wantEnd:   true,
			wantFail:  false,
		},
		{
			name: "begin only (b)",
			scriptLines: []string{
				"#!/bin/bash",
				"#PBS -m b",
			},
			wantBegin: true,
			wantEnd:   false,
			wantFail:  false,
		},
		{
			name: "none disables all (n)",
			scriptLines: []string{
				"#!/bin/bash",
				"#PBS -m n",
			},
			wantBegin: false,
			wantEnd:   false,
			wantFail:  false,
		},
		{
			name: "no email directives",
			scriptLines: []string{
				"#!/bin/bash",
				"#PBS -N testjob",
			},
			wantBegin: false,
			wantEnd:   false,
			wantFail:  false,
		},
		{
			name: "mail user without mail options",
			scriptLines: []string{
				"#!/bin/bash",
				"#PBS -M admin@company.com",
			},
			wantBegin:    false,
			wantEnd:      false,
			wantFail:     false,
			wantMailUser: "admin@company.com",
		},
	}

	for _, tt := range tests {
		t.Run(tt.name, func(t *testing.T) {
			tmpDir := t.TempDir()
			scriptPath := filepath.Join(tmpDir, "test.sh")
			content := strings.Join(tt.scriptLines, "\n")
			if err := os.WriteFile(scriptPath, []byte(content), 0644); err != nil {
				t.Fatalf("Failed to create test script: %v", err)
			}

			pbs := newTestPbsScheduler()
			specs, err := pbs.ReadScriptSpecs(scriptPath)
			if err != nil {
				t.Fatalf("Failed to parse script: %v", err)
			}

			if specs.Control.EmailOnBegin != tt.wantBegin {
				t.Errorf("EmailOnBegin = %v; want %v", specs.Control.EmailOnBegin, tt.wantBegin)
			}
			if specs.Control.EmailOnEnd != tt.wantEnd {
				t.Errorf("EmailOnEnd = %v; want %v", specs.Control.EmailOnEnd, tt.wantEnd)
			}
			if specs.Control.EmailOnFail != tt.wantFail {
				t.Errorf("EmailOnFail = %v; want %v", specs.Control.EmailOnFail, tt.wantFail)
			}
			if specs.Control.MailUser != tt.wantMailUser {
				t.Errorf("MailUser = %q; want %q", specs.Control.MailUser, tt.wantMailUser)
			}
		})
	}
}

func TestPbsEmailScriptGeneration(t *testing.T) {
	tests := []struct {
		name         string
		emailOnBegin bool
		emailOnEnd   bool
		emailOnFail  bool
		mailUser     string
		wantMailOpts string // expected -m value (empty = not present)
		wantMailUser string // expected -M value (empty = not present)
	}{
		{
			name:         "all three notifications",
			emailOnBegin: true,
			emailOnEnd:   true,
			emailOnFail:  true,
			mailUser:     "test@example.com",
			wantMailOpts: "bea",
			wantMailUser: "test@example.com",
		},
		{
			name:         "begin and end only",
			emailOnBegin: true,
			emailOnEnd:   true,
			emailOnFail:  false,
			mailUser:     "user@domain.org",
			wantMailOpts: "be",
			wantMailUser: "user@domain.org",
		},
		{
			name:         "end only",
			emailOnBegin: false,
			emailOnEnd:   true,
			emailOnFail:  false,
			wantMailOpts: "e",
			wantMailUser: "",
		},
		{
			name:         "fail only",
			emailOnBegin: false,
			emailOnEnd:   false,
			emailOnFail:  true,
			mailUser:     "admin@company.com",
			wantMailOpts: "a",
			wantMailUser: "admin@company.com",
		},
		{
			name:         "no email notifications",
			emailOnBegin: false,
			emailOnEnd:   false,
			emailOnFail:  false,
			wantMailOpts: "",
			wantMailUser: "",
		},
	}

	for _, tt := range tests {
		t.Run(tt.name, func(t *testing.T) {
			tmpDir := t.TempDir()
			pbs := newTestPbsScheduler()

			jobSpec := &JobSpec{
				Name:    "test_job",
				Command: "echo 'test'",
				Specs: &ScriptSpecs{
					Control: RuntimeConfig{
						JobName:      "test_job",
						EmailOnBegin: tt.emailOnBegin,
						EmailOnEnd:   tt.emailOnEnd,
						EmailOnFail:  tt.emailOnFail,
						MailUser:     tt.mailUser,
					},
					Spec: &ResourceSpec{
						CpusPerTask:  4,
						MemPerNodeMB: 8000,
						Time:         time.Hour,
					},
					RemainingFlags: []string{},
				},
			}

			scriptPath, err := pbs.CreateScriptWithSpec(jobSpec, tmpDir)
			if err != nil {
				t.Fatalf("Failed to create script: %v", err)
			}

			content, err := os.ReadFile(scriptPath)
			if err != nil {
				t.Fatalf("Failed to read generated script: %v", err)
			}
			scriptContent := string(content)

			// Check -m directive
			if tt.wantMailOpts != "" {
				expectedLine := "#PBS -m " + tt.wantMailOpts
				if !strings.Contains(scriptContent, expectedLine) {
					t.Errorf("Script missing expected line: %q\nScript content:\n%s", expectedLine, scriptContent)
				}
			} else {
				if strings.Contains(scriptContent, "#PBS -m ") {
					t.Errorf("Script should not contain -m directive\nScript content:\n%s", scriptContent)
				}
			}

			// Check -M directive
			if tt.wantMailUser != "" {
				expectedLine := "#PBS -M " + tt.wantMailUser
				if !strings.Contains(scriptContent, expectedLine) {
					t.Errorf("Script missing expected line: %q\nScript content:\n%s", expectedLine, scriptContent)
				}
			} else {
				if strings.Contains(scriptContent, "#PBS -M ") {
					t.Errorf("Script should not contain -M directive\nScript content:\n%s", scriptContent)
				}
			}
		})
	}
}

func TestPbsEmailRoundTrip(t *testing.T) {
	tmpDir := t.TempDir()

	originalScript := `#!/bin/bash
#PBS -N email_test
#PBS -l select=1:ncpus=8:mem=16gb
#PBS -l walltime=02:00:00
#PBS -m abe
#PBS -M roundtrip@example.com

echo "Running job"
`
	scriptPath := filepath.Join(tmpDir, "original.sh")
	if err := os.WriteFile(scriptPath, []byte(originalScript), 0644); err != nil {
		t.Fatalf("Failed to create test script: %v", err)
	}

	pbs := newTestPbsScheduler()

	// Parse the original script
	specs, err := pbs.ReadScriptSpecs(scriptPath)
	if err != nil {
		t.Fatalf("Failed to parse script: %v", err)
	}

	// Verify parsing
	if !specs.Control.EmailOnBegin || !specs.Control.EmailOnEnd || !specs.Control.EmailOnFail {
		t.Errorf("Parsing failed: EmailOnBegin=%v, EmailOnEnd=%v, EmailOnFail=%v",
			specs.Control.EmailOnBegin, specs.Control.EmailOnEnd, specs.Control.EmailOnFail)
	}
	if specs.Control.MailUser != "roundtrip@example.com" {
		t.Errorf("MailUser = %q; want %q", specs.Control.MailUser, "roundtrip@example.com")
	}
	if specs.Spec == nil {
		t.Fatal("Spec is nil; want non-nil")
	}
	if specs.Spec.Time != 2*time.Hour {
		t.Errorf("Time = %v; want 2h", specs.Spec.Time)
	}
	if specs.Spec.CpusPerTask != 8 {
		t.Errorf("CpusPerTask = %d; want 8", specs.Spec.CpusPerTask)
	}
	if specs.Spec.MemPerNodeMB != 16*1024 {
		t.Errorf("MemPerNodeMB = %d; want %d", specs.Spec.MemPerNodeMB, 16*1024)
	}

	// Generate a new script from the parsed specs
	jobSpec := &JobSpec{
		Name:    "email_test",
		Command: "echo 'Running job'",
		Specs:   specs,
	}

	newScriptPath, err := pbs.CreateScriptWithSpec(jobSpec, tmpDir)
	if err != nil {
		t.Fatalf("Failed to create script: %v", err)
	}

	content, err := os.ReadFile(newScriptPath)
	if err != nil {
		t.Fatalf("Failed to read generated script: %v", err)
	}
	scriptContent := string(content)

	// Verify the generated script contains the email directives
	if !strings.Contains(scriptContent, "#PBS -m bea") {
		t.Errorf("Generated script missing -m directive:\n%s", scriptContent)
	}
	if !strings.Contains(scriptContent, "#PBS -M roundtrip@example.com") {
		t.Errorf("Generated script missing -M directive:\n%s", scriptContent)
	}

	// Verify NO duplicates
	mailOptsCount := strings.Count(scriptContent, "#PBS -m ")
	mailUserCount := strings.Count(scriptContent, "#PBS -M ")
	if mailOptsCount != 1 {
		t.Errorf("Expected exactly 1 -m directive, found %d:\n%s", mailOptsCount, scriptContent)
	}
	if mailUserCount != 1 {
		t.Errorf("Expected exactly 1 -M directive, found %d:\n%s", mailUserCount, scriptContent)
	}
}

func TestPbsResourceParsing(t *testing.T) {
	tests := []struct {
		name            string
		lines           []string
		wantCpusPerTask int
		wantMemMB       int64
		wantTime        time.Duration
		wantGpu         *GpuSpec
		wantPassthrough bool
	}{
		{
			name: "select format with ncpus, mem",
			lines: []string{
				"#!/bin/bash",
				"#PBS -l select=1:ncpus=4:mem=8gb",
			},
			wantCpusPerTask: 4,
			wantMemMB:       8 * 1024,
		},
		{
			name: "walltime only",
			lines: []string{
				"#!/bin/bash",
				"#PBS -l walltime=02:30:00",
			},
			wantCpusPerTask: 1, // default
			wantTime:        2*time.Hour + 30*time.Minute,
		},
		{
			name: "nodes/ppn format",
			lines: []string{
				"#!/bin/bash",
				"#PBS -l nodes=1:ppn=8",
			},
			wantPassthrough: true,
		},
		{
			name: "ngpus resource",
			lines: []string{
				"#!/bin/bash",
				"#PBS -l select=1:ncpus=4:mem=16gb:ngpus=2",
			},
			wantCpusPerTask: 4,
			wantMemMB:       16 * 1024,
			wantGpu:         &GpuSpec{Type: "gpu", Count: 2, Raw: "ngpus=2"},
		},
		{
			// #PBS -l gpus=... is old Torque syntax; PBS Pro requires ngpus= inside select=.
			name: "gpus standalone -l is old Torque syntax: passthrough",
			lines: []string{
				"#!/bin/bash",
				"#PBS -l select=1:ncpus=4",
				"#PBS -l gpus=a100:2",
			},
			wantPassthrough: true,
		},
		{
			name: "multiple resource lines",
			lines: []string{
				"#!/bin/bash",
				"#PBS -l select=1:ncpus=16:mem=32gb",
				"#PBS -l walltime=04:00:00",
			},
			wantCpusPerTask: 16,
			wantMemMB:       32 * 1024,
			wantTime:        4 * time.Hour,
		},
		{
			name: "select with all resources",
			lines: []string{
				"#!/bin/bash",
				"#PBS -N fulltest",
				"#PBS -l select=1:ncpus=8:mem=64gb:ngpus=4",
				"#PBS -l walltime=12:00:00",
				"#PBS -m abe",
				"#PBS -M user@example.com",
			},
			wantCpusPerTask: 8,
			wantMemMB:       64 * 1024,
			wantTime:        12 * time.Hour,
			wantGpu:         &GpuSpec{Type: "gpu", Count: 4, Raw: "ngpus=4"},
		},
		{
			name: "memory in mb",
			lines: []string{
				"#!/bin/bash",
				"#PBS -l select=1:ncpus=2:mem=4096mb",
			},
			wantCpusPerTask: 2,
			wantMemMB:       4096,
		},
	}

	for _, tt := range tests {
		t.Run(tt.name, func(t *testing.T) {
			tmpDir := t.TempDir()
			scriptPath := filepath.Join(tmpDir, "test.sh")
			content := strings.Join(tt.lines, "\n")
			if err := os.WriteFile(scriptPath, []byte(content), 0644); err != nil {
				t.Fatalf("Failed to create test script: %v", err)
			}

			pbs := newTestPbsScheduler()
			specs, err := pbs.ReadScriptSpecs(scriptPath)
			if err != nil {
				t.Fatalf("Failed to parse script: %v", err)
			}

			if tt.wantPassthrough {
				if specs.Spec != nil {
					t.Errorf("Spec = %+v; want nil (passthrough)", specs.Spec)
				}
				return
			}

			if specs.Spec == nil {
				t.Fatal("Spec is nil; want non-nil")
			}

			if specs.Spec.CpusPerTask != tt.wantCpusPerTask {
				t.Errorf("CpusPerTask = %d; want %d", specs.Spec.CpusPerTask, tt.wantCpusPerTask)
			}
			if tt.wantMemMB > 0 && specs.Spec.MemPerNodeMB != tt.wantMemMB {
				t.Errorf("MemPerNodeMB = %d; want %d", specs.Spec.MemPerNodeMB, tt.wantMemMB)
			}
			if tt.wantTime > 0 && specs.Spec.Time != tt.wantTime {
				t.Errorf("Time = %v; want %v", specs.Spec.Time, tt.wantTime)
			}
			if tt.wantGpu != nil {
				if specs.Spec.Gpu == nil {
					t.Fatal("Gpu is nil; want non-nil")
				}
				if specs.Spec.Gpu.Type != tt.wantGpu.Type {
					t.Errorf("Gpu.Type = %q; want %q", specs.Spec.Gpu.Type, tt.wantGpu.Type)
				}
				if specs.Spec.Gpu.Count != tt.wantGpu.Count {
					t.Errorf("Gpu.Count = %d; want %d", specs.Spec.Gpu.Count, tt.wantGpu.Count)
				}
			} else if specs.Spec.Gpu != nil {
				t.Errorf("Gpu = %+v; want nil", specs.Spec.Gpu)
			}
		})
	}
}

func TestPbsMemoryParsing(t *testing.T) {
	t.Run("parsePbsMemory (returns KB)", func(t *testing.T) {
		tests := []struct {
			input  string
			wantKB int64
		}{
			{"8gb", 8 * 1024 * 1024},
			{"1024mb", 1024 * 1024},
			{"4096kb", 4096},
			{"1tb", 1024 * 1024 * 1024},
			{"1048576b", 1024},
		}

		for _, tt := range tests {
			kb, err := parsePbsMemory(tt.input)
			if err != nil {
				t.Errorf("parsePbsMemory(%q) error: %v", tt.input, err)
				continue
			}
			if kb != tt.wantKB {
				t.Errorf("parsePbsMemory(%q) = %d KB; want %d KB", tt.input, kb, tt.wantKB)
			}
		}
	})

}

func TestPbsIsAvailable(t *testing.T) {
	t.Run("with binary", func(t *testing.T) {
		pbs := newTestPbsScheduler()
		if !pbs.IsAvailable() {
			t.Error("Expected IsAvailable to return true when binary is set")
		}
	})

	t.Run("no binary", func(t *testing.T) {
		pbs := &PbsScheduler{}
		if pbs.IsAvailable() {
			t.Error("Expected IsAvailable to return false when no binary is set")
		}
	})
}

func TestPbsIsInsideJob(t *testing.T) {
	t.Run("inside job", func(t *testing.T) {
		clearJobEnvVars(t)
		t.Setenv("PBS_JOBID", "67890.pbs-server")
		pbs := newTestPbsScheduler()
		if !pbs.IsInsideJob() {
			t.Error("Expected IsInsideJob to return true when PBS_JOBID is set")
		}
	})

	t.Run("not in job", func(t *testing.T) {
		clearJobEnvVars(t)
		pbs := newTestPbsScheduler()
		if pbs.IsInsideJob() {
			t.Error("Expected IsInsideJob to return false when PBS_JOBID is not set")
		}
	})
}

func TestPbsGetInfo(t *testing.T) {
	clearJobEnvVars(t)
	pbs := newTestPbsScheduler()

	info := pbs.GetInfo()
	if info.Type != "PBS" {
		t.Errorf("Type = %q; want %q", info.Type, "PBS")
	}
	if info.Binary != "/usr/bin/qsub" {
		t.Errorf("Binary = %q; want %q", info.Binary, "/usr/bin/qsub")
	}
	if info.InJob {
		t.Error("InJob should be false")
	}
	if !info.Available {
		t.Error("Available should be true")
	}
}

func TestTryParsePbsScript(t *testing.T) {
	tmpDir := t.TempDir()

	script := `#!/bin/bash
#PBS -N testjob
#PBS -l select=1:ncpus=4:mem=8gb:ngpus=1
#PBS -l walltime=02:00:00
#PBS -m e
#PBS -M user@example.com

echo "Running job"
`
	scriptPath := filepath.Join(tmpDir, "test.sh")
	if err := os.WriteFile(scriptPath, []byte(script), 0644); err != nil {
		t.Fatalf("Failed to create test script: %v", err)
	}

	specs, err := TryParsePbsScript(scriptPath)
	if err != nil {
		t.Fatalf("TryParsePbsScript failed: %v", err)
	}

	if specs.Spec == nil {
		t.Fatal("Spec is nil; want non-nil")
	}
	if specs.Spec.CpusPerTask != 4 {
		t.Errorf("CpusPerTask = %d; want 4", specs.Spec.CpusPerTask)
	}
	if specs.Spec.MemPerNodeMB != 8*1024 {
		t.Errorf("MemPerNodeMB = %d; want %d", specs.Spec.MemPerNodeMB, 8*1024)
	}
	if specs.Spec.Gpu == nil || specs.Spec.Gpu.Count != 1 {
		t.Errorf("Gpu = %+v; want count 1", specs.Spec.Gpu)
	}
	if specs.Spec.Time != 2*time.Hour {
		t.Errorf("Time = %v; want 2h", specs.Spec.Time)
	}
	if !specs.Control.EmailOnEnd {
		t.Error("EmailOnEnd should be true")
	}
	if specs.Control.MailUser != "user@example.com" {
		t.Errorf("MailUser = %q; want %q", specs.Control.MailUser, "user@example.com")
	}
	// RemainingFlags should be empty - all flags above are recognized and parsed into typed fields
	if len(specs.RemainingFlags) != 0 {
		t.Errorf("RemainingFlags count = %d; want 0 (all flags were recognized)", len(specs.RemainingFlags))
	}
}

func TestPbsNodeTaskParsing(t *testing.T) {
	tests := []struct {
		name             string
		lines            []string
		wantNodes        int
		wantTasksPerNode int
		wantCpusPerTask  int
		wantPassthrough  bool
	}{
		{
			name: "select with mpiprocs",
			lines: []string{
				"#!/bin/bash",
				"#PBS -l select=4:ncpus=8:mpiprocs=2:mem=8gb",
			},
			wantNodes:        4,
			wantTasksPerNode: 2,
			wantCpusPerTask:  4, // ncpus(8) / mpiprocs(2) = 4 threads per task
		},
		{
			name: "select without mpiprocs",
			lines: []string{
				"#!/bin/bash",
				"#PBS -l select=2:ncpus=4",
			},
			wantNodes:        2,
			wantTasksPerNode: 1, // default, no mpiprocs
			wantCpusPerTask:  4,
		},
		{
			name: "nodes format",
			lines: []string{
				"#!/bin/bash",
				"#PBS -l nodes=3:ppn=8",
			},
			wantPassthrough: true,
		},
		{
			name: "select=1 single node",
			lines: []string{
				"#!/bin/bash",
				"#PBS -l select=1:ncpus=4:mem=16gb",
			},
			wantNodes:        1,
			wantTasksPerNode: 1,
			wantCpusPerTask:  4,
		},
		{
			name: "defaults no resource directives",
			lines: []string{
				"#!/bin/bash",
				"#PBS -N testjob",
			},
			wantNodes:        1,
			wantTasksPerNode: 0,
			wantCpusPerTask:  1, // default
		},
		{
			name: "select with mpiprocs and ngpus",
			lines: []string{
				"#!/bin/bash",
				"#PBS -l select=2:ncpus=16:mpiprocs=4:mem=32gb:ngpus=1",
			},
			wantNodes:        2,
			wantTasksPerNode: 4,
			wantCpusPerTask:  4, // ncpus(16) / mpiprocs(4) = 4 threads per task
		},
		{
			name: "ncpus not divisible by mpiprocs: passthrough",
			lines: []string{
				"#!/bin/bash",
				"#PBS -l select=4:ncpus=7:mpiprocs=2:mem=8gb",
			},
			wantPassthrough: true, // 7 % 2 != 0 → passthrough
		},
	}

	for _, tt := range tests {
		t.Run(tt.name, func(t *testing.T) {
			tmpDir := t.TempDir()
			scriptPath := filepath.Join(tmpDir, "test.sh")
			content := strings.Join(tt.lines, "\n")
			if err := os.WriteFile(scriptPath, []byte(content), 0644); err != nil {
				t.Fatalf("Failed to create test script: %v", err)
			}

			pbs := newTestPbsScheduler()
			specs, err := pbs.ReadScriptSpecs(scriptPath)
			if err != nil {
				t.Fatalf("Failed to parse script: %v", err)
			}

			if tt.wantPassthrough {
				if specs.Spec != nil {
					t.Errorf("Spec = %+v; want nil (passthrough)", specs.Spec)
				}
				return
			}

			if specs.Spec == nil {
				t.Fatal("Spec is nil; want non-nil")
			}

			if specs.Spec.Nodes != tt.wantNodes {
				t.Errorf("Nodes = %d; want %d", specs.Spec.Nodes, tt.wantNodes)
			}
			if specs.Spec.TasksPerNode != tt.wantTasksPerNode {
				t.Errorf("TasksPerNode = %d; want %d", specs.Spec.TasksPerNode, tt.wantTasksPerNode)
			}
			if specs.Spec.CpusPerTask != tt.wantCpusPerTask {
				t.Errorf("CpusPerTask = %d; want %d", specs.Spec.CpusPerTask, tt.wantCpusPerTask)
			}
		})
	}
}

func TestPbsMultiChunkParsing(t *testing.T) {
	tests := []struct {
		name             string
		lines            []string
		wantCpusPerTask  int
		wantMemPerCpuMB  int64
		wantNodes        int
		wantNtasks       int
		wantTasksPerNode int // Ceiling division for cross-scheduler translation
		wantPassthrough  bool
	}{
		{
			name: "uniform mpiprocs uniform mem",
			lines: []string{
				"#!/bin/bash",
				"#PBS -l select=2:ncpus=8:mpiprocs=2:mem=16gb+2:ncpus=8:mpiprocs=2:mem=16gb",
			},
			// CpusPerTask: 8/2=4; MemPerCpuMB: 16*1024/8=2048; Nodes: 2+2=4; Ntasks: 4+4=8
			// TasksPerNode: (8+4-1)/4 = 2 (ceiling = exact for uniform)
			wantCpusPerTask:  4,
			wantMemPerCpuMB:  2048,
			wantNodes:        4,
			wantNtasks:       8,
			wantTasksPerNode: 2,
		},
		{
			name: "uniform mpiprocs no mem",
			lines: []string{
				"#!/bin/bash",
				"#PBS -l select=2:ncpus=8:mpiprocs=2+2:ncpus=8:mpiprocs=2",
			},
			// TasksPerNode: (8+4-1)/4 = 2
			wantCpusPerTask:  4,
			wantNodes:        4,
			wantNtasks:       8,
			wantTasksPerNode: 2,
		},
		{
			name: "non-uniform mpiprocs same CpusPerTask no mem",
			lines: []string{
				"#!/bin/bash",
				"#PBS -l select=2:ncpus=4:mpiprocs=2+1:ncpus=8:mpiprocs=4",
			},
			// CpusPerTask: 4/2=2 and 8/4=2 (uniform); Nodes: 2+1=3; Ntasks: 4+4=8
			// TasksPerNode: (8+3-1)/3 = 3 (ceiling, non-uniform)
			wantCpusPerTask:  2,
			wantNodes:        3,
			wantNtasks:       8,
			wantTasksPerNode: 3,
		},
		{
			name: "non-uniform mem per cpu: passthrough",
			lines: []string{
				"#!/bin/bash",
				"#PBS -l select=2:ncpus=8:mpiprocs=2:mem=16gb+2:ncpus=8:mpiprocs=2:mem=32gb",
			},
			// MemPerCpu: 16*1024/8=2048 vs 32*1024/8=4096 → non-uniform → passthrough
			wantPassthrough: true,
		},
		{
			name: "non-uniform CpusPerTask: passthrough",
			lines: []string{
				"#!/bin/bash",
				"#PBS -l select=1:ncpus=8:mpiprocs=2+1:ncpus=2:mpiprocs=1",
			},
			// CpusPerTask: 8/2=4 vs 2/1=2 (non-uniform) → passthrough
			wantPassthrough: true,
		},
		{
			name: "GPU in multi-chunk: passthrough",
			lines: []string{
				"#!/bin/bash",
				"#PBS -l select=1:ncpus=4:ngpus=2+1:ncpus=4:ngpus=2",
			},
			// ngpus= inside multi-chunk select is not supported → passthrough
			wantPassthrough: true,
		},
		{
			name: "three-chunk uniform",
			lines: []string{
				"#!/bin/bash",
				"#PBS -l select=1:ncpus=8:mpiprocs=2+2:ncpus=8:mpiprocs=2+1:ncpus=8:mpiprocs=2",
			},
			// Nodes: 1+2+1=4; Ntasks: 2+4+2=8
			// TasksPerNode: (8+4-1)/4 = 2 (ceiling = exact for uniform)
			wantCpusPerTask:  4,
			wantNodes:        4,
			wantNtasks:       8,
			wantTasksPerNode: 2,
		},
	}

	for _, tt := range tests {
		t.Run(tt.name, func(t *testing.T) {
			tmpDir := t.TempDir()
			scriptPath := filepath.Join(tmpDir, "test.sh")
			content := strings.Join(tt.lines, "\n")
			if err := os.WriteFile(scriptPath, []byte(content), 0644); err != nil {
				t.Fatalf("Failed to create test script: %v", err)
			}

			pbs := newTestPbsScheduler()
			specs, err := pbs.ReadScriptSpecs(scriptPath)
			if err != nil {
				t.Fatalf("Failed to parse script: %v", err)
			}

			if tt.wantPassthrough {
				if specs.Spec != nil {
					t.Errorf("Spec = %+v; want nil (passthrough)", specs.Spec)
				}
				return
			}

			if specs.Spec == nil {
				t.Fatal("Spec is nil; want non-nil")
			}

			if specs.Spec.CpusPerTask != tt.wantCpusPerTask {
				t.Errorf("CpusPerTask = %d; want %d", specs.Spec.CpusPerTask, tt.wantCpusPerTask)
			}
			if specs.Spec.MemPerCpuMB != tt.wantMemPerCpuMB {
				t.Errorf("MemPerCpuMB = %d; want %d", specs.Spec.MemPerCpuMB, tt.wantMemPerCpuMB)
			}
			if specs.Spec.Nodes != tt.wantNodes {
				t.Errorf("Nodes = %d; want %d", specs.Spec.Nodes, tt.wantNodes)
			}
			if specs.Spec.Ntasks != tt.wantNtasks {
				t.Errorf("Ntasks = %d; want %d", specs.Spec.Ntasks, tt.wantNtasks)
			}
			if specs.Spec.TasksPerNode != tt.wantTasksPerNode {
				t.Errorf("TasksPerNode = %d; want %d (ceiling for cross-scheduler translation)", specs.Spec.TasksPerNode, tt.wantTasksPerNode)
			}
		})
	}
}

func TestPbsGetJobResources(t *testing.T) {
	sched := &PbsScheduler{}

	t.Run("not in job", func(t *testing.T) {
		clearJobEnvVars(t)
		if res := sched.GetJobResources(); res != nil {
			t.Fatalf("expected nil, got %+v", res)
		}
	})

	t.Run("cpus and gpus", func(t *testing.T) {
		clearJobEnvVars(t)
		t.Setenv("PBS_JOBID", "67890.pbs-server")
		t.Setenv("NCPUS", "8")
		t.Setenv("CUDA_VISIBLE_DEVICES", "0,1,2")

		res := sched.GetJobResources()
		if res == nil {
			t.Fatal("expected non-nil")
		}
		if res.CpusPerTask != 8 {
			t.Errorf("CpusPerTask = %d; want 8", res.CpusPerTask)
		}
		if res.Gpu == nil || res.Gpu.Count != 3 {
			t.Errorf("Gpu.Count = %v; want 3", res.Gpu)
		}
	})

	t.Run("nodes and tasks", func(t *testing.T) {
		clearJobEnvVars(t)
		t.Setenv("PBS_JOBID", "67890.pbs-server")
		t.Setenv("NCPUS", "2") // cpus per task (ompthreads per chunk)
		t.Setenv("PBS_NUM_NODES", "4")
		t.Setenv("PBS_NP", "16")

		res := sched.GetJobResources()
		if res == nil {
			t.Fatal("expected non-nil")
		}
		if res.Nodes != 4 {
			t.Errorf("Nodes = %d; want 4", res.Nodes)
		}
		if res.Ntasks != 16 {
			t.Errorf("Ntasks = %d; want 16", res.Ntasks)
		}
		// TasksPerNode = PBS_NP / PBS_NUM_NODES = 16 / 4 = 4
		if res.TasksPerNode != 4 {
			t.Errorf("TasksPerNode = %d; want 4 (derived from PBS_NP/PBS_NUM_NODES)", res.TasksPerNode)
		}
		// CpusPerTask = NCPUS = 2
		if res.CpusPerTask != 2 {
			t.Errorf("CpusPerTask = %d; want 2", res.CpusPerTask)
		}
	})

	t.Run("PBS_TASKNUM fallback", func(t *testing.T) {
		clearJobEnvVars(t)
		t.Setenv("PBS_JOBID", "67890")
		t.Setenv("PBS_NUM_NODES", "2")
		t.Setenv("PBS_TASKNUM", "8")

		res := sched.GetJobResources()
		if res == nil {
			t.Fatal("expected non-nil")
		}
		if res.Ntasks != 8 {
			t.Errorf("Ntasks = %d; want 8", res.Ntasks)
		}
		// TasksPerNode = PBS_TASKNUM / PBS_NUM_NODES = 8 / 2 = 4
		if res.TasksPerNode != 4 {
			t.Errorf("TasksPerNode = %d; want 4 (derived from PBS_TASKNUM/PBS_NUM_NODES)", res.TasksPerNode)
		}
	})

	t.Run("NCPUS direct", func(t *testing.T) {
		clearJobEnvVars(t)
		t.Setenv("PBS_JOBID", "67890")
		t.Setenv("NCPUS", "12")

		res := sched.GetJobResources()
		if res == nil {
			t.Fatal("expected non-nil")
		}
		if res.CpusPerTask != 12 {
			t.Errorf("CpusPerTask = %d; want 12", res.CpusPerTask)
		}
	})

	t.Run("memory from PBS_VMEM", func(t *testing.T) {
		clearJobEnvVars(t)
		t.Setenv("PBS_JOBID", "67890")
		t.Setenv("PBS_VMEM", "8589934592") // 8 GB in bytes

		res := sched.GetJobResources()
		if res == nil {
			t.Fatal("expected non-nil")
		}
		if res.MemPerNodeMB != 8192 {
			t.Errorf("MemPerNodeMB = %d; want 8192", res.MemPerNodeMB)
		}
	})

	t.Run("partial data", func(t *testing.T) {
		clearJobEnvVars(t)
		t.Setenv("PBS_JOBID", "67890")
		t.Setenv("NCPUS", "4")

		res := sched.GetJobResources()
		if res == nil {
			t.Fatal("expected non-nil")
		}
		if res.CpusPerTask != 4 {
			t.Errorf("CpusPerTask = %d; want 4", res.CpusPerTask)
		}
		if res.MemPerNodeMB != 0 {
			t.Errorf("MemPerNodeMB should be 0 (not set), got %d", res.MemPerNodeMB)
		}
		if res.Gpu != nil {
			t.Errorf("Gpu should be nil, got %+v", res.Gpu)
		}
	})

	t.Run("invalid values", func(t *testing.T) {
		clearJobEnvVars(t)
		t.Setenv("PBS_JOBID", "67890")
		t.Setenv("NCPUS", "not-a-number")
		t.Setenv("PBS_VMEM", "-100")

		res := sched.GetJobResources()
		if res == nil {
			t.Fatal("expected non-nil")
		}
		if res.CpusPerTask != 0 {
			t.Errorf("CpusPerTask should be 0 for invalid value, got %d", res.CpusPerTask)
		}
		if res.MemPerNodeMB != 0 {
			t.Errorf("MemPerNodeMB should be 0 for negative value, got %d", res.MemPerNodeMB)
		}
	})
}

// TestPbsOldStyleResourcePassthrough verifies that Torque-era standalone -l
// resource directives (ncpus=, mem=, ngpus=, gpus=) are rejected and trigger
// full passthrough (Spec == nil), consistent with how nodes= is handled.
func TestPbsOldStyleResourcePassthrough(t *testing.T) {
	tests := []struct {
		name  string
		lines []string
	}{
		{
			name: "standalone ncpus",
			lines: []string{
				"#!/bin/bash",
				"#PBS -l ncpus=8",
			},
		},
		{
			name: "standalone mem",
			lines: []string{
				"#!/bin/bash",
				"#PBS -l mem=16gb",
			},
		},
		{
			name: "standalone ngpus",
			lines: []string{
				"#!/bin/bash",
				"#PBS -l ngpus=2",
			},
		},
		{
			name: "standalone gpus",
			lines: []string{
				"#!/bin/bash",
				"#PBS -l gpus=2",
			},
		},
		{
			name: "mixed: valid select then old-style mem",
			lines: []string{
				"#!/bin/bash",
				"#PBS -l select=1:ncpus=4",
				"#PBS -l mem=8gb",
			},
		},
	}

	for _, tt := range tests {
		t.Run(tt.name, func(t *testing.T) {
			tmpDir := t.TempDir()
			scriptPath := filepath.Join(tmpDir, "test.sh")
			content := strings.Join(tt.lines, "\n")
			if err := os.WriteFile(scriptPath, []byte(content), 0644); err != nil {
				t.Fatalf("Failed to create test script: %v", err)
			}

			pbs := newTestPbsScheduler()
			specs, err := pbs.ReadScriptSpecs(scriptPath)
			if err != nil {
				t.Fatalf("ReadScriptSpecs failed: %v", err)
			}
			if specs.Spec != nil {
				t.Errorf("Spec = %+v; want nil (expected passthrough for old Torque syntax)", specs.Spec)
			}
		})
	}
}

// TestPbsNonUniformMpiScriptGen tests CreateScriptWithSpec for non-uniform MPI distributions
// where Ntasks % Nodes != 0. The scheduler should emit two select= chunks.
func TestPbsNonUniformMpiScriptGen(t *testing.T) {
	tests := []struct {
		name       string
		spec       *ResourceSpec
		wantSelect string // expected substring in #PBS -l line
		wantErr    bool
	}{
		{
			name: "non-uniform with MemPerNodeMB: two chunks same mem",
			spec: &ResourceSpec{
				Nodes:        4,
				Ntasks:       14,
				CpusPerTask:  2,
				MemPerNodeMB: 16384,
			},
			// ceil(14/4)=4 → nFull=3, rem=2, remNodes=1
			// chunk1: 3:ncpus=8:mpiprocs=4:mem=16384mb
			// chunk2: 1:ncpus=4:mpiprocs=2:mem=16384mb
			wantSelect: "select=3:ncpus=8:mpiprocs=4:mem=16384mb+1:ncpus=4:mpiprocs=2:mem=16384mb",
		},
		{
			name: "non-uniform with MemPerCpuMB: derived mem differs per chunk",
			spec: &ResourceSpec{
				Nodes:       4,
				Ntasks:      14,
				CpusPerTask: 2,
				MemPerCpuMB: 2048,
			},
			// full chunk: mem=2048*8=16384mb; remainder: mem=2048*4=8192mb
			wantSelect: "select=3:ncpus=8:mpiprocs=4:mem=16384mb+1:ncpus=4:mpiprocs=2:mem=8192mb",
		},
		{
			name: "non-uniform no mem: two chunks without mem",
			spec: &ResourceSpec{
				Nodes:       4,
				Ntasks:      14,
				CpusPerTask: 2,
			},
			wantSelect: "select=3:ncpus=8:mpiprocs=4+1:ncpus=4:mpiprocs=2",
		},
		{
			name: "uniform Ntasks divisible by Nodes: single chunk (inferred tasksPerNode)",
			spec: &ResourceSpec{
				Nodes:        4,
				Ntasks:       12,
				CpusPerTask:  2,
				MemPerNodeMB: 8192,
			},
			// 12%4==0 → tasksPerNode=3 inferred
			// single chunk: select=4:ncpus=6:mpiprocs=3:mem=8192mb
			wantSelect: "select=4:ncpus=6:mpiprocs=3:mem=8192mb",
		},
		{
			name: "GPU with non-uniform Ntasks: error",
			spec: &ResourceSpec{
				Nodes:       4,
				Ntasks:      14,
				CpusPerTask: 2,
				Gpu:         &GpuSpec{Type: "gpu", Count: 1},
			},
			wantErr: true,
		},
	}

	for _, tt := range tests {
		t.Run(tt.name, func(t *testing.T) {
			tmpDir := t.TempDir()
			pbs := newTestPbsScheduler()

			jobSpec := &JobSpec{
				Name:    "mpi_job",
				Command: "mpirun ./app",
				Specs: &ScriptSpecs{
					Spec:           tt.spec,
					RemainingFlags: []string{},
				},
				Metadata: map[string]string{},
			}

			scriptPath, err := pbs.CreateScriptWithSpec(jobSpec, tmpDir)
			if tt.wantErr {
				if err == nil {
					t.Errorf("CreateScriptWithSpec succeeded; want error")
				}
				return
			}
			if err != nil {
				t.Fatalf("CreateScriptWithSpec failed: %v", err)
			}

			content, err := os.ReadFile(scriptPath)
			if err != nil {
				t.Fatalf("Failed to read generated script: %v", err)
			}
			scriptContent := string(content)

			if !strings.Contains(scriptContent, tt.wantSelect) {
				t.Errorf("Script missing expected select directive:\n  want: %q\n  script:\n%s", tt.wantSelect, scriptContent)
			}
		})
	}
}
