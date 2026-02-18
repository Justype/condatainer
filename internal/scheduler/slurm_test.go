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

// newTestSlurmScheduler creates a SLURM scheduler instance for testing
// without requiring sbatch to be installed
func newTestSlurmScheduler() *SlurmScheduler {
	return &SlurmScheduler{
		sbatchBin:       "/usr/bin/sbatch",   // fake path for testing
		sinfoCommand:    "/usr/bin/sinfo",    // fake path for testing
		scontrolCommand: "/usr/bin/scontrol", // fake path for testing
		directiveRe:     regexp.MustCompile(`^#SBATCH\s+(.+)`),
		jobIDRe:         regexp.MustCompile(`Submitted batch job (\d+)`),
	}
}

func TestSlurmEmailParsing(t *testing.T) {
	tests := []struct {
		name         string
		scriptLines  []string
		wantBegin    bool
		wantEnd      bool
		wantFail     bool
		wantMailUser string
	}{
		{
			name: "BEGIN,END,FAIL",
			scriptLines: []string{
				"#!/bin/bash",
				"#SBATCH --mail-type=BEGIN,END,FAIL",
				"#SBATCH --mail-user=test@example.com",
			},
			wantBegin:    true,
			wantEnd:      true,
			wantFail:     true,
			wantMailUser: "test@example.com",
		},
		{
			name: "ALL type",
			scriptLines: []string{
				"#!/bin/bash",
				"#SBATCH --mail-type=ALL",
				"#SBATCH --mail-user=user@domain.org",
			},
			wantBegin:    true,
			wantEnd:      true,
			wantFail:     true,
			wantMailUser: "user@domain.org",
		},
		{
			name: "Only BEGIN",
			scriptLines: []string{
				"#!/bin/bash",
				"#SBATCH --mail-type=BEGIN",
			},
			wantBegin:    true,
			wantEnd:      false,
			wantFail:     false,
			wantMailUser: "",
		},
		{
			name: "Only END",
			scriptLines: []string{
				"#!/bin/bash",
				"#SBATCH --mail-type=END",
			},
			wantBegin:    false,
			wantEnd:      true,
			wantFail:     false,
			wantMailUser: "",
		},
		{
			name: "FAIL and REQUEUE",
			scriptLines: []string{
				"#!/bin/bash",
				"#SBATCH --mail-type=FAIL,REQUEUE",
			},
			wantBegin: false,
			wantEnd:   false,
			wantFail:  true,
		},
		{
			name: "NONE disables all",
			scriptLines: []string{
				"#!/bin/bash",
				"#SBATCH --mail-type=NONE",
			},
			wantBegin: false,
			wantEnd:   false,
			wantFail:  false,
		},
		{
			name: "Case insensitive",
			scriptLines: []string{
				"#!/bin/bash",
				"#SBATCH --mail-type=begin,end",
				"#SBATCH --mail-user=Test@Example.COM",
			},
			wantBegin:    true,
			wantEnd:      true,
			wantFail:     false,
			wantMailUser: "Test@Example.COM",
		},
		{
			name: "Advanced types ignored (TIME_LIMIT, ARRAY_TASKS)",
			scriptLines: []string{
				"#!/bin/bash",
				"#SBATCH --mail-type=TIME_LIMIT_80,ARRAY_TASKS,END",
			},
			wantBegin: false,
			wantEnd:   true,
			wantFail:  false,
		},
		{
			name: "No email directives",
			scriptLines: []string{
				"#!/bin/bash",
				"#SBATCH --job-name=test",
			},
			wantBegin: false,
			wantEnd:   false,
			wantFail:  false,
		},
	}

	for _, tt := range tests {
		t.Run(tt.name, func(t *testing.T) {
			// Create temporary test script
			tmpDir := t.TempDir()
			scriptPath := filepath.Join(tmpDir, "test.sh")
			content := strings.Join(tt.scriptLines, "\n")
			if err := os.WriteFile(scriptPath, []byte(content), 0644); err != nil {
				t.Fatalf("Failed to create test script: %v", err)
			}

			// Parse the script
			slurm := newTestSlurmScheduler()

			specs, err := slurm.ReadScriptSpecs(scriptPath)
			if err != nil {
				t.Fatalf("Failed to parse script: %v", err)
			}

			// Verify results
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

func TestSlurmEmailScriptGeneration(t *testing.T) {
	tests := []struct {
		name         string
		emailOnBegin bool
		emailOnEnd   bool
		emailOnFail  bool
		mailUser     string
		wantMailType string // expected mail-type value (empty = not present)
		wantMailUser string // expected mail-user value (empty = not present)
	}{
		{
			name:         "All three notifications",
			emailOnBegin: true,
			emailOnEnd:   true,
			emailOnFail:  true,
			mailUser:     "test@example.com",
			wantMailType: "BEGIN,END,FAIL",
			wantMailUser: "test@example.com",
		},
		{
			name:         "Only BEGIN and END",
			emailOnBegin: true,
			emailOnEnd:   true,
			emailOnFail:  false,
			mailUser:     "user@domain.org",
			wantMailType: "BEGIN,END",
			wantMailUser: "user@domain.org",
		},
		{
			name:         "Only END",
			emailOnBegin: false,
			emailOnEnd:   true,
			emailOnFail:  false,
			wantMailType: "END",
			wantMailUser: "",
		},
		{
			name:         "Only FAIL",
			emailOnBegin: false,
			emailOnEnd:   false,
			emailOnFail:  true,
			mailUser:     "admin@company.com",
			wantMailType: "FAIL",
			wantMailUser: "admin@company.com",
		},
		{
			name:         "No email notifications",
			emailOnBegin: false,
			emailOnEnd:   false,
			emailOnFail:  false,
			wantMailType: "",
			wantMailUser: "",
		},
	}

	for _, tt := range tests {
		t.Run(tt.name, func(t *testing.T) {
			tmpDir := t.TempDir()

			slurm := newTestSlurmScheduler()

			jobSpec := &JobSpec{
				Name:    "test_job",
				Command: "echo 'test'",
				Specs: &ScriptSpecs{
					Spec: &ResourceSpec{
						CpusPerTask:  4,
						MemPerNodeMB: 8000,
						Time:         time.Hour,
					},
					Control: RuntimeConfig{
						JobName:      "test_job",
						EmailOnBegin: tt.emailOnBegin,
						EmailOnEnd:   tt.emailOnEnd,
						EmailOnFail:  tt.emailOnFail,
						MailUser:     tt.mailUser,
					},
				},
			}

			scriptPath, err := slurm.CreateScriptWithSpec(jobSpec, tmpDir)
			if err != nil {
				t.Fatalf("Failed to create script: %v", err)
			}

			// Read the generated script
			content, err := os.ReadFile(scriptPath)
			if err != nil {
				t.Fatalf("Failed to read generated script: %v", err)
			}
			scriptContent := string(content)

			// Check mail-type directive
			if tt.wantMailType != "" {
				expectedLine := "#SBATCH --mail-type=" + tt.wantMailType
				if !strings.Contains(scriptContent, expectedLine) {
					t.Errorf("Script missing expected line: %q\nScript content:\n%s", expectedLine, scriptContent)
				}
			} else {
				// Should not contain mail-type
				if strings.Contains(scriptContent, "--mail-type=") {
					t.Errorf("Script should not contain --mail-type directive\nScript content:\n%s", scriptContent)
				}
			}

			// Check mail-user directive
			if tt.wantMailUser != "" {
				expectedLine := "#SBATCH --mail-user=" + tt.wantMailUser
				if !strings.Contains(scriptContent, expectedLine) {
					t.Errorf("Script missing expected line: %q\nScript content:\n%s", expectedLine, scriptContent)
				}
			} else {
				// Should not contain mail-user
				if strings.Contains(scriptContent, "--mail-user=") {
					t.Errorf("Script should not contain --mail-user directive\nScript content:\n%s", scriptContent)
				}
			}
		})
	}
}

func TestSlurmEmailRoundTrip(t *testing.T) {
	// Test that parsing and generating produces consistent results
	tmpDir := t.TempDir()

	// Create a script with email directives
	originalScript := `#!/bin/bash
#SBATCH --job-name=email_test
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --time=02:00:00
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=roundtrip@example.com

echo "Running job"
`
	scriptPath := filepath.Join(tmpDir, "original.sh")
	if err := os.WriteFile(scriptPath, []byte(originalScript), 0644); err != nil {
		t.Fatalf("Failed to create test script: %v", err)
	}

	slurm := newTestSlurmScheduler()

	// Parse the original script
	specs, err := slurm.ReadScriptSpecs(scriptPath)
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

	// Generate a new script from the parsed specs
	jobSpec := &JobSpec{
		Name:    "email_test",
		Command: "echo 'Running job'",
		Specs:   specs,
	}

	newScriptPath, err := slurm.CreateScriptWithSpec(jobSpec, tmpDir)
	if err != nil {
		t.Fatalf("Failed to create script: %v", err)
	}

	// Read the generated script
	content, err := os.ReadFile(newScriptPath)
	if err != nil {
		t.Fatalf("Failed to read generated script: %v", err)
	}
	scriptContent := string(content)

	// Verify the generated script contains the email directives
	if !strings.Contains(scriptContent, "#SBATCH --mail-type=BEGIN,END,FAIL") {
		t.Errorf("Generated script missing mail-type directive:\n%s", scriptContent)
	}
	if !strings.Contains(scriptContent, "#SBATCH --mail-user=roundtrip@example.com") {
		t.Errorf("Generated script missing mail-user directive:\n%s", scriptContent)
	}
	if !strings.Contains(scriptContent, "#SBATCH --time=02:00:00") {
		t.Errorf("Generated script missing time directive:\n%s", scriptContent)
	}

	// Verify NO duplicates - count occurrences
	mailTypeCount := strings.Count(scriptContent, "--mail-type=")
	mailUserCount := strings.Count(scriptContent, "--mail-user=")
	timeCount := strings.Count(scriptContent, "--time=")
	if mailTypeCount != 1 {
		t.Errorf("Expected exactly 1 --mail-type directive, found %d:\n%s", mailTypeCount, scriptContent)
	}
	if mailUserCount != 1 {
		t.Errorf("Expected exactly 1 --mail-user directive, found %d:\n%s", mailUserCount, scriptContent)
	}
	if timeCount != 1 {
		t.Errorf("Expected exactly 1 --time directive, found %d:\n%s", timeCount, scriptContent)
	}
}

// Test that scheduler uses outputDir for log files
func TestSlurmCreateScriptUsesOutputDir(t *testing.T) {
	tmpOutputDir := t.TempDir()

	slurm := newTestSlurmScheduler()

	jobSpec := &JobSpec{
		Name:    "test/job",
		Command: "echo 'hello'",
		Specs: &ScriptSpecs{
			RemainingFlags: []string{},
			Spec: &ResourceSpec{
				CpusPerTask: 1,
			},
		},
	}

	scriptPath, err := slurm.CreateScriptWithSpec(jobSpec, tmpOutputDir)
	if err != nil {
		t.Fatalf("CreateScriptWithSpec failed: %v", err)
	}

	// Ensure the generated script uses outputDir for logs
	logPath := filepath.Join(tmpOutputDir, "test--job.log")
	content, err := os.ReadFile(scriptPath)
	if err != nil {
		t.Fatalf("Failed to read generated script: %v", err)
	}
	if !strings.Contains(string(content), fmt.Sprintf("--output=%s", logPath)) {
		t.Errorf("Generated script does not contain expected output path %s\nScript:\n%s", logPath, string(content))
	}
}

func TestSlurmResourceParsing(t *testing.T) {
	tests := []struct {
		name      string
		lines     []string
		wantCpus  int
		wantMemMB int64
		wantTime  time.Duration
		wantGpu   *GpuSpec
	}{
		{
			name: "basic resources",
			lines: []string{
				"#!/bin/bash",
				"#SBATCH --job-name=testjob",
				"#SBATCH --cpus-per-task=8",
				"#SBATCH --mem=16G",
				"#SBATCH --time=02:00:00",
			},
			wantCpus:  8,
			wantMemMB: 16 * 1024,
			wantTime:  2 * time.Hour,
		},
		{
			name: "short flags",
			lines: []string{
				"#!/bin/bash",
				"#SBATCH -J testjob",
				"#SBATCH -c 4",
				"#SBATCH -t 01:30:00",
			},
			wantCpus: 4,
			wantTime: time.Hour + 30*time.Minute,
		},
		{
			name: "memory in MB",
			lines: []string{
				"#!/bin/bash",
				"#SBATCH --mem=4096M",
			},
			wantCpus:  2, // default
			wantMemMB: 4096,
		},
		{
			name: "GPU gres",
			lines: []string{
				"#!/bin/bash",
				"#SBATCH --cpus-per-task=8",
				"#SBATCH --gres=gpu:a100:2",
			},
			wantCpus: 8,
			wantGpu:  &GpuSpec{Type: "a100", Count: 2},
		},
		{
			name: "GPU gpus flag",
			lines: []string{
				"#!/bin/bash",
				"#SBATCH --gpus=v100:4",
			},
			wantCpus: 2, // default
			wantGpu:  &GpuSpec{Type: "v100", Count: 4},
		},
		{
			name: "GPU gpus-per-node",
			lines: []string{
				"#!/bin/bash",
				"#SBATCH --gpus-per-node=2",
			},
			wantCpus: 2, // default
			wantGpu:  &GpuSpec{Type: "gpu", Count: 2},
		},
		{
			name: "time with days",
			lines: []string{
				"#!/bin/bash",
				"#SBATCH --time=1-12:00:00",
			},
			wantCpus: 2, // default
			wantTime: 36 * time.Hour,
		},
		{
			name: "full job script",
			lines: []string{
				"#!/bin/bash",
				"#SBATCH --job-name=myjob",
				"#SBATCH --cpus-per-task=16",
				"#SBATCH --mem=64G",
				"#SBATCH --time=12:00:00",
				"#SBATCH --gres=gpu:h100:4",
				"#SBATCH --mail-type=ALL",
				"#SBATCH --mail-user=user@example.com",
			},
			wantCpus:  16,
			wantMemMB: 64 * 1024,
			wantTime:  12 * time.Hour,
			wantGpu:   &GpuSpec{Type: "h100", Count: 4},
		},
		{
			name: "no resource directives (defaults)",
			lines: []string{
				"#!/bin/bash",
				"echo hello",
			},
			wantCpus: 2, // default
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

			slurm := newTestSlurmScheduler()
			specs, err := slurm.ReadScriptSpecs(scriptPath)
			if err != nil {
				t.Fatalf("Failed to parse script: %v", err)
			}

			if specs.Spec == nil {
				t.Fatal("Spec is nil; want non-nil")
			}

			if specs.Spec.CpusPerTask != tt.wantCpus {
				t.Errorf("CpusPerTask = %d; want %d", specs.Spec.CpusPerTask, tt.wantCpus)
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

func TestSlurmTimeParsing(t *testing.T) {
	tests := []struct {
		name    string
		input   string
		wantDur time.Duration
		wantErr bool
	}{
		{"HH:MM:SS", "02:30:00", 2*time.Hour + 30*time.Minute, false},
		{"HH:MM", "10:30", 10*time.Hour + 30*time.Minute, false},
		{"minutes only", "90", 90 * time.Minute, false},
		{"with days", "1-12:00:00", 36 * time.Hour, false},
		{"with seconds", "01:00:30", time.Hour + 30*time.Second, false},
		{"empty string", "", 0, false},
		{"large walltime", "7-00:00:00", 168 * time.Hour, false},
	}

	for _, tt := range tests {
		t.Run(tt.name, func(t *testing.T) {
			dur, err := parseSlurmTimeSpec(tt.input)
			if tt.wantErr && err == nil {
				t.Error("Expected error, got nil")
			}
			if !tt.wantErr && err != nil {
				t.Errorf("Unexpected error: %v", err)
			}
			if dur != tt.wantDur {
				t.Errorf("parseSlurmTimeSpec(%q) = %v; want %v", tt.input, dur, tt.wantDur)
			}
		})
	}
}

func TestSlurmMemoryParsing(t *testing.T) {
	tests := []struct {
		input  string
		wantMB int64
	}{
		{"8G", 8 * 1024},
		{"8GB", 8 * 1024},
		{"1024M", 1024},
		{"1024MB", 1024},
		{"4096K", 4},
		{"4096KB", 4},
		{"1T", 1024 * 1024},
		{"1TB", 1024 * 1024},
	}

	for _, tt := range tests {
		t.Run(tt.input, func(t *testing.T) {
			mb, err := parseMemory(tt.input)
			if err != nil {
				t.Errorf("parseMemory(%q) error: %v", tt.input, err)
				return
			}
			if mb != tt.wantMB {
				t.Errorf("parseMemory(%q) = %d MB; want %d MB", tt.input, mb, tt.wantMB)
			}
		})
	}
}

func TestSlurmIsAvailable(t *testing.T) {
	t.Run("not in job", func(t *testing.T) {
		clearJobEnvVars(t)
		slurm := newTestSlurmScheduler()
		if !slurm.IsAvailable() {
			t.Error("Expected IsAvailable to return true when not in a job")
		}
	})

	t.Run("inside job", func(t *testing.T) {
		clearJobEnvVars(t)
		t.Setenv("SLURM_JOB_ID", "12345")
		slurm := newTestSlurmScheduler()
		if slurm.IsAvailable() {
			t.Error("Expected IsAvailable to return false when inside a job")
		}
	})

	t.Run("no binary", func(t *testing.T) {
		clearJobEnvVars(t)
		slurm := &SlurmScheduler{}
		if slurm.IsAvailable() {
			t.Error("Expected IsAvailable to return false when no binary is set")
		}
	})
}

func TestSlurmGetInfo(t *testing.T) {
	clearJobEnvVars(t)
	slurm := newTestSlurmScheduler()

	info := slurm.GetInfo()
	if info.Type != "SLURM" {
		t.Errorf("Type = %q; want %q", info.Type, "SLURM")
	}
	if info.Binary != "/usr/bin/sbatch" {
		t.Errorf("Binary = %q; want %q", info.Binary, "/usr/bin/sbatch")
	}
	if info.InJob {
		t.Error("InJob should be false")
	}
	if !info.Available {
		t.Error("Available should be true")
	}
}

func TestTryParseSlurmScript(t *testing.T) {
	tmpDir := t.TempDir()

	script := `#!/bin/bash
#SBATCH --job-name=testjob
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --time=02:00:00
#SBATCH --gres=gpu:a100:1
#SBATCH --mail-type=END
#SBATCH --mail-user=user@example.com

echo "Running job"
`
	scriptPath := filepath.Join(tmpDir, "test.sh")
	if err := os.WriteFile(scriptPath, []byte(script), 0644); err != nil {
		t.Fatalf("Failed to create test script: %v", err)
	}

	specs, err := TryParseSlurmScript(scriptPath)
	if err != nil {
		t.Fatalf("TryParseSlurmScript failed: %v", err)
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
	if specs.Spec.Gpu == nil || specs.Spec.Gpu.Count != 1 || specs.Spec.Gpu.Type != "a100" {
		t.Errorf("Gpu = %+v; want a100 x 1", specs.Spec.Gpu)
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

// clearJobEnvVars ensures no scheduler job env vars are set for test isolation.
// Shared across all scheduler test files (same package).
func clearJobEnvVars(t *testing.T) {
	t.Helper()
	for _, key := range []string{
		"SLURM_JOB_ID", "SLURM_CPUS_PER_TASK", "SLURM_MEM_PER_NODE",
		"SLURM_NTASKS", "SLURM_JOB_NUM_NODES",
		"PBS_JOBID", "PBS_NCPUS", "NCPUS", "PBS_VMEM",
		"PBS_NUM_NODES", "PBS_NP", "PBS_TASKNUM",
		"LSB_JOBID", "LSB_DJOB_NUMPROC", "LSB_MAX_NUM_PROCESSORS", "LSB_MAX_MEM_RUSAGE",
		"_CONDOR_JOB_AD", "_CONDOR_REQUEST_CPUS", "_CONDOR_REQUEST_MEMORY",
		"CUDA_VISIBLE_DEVICES",
	} {
		os.Unsetenv(key)
		t.Setenv(key, "")
		os.Unsetenv(key)
	}
}

func TestSlurmGetJobResources(t *testing.T) {
	sched := &SlurmScheduler{}

	t.Run("not in job", func(t *testing.T) {
		clearJobEnvVars(t)
		if res := sched.GetJobResources(); res != nil {
			t.Fatalf("expected nil, got %+v", res)
		}
	})

	t.Run("full resources", func(t *testing.T) {
		clearJobEnvVars(t)
		t.Setenv("SLURM_JOB_ID", "12345")
		t.Setenv("SLURM_CPUS_PER_TASK", "16")
		t.Setenv("SLURM_NTASKS", "4")
		t.Setenv("SLURM_JOB_NUM_NODES", "2")
		t.Setenv("SLURM_MEM_PER_NODE", "8192")
		t.Setenv("CUDA_VISIBLE_DEVICES", "0,1")

		res := sched.GetJobResources()
		if res == nil {
			t.Fatal("expected non-nil")
		}
		if res.Ncpus == nil || *res.Ncpus != 16 {
			t.Errorf("Ncpus = %v; want 16", res.Ncpus)
		}
		if res.Ntasks == nil || *res.Ntasks != 4 {
			t.Errorf("Ntasks = %v; want 4", res.Ntasks)
		}
		if res.Nodes == nil || *res.Nodes != 2 {
			t.Errorf("Nodes = %v; want 2", res.Nodes)
		}
		if res.MemMB == nil || *res.MemMB != 8192 {
			t.Errorf("MemMB = %v; want 8192", res.MemMB)
		}
		if res.Ngpus == nil || *res.Ngpus != 2 {
			t.Errorf("Ngpus = %v; want 2", res.Ngpus)
		}
	})

	t.Run("partial data", func(t *testing.T) {
		clearJobEnvVars(t)
		t.Setenv("SLURM_JOB_ID", "12345")
		t.Setenv("SLURM_CPUS_PER_TASK", "4")

		res := sched.GetJobResources()
		if res == nil {
			t.Fatal("expected non-nil")
		}
		if res.Ncpus == nil || *res.Ncpus != 4 {
			t.Errorf("Ncpus = %v; want 4", res.Ncpus)
		}
		if res.MemMB != nil {
			t.Errorf("MemMB should be nil, got %d", *res.MemMB)
		}
		if res.Ngpus != nil {
			t.Errorf("Ngpus should be nil, got %d", *res.Ngpus)
		}
	})

	t.Run("invalid values", func(t *testing.T) {
		clearJobEnvVars(t)
		t.Setenv("SLURM_JOB_ID", "12345")
		t.Setenv("SLURM_CPUS_PER_TASK", "not-a-number")
		t.Setenv("SLURM_MEM_PER_NODE", "-100")

		res := sched.GetJobResources()
		if res == nil {
			t.Fatal("expected non-nil")
		}
		if res.Ncpus != nil {
			t.Errorf("Ncpus should be nil for invalid value, got %d", *res.Ncpus)
		}
		if res.MemMB != nil {
			t.Errorf("MemMB should be nil for negative value, got %d", *res.MemMB)
		}
	})
}

func TestSlurmNodeTaskParsing(t *testing.T) {
	sched := newTestSlurmScheduler()

	tests := []struct {
		name             string
		lines            []string
		wantNodes        int
		wantTasksPerNode int
		wantCpus         int
	}{
		{
			name: "all specified",
			lines: []string{
				"#!/bin/bash",
				"#SBATCH --nodes=4",
				"#SBATCH --ntasks=16",
				"#SBATCH --cpus-per-task=4",
			},
			// ntasks=16, nodes=4 => TasksPerNode = 16/4 = 4
			wantNodes: 4, wantTasksPerNode: 4, wantCpus: 4,
		},
		{
			name: "short flags",
			lines: []string{
				"#!/bin/bash",
				"#SBATCH -N 3",
				"#SBATCH -n 12",
				"#SBATCH -c 2",
			},
			// ntasks=12, nodes=3 => TasksPerNode = 12/3 = 4
			wantNodes: 3, wantTasksPerNode: 4, wantCpus: 2,
		},
		{
			name: "ntasks-per-node",
			lines: []string{
				"#!/bin/bash",
				"#SBATCH --nodes=2",
				"#SBATCH --ntasks-per-node=8",
				"#SBATCH --cpus-per-task=4",
			},
			// ntasks-per-node=8 => TasksPerNode = 8
			wantNodes: 2, wantTasksPerNode: 8, wantCpus: 4,
		},
		{
			name: "ntasks overrides ntasks-per-node",
			lines: []string{
				"#!/bin/bash",
				"#SBATCH --nodes=2",
				"#SBATCH --ntasks=10",
				"#SBATCH --ntasks-per-node=8",
			},
			// ntasks=10 overrides, nodes=2 => TasksPerNode = 10/2 = 5
			wantNodes: 2, wantTasksPerNode: 5, wantCpus: 2,
		},
		{
			name: "ntasks only",
			lines: []string{
				"#!/bin/bash",
				"#SBATCH --ntasks=8",
			},
			// ntasks=8, nodes defaults to 1 => TasksPerNode = 8/1 = 8
			wantNodes: 1, wantTasksPerNode: 8, wantCpus: 2,
		},
		{
			name: "defaults when no node/task flags",
			lines: []string{
				"#!/bin/bash",
				"#SBATCH --mem=8G",
			},
			// defaults: nodes=1, tasks=1 => TasksPerNode = 1
			wantNodes: 1, wantTasksPerNode: 1, wantCpus: 2,
		},
	}

	for _, tt := range tests {
		t.Run(tt.name, func(t *testing.T) {
			tmpFile := filepath.Join(t.TempDir(), "test.sh")
			os.WriteFile(tmpFile, []byte(strings.Join(tt.lines, "\n")), 0644)
			specs, err := sched.ReadScriptSpecs(tmpFile)
			if err != nil {
				t.Fatalf("unexpected error: %v", err)
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
			if specs.Spec.CpusPerTask != tt.wantCpus {
				t.Errorf("CpusPerTask = %d; want %d", specs.Spec.CpusPerTask, tt.wantCpus)
			}
		})
	}
}

func TestSlurmNodeTaskScriptGeneration(t *testing.T) {
	sched := newTestSlurmScheduler()
	outputDir := t.TempDir()

	jobSpec := &JobSpec{
		Name:    "test-job",
		Command: "echo hello",
		Specs: &ScriptSpecs{
			Spec: &ResourceSpec{
				CpusPerTask:  4,
				TasksPerNode: 4, // 16 total tasks / 4 nodes = 4 tasks per node
				Nodes:        4,
				MemPerNodeMB: 8192,
				Time:         2 * time.Hour,
			},
			Control: RuntimeConfig{
				JobName: "test-job",
			},
		},
		Metadata: map[string]string{},
	}

	scriptPath, err := sched.CreateScriptWithSpec(jobSpec, outputDir)
	if err != nil {
		t.Fatalf("CreateScriptWithSpec failed: %v", err)
	}

	content, err := os.ReadFile(scriptPath)
	if err != nil {
		t.Fatalf("failed to read script: %v", err)
	}

	script := string(content)
	if !strings.Contains(script, "#SBATCH --nodes=4") {
		t.Error("script should contain --nodes=4")
	}
	if !strings.Contains(script, "#SBATCH --cpus-per-task=4") {
		t.Error("script should contain --cpus-per-task=4")
	}
	// Should not contain original raw flags duplicated
	count := strings.Count(script, "--nodes=")
	if count != 1 {
		t.Errorf("--nodes= appears %d times; want 1", count)
	}
}

func TestSlurmNodeTaskRoundTrip(t *testing.T) {
	sched := newTestSlurmScheduler()

	// Create a script with node/task directives
	// nodes=4, ntasks=16 => TasksPerNode = 16/4 = 4
	lines := []string{
		"#!/bin/bash",
		"#SBATCH --nodes=4",
		"#SBATCH --ntasks=16",
		"#SBATCH --cpus-per-task=4",
		"#SBATCH --mem=8G",
		"#SBATCH --time=2:00:00",
		"echo hello",
	}
	tmpFile := filepath.Join(t.TempDir(), "input.sh")
	os.WriteFile(tmpFile, []byte(strings.Join(lines, "\n")), 0644)

	// Parse
	specs, err := sched.ReadScriptSpecs(tmpFile)
	if err != nil {
		t.Fatalf("parse failed: %v", err)
	}
	if specs.Spec == nil {
		t.Fatal("Spec is nil; want non-nil")
	}
	if specs.Spec.Nodes != 4 || specs.Spec.TasksPerNode != 4 || specs.Spec.CpusPerTask != 4 {
		t.Fatalf("parse mismatch: Nodes=%d TasksPerNode=%d CpusPerTask=%d",
			specs.Spec.Nodes, specs.Spec.TasksPerNode, specs.Spec.CpusPerTask)
	}

	// Generate
	outputDir := t.TempDir()
	jobSpec := &JobSpec{
		Name:     "roundtrip",
		Command:  "echo hello",
		Specs:    specs,
		Metadata: map[string]string{},
	}
	scriptPath, err := sched.CreateScriptWithSpec(jobSpec, outputDir)
	if err != nil {
		t.Fatalf("generate failed: %v", err)
	}

	// Re-parse generated script
	specs2, err := sched.ReadScriptSpecs(scriptPath)
	if err != nil {
		t.Fatalf("re-parse failed: %v", err)
	}
	if specs2.Spec == nil {
		t.Fatal("re-parsed Spec is nil; want non-nil")
	}
	if specs2.Spec.Nodes != 4 {
		t.Errorf("round-trip Nodes = %d; want 4", specs2.Spec.Nodes)
	}
	if specs2.Spec.TasksPerNode != 4 {
		t.Errorf("round-trip TasksPerNode = %d; want 4", specs2.Spec.TasksPerNode)
	}
	if specs2.Spec.CpusPerTask != 4 {
		t.Errorf("round-trip CpusPerTask = %d; want 4", specs2.Spec.CpusPerTask)
	}

	// Check no duplicates
	content, _ := os.ReadFile(scriptPath)
	script := string(content)
	if strings.Count(script, "--nodes=") != 1 {
		t.Error("--nodes= should appear exactly once")
	}
	if strings.Count(script, "--cpus-per-task=") != 1 {
		t.Error("--cpus-per-task= should appear exactly once")
	}
}

func TestSlurmDirectivesWithComments(t *testing.T) {
	sched := newTestSlurmScheduler()
	tmpDir := t.TempDir()
	tmpFile := filepath.Join(tmpDir, "test.sh")

	script := `#!/bin/bash
#SBATCH --job-name=test_job  # This is a comment
#SBATCH --cpus-per-task=8 # Number of CPUs
#SBATCH --mem=16G#No space before comment
#SBATCH --time=02:00:00   #   Time limit with spaces
#SBATCH --output=out.log # Output file

echo "hello"
`
	if err := os.WriteFile(tmpFile, []byte(script), 0644); err != nil {
		t.Fatalf("failed to create test script: %v", err)
	}

	specs, err := sched.ReadScriptSpecs(tmpFile)
	if err != nil {
		t.Fatalf("parse failed: %v", err)
	}

	// Verify that comments were stripped and values parsed correctly
	if specs.Control.JobName != "test_job" {
		t.Errorf("JobName = %q; want %q", specs.Control.JobName, "test_job")
	}
	if specs.Spec == nil {
		t.Fatal("Spec is nil; want non-nil")
	}
	if specs.Spec.CpusPerTask != 8 {
		t.Errorf("CpusPerTask = %d; want 8", specs.Spec.CpusPerTask)
	}
	if specs.Spec.MemPerNodeMB != 16*1024 {
		t.Errorf("MemPerNodeMB = %d; want %d", specs.Spec.MemPerNodeMB, 16*1024)
	}
	if specs.Spec.Time != 2*time.Hour {
		t.Errorf("Time = %v; want %v", specs.Spec.Time, 2*time.Hour)
	}
	if specs.Control.Stdout != "out.log" {
		t.Errorf("Stdout = %q; want %q", specs.Control.Stdout, "out.log")
	}
}
