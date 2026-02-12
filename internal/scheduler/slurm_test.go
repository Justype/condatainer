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
			if specs.EmailOnBegin != tt.wantBegin {
				t.Errorf("EmailOnBegin = %v; want %v", specs.EmailOnBegin, tt.wantBegin)
			}
			if specs.EmailOnEnd != tt.wantEnd {
				t.Errorf("EmailOnEnd = %v; want %v", specs.EmailOnEnd, tt.wantEnd)
			}
			if specs.EmailOnFail != tt.wantFail {
				t.Errorf("EmailOnFail = %v; want %v", specs.EmailOnFail, tt.wantFail)
			}
			if specs.MailUser != tt.wantMailUser {
				t.Errorf("MailUser = %q; want %q", specs.MailUser, tt.wantMailUser)
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
					JobName:      "test_job",
					Ncpus:        4,
					MemMB:        8000,
					Time:         time.Hour,
					EmailOnBegin: tt.emailOnBegin,
					EmailOnEnd:   tt.emailOnEnd,
					EmailOnFail:  tt.emailOnFail,
					MailUser:     tt.mailUser,
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
	if !specs.EmailOnBegin || !specs.EmailOnEnd || !specs.EmailOnFail {
		t.Errorf("Parsing failed: EmailOnBegin=%v, EmailOnEnd=%v, EmailOnFail=%v",
			specs.EmailOnBegin, specs.EmailOnEnd, specs.EmailOnFail)
	}
	if specs.MailUser != "roundtrip@example.com" {
		t.Errorf("MailUser = %q; want %q", specs.MailUser, "roundtrip@example.com")
	}
	if specs.Time != 2*time.Hour {
		t.Errorf("Time = %v; want 2h", specs.Time)
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
			RawFlags: []string{},
			Ncpus:    1,
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

// clearJobEnvVars ensures no scheduler job env vars are set for test isolation.
// Shared across all scheduler test files (same package).
func clearJobEnvVars(t *testing.T) {
	t.Helper()
	for _, key := range []string{
		"SLURM_JOB_ID", "SLURM_CPUS_PER_TASK", "SLURM_MEM_PER_NODE",
		"PBS_JOBID", "PBS_NCPUS", "NCPUS", "PBS_VMEM",
		"LSB_JOBID", "LSB_DJOB_NUMPROC", "LSB_MAX_NUM_PROCESSORS", "LSB_MAX_MEM_RUSAGE",
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
		t.Setenv("SLURM_MEM_PER_NODE", "8192")
		t.Setenv("CUDA_VISIBLE_DEVICES", "0,1")

		res := sched.GetJobResources()
		if res == nil {
			t.Fatal("expected non-nil")
		}
		if res.Ncpus == nil || *res.Ncpus != 16 {
			t.Errorf("Ncpus = %v; want 16", res.Ncpus)
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
