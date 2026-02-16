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
			RawFlags: []string{},
			Ncpus:    1,
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
					JobName:      "test_job",
					Ncpus:        4,
					MemMB:        8000,
					Time:         time.Hour,
					EmailOnBegin: tt.emailOnBegin,
					EmailOnEnd:   tt.emailOnEnd,
					EmailOnFail:  tt.emailOnFail,
					MailUser:     tt.mailUser,
					RawFlags:     []string{},
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
	if specs.Ncpus != 8 {
		t.Errorf("Ncpus = %d; want 8", specs.Ncpus)
	}
	if specs.MemMB != 16*1024 {
		t.Errorf("MemMB = %d; want %d", specs.MemMB, 16*1024)
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
		name      string
		lines     []string
		wantNcpus int
		wantMemMB int64
		wantTime  time.Duration
		wantGpu   *GpuSpec
	}{
		{
			name: "select format with ncpus, mem",
			lines: []string{
				"#!/bin/bash",
				"#PBS -l select=1:ncpus=4:mem=8gb",
			},
			wantNcpus: 4,
			wantMemMB: 8 * 1024,
		},
		{
			name: "walltime only",
			lines: []string{
				"#!/bin/bash",
				"#PBS -l walltime=02:30:00",
			},
			wantNcpus: 4, // default
			wantTime:  2*time.Hour + 30*time.Minute,
		},
		{
			name: "nodes/ppn format",
			lines: []string{
				"#!/bin/bash",
				"#PBS -l nodes=1:ppn=8",
			},
			wantNcpus: 8,
		},
		{
			name: "ngpus resource",
			lines: []string{
				"#!/bin/bash",
				"#PBS -l select=1:ncpus=4:mem=16gb:ngpus=2",
			},
			wantNcpus: 4,
			wantMemMB: 16 * 1024,
			wantGpu:   &GpuSpec{Type: "gpu", Count: 2, Raw: "ngpus=2"},
		},
		{
			name: "gpus with type on separate line",
			lines: []string{
				"#!/bin/bash",
				"#PBS -l select=1:ncpus=4",
				"#PBS -l gpus=a100:2",
			},
			wantNcpus: 4,
			wantGpu:   &GpuSpec{Type: "a100", Count: 2, Raw: "a100:2"},
		},
		{
			name: "multiple resource lines",
			lines: []string{
				"#!/bin/bash",
				"#PBS -l select=1:ncpus=16:mem=32gb",
				"#PBS -l walltime=04:00:00",
			},
			wantNcpus: 16,
			wantMemMB: 32 * 1024,
			wantTime:  4 * time.Hour,
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
			wantNcpus: 8,
			wantMemMB: 64 * 1024,
			wantTime:  12 * time.Hour,
			wantGpu:   &GpuSpec{Type: "gpu", Count: 4, Raw: "ngpus=4"},
		},
		{
			name: "memory in mb",
			lines: []string{
				"#!/bin/bash",
				"#PBS -l select=1:ncpus=2:mem=4096mb",
			},
			wantNcpus: 2,
			wantMemMB: 4096,
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

			if specs.Ncpus != tt.wantNcpus {
				t.Errorf("Ncpus = %d; want %d", specs.Ncpus, tt.wantNcpus)
			}
			if tt.wantMemMB > 0 && specs.MemMB != tt.wantMemMB {
				t.Errorf("MemMB = %d; want %d", specs.MemMB, tt.wantMemMB)
			}
			if tt.wantTime > 0 && specs.Time != tt.wantTime {
				t.Errorf("Time = %v; want %v", specs.Time, tt.wantTime)
			}
			if tt.wantGpu != nil {
				if specs.Gpu == nil {
					t.Fatal("Gpu is nil; want non-nil")
				}
				if specs.Gpu.Type != tt.wantGpu.Type {
					t.Errorf("Gpu.Type = %q; want %q", specs.Gpu.Type, tt.wantGpu.Type)
				}
				if specs.Gpu.Count != tt.wantGpu.Count {
					t.Errorf("Gpu.Count = %d; want %d", specs.Gpu.Count, tt.wantGpu.Count)
				}
			} else if specs.Gpu != nil {
				t.Errorf("Gpu = %+v; want nil", specs.Gpu)
			}
		})
	}
}

func TestPbsTimeParsing(t *testing.T) {
	tests := []struct {
		name    string
		input   string
		wantDur time.Duration
		wantErr bool
	}{
		{
			name:    "HH:MM:SS",
			input:   "02:30:00",
			wantDur: 2*time.Hour + 30*time.Minute,
		},
		{
			name:    "HH:MM",
			input:   "10:30",
			wantDur: 10*time.Hour + 30*time.Minute,
		},
		{
			name:    "minutes only",
			input:   "90",
			wantDur: 90 * time.Minute,
		},
		{
			name:    "with seconds",
			input:   "01:00:30",
			wantDur: time.Hour + 30*time.Second,
		},
		{
			name:    "empty string",
			input:   "",
			wantDur: 0,
		},
		{
			name:    "large walltime",
			input:   "168:00:00",
			wantDur: 168 * time.Hour,
		},
	}

	for _, tt := range tests {
		t.Run(tt.name, func(t *testing.T) {
			dur, err := parsePbsTime(tt.input)
			if tt.wantErr && err == nil {
				t.Error("Expected error, got nil")
			}
			if !tt.wantErr && err != nil {
				t.Errorf("Unexpected error: %v", err)
			}
			if dur != tt.wantDur {
				t.Errorf("parsePbsTime(%q) = %v; want %v", tt.input, dur, tt.wantDur)
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

	t.Run("parseMemoryString (returns MB)", func(t *testing.T) {
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
			mb, err := parseMemoryString(tt.input)
			if err != nil {
				t.Errorf("parseMemoryString(%q) error: %v", tt.input, err)
				continue
			}
			if mb != tt.wantMB {
				t.Errorf("parseMemoryString(%q) = %d MB; want %d MB", tt.input, mb, tt.wantMB)
			}
		}
	})
}

func TestPbsGpuParsing(t *testing.T) {
	tests := []struct {
		name      string
		input     string
		wantType  string
		wantCount int
	}{
		{
			name:      "count only",
			input:     "2",
			wantType:  "gpu",
			wantCount: 2,
		},
		{
			name:      "type only",
			input:     "a100",
			wantType:  "a100",
			wantCount: 1,
		},
		{
			name:      "type:count",
			input:     "a100:2",
			wantType:  "a100",
			wantCount: 2,
		},
		{
			name:      "gpu:type:count",
			input:     "gpu:a100:4",
			wantType:  "a100",
			wantCount: 4,
		},
		{
			name:      "MIG profile without count",
			input:     "nvidia_h100_80gb_hbm3_1g.10gb",
			wantType:  "nvidia_h100_80gb_hbm3_1g.10gb",
			wantCount: 1,
		},
		{
			name:      "MIG profile with count",
			input:     "nvidia_h100_80gb_hbm3_1g.10gb:2",
			wantType:  "nvidia_h100_80gb_hbm3_1g.10gb",
			wantCount: 2,
		},
		{
			name:      "v100 single",
			input:     "v100:1",
			wantType:  "v100",
			wantCount: 1,
		},
	}

	for _, tt := range tests {
		t.Run(tt.name, func(t *testing.T) {
			spec, err := parseGpuString(tt.input)
			if err != nil {
				t.Fatalf("parseGpuString(%q) error: %v", tt.input, err)
			}
			if spec == nil {
				t.Fatal("parseGpuString returned nil")
			}
			if spec.Type != tt.wantType {
				t.Errorf("Type = %q; want %q", spec.Type, tt.wantType)
			}
			if spec.Count != tt.wantCount {
				t.Errorf("Count = %d; want %d", spec.Count, tt.wantCount)
			}
		})
	}

	// Test empty string returns nil
	t.Run("empty string", func(t *testing.T) {
		spec, err := parseGpuString("")
		if err != nil {
			t.Fatalf("Unexpected error: %v", err)
		}
		if spec != nil {
			t.Errorf("Expected nil for empty string, got %+v", spec)
		}
	})
}

func TestPbsIsAvailable(t *testing.T) {
	t.Run("not in job", func(t *testing.T) {
		clearJobEnvVars(t)
		pbs := newTestPbsScheduler()
		if !pbs.IsAvailable() {
			t.Error("Expected IsAvailable to return true when not in a job")
		}
	})

	t.Run("inside job", func(t *testing.T) {
		clearJobEnvVars(t)
		t.Setenv("PBS_JOBID", "67890.pbs-server")
		pbs := newTestPbsScheduler()
		if pbs.IsAvailable() {
			t.Error("Expected IsAvailable to return false when inside a job")
		}
	})

	t.Run("no binary", func(t *testing.T) {
		clearJobEnvVars(t)
		pbs := &PbsScheduler{}
		if pbs.IsAvailable() {
			t.Error("Expected IsAvailable to return false when no binary is set")
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

	if specs.Ncpus != 4 {
		t.Errorf("Ncpus = %d; want 4", specs.Ncpus)
	}
	if specs.MemMB != 8*1024 {
		t.Errorf("MemMB = %d; want %d", specs.MemMB, 8*1024)
	}
	if specs.Gpu == nil || specs.Gpu.Count != 1 {
		t.Errorf("Gpu = %+v; want count 1", specs.Gpu)
	}
	if specs.Time != 2*time.Hour {
		t.Errorf("Time = %v; want 2h", specs.Time)
	}
	if !specs.EmailOnEnd {
		t.Error("EmailOnEnd should be true")
	}
	if specs.MailUser != "user@example.com" {
		t.Errorf("MailUser = %q; want %q", specs.MailUser, "user@example.com")
	}
	if len(specs.RawFlags) != 5 {
		t.Errorf("RawFlags count = %d; want 5", len(specs.RawFlags))
	}
}

func TestPbsNodeTaskParsing(t *testing.T) {
	tests := []struct {
		name       string
		lines      []string
		wantNodes  int
		wantNtasks int
		wantNcpus  int
	}{
		{
			name: "select with mpiprocs",
			lines: []string{
				"#!/bin/bash",
				"#PBS -l select=4:ncpus=8:mpiprocs=2:mem=8gb",
			},
			wantNodes:  4,
			wantNtasks: 8, // 4 nodes * 2 mpiprocs
			wantNcpus:  8,
		},
		{
			name: "select without mpiprocs",
			lines: []string{
				"#!/bin/bash",
				"#PBS -l select=2:ncpus=4",
			},
			wantNodes:  2,
			wantNtasks: 1, // default, no mpiprocs
			wantNcpus:  4,
		},
		{
			name: "nodes format",
			lines: []string{
				"#!/bin/bash",
				"#PBS -l nodes=3:ppn=8",
			},
			wantNodes:  3,
			wantNtasks: 1, // default, nodes= doesn't set ntasks
			wantNcpus:  8,
		},
		{
			name: "select=1 single node",
			lines: []string{
				"#!/bin/bash",
				"#PBS -l select=1:ncpus=4:mem=16gb",
			},
			wantNodes:  1,
			wantNtasks: 1,
			wantNcpus:  4,
		},
		{
			name: "defaults no resource directives",
			lines: []string{
				"#!/bin/bash",
				"#PBS -N testjob",
			},
			wantNodes:  1,
			wantNtasks: 1,
			wantNcpus:  4, // default
		},
		{
			name: "select with mpiprocs and ngpus",
			lines: []string{
				"#!/bin/bash",
				"#PBS -l select=2:ncpus=16:mpiprocs=4:mem=32gb:ngpus=1",
			},
			wantNodes:  2,
			wantNtasks: 8, // 2 nodes * 4 mpiprocs
			wantNcpus:  16,
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

			if specs.Nodes != tt.wantNodes {
				t.Errorf("Nodes = %d; want %d", specs.Nodes, tt.wantNodes)
			}
			if specs.Ntasks != tt.wantNtasks {
				t.Errorf("Ntasks = %d; want %d", specs.Ntasks, tt.wantNtasks)
			}
			if specs.Ncpus != tt.wantNcpus {
				t.Errorf("Ncpus = %d; want %d", specs.Ncpus, tt.wantNcpus)
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
		t.Setenv("PBS_NCPUS", "8")
		t.Setenv("CUDA_VISIBLE_DEVICES", "0,1,2")

		res := sched.GetJobResources()
		if res == nil {
			t.Fatal("expected non-nil")
		}
		if res.Ncpus == nil || *res.Ncpus != 8 {
			t.Errorf("Ncpus = %v; want 8", res.Ncpus)
		}
		if res.Ngpus == nil || *res.Ngpus != 3 {
			t.Errorf("Ngpus = %v; want 3", res.Ngpus)
		}
	})

	t.Run("nodes and tasks", func(t *testing.T) {
		clearJobEnvVars(t)
		t.Setenv("PBS_JOBID", "67890.pbs-server")
		t.Setenv("PBS_NCPUS", "8")
		t.Setenv("PBS_NUM_NODES", "4")
		t.Setenv("PBS_NP", "16")

		res := sched.GetJobResources()
		if res == nil {
			t.Fatal("expected non-nil")
		}
		if res.Nodes == nil || *res.Nodes != 4 {
			t.Errorf("Nodes = %v; want 4", res.Nodes)
		}
		if res.Ntasks == nil || *res.Ntasks != 16 {
			t.Errorf("Ntasks = %v; want 16", res.Ntasks)
		}
	})

	t.Run("PBS_TASKNUM fallback", func(t *testing.T) {
		clearJobEnvVars(t)
		t.Setenv("PBS_JOBID", "67890")
		t.Setenv("PBS_TASKNUM", "8")

		res := sched.GetJobResources()
		if res == nil {
			t.Fatal("expected non-nil")
		}
		if res.Ntasks == nil || *res.Ntasks != 8 {
			t.Errorf("Ntasks = %v; want 8", res.Ntasks)
		}
	})

	t.Run("NCPUS fallback", func(t *testing.T) {
		clearJobEnvVars(t)
		t.Setenv("PBS_JOBID", "67890")
		t.Setenv("NCPUS", "12")

		res := sched.GetJobResources()
		if res == nil {
			t.Fatal("expected non-nil")
		}
		if res.Ncpus == nil || *res.Ncpus != 12 {
			t.Errorf("Ncpus = %v; want 12", res.Ncpus)
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
		if res.MemMB == nil || *res.MemMB != 8192 {
			t.Errorf("MemMB = %v; want 8192", res.MemMB)
		}
	})

	t.Run("partial data", func(t *testing.T) {
		clearJobEnvVars(t)
		t.Setenv("PBS_JOBID", "67890")
		t.Setenv("PBS_NCPUS", "4")

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
		t.Setenv("PBS_JOBID", "67890")
		t.Setenv("PBS_NCPUS", "not-a-number")
		t.Setenv("PBS_VMEM", "-100")

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
