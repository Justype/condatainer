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

// newTestLsfScheduler creates an LSF scheduler instance for testing
// without requiring bsub to be installed
func newTestLsfScheduler() *LsfScheduler {
	return &LsfScheduler{
		bsubBin:     "/usr/bin/bsub",
		bjobsBin:    "/usr/bin/bjobs",
		bhostsBin:   "/usr/bin/bhosts",
		directiveRe: regexp.MustCompile(`^\s*#BSUB\s+(.+)$`),
		jobIDRe:     regexp.MustCompile(`Job <(\d+)> is submitted`),
	}
}

// Test that scheduler uses outputDir for log files
func TestLsfCreateScriptUsesOutputDir(t *testing.T) {
	tmpOutputDir := t.TempDir()

	lsf := newTestLsfScheduler()

	jobSpec := &JobSpec{
		Name:    "lsftest/job",
		Command: "echo 'hello'",
		Specs: &ScriptSpecs{
			RawFlags: []string{},
			Ncpus:    1,
		},
	}

	scriptPath, err := lsf.CreateScriptWithSpec(jobSpec, tmpOutputDir)
	if err != nil {
		t.Fatalf("CreateScriptWithSpec failed: %v", err)
	}

	// Ensure the generated script uses outputDir for logs
	logPath := filepath.Join(tmpOutputDir, "lsftest--job.log")
	content, err := os.ReadFile(scriptPath)
	if err != nil {
		t.Fatalf("Failed to read generated script: %v", err)
	}
	if !strings.Contains(string(content), fmt.Sprintf("-o %s", logPath)) {
		t.Errorf("Generated script does not contain expected output path %s\nScript:\n%s", logPath, string(content))
	}
	// Verify single-node enforcement via span
	if !strings.Contains(string(content), `span[hosts=1]`) {
		t.Errorf("Generated script does not contain span[hosts=1]\nScript:\n%s", string(content))
	}
}

func TestLsfEmailParsing(t *testing.T) {
	tests := []struct {
		name         string
		scriptLines  []string
		wantBegin    bool
		wantEnd      bool
		wantFail     bool
		wantMailUser string
	}{
		{
			name: "begin and end notifications",
			scriptLines: []string{
				"#!/bin/bash",
				"#BSUB -B",
				"#BSUB -N",
				"#BSUB -u test@example.com",
			},
			wantBegin:    true,
			wantEnd:      true,
			wantFail:     false,
			wantMailUser: "test@example.com",
		},
		{
			name: "begin only",
			scriptLines: []string{
				"#!/bin/bash",
				"#BSUB -B",
			},
			wantBegin: true,
			wantEnd:   false,
			wantFail:  false,
		},
		{
			name: "end only",
			scriptLines: []string{
				"#!/bin/bash",
				"#BSUB -N",
			},
			wantBegin: false,
			wantEnd:   true,
			wantFail:  false,
		},
		{
			name: "mail user without notifications",
			scriptLines: []string{
				"#!/bin/bash",
				"#BSUB -u admin@company.com",
			},
			wantBegin:    false,
			wantEnd:      false,
			wantFail:     false,
			wantMailUser: "admin@company.com",
		},
		{
			name: "no email directives",
			scriptLines: []string{
				"#!/bin/bash",
				"#BSUB -J testjob",
			},
			wantBegin: false,
			wantEnd:   false,
			wantFail:  false,
		},
		{
			name: "begin with mail user",
			scriptLines: []string{
				"#!/bin/bash",
				"#BSUB -B",
				"#BSUB -u user@domain.org",
			},
			wantBegin:    true,
			wantEnd:      false,
			wantFail:     false,
			wantMailUser: "user@domain.org",
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

			lsf := newTestLsfScheduler()
			specs, err := lsf.ReadScriptSpecs(scriptPath)
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

func TestLsfEmailScriptGeneration(t *testing.T) {
	tests := []struct {
		name         string
		emailOnBegin bool
		emailOnEnd   bool
		mailUser     string
		wantB        bool   // expect #BSUB -B
		wantN        bool   // expect #BSUB -N
		wantMailUser string // expected -u value (empty = not present)
	}{
		{
			name:         "begin and end with user",
			emailOnBegin: true,
			emailOnEnd:   true,
			mailUser:     "test@example.com",
			wantB:        true,
			wantN:        true,
			wantMailUser: "test@example.com",
		},
		{
			name:         "begin only",
			emailOnBegin: true,
			emailOnEnd:   false,
			wantB:        true,
			wantN:        false,
			wantMailUser: "",
		},
		{
			name:         "end only with user",
			emailOnBegin: false,
			emailOnEnd:   true,
			mailUser:     "user@domain.org",
			wantB:        false,
			wantN:        true,
			wantMailUser: "user@domain.org",
		},
		{
			name:         "no email notifications",
			emailOnBegin: false,
			emailOnEnd:   false,
			wantB:        false,
			wantN:        false,
			wantMailUser: "",
		},
	}

	for _, tt := range tests {
		t.Run(tt.name, func(t *testing.T) {
			tmpDir := t.TempDir()
			lsf := newTestLsfScheduler()

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
					MailUser:     tt.mailUser,
					RawFlags:     []string{},
				},
			}

			scriptPath, err := lsf.CreateScriptWithSpec(jobSpec, tmpDir)
			if err != nil {
				t.Fatalf("Failed to create script: %v", err)
			}

			content, err := os.ReadFile(scriptPath)
			if err != nil {
				t.Fatalf("Failed to read generated script: %v", err)
			}
			scriptContent := string(content)

			// Check -B directive
			if tt.wantB {
				if !strings.Contains(scriptContent, "#BSUB -B") {
					t.Errorf("Script missing #BSUB -B\nScript content:\n%s", scriptContent)
				}
			} else {
				if strings.Contains(scriptContent, "#BSUB -B") {
					t.Errorf("Script should not contain #BSUB -B\nScript content:\n%s", scriptContent)
				}
			}

			// Check -N directive
			if tt.wantN {
				if !strings.Contains(scriptContent, "#BSUB -N") {
					t.Errorf("Script missing #BSUB -N\nScript content:\n%s", scriptContent)
				}
			} else {
				if strings.Contains(scriptContent, "#BSUB -N\n") {
					t.Errorf("Script should not contain #BSUB -N\nScript content:\n%s", scriptContent)
				}
			}

			// Check -u directive
			if tt.wantMailUser != "" {
				expectedLine := "#BSUB -u " + tt.wantMailUser
				if !strings.Contains(scriptContent, expectedLine) {
					t.Errorf("Script missing expected line: %q\nScript content:\n%s", expectedLine, scriptContent)
				}
			} else {
				if strings.Contains(scriptContent, "#BSUB -u ") {
					t.Errorf("Script should not contain -u directive\nScript content:\n%s", scriptContent)
				}
			}
		})
	}
}

func TestLsfEmailRoundTrip(t *testing.T) {
	tmpDir := t.TempDir()

	originalScript := `#!/bin/bash
#BSUB -J email_test
#BSUB -n 8
#BSUB -M 16384
#BSUB -W 02:00
#BSUB -B
#BSUB -N
#BSUB -u roundtrip@example.com

echo "Running job"
`
	scriptPath := filepath.Join(tmpDir, "original.sh")
	if err := os.WriteFile(scriptPath, []byte(originalScript), 0644); err != nil {
		t.Fatalf("Failed to create test script: %v", err)
	}

	lsf := newTestLsfScheduler()

	// Parse the original script
	specs, err := lsf.ReadScriptSpecs(scriptPath)
	if err != nil {
		t.Fatalf("Failed to parse script: %v", err)
	}

	// Verify parsing
	if !specs.EmailOnBegin || !specs.EmailOnEnd {
		t.Errorf("Parsing failed: EmailOnBegin=%v, EmailOnEnd=%v",
			specs.EmailOnBegin, specs.EmailOnEnd)
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

	// Generate a new script from the parsed specs
	jobSpec := &JobSpec{
		Name:    "email_test",
		Command: "echo 'Running job'",
		Specs:   specs,
	}

	newScriptPath, err := lsf.CreateScriptWithSpec(jobSpec, tmpDir)
	if err != nil {
		t.Fatalf("Failed to create script: %v", err)
	}

	content, err := os.ReadFile(newScriptPath)
	if err != nil {
		t.Fatalf("Failed to read generated script: %v", err)
	}
	scriptContent := string(content)

	// Verify the generated script contains the email directives
	if !strings.Contains(scriptContent, "#BSUB -B") {
		t.Errorf("Generated script missing -B directive:\n%s", scriptContent)
	}
	if !strings.Contains(scriptContent, "#BSUB -N") {
		t.Errorf("Generated script missing -N directive:\n%s", scriptContent)
	}
	if !strings.Contains(scriptContent, "#BSUB -u roundtrip@example.com") {
		t.Errorf("Generated script missing -u directive:\n%s", scriptContent)
	}

	// Verify NO duplicates
	bCount := strings.Count(scriptContent, "#BSUB -B")
	nCount := strings.Count(scriptContent, "#BSUB -N\n")
	uCount := strings.Count(scriptContent, "#BSUB -u ")
	if bCount != 1 {
		t.Errorf("Expected exactly 1 -B directive, found %d:\n%s", bCount, scriptContent)
	}
	if nCount != 1 {
		t.Errorf("Expected exactly 1 -N directive, found %d:\n%s", nCount, scriptContent)
	}
	if uCount != 1 {
		t.Errorf("Expected exactly 1 -u directive, found %d:\n%s", uCount, scriptContent)
	}
}

func TestLsfResourceParsing(t *testing.T) {
	tests := []struct {
		name      string
		lines     []string
		wantNcpus int
		wantMemMB int64
		wantTime  time.Duration
		wantGpu   *GpuSpec
	}{
		{
			name: "basic resources",
			lines: []string{
				"#!/bin/bash",
				"#BSUB -J testjob",
				"#BSUB -n 4",
				"#BSUB -M 8192",
				"#BSUB -W 01:30",
			},
			wantNcpus: 4,
			wantMemMB: 8, // 8192 KB = 8 MB
			wantTime:  time.Hour + 30*time.Minute,
		},
		{
			name: "ncpus only",
			lines: []string{
				"#!/bin/bash",
				"#BSUB -n 16",
			},
			wantNcpus: 16,
		},
		{
			name: "walltime only",
			lines: []string{
				"#!/bin/bash",
				"#BSUB -W 04:00",
			},
			wantNcpus: 2, // default
			wantTime:  4 * time.Hour,
		},
		{
			name: "gpu directive",
			lines: []string{
				"#!/bin/bash",
				"#BSUB -n 8",
				`#BSUB -gpu "num=2:type=a100"`,
			},
			wantNcpus: 8,
			wantGpu:   &GpuSpec{Type: "a100", Count: 2},
		},
		{
			name: "gpu num only",
			lines: []string{
				"#!/bin/bash",
				`#BSUB -gpu "num=4"`,
			},
			wantNcpus: 2, // default
			wantGpu:   &GpuSpec{Type: "gpu", Count: 4},
		},
		{
			name: "rusage mem",
			lines: []string{
				"#!/bin/bash",
				"#BSUB -n 4",
				`#BSUB -R "rusage[mem=4096]"`,
			},
			wantNcpus: 4,
			wantMemMB: 4, // 4096 KB = 4 MB
		},
		{
			name: "rusage ngpus_physical",
			lines: []string{
				"#!/bin/bash",
				"#BSUB -n 8",
				`#BSUB -R "rusage[ngpus_physical=2]"`,
			},
			wantNcpus: 8,
			wantGpu:   &GpuSpec{Type: "gpu", Count: 2},
		},
		{
			name: "full job script",
			lines: []string{
				"#!/bin/bash",
				"#BSUB -J myjob",
				"#BSUB -n 32",
				"#BSUB -M 65536",
				"#BSUB -W 12:00",
				"#BSUB -B",
				"#BSUB -N",
				"#BSUB -u user@example.com",
				`#BSUB -gpu "num=4:type=v100"`,
			},
			wantNcpus: 32,
			wantMemMB: 64, // 65536 KB = 64 MB
			wantTime:  12 * time.Hour,
			wantGpu:   &GpuSpec{Type: "v100", Count: 4},
		},
		{
			name: "memory with MB suffix",
			lines: []string{
				"#!/bin/bash",
				"#BSUB -M 4096MB",
			},
			wantNcpus: 2, // default
			wantMemMB: 4096,
		},
		{
			name: "memory with GB suffix",
			lines: []string{
				"#!/bin/bash",
				"#BSUB -M 8GB",
			},
			wantNcpus: 2, // default
			wantMemMB: 8 * 1024,
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

			lsf := newTestLsfScheduler()
			specs, err := lsf.ReadScriptSpecs(scriptPath)
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

func TestLsfTimeParsing(t *testing.T) {
	tests := []struct {
		name    string
		input   string
		wantDur time.Duration
		wantErr bool
	}{
		{
			name:    "HH:MM",
			input:   "02:30",
			wantDur: 2*time.Hour + 30*time.Minute,
		},
		{
			name:    "minutes only",
			input:   "90",
			wantDur: 90 * time.Minute,
		},
		{
			name:    "HH:MM:SS",
			input:   "01:30:45",
			wantDur: time.Hour + 30*time.Minute + 45*time.Second,
		},
		{
			name:    "empty string",
			input:   "",
			wantDur: 0,
		},
		{
			name:    "large walltime",
			input:   "168:00",
			wantDur: 168 * time.Hour,
		},
		{
			name:    "zero",
			input:   "00:00",
			wantDur: 0,
		},
	}

	for _, tt := range tests {
		t.Run(tt.name, func(t *testing.T) {
			dur, err := parseLsfTime(tt.input)
			if tt.wantErr && err == nil {
				t.Error("Expected error, got nil")
			}
			if !tt.wantErr && err != nil {
				t.Errorf("Unexpected error: %v", err)
			}
			if dur != tt.wantDur {
				t.Errorf("parseLsfTime(%q) = %v; want %v", tt.input, dur, tt.wantDur)
			}
		})
	}
}

func TestLsfMemoryParsing(t *testing.T) {
	tests := []struct {
		input  string
		wantMB int64
	}{
		{"8192", 8},          // 8192 KB = 8 MB
		{"1048576", 1024},    // 1048576 KB = 1024 MB
		{"4096MB", 4096},     // 4096 MB
		{"8GB", 8 * 1024},    // 8 GB
		{"1TB", 1024 * 1024}, // 1 TB
		{"512", 0},           // 512 KB = 0 MB (integer division)
		{"1024", 1},          // 1024 KB = 1 MB
	}

	for _, tt := range tests {
		t.Run(tt.input, func(t *testing.T) {
			mb, err := parseLsfMemory(tt.input)
			if err != nil {
				t.Errorf("parseLsfMemory(%q) error: %v", tt.input, err)
				return
			}
			if mb != tt.wantMB {
				t.Errorf("parseLsfMemory(%q) = %d MB; want %d MB", tt.input, mb, tt.wantMB)
			}
		})
	}
}

func TestLsfGpuParsing(t *testing.T) {
	tests := []struct {
		name      string
		input     string
		wantType  string
		wantCount int
	}{
		{
			name:      "num only",
			input:     "num=2",
			wantType:  "gpu",
			wantCount: 2,
		},
		{
			name:      "num and type",
			input:     "num=4:type=a100",
			wantType:  "a100",
			wantCount: 4,
		},
		{
			name:      "num and type reversed order",
			input:     "type=v100:num=2",
			wantType:  "v100",
			wantCount: 2,
		},
		{
			name:      "single gpu",
			input:     "num=1:type=h100",
			wantType:  "h100",
			wantCount: 1,
		},
		{
			name:      "num with extra options",
			input:     "num=2:type=a100:mode=exclusive_process",
			wantType:  "a100",
			wantCount: 2,
		},
	}

	for _, tt := range tests {
		t.Run(tt.name, func(t *testing.T) {
			spec := parseLsfGpuDirective(tt.input)
			if spec == nil {
				t.Fatal("parseLsfGpuDirective returned nil")
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
		spec := parseLsfGpuDirective("")
		if spec != nil {
			t.Errorf("Expected nil for empty string, got %+v", spec)
		}
	})
}

func TestLsfTimeFormatting(t *testing.T) {
	tests := []struct {
		input time.Duration
		want  string
	}{
		{time.Hour, "01:00"},
		{2*time.Hour + 30*time.Minute, "02:30"},
		{30 * time.Minute, "00:30"},
		{24 * time.Hour, "24:00"},
		{168 * time.Hour, "168:00"},
		{0, ""},
	}

	for _, tt := range tests {
		t.Run(tt.want, func(t *testing.T) {
			got := formatLsfTime(tt.input)
			if got != tt.want {
				t.Errorf("formatLsfTime(%v) = %q; want %q", tt.input, got, tt.want)
			}
		})
	}
}

func TestLsfRusageParsing(t *testing.T) {
	tests := []struct {
		name      string
		lines     []string
		wantMemMB int64
		wantGpu   *GpuSpec
	}{
		{
			name: "rusage mem and ngpus",
			lines: []string{
				"#!/bin/bash",
				`#BSUB -R "rusage[mem=8192,ngpus_physical=2]"`,
			},
			wantMemMB: 8,
			wantGpu:   &GpuSpec{Type: "gpu", Count: 2},
		},
		{
			name: "rusage mem colon separated",
			lines: []string{
				"#!/bin/bash",
				`#BSUB -R "rusage[mem=4096:ngpus_physical=1]"`,
			},
			wantMemMB: 4,
			wantGpu:   &GpuSpec{Type: "gpu", Count: 1},
		},
		{
			name: "rusage with span",
			lines: []string{
				"#!/bin/bash",
				`#BSUB -R "span[hosts=1] rusage[mem=16384]"`,
			},
			wantMemMB: 16,
		},
		{
			name: "rusage ngpus shorthand",
			lines: []string{
				"#!/bin/bash",
				`#BSUB -R "rusage[ngpus=4]"`,
			},
			wantGpu: &GpuSpec{Type: "gpu", Count: 4},
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

			lsf := newTestLsfScheduler()
			specs, err := lsf.ReadScriptSpecs(scriptPath)
			if err != nil {
				t.Fatalf("Failed to parse script: %v", err)
			}

			if tt.wantMemMB > 0 && specs.MemMB != tt.wantMemMB {
				t.Errorf("MemMB = %d; want %d", specs.MemMB, tt.wantMemMB)
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

func TestLsfSpanHostsParsing(t *testing.T) {
	tests := []struct {
		name      string
		lines     []string
		wantNodes int
	}{
		{
			name: "span hosts=2 with rusage",
			lines: []string{
				"#!/bin/bash",
				"#BSUB -n 8",
				`#BSUB -R "span[hosts=2] rusage[mem=8192]"`,
			},
			wantNodes: 2,
		},
		{
			name: "span hosts=1",
			lines: []string{
				"#!/bin/bash",
				"#BSUB -n 4",
				`#BSUB -R "span[hosts=1]"`,
			},
			wantNodes: 1,
		},
		{
			name: "no span directive",
			lines: []string{
				"#!/bin/bash",
				"#BSUB -n 16",
				`#BSUB -R "rusage[mem=4096]"`,
			},
			wantNodes: 1, // default
		},
		{
			name: "span hosts=4 on separate line",
			lines: []string{
				"#!/bin/bash",
				"#BSUB -n 32",
				`#BSUB -R "span[hosts=4]"`,
				`#BSUB -R "rusage[mem=16384]"`,
			},
			wantNodes: 4,
		},
		{
			name: "no resource directives at all",
			lines: []string{
				"#!/bin/bash",
				"#BSUB -J testjob",
			},
			wantNodes: 1, // default
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

			lsf := newTestLsfScheduler()
			specs, err := lsf.ReadScriptSpecs(scriptPath)
			if err != nil {
				t.Fatalf("Failed to parse script: %v", err)
			}

			if specs.Nodes != tt.wantNodes {
				t.Errorf("Nodes = %d; want %d", specs.Nodes, tt.wantNodes)
			}
			// Ntasks should always be 1 (default) for LSF
			if specs.Ntasks != 1 {
				t.Errorf("Ntasks = %d; want 1", specs.Ntasks)
			}
		})
	}
}

func TestLsfScriptGenerationNodesCpus(t *testing.T) {
	tests := []struct {
		name         string
		specs        *ScriptSpecs
		wantSpanHost string // expected span[hosts=N]
		wantNcpus    string // expected -n N
	}{
		{
			name: "default single node",
			specs: &ScriptSpecs{
				Ncpus:    4,
				RawFlags: []string{},
			},
			wantSpanHost: `span[hosts=1]`,
			wantNcpus:    "#BSUB -n 4",
		},
		{
			name: "multi-node",
			specs: &ScriptSpecs{
				Ncpus:    8,
				Nodes:    2,
				RawFlags: []string{},
			},
			wantSpanHost: `span[hosts=2]`,
			wantNcpus:    "#BSUB -n 8",
		},
		{
			name: "raw span no longer in RawFlags - only rusage preserved",
			specs: &ScriptSpecs{
				Ncpus:    16,
				Nodes:    3,
				RawFlags: []string{`-R "rusage[mem=4096]"`}, // span parsed out, only rusage remains
			},
			wantSpanHost: `span[hosts=3]`,
			wantNcpus:    "#BSUB -n 16",
		},
		{
			name: "raw -n no longer in RawFlags - comes from Ncpus field",
			specs: &ScriptSpecs{
				Ncpus:    8,
				Nodes:    1,
				RawFlags: []string{}, // -n not in RawFlags, comes from Ncpus
			},
			wantSpanHost: `span[hosts=1]`,
			wantNcpus:    "#BSUB -n 8",
		},
	}

	for _, tt := range tests {
		t.Run(tt.name, func(t *testing.T) {
			tmpDir := t.TempDir()
			lsf := newTestLsfScheduler()

			jobSpec := &JobSpec{
				Name:    "test_job",
				Command: "echo 'test'",
				Specs:   tt.specs,
			}

			scriptPath, err := lsf.CreateScriptWithSpec(jobSpec, tmpDir)
			if err != nil {
				t.Fatalf("Failed to create script: %v", err)
			}

			content, err := os.ReadFile(scriptPath)
			if err != nil {
				t.Fatalf("Failed to read generated script: %v", err)
			}
			scriptContent := string(content)

			// Check span directive
			if !strings.Contains(scriptContent, tt.wantSpanHost) {
				t.Errorf("Script missing expected span: %q\nScript:\n%s", tt.wantSpanHost, scriptContent)
			}

			// Check -n directive
			if !strings.Contains(scriptContent, tt.wantNcpus) {
				t.Errorf("Script missing expected -n: %q\nScript:\n%s", tt.wantNcpus, scriptContent)
			}

			// Verify no duplicate span directives
			spanCount := strings.Count(scriptContent, "span[hosts=")
			if spanCount != 1 {
				t.Errorf("span[hosts=] appears %d times; want 1\nScript:\n%s", spanCount, scriptContent)
			}

			// Verify no duplicate -n directives
			nCount := strings.Count(scriptContent, "#BSUB -n ")
			if nCount != 1 {
				t.Errorf("#BSUB -n appears %d times; want 1\nScript:\n%s", nCount, scriptContent)
			}
		})
	}
}

func TestLsfScriptPreservesRusageWithSpan(t *testing.T) {
	tmpDir := t.TempDir()
	lsf := newTestLsfScheduler()

	jobSpec := &JobSpec{
		Name:    "test_job",
		Command: "echo 'test'",
		Specs: &ScriptSpecs{
			Ncpus:    8,
			Nodes:    2,
			MemMB:    4, // 4096 KB = 4 MB parsed from rusage
			RawFlags: []string{`-R "span[hosts=2] rusage[mem=4096]"`},
		},
	}

	scriptPath, err := lsf.CreateScriptWithSpec(jobSpec, tmpDir)
	if err != nil {
		t.Fatalf("Failed to create script: %v", err)
	}

	content, err := os.ReadFile(scriptPath)
	if err != nil {
		t.Fatalf("Failed to read script: %v", err)
	}
	scriptContent := string(content)

	// Rusage should be preserved (span stripped, rusage kept)
	if !strings.Contains(scriptContent, "rusage[mem=4096]") {
		t.Errorf("Script should preserve rusage after stripping span\nScript:\n%s", scriptContent)
	}

	// Span should be regenerated
	if !strings.Contains(scriptContent, `span[hosts=2]`) {
		t.Errorf("Script should contain regenerated span[hosts=2]\nScript:\n%s", scriptContent)
	}
}

func TestLsfNodeTaskRoundTrip(t *testing.T) {
	lsf := newTestLsfScheduler()

	// Create a script with span and -n directives
	lines := []string{
		"#!/bin/bash",
		"#BSUB -J roundtrip_test",
		"#BSUB -n 8",
		`#BSUB -R "span[hosts=2]"`,
		"#BSUB -W 02:00",
		"echo hello",
	}
	tmpFile := filepath.Join(t.TempDir(), "input.sh")
	os.WriteFile(tmpFile, []byte(strings.Join(lines, "\n")), 0644)

	// Parse
	specs, err := lsf.ReadScriptSpecs(tmpFile)
	if err != nil {
		t.Fatalf("parse failed: %v", err)
	}
	if specs.Nodes != 2 {
		t.Fatalf("parse Nodes = %d; want 2", specs.Nodes)
	}
	if specs.Ncpus != 8 {
		t.Fatalf("parse Ncpus = %d; want 8", specs.Ncpus)
	}

	// Generate
	outputDir := t.TempDir()
	jobSpec := &JobSpec{
		Name:    "roundtrip_test",
		Command: "echo hello",
		Specs:   specs,
	}
	scriptPath, err := lsf.CreateScriptWithSpec(jobSpec, outputDir)
	if err != nil {
		t.Fatalf("generate failed: %v", err)
	}

	// Re-parse generated script
	specs2, err := lsf.ReadScriptSpecs(scriptPath)
	if err != nil {
		t.Fatalf("re-parse failed: %v", err)
	}
	if specs2.Nodes != 2 {
		t.Errorf("round-trip Nodes = %d; want 2", specs2.Nodes)
	}
	if specs2.Ncpus != 8 {
		t.Errorf("round-trip Ncpus = %d; want 8", specs2.Ncpus)
	}

	// Check no duplicates in generated script
	content, _ := os.ReadFile(scriptPath)
	script := string(content)
	if strings.Count(script, "span[hosts=") != 1 {
		t.Error("span[hosts=] should appear exactly once")
	}
	if strings.Count(script, "#BSUB -n ") != 1 {
		t.Error("#BSUB -n should appear exactly once")
	}
}

func TestLsfIsAvailable(t *testing.T) {
	t.Run("not in job", func(t *testing.T) {
		clearJobEnvVars(t)
		lsf := newTestLsfScheduler()
		if !lsf.IsAvailable() {
			t.Error("Expected IsAvailable to return true when not in a job")
		}
	})

	t.Run("inside job", func(t *testing.T) {
		clearJobEnvVars(t)
		t.Setenv("LSB_JOBID", "99999")
		lsf := newTestLsfScheduler()
		if lsf.IsAvailable() {
			t.Error("Expected IsAvailable to return false when inside a job")
		}
	})

	t.Run("no binary", func(t *testing.T) {
		clearJobEnvVars(t)
		lsf := &LsfScheduler{}
		if lsf.IsAvailable() {
			t.Error("Expected IsAvailable to return false when no binary is set")
		}
	})
}

func TestLsfGetInfo(t *testing.T) {
	clearJobEnvVars(t)
	lsf := newTestLsfScheduler()

	info := lsf.GetInfo()
	if info.Type != "LSF" {
		t.Errorf("Type = %q; want %q", info.Type, "LSF")
	}
	if info.Binary != "/usr/bin/bsub" {
		t.Errorf("Binary = %q; want %q", info.Binary, "/usr/bin/bsub")
	}
	if info.InJob {
		t.Error("InJob should be false")
	}
	if !info.Available {
		t.Error("Available should be true")
	}
}

func TestTryParseLsfScript(t *testing.T) {
	tmpDir := t.TempDir()

	script := `#!/bin/bash
#BSUB -J testjob
#BSUB -n 4
#BSUB -M 8GB
#BSUB -W 02:00
#BSUB -B
#BSUB -N
#BSUB -u user@example.com

echo "Running job"
`
	scriptPath := filepath.Join(tmpDir, "test.sh")
	if err := os.WriteFile(scriptPath, []byte(script), 0644); err != nil {
		t.Fatalf("Failed to create test script: %v", err)
	}

	specs, err := TryParseLsfScript(scriptPath)
	if err != nil {
		t.Fatalf("TryParseLsfScript failed: %v", err)
	}

	if specs.Ncpus != 4 {
		t.Errorf("Ncpus = %d; want 4", specs.Ncpus)
	}
	if specs.MemMB != 8*1024 {
		t.Errorf("MemMB = %d; want %d", specs.MemMB, 8*1024)
	}
	if specs.Time != 2*time.Hour {
		t.Errorf("Time = %v; want 2h", specs.Time)
	}
	if !specs.EmailOnBegin {
		t.Error("EmailOnBegin should be true")
	}
	if !specs.EmailOnEnd {
		t.Error("EmailOnEnd should be true")
	}
	if specs.MailUser != "user@example.com" {
		t.Errorf("MailUser = %q; want %q", specs.MailUser, "user@example.com")
	}
	// RawFlags should be empty - all flags above are recognized and parsed into typed fields
	if len(specs.RawFlags) != 0 {
		t.Errorf("RawFlags count = %d; want 0 (all flags were recognized)", len(specs.RawFlags))
	}
}

func TestLsfGetJobResources(t *testing.T) {
	sched := &LsfScheduler{}

	t.Run("not in job", func(t *testing.T) {
		clearJobEnvVars(t)
		if res := sched.GetJobResources(); res != nil {
			t.Fatalf("expected nil, got %+v", res)
		}
	})

	t.Run("full resources", func(t *testing.T) {
		clearJobEnvVars(t)
		t.Setenv("LSB_JOBID", "99999")
		t.Setenv("LSB_DJOB_NUMPROC", "32")
		t.Setenv("LSB_MAX_MEM_RUSAGE", "8388608") // 8 GB in KB
		t.Setenv("CUDA_VISIBLE_DEVICES", "0")

		res := sched.GetJobResources()
		if res == nil {
			t.Fatal("expected non-nil")
		}
		if res.Ncpus == nil || *res.Ncpus != 32 {
			t.Errorf("Ncpus = %v; want 32", res.Ncpus)
		}
		if res.MemMB == nil || *res.MemMB != 8192 {
			t.Errorf("MemMB = %v; want 8192", res.MemMB)
		}
		if res.Ngpus == nil || *res.Ngpus != 1 {
			t.Errorf("Ngpus = %v; want 1", res.Ngpus)
		}
	})

	t.Run("LSB_MAX_NUM_PROCESSORS fallback", func(t *testing.T) {
		clearJobEnvVars(t)
		t.Setenv("LSB_JOBID", "99999")
		t.Setenv("LSB_MAX_NUM_PROCESSORS", "64")

		res := sched.GetJobResources()
		if res == nil {
			t.Fatal("expected non-nil")
		}
		if res.Ncpus == nil || *res.Ncpus != 64 {
			t.Errorf("Ncpus = %v; want 64", res.Ncpus)
		}
	})

	t.Run("partial data", func(t *testing.T) {
		clearJobEnvVars(t)
		t.Setenv("LSB_JOBID", "99999")
		t.Setenv("LSB_DJOB_NUMPROC", "4")

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
		t.Setenv("LSB_JOBID", "99999")
		t.Setenv("LSB_DJOB_NUMPROC", "not-a-number")
		t.Setenv("LSB_MAX_MEM_RUSAGE", "-100")

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
