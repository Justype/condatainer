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
			RemainingFlags: []string{},
			Spec:           &ResourceSpec{CpusPerTask: 1, Nodes: 1},
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
					Control: RuntimeConfig{
						JobName:      "test_job",
						EmailOnBegin: tt.emailOnBegin,
						EmailOnEnd:   tt.emailOnEnd,
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
	if !specs.Control.EmailOnBegin || !specs.Control.EmailOnEnd {
		t.Errorf("Parsing failed: EmailOnBegin=%v, EmailOnEnd=%v",
			specs.Control.EmailOnBegin, specs.Control.EmailOnEnd)
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
	// No span: -n 8 → MPI free-dist (Ntasks=8, CpusPerTask=1)
	if specs.Spec.Ntasks != 8 {
		t.Errorf("Ntasks = %d; want 8 (MPI free-dist)", specs.Spec.Ntasks)
	}
	if specs.Spec.CpusPerTask != 1 {
		t.Errorf("CpusPerTask = %d; want 1 (MPI free-dist)", specs.Spec.CpusPerTask)
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
			// no-span: -n=4 → Ntasks=4 (MPI free-dist), CpusPerTask=1, MemPerCpuMB=8MB
			wantNcpus: 1,
			wantMemMB: 0, // per-node undetermined without TasksPerNode; MemPerCpuMB=8MB
			wantTime:  time.Hour + 30*time.Minute,
		},
		{
			name: "ncpus only",
			lines: []string{
				"#!/bin/bash",
				"#BSUB -n 16",
			},
			// no-span: -n=16 → Ntasks=16 (MPI free-dist), CpusPerTask=1
			wantNcpus: 1,
		},
		{
			name: "walltime only",
			lines: []string{
				"#!/bin/bash",
				"#BSUB -W 04:00",
			},
			wantNcpus: 1, // default
			wantTime:  4 * time.Hour,
		},
		{
			name: "gpu directive",
			lines: []string{
				"#!/bin/bash",
				"#BSUB -n 8",
				`#BSUB -gpu "num=2:type=a100"`,
			},
			// no-span: -n=8 → Ntasks=8 (MPI free-dist), CpusPerTask=1
			wantNcpus: 1,
			wantGpu:   &GpuSpec{Type: "a100", Count: 2},
		},
		{
			name: "gpu num only",
			lines: []string{
				"#!/bin/bash",
				`#BSUB -gpu "num=4"`,
			},
			wantNcpus: 1, // default
			wantGpu:   &GpuSpec{Type: "gpu", Count: 4},
		},
		{
			name: "rusage mem",
			lines: []string{
				"#!/bin/bash",
				"#BSUB -n 4",
				`#BSUB -R "rusage[mem=4096]"`,
			},
			// no-span: -n=4 → Ntasks=4 (MPI free-dist), CpusPerTask=1, MemPerCpuMB=4096MB
			wantNcpus: 1,
			wantMemMB: 0, // per-node undetermined without TasksPerNode; MemPerCpuMB=4096MB
		},
		{
			name: "rusage ngpus_physical",
			lines: []string{
				"#!/bin/bash",
				"#BSUB -n 8",
				`#BSUB -R "rusage[ngpus_physical=2]"`,
			},
			// no-span: -n=8 → Ntasks=8 (MPI free-dist), CpusPerTask=1
			wantNcpus: 1,
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
			// no-span: -n=32 → Ntasks=32 (MPI free-dist), CpusPerTask=1, MemPerCpuMB=64MB
			wantNcpus: 1,
			wantMemMB: 0, // per-node undetermined without TasksPerNode; MemPerCpuMB=64MB
			wantTime:  12 * time.Hour,
			wantGpu:   &GpuSpec{Type: "v100", Count: 4},
		},
		{
			name: "memory with MB suffix",
			lines: []string{
				"#!/bin/bash",
				"#BSUB -M 4096MB",
			},
			wantNcpus: 1, // default
			wantMemMB: 4096,
		},
		{
			name: "memory with GB suffix",
			lines: []string{
				"#!/bin/bash",
				"#BSUB -M 8GB",
			},
			wantNcpus: 1, // default
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

			if specs.Spec == nil {
				t.Fatal("Spec is nil; want non-nil")
			}
			if specs.Spec.CpusPerTask != tt.wantNcpus {
				t.Errorf("CpusPerTask = %d; want %d", specs.Spec.CpusPerTask, tt.wantNcpus)
			}
			if tt.wantMemMB > 0 && specs.Spec.GetMemPerNodeMB() != tt.wantMemMB {
				t.Errorf("GetMemPerNodeMB() = %d; want %d (MemPerCpuMB=%d MemPerNodeMB=%d)", specs.Spec.GetMemPerNodeMB(), tt.wantMemMB, specs.Spec.MemPerCpuMB, specs.Spec.MemPerNodeMB)
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
		{
			name:      "gmodel key (LSF native format)",
			input:     "num=4:gmodel=a100",
			wantType:  "a100",
			wantCount: 4,
		},
		{
			name:      "gmodel with extra options",
			input:     "num=2:gmodel=h100:mode=exclusive_process",
			wantType:  "h100",
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
			wantMemMB: 8192, // 8192 MB (bare number in rusage[] defaults to MB per LSF docs)
			wantGpu:   &GpuSpec{Type: "gpu", Count: 2},
		},
		{
			name: "rusage mem colon separated",
			lines: []string{
				"#!/bin/bash",
				`#BSUB -R "rusage[mem=4096:ngpus_physical=1]"`,
			},
			wantMemMB: 4096, // 4096 MB
			wantGpu:   &GpuSpec{Type: "gpu", Count: 1},
		},
		{
			name: "rusage with span",
			lines: []string{
				"#!/bin/bash",
				`#BSUB -R "span[hosts=1] rusage[mem=16384]"`,
			},
			wantMemMB: 16384, // 16384 MB
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

			if specs.Spec == nil {
				t.Fatal("Spec is nil; want non-nil")
			}
			if tt.wantMemMB > 0 && specs.Spec.GetMemPerNodeMB() != tt.wantMemMB {
				t.Errorf("GetMemPerNodeMB() = %d; want %d (MemPerCpuMB=%d MemPerNodeMB=%d)", specs.Spec.GetMemPerNodeMB(), tt.wantMemMB, specs.Spec.MemPerCpuMB, specs.Spec.MemPerNodeMB)
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

func TestLsfSpanHostsParsing(t *testing.T) {
	tests := []struct {
		name             string
		lines            []string
		wantPassthrough  bool // true = expect Spec == nil (passthrough mode)
		wantNodes        int  // 0 = don't check
		wantNodesZero    bool // true = assert Nodes == 0 (no-span: no span emitted)
		wantTasksPerNode int  // 0 = don't check
		wantNtasks       int  // 0 = don't check
		wantCpusPerTask  int  // 0 = don't check
	}{
		{
			name: "span hosts=1: single-node OpenMP",
			lines: []string{
				"#!/bin/bash",
				"#BSUB -n 4",
				`#BSUB -R "span[hosts=1]"`,
			},
			wantNodes:        1,
			wantTasksPerNode: 1,
			wantNtasks:       1, // OpenMP: 1 MPI rank
			wantCpusPerTask:  4,
		},
		{
			name: "span ptile=4: pure MPI, 16 total tasks → 4 nodes",
			lines: []string{
				"#!/bin/bash",
				"#BSUB -n 16",
				`#BSUB -R "span[ptile=4]"`,
			},
			wantNodes:        4,
			wantTasksPerNode: 4,
			wantNtasks:       16,
			wantCpusPerTask:  1, // pure MPI: single-threaded tasks
		},
		{
			name: "span ptile=4 affinity[cores(8)]: MPI+OpenMP",
			lines: []string{
				"#!/bin/bash",
				"#BSUB -n 16",
				`#BSUB -R "span[ptile=4] affinity[cores(8)]"`,
			},
			wantNodes:        4,
			wantTasksPerNode: 4,
			wantNtasks:       16,
			wantCpusPerTask:  8,
		},
		{
			name: "span ptile non-divisible: passthrough (even distribution required)",
			lines: []string{
				"#!/bin/bash",
				"#BSUB -n 10",
				`#BSUB -R "span[ptile=3]"`,
			},
			wantPassthrough: true, // 10 % 3 != 0 → passthrough
		},
		{
			name: "span ptile with rusage",
			lines: []string{
				"#!/bin/bash",
				"#BSUB -n 8",
				`#BSUB -R "span[ptile=4] rusage[mem=8192]"`,
			},
			wantNodes:        2,
			wantTasksPerNode: 4,
			wantNtasks:       8,
			wantCpusPerTask:  1, // pure MPI (no affinity)
		},
		{
			name: "span hosts=1 affinity[cores(T)]: hybrid single-node",
			lines: []string{
				"#!/bin/bash",
				"#BSUB -n 4",
				`#BSUB -R "span[hosts=1] affinity[cores(8)]"`,
			},
			wantNodes:        1,
			wantTasksPerNode: 4, // 4 processes on 1 node
			wantNtasks:       4, // -n 4 = 4 processes
			wantCpusPerTask:  8, // affinity[cores(8)] = 8 cores per process
		},
		{
			name: "span hosts=1 on separate -R line",
			lines: []string{
				"#!/bin/bash",
				"#BSUB -n 32",
				`#BSUB -R "span[hosts=1]"`,
				`#BSUB -R "rusage[mem=16384]"`,
			},
			wantNodes:        1,
			wantTasksPerNode: 1,
			wantNtasks:       1, // OpenMP: 1 MPI rank
			wantCpusPerTask:  32,
		},
		{
			name: "span hosts>1 is ignored (warns): falls back to no-span MPI",
			lines: []string{
				"#!/bin/bash",
				"#BSUB -n 8",
				`#BSUB -R "span[hosts=4]"`,
			},
			wantNodesZero:   true, // no-span: Nodes=0, free distribution
			wantNtasks:      8,    // -n = total tasks
			wantCpusPerTask: 1,    // pure MPI: one thread per task
		},
		{
			name: "no span directive: free-distribution MPI",
			lines: []string{
				"#!/bin/bash",
				"#BSUB -n 16",
				`#BSUB -R "rusage[mem=4096]"`,
			},
			wantNodesZero:   true, // no-span: Nodes=0, free distribution
			wantNtasks:      16,   // -n = total tasks
			wantCpusPerTask: 1,    // pure MPI: one thread per task
		},
		{
			name: "no resource directives at all",
			lines: []string{
				"#!/bin/bash",
				"#BSUB -J testjob",
			},
			wantNodes:        0, // passthrough
			wantTasksPerNode: 0, // passthrough
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

			if tt.wantPassthrough {
				if specs.Spec != nil {
					t.Errorf("Spec = %+v; want nil (passthrough)", specs.Spec)
				}
				return
			}

			if specs.Spec == nil {
				t.Fatal("Spec is nil; want non-nil")
			}
			if tt.wantNodesZero {
				if specs.Spec.Nodes != 0 {
					t.Errorf("Nodes = %d; want 0 (no-span)", specs.Spec.Nodes)
				}
			} else if tt.wantNodes > 0 {
				if specs.Spec.Nodes != tt.wantNodes {
					t.Errorf("Nodes = %d; want %d", specs.Spec.Nodes, tt.wantNodes)
				}
			}
			if tt.wantTasksPerNode > 0 && specs.Spec.TasksPerNode != tt.wantTasksPerNode {
				t.Errorf("TasksPerNode = %d; want %d", specs.Spec.TasksPerNode, tt.wantTasksPerNode)
			}
			if tt.wantNtasks > 0 && specs.Spec.Ntasks != tt.wantNtasks {
				t.Errorf("Ntasks = %d; want %d", specs.Spec.Ntasks, tt.wantNtasks)
			}
			if tt.wantCpusPerTask > 0 && specs.Spec.CpusPerTask != tt.wantCpusPerTask {
				t.Errorf("CpusPerTask = %d; want %d", specs.Spec.CpusPerTask, tt.wantCpusPerTask)
			}
		})
	}
}

func TestLsfScriptGenerationNodesCpus(t *testing.T) {
	tests := []struct {
		name      string
		specs     *ScriptSpecs
		wantSpan  string // expected span[...] content
		wantNcpus string // expected -n N line
	}{
		{
			name: "single-node OpenMP",
			specs: &ScriptSpecs{
				Spec:           &ResourceSpec{CpusPerTask: 4, Nodes: 1, TasksPerNode: 1},
				RemainingFlags: []string{},
			},
			wantSpan:  `span[hosts=1]`,
			wantNcpus: "#BSUB -n 4",
		},
		{
			name: "pure MPI: 4 nodes × 8 tasks/node = 32 total",
			specs: &ScriptSpecs{
				Spec:           &ResourceSpec{Nodes: 4, TasksPerNode: 8, CpusPerTask: 1},
				RemainingFlags: []string{},
			},
			wantSpan:  `span[ptile=8]`,
			wantNcpus: "#BSUB -n 32",
		},
		{
			name: "MPI+OpenMP: 4 nodes × 4 tasks/node × 8 threads",
			specs: &ScriptSpecs{
				Spec:           &ResourceSpec{Nodes: 4, TasksPerNode: 4, CpusPerTask: 8},
				RemainingFlags: []string{},
			},
			wantSpan:  `span[ptile=4] affinity[cores(8)]`,
			wantNcpus: "#BSUB -n 16",
		},
		{
			name: "explicit Ntasks used in -n output",
			specs: &ScriptSpecs{
				Spec:           &ResourceSpec{Ntasks: 16, Nodes: 4, TasksPerNode: 4},
				RemainingFlags: []string{},
			},
			wantSpan:  `span[ptile=4]`,
			wantNcpus: "#BSUB -n 16",
		},
		{
			name: "single-node from Ncpus field only",
			specs: &ScriptSpecs{
				Spec:           &ResourceSpec{CpusPerTask: 8, Nodes: 1},
				RemainingFlags: []string{},
			},
			wantSpan:  `span[hosts=1]`,
			wantNcpus: "#BSUB -n 8",
		},
		{
			name: "rusage preserved alongside new span",
			specs: &ScriptSpecs{
				Spec:           &ResourceSpec{Nodes: 3, TasksPerNode: 4, CpusPerTask: 1},
				RemainingFlags: []string{`-R "rusage[mem=4096]"`},
			},
			wantSpan:  `span[ptile=4]`,
			wantNcpus: "#BSUB -n 12",
		},
		{
			name: "free-dist pure MPI: Ntasks set, TasksPerNode=0",
			specs: &ScriptSpecs{
				Spec:           &ResourceSpec{Ntasks: 14, CpusPerTask: 1},
				RemainingFlags: []string{},
			},
			wantSpan:  "", // no span[ptile=] emitted
			wantNcpus: "#BSUB -n 14",
		},
		{
			name: "free-dist hybrid: Ntasks set, CpusPerTask>1, TasksPerNode=0",
			specs: &ScriptSpecs{
				Spec:           &ResourceSpec{Ntasks: 14, CpusPerTask: 4},
				RemainingFlags: []string{},
			},
			wantSpan:  "affinity[cores(4)]",
			wantNcpus: "#BSUB -n 14",
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

			// Check span directive (empty wantSpan means assert span[ptile= is absent)
			if tt.wantSpan == "" {
				if strings.Contains(scriptContent, "span[ptile=") {
					t.Errorf("Script contains unexpected span[ptile=]\nScript:\n%s", scriptContent)
				}
			} else if !strings.Contains(scriptContent, tt.wantSpan) {
				t.Errorf("Script missing expected span: %q\nScript:\n%s", tt.wantSpan, scriptContent)
			}

			// Check -n directive
			if !strings.Contains(scriptContent, tt.wantNcpus) {
				t.Errorf("Script missing expected -n: %q\nScript:\n%s", tt.wantNcpus, scriptContent)
			}

			// Verify no duplicate -n directives
			nCount := strings.Count(scriptContent, "#BSUB -n ")
			if nCount != 1 {
				t.Errorf("#BSUB -n appears %d times; want 1\nScript:\n%s", nCount, scriptContent)
			}
		})
	}
}

func TestLsfScriptGenerationGpu(t *testing.T) {
	tests := []struct {
		name       string
		gpu        *GpuSpec
		wantLine   string // expected in output; empty means assert absent
		wantAbsent string // must NOT appear in output
		wantCount  int    // expected occurrence count of wantLine (0 = use absent check)
	}{
		{
			name:       "no GPU",
			gpu:        nil,
			wantLine:   "",
			wantAbsent: "#BSUB -gpu",
		},
		{
			name:       "GPU count only (generic type)",
			gpu:        &GpuSpec{Type: "gpu", Count: 2},
			wantLine:   `#BSUB -gpu "num=2"`,
			wantAbsent: "gmodel=",
			wantCount:  1,
		},
		{
			name:      "GPU with model type",
			gpu:       &GpuSpec{Type: "a100", Count: 4},
			wantLine:  `#BSUB -gpu "num=4:gmodel=a100"`,
			wantCount: 1,
		},
		{
			name:      "GPU with type — no duplicate",
			gpu:       &GpuSpec{Type: "h100", Count: 1},
			wantLine:  `#BSUB -gpu "num=1:gmodel=h100"`,
			wantCount: 1,
		},
	}

	for _, tt := range tests {
		t.Run(tt.name, func(t *testing.T) {
			tmpDir := t.TempDir()
			lsf := newTestLsfScheduler()

			jobSpec := &JobSpec{
				Name:    "test_gpu",
				Command: "echo 'test'",
				Specs: &ScriptSpecs{
					Spec:           &ResourceSpec{CpusPerTask: 4, Gpu: tt.gpu},
					RemainingFlags: []string{},
				},
			}

			scriptPath, err := lsf.CreateScriptWithSpec(jobSpec, tmpDir)
			if err != nil {
				t.Fatalf("CreateScriptWithSpec failed: %v", err)
			}
			content, err := os.ReadFile(scriptPath)
			if err != nil {
				t.Fatalf("Failed to read script: %v", err)
			}
			scriptContent := string(content)

			if tt.wantAbsent != "" && strings.Contains(scriptContent, tt.wantAbsent) {
				t.Errorf("Script contains unexpected %q\nScript:\n%s", tt.wantAbsent, scriptContent)
			}
			if tt.wantLine != "" {
				if !strings.Contains(scriptContent, tt.wantLine) {
					t.Errorf("Script missing %q\nScript:\n%s", tt.wantLine, scriptContent)
				}
				if tt.wantCount > 0 {
					n := strings.Count(scriptContent, tt.wantLine)
					if n != tt.wantCount {
						t.Errorf("%q appears %d times; want %d\nScript:\n%s", tt.wantLine, n, tt.wantCount, scriptContent)
					}
				}
			}
		})
	}
}

func TestLsfScriptGenerationMem(t *testing.T) {
	tests := []struct {
		name       string
		spec       *ResourceSpec
		wantLine   string
		wantAbsent string
	}{
		{
			name:       "no memory",
			spec:       &ResourceSpec{Nodes: 1, CpusPerTask: 4},
			wantAbsent: "rusage[mem=",
		},
		{
			name: "MemPerNodeMB only — OpenMP",
			spec: &ResourceSpec{Nodes: 1, CpusPerTask: 8, MemPerNodeMB: 16384},
			// slotsPerNode = CpusPerTask = 8; per-slot = 16384/8 = 2048
			wantLine: `#BSUB -R "span[hosts=1] rusage[mem=2048MB]"`,
		},
		{
			name: "MemPerCpuMB only — OpenMP (single node)",
			spec: &ResourceSpec{Nodes: 1, CpusPerTask: 8, MemPerCpuMB: 2048},
			// cpuPerSlot = 1 (not hybrid); rusage_mem = 2048
			wantLine: `#BSUB -R "span[hosts=1] rusage[mem=2048MB]"`,
		},
		{
			name: "MemPerCpuMB only — pure MPI",
			spec: &ResourceSpec{Nodes: 2, TasksPerNode: 4, CpusPerTask: 1, MemPerCpuMB: 4096},
			// cpuPerSlot = 1 (CpusPerTask <= 1); rusage_mem = 4096
			wantLine: `#BSUB -R "span[ptile=4] rusage[mem=4096MB]"`,
		},
		{
			name: "MemPerCpuMB only — hybrid MPI+OpenMP",
			spec: &ResourceSpec{Nodes: 2, TasksPerNode: 4, CpusPerTask: 4, MemPerCpuMB: 1024},
			// cpuPerSlot = CpusPerTask = 4; rusage_mem = 1024*4 = 4096
			wantLine: `#BSUB -R "span[ptile=4] affinity[cores(4)] rusage[mem=4096MB]"`,
		},
		{
			name: "MemPerCpuMB — hybrid single-node (Ntasks=TasksPerNode=Nodes=1 forced via Nodes=1)",
			// Nodes=1, TasksPerNode=4, CpusPerTask=8 → isMPI=true, TasksPerNode>0 → ptile path
			spec: &ResourceSpec{Nodes: 1, TasksPerNode: 4, CpusPerTask: 8, MemPerCpuMB: 512},
			// cpuPerSlot = 8; rusage_mem = 512*8 = 4096 MB per slot
			wantLine: `#BSUB -R "span[ptile=4] affinity[cores(8)] rusage[mem=4096MB]"`,
		},
	}

	for _, tt := range tests {
		t.Run(tt.name, func(t *testing.T) {
			tmpDir := t.TempDir()
			lsf := newTestLsfScheduler()

			jobSpec := &JobSpec{
				Name:    "test_mem",
				Command: "echo 'test'",
				Specs: &ScriptSpecs{
					Spec:           tt.spec,
					RemainingFlags: []string{},
				},
			}

			scriptPath, err := lsf.CreateScriptWithSpec(jobSpec, tmpDir)
			if err != nil {
				t.Fatalf("CreateScriptWithSpec failed: %v", err)
			}
			content, err := os.ReadFile(scriptPath)
			if err != nil {
				t.Fatalf("Failed to read script: %v", err)
			}
			scriptContent := string(content)

			if tt.wantAbsent != "" && strings.Contains(scriptContent, tt.wantAbsent) {
				t.Errorf("Script contains unexpected %q\nScript:\n%s", tt.wantAbsent, scriptContent)
			}
			if tt.wantLine != "" && !strings.Contains(scriptContent, tt.wantLine) {
				t.Errorf("Script missing %q\nScript:\n%s", tt.wantLine, scriptContent)
			}
		})
	}
}

func TestLsfMemPerCpuRoundTrip(t *testing.T) {
	tests := []struct {
		name     string
		script   string
		wantLine string
	}{
		{
			name:   "OpenMP: rusage mem round-trips via MemPerCpuMB",
			script: "#!/bin/bash\n#BSUB -J rt_openmp\n#BSUB -n 8\n#BSUB -R \"span[hosts=1] rusage[mem=4096MB]\"\necho hello\n",
			// parsed: singleNode+no-affinity → CpusPerTask=8, Ntasks=1; MemPerCpuMB=4096 (slot=core)
			// generated: cpuPerSlot=1 (isMPI=false) → rusage[mem=4096MB]
			wantLine: `#BSUB -R "span[hosts=1] rusage[mem=4096MB]"`,
		},
		{
			name:   "Hybrid single-node: rusage per slot (=per process) round-trips",
			script: "#!/bin/bash\n#BSUB -J rt_hybrid_sn\n#BSUB -n 4\n#BSUB -R \"span[hosts=1] affinity[cores(8)] rusage[mem=8192MB]\"\necho hello\n",
			// parsed: singleNode+affinity(8) → Ntasks=4, TasksPerNode=4, CpusPerTask=8
			//         MemPerCpuMB = 8192/8 = 1024 (per core)
			// generated: isMPI=true, TasksPerNode=4 → span[ptile=4]; cpuPerSlot=8 → rusage[mem=8192MB]
			wantLine: `#BSUB -R "span[ptile=4] affinity[cores(8)] rusage[mem=8192MB]"`,
		},
		{
			name:   "Pure MPI: rusage mem round-trips via MemPerCpuMB",
			script: "#!/bin/bash\n#BSUB -J rt_mpi\n#BSUB -n 16\n#BSUB -R \"span[ptile=4] rusage[mem=2048MB]\"\necho hello\n",
			// parsed: MemPerCpuMB=2048 (rusage/cpt=1), Ntasks=16, TasksPerNode=4, Nodes=4
			// generated: cpuPerSlot=1 → rusage[mem=2048MB]
			wantLine: `#BSUB -R "span[ptile=4] rusage[mem=2048MB]"`,
		},
		{
			name:   "Hybrid MPI+OpenMP: rusage mem round-trips via MemPerCpuMB",
			script: "#!/bin/bash\n#BSUB -J rt_hybrid\n#BSUB -n 8\n#BSUB -R \"span[ptile=4] affinity[cores(4)] rusage[mem=4096MB]\"\necho hello\n",
			// parsed: MemPerCpuMB=4096/4=1024, CpusPerTask=4, TasksPerNode=4, Ntasks=8, Nodes=2
			// generated: cpuPerSlot=4 → rusage_mem=1024*4=4096MB
			wantLine: `#BSUB -R "span[ptile=4] affinity[cores(4)] rusage[mem=4096MB]"`,
		},
	}

	for _, tt := range tests {
		t.Run(tt.name, func(t *testing.T) {
			tmpDir := t.TempDir()
			lsf := newTestLsfScheduler()

			scriptPath := filepath.Join(tmpDir, "original.sh")
			if err := os.WriteFile(scriptPath, []byte(tt.script), 0644); err != nil {
				t.Fatalf("Write original script: %v", err)
			}

			specs, err := lsf.ReadScriptSpecs(scriptPath)
			if err != nil {
				t.Fatalf("ReadScriptSpecs: %v", err)
			}
			if specs.Spec == nil {
				t.Fatal("Spec is nil")
			}
			if specs.Spec.MemPerCpuMB == 0 {
				t.Errorf("Expected MemPerCpuMB > 0 after parsing, got 0 (MemPerNodeMB=%d)",
					specs.Spec.MemPerNodeMB)
			}

			jobSpec := &JobSpec{
				Name:    "rt_test",
				Command: "echo hello",
				Specs:   specs,
			}
			newPath, err := lsf.CreateScriptWithSpec(jobSpec, tmpDir)
			if err != nil {
				t.Fatalf("CreateScriptWithSpec: %v", err)
			}
			content, err := os.ReadFile(newPath)
			if err != nil {
				t.Fatalf("Read generated script: %v", err)
			}
			scriptContent := string(content)

			if !strings.Contains(scriptContent, tt.wantLine) {
				t.Errorf("Generated script missing %q\nScript:\n%s", tt.wantLine, scriptContent)
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
			Spec:           &ResourceSpec{TasksPerNode: 4, Nodes: 2},
			RemainingFlags: []string{`-R "rusage[mem=4096]"`},
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

	// Rusage should be preserved in the RemainingFlags passthrough
	if !strings.Contains(scriptContent, "rusage[mem=4096]") {
		t.Errorf("Script should preserve rusage\nScript:\n%s", scriptContent)
	}

	// Span should be regenerated as ptile (multi-node MPI)
	if !strings.Contains(scriptContent, `span[ptile=4]`) {
		t.Errorf("Script should contain regenerated span[ptile=4]\nScript:\n%s", scriptContent)
	}
}

func TestLsfNodeTaskRoundTrip(t *testing.T) {
	lsf := newTestLsfScheduler()

	// Create a script with span[ptile=4] and -n 16 → 4 nodes × 4 tasks/node
	lines := []string{
		"#!/bin/bash",
		"#BSUB -J roundtrip_test",
		"#BSUB -n 16",
		`#BSUB -R "span[ptile=4]"`,
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
	if specs.Spec == nil {
		t.Fatal("Spec is nil; want non-nil")
	}
	if specs.Spec.Nodes != 4 {
		t.Fatalf("parse Nodes = %d; want 4", specs.Spec.Nodes)
	}
	if specs.Spec.TasksPerNode != 4 {
		t.Fatalf("parse TasksPerNode = %d; want 4", specs.Spec.TasksPerNode)
	}
	if specs.Spec.Ntasks != 16 {
		t.Fatalf("parse Ntasks = %d; want 16", specs.Spec.Ntasks)
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
	if specs2.Spec == nil {
		t.Fatal("specs2.Spec is nil; want non-nil")
	}
	if specs2.Spec.Nodes != 4 {
		t.Errorf("round-trip Nodes = %d; want 4", specs2.Spec.Nodes)
	}
	if specs2.Spec.TasksPerNode != 4 {
		t.Errorf("round-trip TasksPerNode = %d; want 4", specs2.Spec.TasksPerNode)
	}
	if specs2.Spec.Ntasks != 16 {
		t.Errorf("round-trip Ntasks = %d; want 16", specs2.Spec.Ntasks)
	}

	// Check no duplicates in generated script
	content, _ := os.ReadFile(scriptPath)
	script := string(content)
	if strings.Count(script, "span[ptile=") != 1 {
		t.Error("span[ptile=] should appear exactly once")
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

	if specs.Spec == nil {
		t.Fatal("Spec is nil; want non-nil")
	}
	// No span: -n 4 → MPI free-dist (Ntasks=4, CpusPerTask=1)
	if specs.Spec.Ntasks != 4 {
		t.Errorf("Ntasks = %d; want 4 (MPI free-dist)", specs.Spec.Ntasks)
	}
	if specs.Spec.CpusPerTask != 1 {
		t.Errorf("CpusPerTask = %d; want 1 (MPI free-dist)", specs.Spec.CpusPerTask)
	}
	// -M 8GB → MemPerCpuMB=8192 MB (per task); GetMemPerNodeMB()=0 since TasksPerNode=0 (free-dist)
	if specs.Spec.MemPerCpuMB != 8192 {
		t.Errorf("MemPerCpuMB = %d; want 8192 (8 GB per task)", specs.Spec.MemPerCpuMB)
	}
	if specs.Spec.Time != 2*time.Hour {
		t.Errorf("Time = %v; want 2h", specs.Spec.Time)
	}
	if !specs.Control.EmailOnBegin {
		t.Error("EmailOnBegin should be true")
	}
	if !specs.Control.EmailOnEnd {
		t.Error("EmailOnEnd should be true")
	}
	if specs.Control.MailUser != "user@example.com" {
		t.Errorf("MailUser = %q; want %q", specs.Control.MailUser, "user@example.com")
	}
	// -M passes through to RemainingFlags (it is a limit/ulimit, not an allocation directive)
	if len(specs.RemainingFlags) != 1 || specs.RemainingFlags[0] != "-M 8GB" {
		t.Errorf("RemainingFlags = %v; want [\"-M 8GB\"]", specs.RemainingFlags)
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

	// Pure MPI: LSB_DJOB_NUMPROC = total tasks; no NCPUS injected.
	t.Run("pure MPI via LSB_DJOB_NUMPROC", func(t *testing.T) {
		clearJobEnvVars(t)
		t.Setenv("LSB_JOBID", "99999")
		t.Setenv("LSB_DJOB_NUMPROC", "32")
		t.Setenv("CUDA_VISIBLE_DEVICES", "0")

		res := sched.GetJobResources()
		if res == nil {
			t.Fatal("expected non-nil")
		}
		if res.Ntasks != 32 {
			t.Errorf("Ntasks = %d; want 32", res.Ntasks)
		}
		if res.CpusPerTask != 1 {
			t.Errorf("CpusPerTask = %d; want 1", res.CpusPerTask)
		}
		if res.Gpu == nil || res.Gpu.Count != 1 {
			t.Errorf("Gpu.Count = %v; want 1", res.Gpu)
		}
	})

	// Hybrid: injected NCPUS splits the total slot count correctly.
	t.Run("hybrid MPI+OpenMP via injected NCPUS", func(t *testing.T) {
		clearJobEnvVars(t)
		t.Setenv("LSB_JOBID", "99999")
		t.Setenv("LSB_DJOB_NUMPROC", "32") // 8 tasks × 4 threads
		t.Setenv("NCPUS", "4")

		res := sched.GetJobResources()
		if res == nil {
			t.Fatal("expected non-nil")
		}
		if res.Ntasks != 8 {
			t.Errorf("Ntasks = %d; want 8 (32 slots / 4 cpus-per-task)", res.Ntasks)
		}
		if res.CpusPerTask != 4 {
			t.Errorf("CpusPerTask = %d; want 4", res.CpusPerTask)
		}
	})

	// Full geometry from injected env vars wins over LSF native vars.
	t.Run("full geometry from injected env vars", func(t *testing.T) {
		clearJobEnvVars(t)
		t.Setenv("LSB_JOBID", "99999")
		t.Setenv("LSB_DJOB_NUMPROC", "64") // ignored when NTASKS is already set
		t.Setenv("NNODES", "4")
		t.Setenv("NTASKS", "16")
		t.Setenv("NTASKS_PER_NODE", "4")
		t.Setenv("NCPUS", "2")

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
		if res.TasksPerNode != 4 {
			t.Errorf("TasksPerNode = %d; want 4", res.TasksPerNode)
		}
		if res.CpusPerTask != 2 {
			t.Errorf("CpusPerTask = %d; want 2", res.CpusPerTask)
		}
	})

	// LSB_MAX_NUM_PROCESSORS is the fallback for LSB_DJOB_NUMPROC.
	t.Run("LSB_MAX_NUM_PROCESSORS fallback", func(t *testing.T) {
		clearJobEnvVars(t)
		t.Setenv("LSB_JOBID", "99999")
		t.Setenv("LSB_MAX_NUM_PROCESSORS", "64")

		res := sched.GetJobResources()
		if res == nil {
			t.Fatal("expected non-nil")
		}
		if res.Ntasks != 64 {
			t.Errorf("Ntasks = %d; want 64", res.Ntasks)
		}
	})

	// LSB_MAX_MEM_RUSAGE is per-slot (per-task) in KB → MemPerCpuMB.
	t.Run("memory from LSB_MAX_MEM_RUSAGE", func(t *testing.T) {
		clearJobEnvVars(t)
		t.Setenv("LSB_JOBID", "99999")
		t.Setenv("LSB_MAX_MEM_RUSAGE", "8388608") // 8 GB in KB

		res := sched.GetJobResources()
		if res == nil {
			t.Fatal("expected non-nil")
		}
		if res.MemPerCpuMB != 8192 {
			t.Errorf("MemPerCpuMB = %d; want 8192", res.MemPerCpuMB)
		}
		if res.MemPerNodeMB != 0 {
			t.Errorf("MemPerNodeMB should be 0 (per-slot memory goes to MemPerCpuMB), got %d", res.MemPerNodeMB)
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
		if res.Ntasks != 4 {
			t.Errorf("Ntasks = %d; want 4", res.Ntasks)
		}
		if res.MemPerCpuMB != 0 {
			t.Errorf("MemPerCpuMB should be 0 (not set), got %d", res.MemPerCpuMB)
		}
		if res.Gpu != nil {
			t.Errorf("Gpu should be nil, got %+v", res.Gpu)
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
		if res.Ntasks != 0 {
			t.Errorf("Ntasks should be 0 for invalid value, got %d", res.Ntasks)
		}
		if res.MemPerCpuMB != 0 {
			t.Errorf("MemPerCpuMB should be 0 for negative value, got %d", res.MemPerCpuMB)
		}
	})
}
