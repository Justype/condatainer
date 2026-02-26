package utils

import (
	"os"
	"testing"
	"time"
)

func TestParseMemoryMB(t *testing.T) {
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
			mb, err := ParseMemoryMB(tt.input)
			if err != nil {
				t.Errorf("ParseMemoryMB(%q) error: %v", tt.input, err)
				return
			}
			if mb != tt.wantMB {
				t.Errorf("ParseMemoryMB(%q) = %d MB; want %d MB", tt.input, mb, tt.wantMB)
			}
		})
	}
}

func TestParseWalltime(t *testing.T) {
	hour := time.Hour
	min := time.Minute
	sec := time.Second
	day := 24 * time.Hour

	tests := []struct {
		input   string
		want    time.Duration
		wantErr bool
	}{
		// Compound: Go-style with optional integer days
		{"4d12h", 4*day + 12*hour, false},
		{"2h30m", 2*hour + 30*min, false},
		{"3h", 3 * hour, false},
		{"3H", 3 * hour, false}, // case-insensitive
		{"90m", 90 * min, false},
		{"1.5h", 90 * min, false},
		{"1d2h30m45s", day + 2*hour + 30*min + 45*sec, false},
		{"4d", 4 * day, false},
		{"", 0, false},
		// Colon-separated
		{"01:30:00", hour + 30*min, false},
		{"1:30", hour + 30*min, false},
		{"02:30:00", 2*hour + 30*min, false},
		{"90", 90 * min, false}, // minutes only
		// D-HH:MM:SS
		{"1-12:00:00", day + 12*hour, false},
		{"2-06:00:00", 2*day + 6*hour, false},
		// Errors
		{"abc", 0, true},      // no valid unit letters
		{"1.5d", 0, true},     // fractional days not supported (integer only)
		{"bad:time", 0, true}, // letters in colon-separated path
	}
	for _, tt := range tests {
		t.Run(tt.input, func(t *testing.T) {
			dur, err := ParseWalltime(tt.input)
			if tt.wantErr {
				if err == nil {
					t.Errorf("ParseWalltime(%q): expected error, got %v", tt.input, dur)
				}
				return
			}
			if err != nil {
				t.Fatalf("ParseWalltime(%q) unexpected error: %v", tt.input, err)
			}
			if dur != tt.want {
				t.Errorf("ParseWalltime(%q) = %v; want %v", tt.input, dur, tt.want)
			}
		})
	}
}

func TestStripInlineComment(t *testing.T) {
	tests := []struct {
		name  string
		input string
		want  string
	}{
		{"No comment", "--cpus-per-task=8", "--cpus-per-task=8"},
		{"With comment", "--cpus-per-task=8  # This is a comment", "--cpus-per-task=8"},
		{"Comment only", "# Just a comment", ""},
		{"Multiple hashes", "--mem=16G # First # Second", "--mem=16G"},
		{"Hash in value needs escaping", "foo=bar#baz", "foo=bar"},
		{"Whitespace around comment", "--time=02:00:00   #   Time limit  ", "--time=02:00:00"},
		{"Empty after hash", "value #", "value"},
	}

	for _, tt := range tests {
		t.Run(tt.name, func(t *testing.T) {
			got := StripInlineComment(tt.input)
			if got != tt.want {
				t.Errorf("StripInlineComment(%q) = %q, want %q", tt.input, got, tt.want)
			}
		})
	}
}

func TestGetInteractivePromptsFromScript(t *testing.T) {
	tmp, err := os.CreateTemp("", "script-*.sh")
	if err != nil {
		t.Fatalf("failed to create temp file: %v", err)
	}
	defer os.Remove(tmp.Name())

	content := `#!/bin/bash
#DEP:foo/1.0
#INTERACTIVE:Please enter the license key
#INTERACTIVE:Second prompt text
echo hello
`
	if _, err := tmp.WriteString(content); err != nil {
		t.Fatalf("failed to write to temp file: %v", err)
	}
	tmp.Close()

	prompts, err := GetInteractivePromptsFromScript(tmp.Name())
	if err != nil {
		t.Fatalf("unexpected error: %v", err)
	}
	if len(prompts) != 2 {
		t.Fatalf("expected 2 prompts, got %d", len(prompts))
	}
	if prompts[0] != "Please enter the license key" {
		t.Fatalf("unexpected first prompt: %s", prompts[0])
	}
	if prompts[1] != "Second prompt text" {
		t.Fatalf("unexpected second prompt: %s", prompts[1])
	}
}

func TestGetInteractivePromptsFromScript_Missing(t *testing.T) {
	_, err := GetInteractivePromptsFromScript("/nonexistent/path/to/script.sh")
	if err == nil {
		t.Fatalf("expected error for missing file")
	}
}

func TestGetDependenciesFromScript_MLAndModuleLoad(t *testing.T) {
	tmp, err := os.CreateTemp("", "script-deps-*.sh")
	if err != nil {
		t.Fatalf("failed to create temp file: %v", err)
	}
	defer os.Remove(tmp.Name())

	content := `#!/bin/bash
#DEP:foo/1.0
module load alpha/1.0 beta/2.0
ml load gamma/3.0 delta/4.0
ml purge
ml foo/5.0
`
	if _, err := tmp.WriteString(content); err != nil {
		t.Fatalf("failed to write to temp file: %v", err)
	}
	tmp.Close()

	deps, err := GetDependenciesFromScript(tmp.Name(), true)
	if err != nil {
		t.Fatalf("unexpected error: %v", err)
	}

	expected := []string{"foo/1.0", "alpha/1.0", "beta/2.0", "gamma/3.0", "delta/4.0", "foo/5.0"}
	if len(deps) != len(expected) {
		t.Fatalf("expected %d deps, got %d: %v", len(expected), len(deps), deps)
	}
	for i := range expected {
		if deps[i] != expected[i] {
			t.Fatalf("dep[%d] mismatch: expected %s, got %s", i, expected[i], deps[i])
		}
	}
}

func TestGetDependenciesFromScript_WithComments(t *testing.T) {
	tmp, err := os.CreateTemp("", "script-deps-comments-*.sh")
	if err != nil {
		t.Fatalf("failed to create temp file: %v", err)
	}
	defer os.Remove(tmp.Name())

	content := `#!/bin/bash
#DEP:foo/1.0  # This is a dependency comment
#DEP:bar/2.0 # Another comment
#DEP:baz/3.0#No space before hash
module load alpha/1.0 beta/2.0
`
	if _, err := tmp.WriteString(content); err != nil {
		t.Fatalf("failed to write to temp file: %v", err)
	}
	tmp.Close()

	deps, err := GetDependenciesFromScript(tmp.Name(), true)
	if err != nil {
		t.Fatalf("unexpected error: %v", err)
	}

	expected := []string{"foo/1.0", "bar/2.0", "baz/3.0", "alpha/1.0", "beta/2.0"}
	if len(deps) != len(expected) {
		t.Fatalf("expected %d deps, got %d: %v", len(expected), len(deps), deps)
	}
	for i := range expected {
		if deps[i] != expected[i] {
			t.Fatalf("dep[%d] mismatch: expected %s, got %s", i, expected[i], deps[i])
		}
	}
}
