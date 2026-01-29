package utils

import (
	"os"
	"testing"
)

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

	deps, err := GetDependenciesFromScript(tmp.Name())
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
