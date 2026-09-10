package conda

import (
	"os"
	"path/filepath"
	"strings"
	"testing"
)

func TestReactivateScriptOrdersDeactivateReverseActivateForward(t *testing.T) {
	root := t.TempDir()
	deactivateDir := filepath.Join(root, "etc", "conda", "deactivate.d")
	activateDir := filepath.Join(root, "etc", "conda", "activate.d")
	if err := os.MkdirAll(deactivateDir, 0o755); err != nil {
		t.Fatal(err)
	}
	if err := os.MkdirAll(activateDir, 0o755); err != nil {
		t.Fatal(err)
	}
	for _, name := range []string{"a_deactivate.sh", "b_deactivate.sh"} {
		if err := os.WriteFile(filepath.Join(deactivateDir, name), nil, 0o644); err != nil {
			t.Fatal(err)
		}
	}
	for _, name := range []string{"a_activate.sh", "b_activate.sh"} {
		if err := os.WriteFile(filepath.Join(activateDir, name), nil, 0o644); err != nil {
			t.Fatal(err)
		}
	}

	got := ReactivateScript(root, ShellBash)

	bIdx := strings.Index(got, "b_deactivate.sh")
	aIdx := strings.Index(got, "a_deactivate.sh")
	if bIdx == -1 || aIdx == -1 || bIdx > aIdx {
		t.Errorf("deactivate scripts must run in reverse filename order, got:\n%s", got)
	}
	aActIdx := strings.Index(got, "a_activate.sh")
	bActIdx := strings.Index(got, "b_activate.sh")
	if aActIdx == -1 || bActIdx == -1 || aActIdx > bActIdx {
		t.Errorf("activate scripts must run in forward filename order, got:\n%s", got)
	}
	if aIdx > aActIdx {
		t.Errorf("deactivate.d must run before activate.d, got:\n%s", got)
	}
	if !strings.Contains(got, "CONDA_PREFIX="+shellQuote(root)) {
		t.Errorf("activate lines must scope CONDA_PREFIX to root, got:\n%s", got)
	}
}

func TestReactivateScriptEmptyWithNoHooks(t *testing.T) {
	root := t.TempDir()
	if got := ReactivateScript(root, ShellBash); got != "" {
		t.Errorf("ReactivateScript(%q, ShellBash) = %q, want empty", root, got)
	}
}

func TestReactivateScriptFishSkipsShOnlyHooks(t *testing.T) {
	root := t.TempDir()
	activateDir := filepath.Join(root, "etc", "conda", "activate.d")
	if err := os.MkdirAll(activateDir, 0o755); err != nil {
		t.Fatal(err)
	}
	if err := os.WriteFile(filepath.Join(activateDir, "a_activate.sh"), nil, 0o644); err != nil {
		t.Fatal(err)
	}

	if got := ReactivateScript(root, ShellFish); got != "" {
		t.Errorf("ReactivateScript(%q, ShellFish) = %q, want empty: a .sh-only hook must not run under fish", root, got)
	}
}

func TestReactivateScriptFishSourcesFishHooks(t *testing.T) {
	root := t.TempDir()
	activateDir := filepath.Join(root, "etc", "conda", "activate.d")
	if err := os.MkdirAll(activateDir, 0o755); err != nil {
		t.Fatal(err)
	}
	if err := os.WriteFile(filepath.Join(activateDir, "a_activate.fish"), nil, 0o644); err != nil {
		t.Fatal(err)
	}

	got := ReactivateScript(root, ShellFish)

	if !strings.Contains(got, "source ") || !strings.Contains(got, "a_activate.fish") {
		t.Errorf("ReactivateScript(%q, ShellFish) = %q, want it to source the .fish hook", root, got)
	}
	if strings.Contains(got, "CONDA_PREFIX=") {
		t.Errorf("ReactivateScript(%q, ShellFish) = %q, want fish syntax, not a POSIX prefix assignment", root, got)
	}
}
