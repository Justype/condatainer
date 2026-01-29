package build

import (
	"os"
	"testing"
)

func TestParseScriptMetadata_RequiresTTY(t *testing.T) {
	// Create a temp script with INTERACTIVE tag
	tmp, err := os.CreateTemp("", "script-*.sh")
	if err != nil {
		t.Fatalf("failed to create temp file: %v", err)
	}
	defer os.Remove(tmp.Name())

	content := "#INTERACTIVE:Please enter something\n"
	if _, err := tmp.WriteString(content); err != nil {
		t.Fatalf("failed to write to temp file: %v", err)
	}
	tmp.Close()

	base := &BaseBuildObject{
		nameVersion:       "foo/bar",
		buildSource:       tmp.Name(),
		ncpus:             1,
		cntDirPath:        "",
		tmpOverlayPath:    "",
		targetOverlayPath: "",
	}

	err = base.parseScriptMetadata()
	if err == nil {
		t.Fatalf("expected error when interactive prompts are present but no TTY is available")
	}
}

func TestParseScriptMetadata_NoInteractive(t *testing.T) {
	// Create a temp script without INTERACTIVE tag
	tmp, err := os.CreateTemp("", "script-*.sh")
	if err != nil {
		t.Fatalf("failed to create temp file: %v", err)
	}
	defer os.Remove(tmp.Name())

	content := "#DEP:foo/1.0\necho hi\n"
	if _, err := tmp.WriteString(content); err != nil {
		t.Fatalf("failed to write to temp file: %v", err)
	}
	tmp.Close()

	base := &BaseBuildObject{
		nameVersion:       "foo/bar",
		buildSource:       tmp.Name(),
		ncpus:             1,
		cntDirPath:        "",
		tmpOverlayPath:    "",
		targetOverlayPath: "",
	}

	err = base.parseScriptMetadata()
	if err != nil {
		t.Fatalf("unexpected error: %v", err)
	}
	if len(base.dependencies) == 0 {
		t.Fatalf("expected dependencies to be parsed")
	}
}
