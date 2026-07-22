package build

import (
	"os"
	"path/filepath"
	"strings"
	"testing"
)

func TestSynthesizeDefFromURI(t *testing.T) {
	tmpDir := t.TempDir()

	t.Run("docker URI maps scheme and image", func(t *testing.T) {
		path, err := synthesizeDefFromURI("docker://ubuntu:22.04", tmpDir)
		if err != nil {
			t.Fatalf("unexpected error: %v", err)
		}
		data, err := os.ReadFile(path)
		if err != nil {
			t.Fatalf("read def: %v", err)
		}
		content := string(data)
		if !strings.Contains(content, "Bootstrap: docker") {
			t.Errorf("missing Bootstrap: docker, got:\n%s", content)
		}
		if !strings.Contains(content, "From: ubuntu:22.04") {
			t.Errorf("missing From: ubuntu:22.04, got:\n%s", content)
		}
		if !strings.Contains(content, "docker://ubuntu:22.04") {
			t.Errorf("header should record the source URI, got:\n%s", content)
		}
	})

	t.Run("rejects non-URI source", func(t *testing.T) {
		if _, err := synthesizeDefFromURI("/path/to/local.def", tmpDir); err == nil {
			t.Error("expected error for non-URI source, got nil")
		}
	})
}

func TestWriteRecordingDef(t *testing.T) {
	tmpDir := t.TempDir()
	cleanPath := filepath.Join(tmpDir, "clean.def")
	// No trailing newline, to exercise the newline-normalizing branch.
	if err := os.WriteFile(cleanPath, []byte("Bootstrap: docker\nFrom: alpine:3.19"), 0o644); err != nil {
		t.Fatalf("write clean def: %v", err)
	}

	recPath, err := writeRecordingDef(cleanPath, tmpDir)
	if err != nil {
		t.Fatalf("unexpected error: %v", err)
	}
	data, err := os.ReadFile(recPath)
	if err != nil {
		t.Fatalf("read recording def: %v", err)
	}
	content := string(data)

	// Original directives preserved.
	if !strings.Contains(content, "Bootstrap: docker") || !strings.Contains(content, "From: alpine:3.19") {
		t.Errorf("recording def dropped original directives:\n%s", content)
	}
	// Appended %files section embeds the clean def at the recorded path.
	absClean, _ := filepath.Abs(cleanPath)
	wantLine := "    " + absClean + " " + recordedDefPath
	if !strings.Contains(content, "%files") || !strings.Contains(content, wantLine) {
		t.Errorf("recording def missing %%files embed line %q:\n%s", wantLine, content)
	}
}
