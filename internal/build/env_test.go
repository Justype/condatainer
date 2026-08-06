package build

import (
	"os"
	"path/filepath"
	"strings"
	"testing"
)

func TestSaveEnvFile(t *testing.T) {
	tmpDir := t.TempDir()
	overlayPath := filepath.Join(tmpDir, "test.sqf")

	envDict := map[string]EnvEntry{
		"STAR_INDEX_DIR": {
			Value: "{prefix}/star",
			Note:  "STAR index directory",
		},
		"GENOME_FASTA": {
			Value: "{prefix}/fasta/genome.fa",
			Note:  "genome fasta file",
		},
		"NO_NOTE": {
			Value: "{prefix}/data",
			Note:  "",
		},
	}

	// Test saving
	err := SaveEnvFile(overlayPath, envDict, "star/2.7.11a", "")
	if err != nil {
		t.Fatalf("SaveEnvFile failed: %v", err)
	}

	// Verify file was created
	envFilePath := overlayPath + ".env"
	if _, err := os.Stat(envFilePath); os.IsNotExist(err) {
		t.Fatalf("ENV file was not created at %s", envFilePath)
	}

	// Read and verify content
	content, err := os.ReadFile(envFilePath)
	if err != nil {
		t.Fatalf("Failed to read ENV file: %v", err)
	}

	contentStr := string(content)

	// Check that {prefix} was replaced
	if contains(contentStr, "{prefix}") {
		t.Error("{prefix} placeholder was not replaced in ENV file")
	}

	// Check for expected content
	expectedLines := []string{
		"STAR_INDEX_DIR=/cnt/star/2.7.11a/star",
		"GENOME_FASTA=/cnt/star/2.7.11a/fasta/genome.fa",
		"NO_NOTE=/cnt/star/2.7.11a/data",
		"#ENVNOTE:STAR_INDEX_DIR=STAR index directory",
		"#ENVNOTE:GENOME_FASTA=genome fasta file",
	}

	for _, expected := range expectedLines {
		if !contains(contentStr, expected) {
			t.Errorf("ENV file does not contain expected line: %s", expected)
		}
	}
}

func TestSaveEnvFileDescription(t *testing.T) {
	tmpDir := t.TempDir()

	t.Run("description only", func(t *testing.T) {
		overlayPath := filepath.Join(tmpDir, "description_only.sqf")
		err := SaveEnvFile(overlayPath, map[string]EnvEntry{}, "", "SAMtools alignment toolkit")
		if err != nil {
			t.Fatalf("SaveEnvFile failed: %v", err)
		}
		content, _ := os.ReadFile(overlayPath + ".env")
		if !contains(string(content), "#DESCRIPTION:SAMtools alignment toolkit") {
			t.Errorf("missing #DESCRIPTION: line, got: %s", string(content))
		}
	})

	t.Run("description with env", func(t *testing.T) {
		overlayPath := filepath.Join(tmpDir, "description_env.sqf")
		envDict := map[string]EnvEntry{"BIN": {Value: "{prefix}/bin", Note: "binary dir"}}
		err := SaveEnvFile(overlayPath, envDict, "tool/1.0", "My Tool")
		if err != nil {
			t.Fatalf("SaveEnvFile failed: %v", err)
		}
		content, _ := os.ReadFile(overlayPath + ".env")
		contentStr := string(content)
		if !contains(contentStr, "#DESCRIPTION:My Tool") {
			t.Errorf("missing #DESCRIPTION: line")
		}
		if !contains(contentStr, "BIN=/cnt/tool/1.0/bin") {
			t.Errorf("missing BIN env line")
		}
		// #DESCRIPTION: must appear before env lines
		descriptionIdx := strings.Index(contentStr, "#DESCRIPTION:")
		binIdx := strings.Index(contentStr, "BIN=")
		if descriptionIdx > binIdx {
			t.Errorf("#DESCRIPTION: should appear before env lines")
		}
	})

	t.Run("empty both", func(t *testing.T) {
		overlayPath := filepath.Join(tmpDir, "empty.sqf")
		err := SaveEnvFile(overlayPath, map[string]EnvEntry{}, "", "")
		if err != nil {
			t.Fatalf("SaveEnvFile failed: %v", err)
		}
		if _, err := os.Stat(overlayPath + ".env"); !os.IsNotExist(err) {
			t.Error("env file should not be created when both description and envDict are empty")
		}
	})
}

func contains(s, substr string) bool {
	return len(s) >= len(substr) && (s == substr || len(s) > len(substr) &&
		(s[:len(substr)] == substr || s[len(s)-len(substr):] == substr ||
			len(s) > len(substr)+1 && containsSubstr(s, substr)))
}

func containsSubstr(s, substr string) bool {
	for i := 0; i <= len(s)-len(substr); i++ {
		if s[i:i+len(substr)] == substr {
			return true
		}
	}
	return false
}
