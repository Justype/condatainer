package build

import (
	"os"
	"path/filepath"
	"testing"
)

func TestGetEnvDictFromBuildScript(t *testing.T) {
	// Create a temporary build script
	tmpDir := t.TempDir()
	scriptPath := filepath.Join(tmpDir, "test.sh")

	scriptContent := `#!/bin/bash
#ENV:STAR_INDEX_DIR=$app_root/star
#ENVNOTE:STAR_INDEX_DIR STAR index directory
#ENV:GENOME_FASTA=$app_root/fasta/genome.fa
#ENVNOTE:GENOME_FASTA genome fasta file
#ENV:NO_NOTE=$app_root/data
#ENV:INVALID_LINE
echo "Build script"
`

	if err := os.WriteFile(scriptPath, []byte(scriptContent), 0o644); err != nil {
		t.Fatalf("Failed to create test script: %v", err)
	}

	// Test parsing
	envDict, err := GetEnvDictFromBuildScript(scriptPath)
	if err != nil {
		t.Fatalf("GetEnvDictFromBuildScript failed: %v", err)
	}

	// Verify results
	expectedCount := 3 // STAR_INDEX_DIR, GENOME_FASTA, NO_NOTE (INVALID_LINE should be skipped)
	if len(envDict) != expectedCount {
		t.Errorf("Expected %d env entries, got %d", expectedCount, len(envDict))
	}

	// Check STAR_INDEX_DIR
	if entry, ok := envDict["STAR_INDEX_DIR"]; !ok {
		t.Error("STAR_INDEX_DIR not found in envDict")
	} else {
		if entry.Value != "$app_root/star" {
			t.Errorf("STAR_INDEX_DIR value = %q, want %q", entry.Value, "$app_root/star")
		}
		if entry.Note != "STAR_INDEX_DIR STAR index directory" {
			t.Errorf("STAR_INDEX_DIR note = %q, want %q", entry.Note, "STAR_INDEX_DIR STAR index directory")
		}
	}

	// Check GENOME_FASTA
	if entry, ok := envDict["GENOME_FASTA"]; !ok {
		t.Error("GENOME_FASTA not found in envDict")
	} else {
		if entry.Value != "$app_root/fasta/genome.fa" {
			t.Errorf("GENOME_FASTA value = %q, want %q", entry.Value, "$app_root/fasta/genome.fa")
		}
		if entry.Note != "GENOME_FASTA genome fasta file" {
			t.Errorf("GENOME_FASTA note = %q, want %q", entry.Note, "GENOME_FASTA genome fasta file")
		}
	}

	// Check NO_NOTE
	if entry, ok := envDict["NO_NOTE"]; !ok {
		t.Error("NO_NOTE not found in envDict")
	} else {
		if entry.Value != "$app_root/data" {
			t.Errorf("NO_NOTE value = %q, want %q", entry.Value, "$app_root/data")
		}
		if entry.Note != "" {
			t.Errorf("NO_NOTE note = %q, want empty string", entry.Note)
		}
	}
}

func TestSaveEnvFile(t *testing.T) {
	tmpDir := t.TempDir()
	overlayPath := filepath.Join(tmpDir, "test.sqf")

	envDict := map[string]EnvEntry{
		"STAR_INDEX_DIR": {
			Value: "$app_root/star",
			Note:  "STAR index directory",
		},
		"GENOME_FASTA": {
			Value: "$app_root/fasta/genome.fa",
			Note:  "genome fasta file",
		},
		"NO_NOTE": {
			Value: "$app_root/data",
			Note:  "",
		},
	}

	// Test saving
	err := SaveEnvFile(overlayPath, envDict, "star/2.7.11a")
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

	// Check that $app_root was replaced
	if contains(contentStr, "$app_root") {
		t.Error("$app_root placeholder was not replaced in ENV file")
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
