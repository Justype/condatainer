package build

import (
	"context"
	"os"
	"path/filepath"
	"strings"
	"testing"

	"github.com/Justype/condatainer/internal/config"
	"github.com/Justype/condatainer/internal/utils"
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
		cntDirPath:        "",
		tmpOverlayPath:    "",
		targetOverlayPath: "",
	}

	err = base.parseScriptMetadata(context.Background())
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
		cntDirPath:        "",
		tmpOverlayPath:    "",
		targetOverlayPath: "",
	}

	err = base.parseScriptMetadata(context.Background())
	if err != nil {
		t.Fatalf("unexpected error: %v", err)
	}
	if len(base.dependencies) == 0 {
		t.Fatalf("expected dependencies to be parsed")
	}
}

func TestParseScriptMetadata_NcpusFromSlurm(t *testing.T) {
	// Create a temp script with SLURM cpus-per-task directive
	tmp, err := os.CreateTemp("", "script-*.sh")
	if err != nil {
		t.Fatalf("failed to create temp file: %v", err)
	}
	defer os.Remove(tmp.Name())

	content := "#!/bin/bash\n#SBATCH --cpus-per-task=8\necho hi\n"
	if _, err := tmp.WriteString(content); err != nil {
		t.Fatalf("failed to write to temp file: %v", err)
	}
	tmp.Close()

	base := &BaseBuildObject{
		nameVersion:       "foo/bar",
		buildSource:       tmp.Name(),
		cntDirPath:        "",
		tmpOverlayPath:    "",
		targetOverlayPath: "",
	}

	err = base.parseScriptMetadata(context.Background())
	if err != nil {
		t.Fatalf("unexpected error: %v", err)
	}
	if base.effectiveNcpus() != 8 {
		t.Fatalf("expected effectiveNcpus=8 from SLURM directive, got %d", base.effectiveNcpus())
	}
}

func TestParseScriptMetadata_NcpusFromPBS(t *testing.T) {
	// Create a temp script with PBS ncpus directive
	tmp, err := os.CreateTemp("", "script-*.sh")
	if err != nil {
		t.Fatalf("failed to create temp file: %v", err)
	}
	defer os.Remove(tmp.Name())

	content := "#!/bin/bash\n#PBS -l ncpus=16\necho hi\n"
	if _, err := tmp.WriteString(content); err != nil {
		t.Fatalf("failed to write to temp file: %v", err)
	}
	tmp.Close()

	base := &BaseBuildObject{
		nameVersion:       "foo/bar",
		buildSource:       tmp.Name(),
		cntDirPath:        "",
		tmpOverlayPath:    "",
		targetOverlayPath: "",
	}

	err = base.parseScriptMetadata(context.Background())
	if err != nil {
		t.Fatalf("unexpected error: %v", err)
	}
	if base.effectiveNcpus() != 16 {
		t.Fatalf("expected effectiveNcpus=16 from PBS directive, got %d", base.effectiveNcpus())
	}
}

func TestNewBuildObject_DoesNotParseInteractiveWhenInstalled(t *testing.T) {
	imagesDir := t.TempDir()
	tmpDir := t.TempDir()

	// Create a temporary base dir and add it to extra base dirs so it is searched first
	baseDir := t.TempDir()
	defer os.RemoveAll(baseDir)
	if err := os.Setenv("CNT_EXTRA_BASE_DIRS", baseDir); err != nil {
		t.Fatalf("failed to set env: %v", err)
	}
	defer os.Unsetenv("CNT_EXTRA_BASE_DIRS")

	// Create a build-scripts entry with an INTERACTIVE prompt to trigger parsing if it were called
	scriptsDir := filepath.Join(baseDir, "build-scripts")
	if err := os.MkdirAll(filepath.Join(scriptsDir, "cellranger"), 0o775); err != nil {
		t.Fatalf("failed to create scripts dir: %v", err)
	}
	scriptPath := filepath.Join(scriptsDir, "cellranger", "8.0.1")
	scriptContent := "#!/bin/bash\n#INTERACTIVE:Please enter the license key\n"
	if err := os.WriteFile(scriptPath, []byte(scriptContent), 0o664); err != nil {
		t.Fatalf("failed to write script: %v", err)
	}

	// Create an overlay file in imagesDir to simulate already-installed overlay
	nameVersion := "cellranger/8.0.1"
	fileName := strings.ReplaceAll(utils.NormalizeNameVersion(nameVersion), "/", "--") + ".sqf"
	overlayPath := filepath.Join(imagesDir, fileName)
	if err := os.WriteFile(overlayPath, []byte{}, 0o664); err != nil {
		t.Fatalf("failed to create overlay file: %v", err)
	}

	// Re-init data paths so the test's extra base dir is picked up
	config.InitDataPaths()

	// Call NewBuildObject: should return without attempting to parse the interactive script
	bo, err := NewBuildObject(context.Background(), nameVersion, false, imagesDir, tmpDir, false)
	if err != nil {
		t.Fatalf("NewBuildObject returned error: %v", err)
	}
	if bo == nil {
		t.Fatalf("expected non-nil BuildObject")
	}
}

func TestNewBuildObject_DoesNotParseInteractiveWhenTmpOverlayExists(t *testing.T) {
	imagesDir := t.TempDir()
	tmpDir := t.TempDir()

	// Create a temporary base dir and add it to extra base dirs so it is searched first
	baseDir := t.TempDir()
	defer os.RemoveAll(baseDir)
	if err := os.Setenv("CNT_EXTRA_BASE_DIRS", baseDir); err != nil {
		t.Fatalf("failed to set env: %v", err)
	}
	defer os.Unsetenv("CNT_EXTRA_BASE_DIRS")

	// Create a build-scripts entry with an INTERACTIVE prompt to trigger parsing if it were called
	scriptsDir := filepath.Join(baseDir, "build-scripts")
	if err := os.MkdirAll(filepath.Join(scriptsDir, "cellranger"), 0o775); err != nil {
		t.Fatalf("failed to create scripts dir: %v", err)
	}
	scriptPath := filepath.Join(scriptsDir, "cellranger", "8.0.1")
	scriptContent := "#!/bin/bash\n#INTERACTIVE:Please enter the license key\n"
	if err := os.WriteFile(scriptPath, []byte(scriptContent), 0o664); err != nil {
		t.Fatalf("failed to write script: %v", err)
	}

	// Create a tmp overlay file in tmpDir to simulate a build in progress
	nameVersion := "cellranger/8.0.1"
	// tmp overlay filename uses -- for / and .img suffix (see getTmpOverlayPath)
	tmpFileName := strings.ReplaceAll(utils.NormalizeNameVersion(nameVersion), "/", "--") + ".img"
	tmpOverlayPath := filepath.Join(tmpDir, tmpFileName)
	if err := os.WriteFile(tmpOverlayPath, []byte{}, 0o664); err != nil {
		t.Fatalf("failed to create tmp overlay file: %v", err)
	}

	// Re-init data paths so the test's extra base dir is picked up
	config.InitDataPaths()

	// Call NewBuildObject: should return without attempting to parse the interactive script
	bo, err := NewBuildObject(context.Background(), nameVersion, false, imagesDir, tmpDir, false)
	if err != nil {
		t.Fatalf("NewBuildObject returned error: %v", err)
	}
	if bo == nil {
		t.Fatalf("expected non-nil BuildObject")
	}
}
