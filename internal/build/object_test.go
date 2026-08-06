package build

import (
	"context"
	"encoding/json"
	"os"
	"path/filepath"
	"strings"
	"testing"
	"time"

	"github.com/Justype/condatainer/catalog"
	"github.com/Justype/condatainer/internal/config"
)

// writeRecipe creates a one-recipe collection and returns its root.
func writeRecipe(t *testing.T, nameVersion, content string) string {
	t.Helper()
	root := t.TempDir()
	path := filepath.Join(root, "recipes", filepath.FromSlash(nameVersion))
	if err := os.MkdirAll(filepath.Dir(path), 0o775); err != nil {
		t.Fatalf("failed to create recipes dir: %v", err)
	}
	if err := os.WriteFile(path, []byte(content), 0o664); err != nil {
		t.Fatalf("failed to write recipe: %v", err)
	}
	return root
}

// setTestSource points the catalog at a single filesystem collection.
func setTestSource(t *testing.T, root string) {
	t.Helper()
	oldSources := config.Global.Sources
	config.Global.Sources = []catalog.Spec{{Name: "test", Base: root}}
	config.ResetCatalog()
	t.Cleanup(func() {
		config.Global.Sources = oldSources
		config.ResetCatalog()
	})
}

func TestParseScriptMetadata_RequiresTTY(t *testing.T) {
	// Create a temp script with an #INPUT: declaration
	tmp, err := os.CreateTemp("", "script-*.sh")
	if err != nil {
		t.Fatalf("failed to create temp file: %v", err)
	}
	defer os.Remove(tmp.Name())

	content := "#INPUT:Please enter something\n"
	if _, err := tmp.WriteString(content); err != nil {
		t.Fatalf("failed to write to temp file: %v", err)
	}
	tmp.Close()

	base := &BuildObject{
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
	// Create a temp script with no #INPUT: declaration
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

	base := &BuildObject{
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

	base := &BuildObject{
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

	content := "#!/bin/bash\n#PBS -l select=1:ncpus=16\necho hi\n"
	if _, err := tmp.WriteString(content); err != nil {
		t.Fatalf("failed to write to temp file: %v", err)
	}
	tmp.Close()

	base := &BuildObject{
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

	// Point the catalog at a collection providing a recipe with an #INPUT:
	// prompt, so parsing it would block on a prompt if it were attempted.
	setTestSource(t, writeRecipe(t, "cellranger/8.0.1",
		"#!/bin/bash\n#INPUT:Please enter the license key\n"))

	// Create an overlay file in imagesDir to simulate already-installed overlay
	nameVersion := "cellranger/8.0.1"
	fileName := strings.ReplaceAll(catalog.Normalize(nameVersion), "/", "--") + ".sqf"
	overlayPath := filepath.Join(imagesDir, fileName)
	if err := os.WriteFile(overlayPath, []byte{}, 0o664); err != nil {
		t.Fatalf("failed to create overlay file: %v", err)
	}

	// Re-init data paths so the test's extra base dir is picked up
	config.InitDataPaths()

	// Call NewBuildObject: should return without attempting to parse the recipe inputs
	bo, err := NewBuildObject(context.Background(), nameVersion, false, imagesDir, tmpDir, false)
	if err != nil {
		t.Fatalf("NewBuildObject returned error: %v", err)
	}
	if bo == nil {
		t.Fatalf("expected non-nil BuildObject")
	}
}

func TestNewBuildObject_ErrorsWhenBuildLockExists(t *testing.T) {
	imagesDir := t.TempDir()
	tmpDir := t.TempDir()

	// Point the catalog at a collection providing a recipe with an #INPUT:
	// prompt, so parsing it would block on a prompt if it were attempted.
	setTestSource(t, writeRecipe(t, "cellranger/8.0.1",
		"#!/bin/bash\n#INPUT:Please enter the license key\n"))

	// Simulate a build in progress: create a live JSON lock file next to the target image.
	// Use the current process PID and hostname so isBuildLockStale() treats it as active.
	// NewBuildObject overrides tmpDir via resolveTmpDirForConda(), so tmp artifacts
	// in the caller-supplied tmpDir are invisible to it. The build-in-progress guard
	// checks base.buildLockPath() = targetOverlayPath + ".lock" (lives in imagesDir).
	nameVersion := "cellranger/8.0.1"
	sqfName := strings.ReplaceAll(catalog.Normalize(nameVersion), "/", "--") + ".sqf"
	lockPath := filepath.Join(imagesDir, sqfName+".lock")
	liveLock := BuildLockInfo{
		Type:      "local",
		Node:      shortHostname(),
		PID:       os.Getpid(), // this process is definitely alive
		CreatedAt: time.Now().Format(time.RFC3339),
	}
	lockData, err := json.Marshal(liveLock)
	if err != nil {
		t.Fatalf("failed to marshal lock: %v", err)
	}
	if err := os.WriteFile(lockPath, lockData, 0o664); err != nil {
		t.Fatalf("failed to create lock file: %v", err)
	}
	defer os.Remove(lockPath)

	// Re-init data paths so the test's extra base dir is picked up
	config.InitDataPaths()

	// Call NewBuildObject: should return an error with lock details.
	_, err = NewBuildObject(context.Background(), nameVersion, false, imagesDir, tmpDir, false)
	if err == nil {
		t.Fatal("expected error when build lock file exists, got nil")
	}
	if !strings.Contains(err.Error(), "build lock found") {
		t.Fatalf("unexpected error message: %v", err)
	}
}
