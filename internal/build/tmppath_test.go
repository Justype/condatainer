package build

import (
	"context"
	"os"
	"path/filepath"
	"testing"

	"github.com/Justype/condatainer/internal/utils"
)

// Test that CreateTmpOverlay creates parent directory when it does not exist
// even if overlay creation fails (system tools may not be present in test env).
func TestCreateTmpOverlay_CreatesParentDir(t *testing.T) {
	base := &BaseBuildObject{
		nameVersion:       "foo/bar",
		tmpOverlayPath:    filepath.Join(t.TempDir(), "nonexistent", "foo.img"),
		targetOverlayPath: "",
		cntDirPath:        "",
	}

	// Ensure parent dir does not exist initially
	parent := filepath.Dir(base.tmpOverlayPath)
	if _, err := os.Stat(parent); !os.IsNotExist(err) {
		// If it exists, remove to ensure test validity
		_ = os.RemoveAll(parent)
	}

	err := base.CreateTmpOverlay(context.Background(), false)
	// We expect an error from overlay creation in CI environments where dd/mke2fs/debugfs might be missing
	if err == nil {
		// If overlay creation succeeded unexpectedly, cleanup the file and return success
		_ = os.Remove(base.tmpOverlayPath)
		return
	}

	// Parent directory should have been created regardless
	if info, err := os.Stat(parent); err != nil {
		t.Fatalf("expected parent dir %s to exist, got error: %v", parent, err)
	} else if !info.IsDir() {
		t.Fatalf("expected %s to be a directory", parent)
	}
}

// TestCreateBuildDirs_CreatesDirs verifies that CreateBuildDirs creates both
// cnt/ and tmp/ subdirectories under the build directory.
func TestCreateBuildDirs_CreatesDirs(t *testing.T) {
	tmpDir := t.TempDir()
	base := &BaseBuildObject{
		nameVersion: "foo/bar",
		cntDirPath:  filepath.Join(tmpDir, "build_foo_bar", "cnt"),
	}

	if err := base.CreateBuildDirs(context.Background(), false); err != nil {
		t.Fatalf("CreateBuildDirs failed: %v", err)
	}

	buildDir := filepath.Dir(base.cntDirPath)
	for _, sub := range []string{"cnt", "tmp"} {
		info, err := os.Stat(filepath.Join(buildDir, sub))
		if err != nil {
			t.Fatalf("expected %s dir to exist under buildDir, got error: %v", sub, err)
		}
		if !info.IsDir() {
			t.Fatalf("expected %s to be a directory", sub)
		}
	}
}

// TestCreateBuildDirs_StaleDir detects existing buildDir as stale.
func TestCreateBuildDirs_StaleDir(t *testing.T) {
	tmpDir := t.TempDir()
	cntDir := filepath.Join(tmpDir, "build_foo_bar", "cnt")
	base := &BaseBuildObject{
		nameVersion: "foo/bar",
		cntDirPath:  cntDir,
	}

	// Pre-create the build dir to simulate a stale build
	_ = os.MkdirAll(filepath.Dir(cntDir), 0o775)

	err := base.CreateBuildDirs(context.Background(), false)
	if err == nil {
		t.Fatal("expected ErrTmpOverlayExists when buildDir already exists, got nil")
	}
}

// TestCreateBuildDirs_ForceRemovesStale verifies force=true cleans stale build dir.
func TestCreateBuildDirs_ForceRemovesStale(t *testing.T) {
	tmpDir := t.TempDir()
	cntDir := filepath.Join(tmpDir, "build_foo_bar", "cnt")
	base := &BaseBuildObject{
		nameVersion: "foo/bar",
		cntDirPath:  cntDir,
	}

	// Pre-create the build dir with a marker file to simulate stale state
	_ = os.MkdirAll(filepath.Dir(cntDir), 0o775)
	markerFile := filepath.Join(filepath.Dir(cntDir), "stale_marker")
	_ = os.WriteFile(markerFile, []byte{}, 0o664)

	if err := base.CreateBuildDirs(context.Background(), true); err != nil {
		t.Fatalf("CreateBuildDirs with force=true failed: %v", err)
	}

	// Marker should be gone (dir was removed and recreated)
	if _, err := os.Stat(markerFile); !os.IsNotExist(err) {
		t.Fatal("expected stale marker to be removed after force cleanup")
	}
}

func TestResolveTmpDirForExternal_DataUsesTargetDir(t *testing.T) {
	targetDir := t.TempDir()
	got := resolveTmpDirForExternal(targetDir, "data")
	if got != targetDir {
		t.Fatalf("resolveTmpDirForExternal(data) = %q, want %q", got, targetDir)
	}
}

func TestResolveTmpDirForExternal_AppUsesScratch(t *testing.T) {
	targetDir := t.TempDir()
	got := resolveTmpDirForExternal(targetDir, "app")
	want := utils.GetTmpDir()
	if got != want {
		t.Fatalf("resolveTmpDirForExternal(app) = %q, want %q", got, want)
	}
}

func TestResolveTmpDirForExternal_EmptyDefaultsToScratch(t *testing.T) {
	targetDir := t.TempDir()
	got := resolveTmpDirForExternal(targetDir, "")
	want := utils.GetTmpDir()
	if got != want {
		t.Fatalf("resolveTmpDirForExternal(empty) = %q, want %q", got, want)
	}
}

func TestResolveTmpDirForExternal_CNTTMPDIROverridesData(t *testing.T) {
	override := t.TempDir()
	if err := os.Setenv("CNT_TMPDIR", override); err != nil {
		t.Fatalf("failed to set CNT_TMPDIR: %v", err)
	}
	defer os.Unsetenv("CNT_TMPDIR")

	targetDir := t.TempDir()
	got := resolveTmpDirForExternal(targetDir, "data")
	want := utils.GetTmpDir()
	if got != want {
		t.Fatalf("resolveTmpDirForExternal(data) with CNT_TMPDIR = %q, want %q", got, want)
	}
}
