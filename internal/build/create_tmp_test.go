package build

import (
	"context"
	"os"
	"path/filepath"
	"testing"
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
