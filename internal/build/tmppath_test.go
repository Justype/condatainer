package build

import (
	"context"
	"os"
	"path/filepath"
	"strings"
	"testing"

	"github.com/Justype/condatainer/catalog"
	"github.com/Justype/condatainer/internal/config"
	"github.com/Justype/condatainer/internal/utils"
)

// Test that CreateTmpOverlay creates parent directory when it does not exist
// even if overlay creation fails (system tools may not be present in test env).
func TestCreateTmpOverlay_CreatesParentDir(t *testing.T) {
	base := &BuildObject{
		spec: Spec{Image: ImageSpec{Name: "foo/bar"}},
		ws:   Workspace{Overlay: filepath.Join(t.TempDir(), "nonexistent", "foo.img")},
	}

	// Ensure parent dir does not exist initially
	parent := filepath.Dir(base.ws.Overlay)
	if _, err := os.Stat(parent); !os.IsNotExist(err) {
		// If it exists, remove to ensure test validity
		_ = os.RemoveAll(parent)
	}

	err := base.CreateTmpOverlay(context.Background(), false)
	// We expect an error from overlay creation in CI environments where dd/mke2fs/debugfs might be missing
	if err == nil {
		// If overlay creation succeeded unexpectedly, cleanup the file and return success
		_ = os.Remove(base.ws.Overlay)
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
	base := &BuildObject{
		spec: Spec{Image: ImageSpec{Name: "foo/bar"}},
		ws:   workspaceFor("foo/bar", tmpDir, ""),
	}

	if err := base.CreateBuildDirs(context.Background(), false); err != nil {
		t.Fatalf("CreateBuildDirs failed: %v", err)
	}

	buildDir := base.ws.BuildDir
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
	base := &BuildObject{
		spec: Spec{Image: ImageSpec{Name: "foo/bar"}},
		ws:   workspaceFor("foo/bar", tmpDir, ""),
	}

	// Pre-create the build dir to simulate a stale build
	_ = os.MkdirAll(base.ws.BuildDir, 0o775)

	err := base.CreateBuildDirs(context.Background(), false)
	if err == nil {
		t.Fatal("expected ErrTmpOverlayExists when buildDir already exists, got nil")
	}
}

// TestCreateBuildDirs_ForceRemovesStale verifies force=true cleans stale build dir.
func TestCreateBuildDirs_ForceRemovesStale(t *testing.T) {
	tmpDir := t.TempDir()
	base := &BuildObject{
		spec: Spec{Image: ImageSpec{Name: "foo/bar"}},
		ws:   workspaceFor("foo/bar", tmpDir, ""),
	}

	// Pre-create the build dir with a marker file to simulate stale state
	_ = os.MkdirAll(base.ws.BuildDir, 0o775)
	markerFile := filepath.Join(base.ws.BuildDir, "stale_marker")
	_ = os.WriteFile(markerFile, []byte{}, 0o664)

	if err := base.CreateBuildDirs(context.Background(), true); err != nil {
		t.Fatalf("CreateBuildDirs with force=true failed: %v", err)
	}

	// Marker should be gone (dir was removed and recreated)
	if _, err := os.Stat(markerFile); !os.IsNotExist(err) {
		t.Fatal("expected stale marker to be removed after force cleanup")
	}
}

// An external build routes by type: an app takes fast local scratch, while data
// and definitions keep their large intermediates beside the target the user chose.
func TestTmpRootForExternal(t *testing.T) {
	targetDir := t.TempDir()
	scratch := utils.GetTmpDir()

	tests := []struct {
		name  string
		typ   catalog.Type
		isDef bool
		want  string
	}{
		{name: "app", typ: catalog.TypeApp, want: scratch},
		{name: "unknown type defaults to scratch", typ: "", want: scratch},
		{name: "data", typ: catalog.TypeData, want: targetDir},
		{name: "definition", typ: catalog.TypeOS, isDef: true, want: targetDir},
	}
	for _, tt := range tests {
		t.Run(tt.name, func(t *testing.T) {
			if got := tmpRootForExternal(targetDir, tt.typ, tt.isDef); got != tt.want {
				t.Errorf("tmpRootForExternal(%q, %v) = %q, want %q", tt.typ, tt.isDef, got, tt.want)
			}
		})
	}
}

// CNT_TMPDIR selects the fast root. It must not pull a build off a root chosen
// for being large and stable — that is what made a data build land on scratch
// the job wipes.
func TestCNTTmpDirMovesOnlyTheFastRoot(t *testing.T) {
	override := t.TempDir()
	t.Setenv("CNT_TMPDIR", override)

	if got := tmpRootForType(catalog.TypeApp); got != utils.GetTmpDir() {
		t.Errorf("app root = %q, want the fast root %q", got, utils.GetTmpDir())
	}
	if !strings.HasPrefix(utils.GetTmpDir(), override) {
		t.Errorf("fast root %q does not sit under CNT_TMPDIR %q", utils.GetTmpDir(), override)
	}

	stable := config.GetWritableTmpDir()
	if strings.HasPrefix(stable, override) {
		t.Errorf("stable root %q followed CNT_TMPDIR", stable)
	}

	targetDir := t.TempDir()
	if got := tmpRootForExternal(targetDir, catalog.TypeData, false); got != targetDir {
		t.Errorf("external data root = %q, want the target dir %q", got, targetDir)
	}
}
