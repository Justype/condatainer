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

// TestMain confines every test's build workspace to a scratch dir removed at
// the end of the run. A constructor test that builds no farther than the
// manifest still resolves its workspace through tmpRootForType/tmpRootForDef,
// which fall back to the real host tmp when nothing overrides them — without
// this, every such test leaves an orphaned build_<name> dir in /tmp/cnt-$USER.
func TestMain(m *testing.M) {
	dir, err := os.MkdirTemp("", "condatainer-build-test-")
	if err != nil {
		panic(err)
	}
	os.Setenv("CNT_TMPDIR", dir)
	code := m.Run()
	os.RemoveAll(dir)
	os.Exit(code)
}

// TestCreateBuildDirs_CreatesDirs verifies that CreateBuildDirs creates both
// cnt/ and tmp/ subdirectories under the build directory.
func TestCreateBuildDirs_CreatesDirs(t *testing.T) {
	tmpDir := t.TempDir()
	base := &BuildObject{
		spec: Spec{Image: ImageSpec{Name: "foo/bar"}},
		ws:   workspaceFor("foo/bar", tmpDir, false),
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
		ws:   workspaceFor("foo/bar", tmpDir, false),
	}

	// Pre-create the build dir to simulate a stale build
	_ = os.MkdirAll(base.ws.BuildDir, 0o775)

	err := base.CreateBuildDirs(context.Background(), false)
	if err == nil {
		t.Fatal("expected ErrBuildDirExists when buildDir already exists, got nil")
	}
}

// TestCreateBuildDirs_ForceRemovesStale verifies force=true cleans stale build dir.
func TestCreateBuildDirs_ForceRemovesStale(t *testing.T) {
	tmpDir := t.TempDir()
	base := &BuildObject{
		spec: Spec{Image: ImageSpec{Name: "foo/bar"}},
		ws:   workspaceFor("foo/bar", tmpDir, false),
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

// An external build routes by type: data keeps its large intermediates beside the
// target the user chose, while an app and a definition take fast local scratch —
// a definition because --fakeroot has to work where its sandbox is written.
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
		{name: "definition", typ: catalog.TypeOS, isDef: true, want: scratch},
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
