package server

import (
	"net/http"
	"net/http/httptest"
	"os"
	"path/filepath"
	"strings"
	"testing"

	"github.com/Justype/condatainer/internal/config"
	"github.com/Justype/condatainer/internal/container"
	"github.com/Justype/condatainer/internal/overlay"
)

// setTestImagesDir points the installed-overlay scan at dir for the test.
func setTestImagesDir(t *testing.T, dir string) {
	t.Helper()
	old := config.GlobalDataPaths.ImagesDirs
	config.GlobalDataPaths.ImagesDirs = []string{dir}
	container.InvalidateInstalledOverlaysCache()
	t.Cleanup(func() {
		config.GlobalDataPaths.ImagesDirs = old
		container.InvalidateInstalledOverlaysCache()
	})
}

func postRemove(path string) *httptest.ResponseRecorder {
	r := httptest.NewRequest(http.MethodPost, "/api/remove",
		strings.NewReader(`{"path":"`+path+`"}`))
	r.Header.Set("Content-Type", "application/json")
	w := httptest.NewRecorder()
	(&srv{}).handleRemove(w, r)
	return w
}

// TestRemoveOverlay checks the happy path: an installed .sqf overlay and its
// .env sidecar are both removed.
func TestRemoveOverlay(t *testing.T) {
	root := t.TempDir()
	ovPath := filepath.Join(root, "foo--1.0.sqf")
	envPath := ovPath + ".env"
	for _, p := range []string{ovPath, envPath} {
		if err := os.WriteFile(p, []byte("x"), 0o664); err != nil {
			t.Fatal(err)
		}
	}
	setTestImagesDir(t, root)

	w := postRemove(ovPath)
	if w.Code != http.StatusOK {
		t.Fatalf("status = %d, want 200 (body: %s)", w.Code, w.Body.String())
	}
	for _, p := range []string{ovPath, envPath} {
		if _, err := os.Stat(p); !os.IsNotExist(err) {
			t.Errorf("%s still exists after remove", p)
		}
	}
}

// TestRemoveOverlayUnknownPath checks that paths outside the installed-overlay
// map are rejected instead of deleted.
func TestRemoveOverlayUnknownPath(t *testing.T) {
	root := t.TempDir()
	setTestImagesDir(t, root)

	stray := filepath.Join(root, "sub", "stray--1.0.sqf")
	if err := os.MkdirAll(filepath.Dir(stray), 0o775); err != nil {
		t.Fatal(err)
	}
	if err := os.WriteFile(stray, []byte("x"), 0o664); err != nil {
		t.Fatal(err)
	}

	w := postRemove(stray)
	if w.Code != http.StatusBadRequest {
		t.Fatalf("status = %d, want 400 (body: %s)", w.Code, w.Body.String())
	}
	if _, err := os.Stat(stray); err != nil {
		t.Errorf("stray file was removed: %v", err)
	}
}

// TestRemoveOverlayInUse checks that an overlay with a shared lock held (as
// exec/run do during execution) is reported as in use with 409.
func TestRemoveOverlayInUse(t *testing.T) {
	root := t.TempDir()
	ovPath := filepath.Join(root, "busy--1.0.sqf")
	if err := os.WriteFile(ovPath, []byte("x"), 0o664); err != nil {
		t.Fatal(err)
	}
	setTestImagesDir(t, root)

	lock, err := overlay.AcquireLock(ovPath, false)
	if err != nil {
		t.Fatal(err)
	}
	defer lock.Close()

	w := postRemove(ovPath)
	if w.Code != http.StatusConflict {
		t.Fatalf("status = %d, want 409 (body: %s)", w.Code, w.Body.String())
	}
	if _, err := os.Stat(ovPath); err != nil {
		t.Errorf("locked overlay was removed: %v", err)
	}
}
