package server

import (
	"net/http"
	"net/http/httptest"
	"net/url"
	"os"
	"path/filepath"
	"strings"
	"testing"
)

// TestMoveIntoItselfRejected checks that moving a directory into itself
// (or its own subtree) is rejected with a clear 400 instead of a kernel
// EINVAL or the EXDEV copy fallback.
func TestMoveIntoItselfRejected(t *testing.T) {
	root := t.TempDir()
	src := filepath.Join(root, "src")
	if err := os.MkdirAll(filepath.Join(src, "sub"), 0o755); err != nil {
		t.Fatal(err)
	}

	for _, dest := range []string{src, filepath.Join(src, "sub")} {
		r := httptest.NewRequest(http.MethodPatch,
			"/api/fs?path="+url.QueryEscape(src),
			strings.NewReader(`{"dest_dir":"`+dest+`"}`))
		r.Header.Set("Content-Type", "application/json")
		w := httptest.NewRecorder()
		(&srv{}).handleFSRename(w, r)
		if w.Code != http.StatusBadRequest {
			t.Fatalf("dest %s: status = %d, want 400 (body: %s)", dest, w.Code, w.Body.String())
		}
	}
}
