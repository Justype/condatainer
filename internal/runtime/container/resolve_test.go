package container

import (
	"path/filepath"
	"strings"
	"testing"
)

// A missing .sqf gets the plain message: it never pairs with anything.
func TestMissingOverlayErrorSqf(t *testing.T) {
	err := missingOverlayError("/proj/tool.sqf")
	if err == nil || !strings.Contains(err.Error(), "/proj/tool.sqf") || strings.Contains(err.Error(), "snapshot") {
		t.Fatalf("err = %v, want a plain not-found message", err)
	}
}

// A missing .img with no paired snapshot also gets the plain message —
// nothing to point the reader at.
func TestMissingOverlayErrorImgNoSnapshot(t *testing.T) {
	dir := t.TempDir()
	path := filepath.Join(dir, "env.img")
	err := missingOverlayError(path)
	if err == nil || strings.Contains(err.Error(), "snapshot") {
		t.Fatalf("err = %v, want a plain not-found message", err)
	}
}

// A missing .img that would have paired with a real snapshot names it and the
// two ways to proceed, rather than silently substituting the read-only .sqf.
func TestMissingOverlayErrorImgWithSnapshot(t *testing.T) {
	requireSquashfsTools(t)
	t.Setenv("USER", "alice")
	dir := t.TempDir()
	snapshot := filepath.Join(dir, "env-alice.sqf")
	packRuntimeSqf(t, snapshot, envRuntime())
	imgPath := filepath.Join(dir, "env-alice.img")

	err := missingOverlayError(imgPath)
	if err == nil {
		t.Fatal("expected an error")
	}
	for _, want := range []string{imgPath, snapshot, "overlay create", "exec -o"} {
		if !strings.Contains(err.Error(), want) {
			t.Errorf("error %q does not mention %q", err, want)
		}
	}
}
