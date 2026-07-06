package server

import (
	"context"
	"net/http"
	"net/http/httptest"
	"net/url"
	"os"
	"path/filepath"
	"strings"
	"testing"
	"time"
)

// TestCopyWithNewName checks that new_name copies the file under a
// different name in the destination directory.
func TestCopyWithNewName(t *testing.T) {
	root := t.TempDir()
	src := filepath.Join(root, "orig.txt")
	if err := os.WriteFile(src, []byte("data"), 0o644); err != nil {
		t.Fatal(err)
	}
	dest := filepath.Join(root, "dest")
	if err := os.MkdirAll(dest, 0o755); err != nil {
		t.Fatal(err)
	}

	s := &srv{ctx: context.Background()}
	r := httptest.NewRequest(http.MethodPost,
		"/api/fs/copy?path="+url.QueryEscape(src),
		strings.NewReader(`{"dest_dir":"`+dest+`","new_name":"renamed.txt"}`))
	r.Header.Set("Content-Type", "application/json")
	w := httptest.NewRecorder()
	s.handleFSCopy(w, r)
	if w.Code != http.StatusOK {
		t.Fatalf("status = %d, body: %s", w.Code, w.Body.String())
	}

	// the copy runs in a background goroutine; poll for the result
	target := filepath.Join(dest, "renamed.txt")
	deadline := time.Now().Add(2 * time.Second)
	for {
		if got, err := os.ReadFile(target); err == nil && string(got) == "data" {
			break
		}
		if time.Now().After(deadline) {
			t.Fatalf("%s not created within deadline", target)
		}
		time.Sleep(10 * time.Millisecond)
	}
}

// TestCopyIntoItselfRejected checks that copying a directory into itself
// (or its own subtree) is rejected instead of starting a runaway copy —
// including when the destination reaches the subtree through a symlink.
func TestCopyIntoItselfRejected(t *testing.T) {
	root := t.TempDir()
	src := filepath.Join(root, "src")
	if err := os.MkdirAll(filepath.Join(src, "sub"), 0o755); err != nil {
		t.Fatal(err)
	}
	link := filepath.Join(root, "link")
	if err := os.Symlink(filepath.Join(src, "sub"), link); err != nil {
		t.Fatal(err)
	}

	for _, dest := range []string{src, filepath.Join(src, "sub"), link} {
		r := httptest.NewRequest(http.MethodPost,
			"/api/fs/copy?path="+url.QueryEscape(src),
			strings.NewReader(`{"dest_dir":`+`"`+dest+`"}`))
		r.Header.Set("Content-Type", "application/json")
		w := httptest.NewRecorder()
		(&srv{}).handleFSCopy(w, r)
		if w.Code != http.StatusBadRequest {
			t.Fatalf("dest %s: status = %d, want 400 (body: %s)", dest, w.Code, w.Body.String())
		}
	}
}
