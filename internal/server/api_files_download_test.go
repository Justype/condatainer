package server

import (
	"bytes"
	"compress/gzip"
	"io"
	"net/http"
	"net/http/httptest"
	"net/url"
	"os"
	"path/filepath"
	"strings"
	"testing"
)

// TestDownloadGz checks that gz=1 streams a single file gzip-compressed,
// named <base>.gz, and that it decompresses back to the original content.
func TestDownloadGz(t *testing.T) {
	dir := t.TempDir()
	content := bytes.Repeat([]byte("condatainer\n"), 1000)
	path := filepath.Join(dir, "notes.txt")
	if err := os.WriteFile(path, content, 0o644); err != nil {
		t.Fatal(err)
	}

	r := httptest.NewRequest(http.MethodGet,
		"/api/fs/download?path="+url.QueryEscape(path)+"&gz=1", nil)
	w := httptest.NewRecorder()
	(&srv{}).handleFSDownload(w, r)

	if w.Code != http.StatusOK {
		t.Fatalf("status = %d, body: %s", w.Code, w.Body.String())
	}
	if cd := w.Header().Get("Content-Disposition"); !strings.Contains(cd, `"notes.txt.gz"`) {
		t.Fatalf("Content-Disposition = %q, want notes.txt.gz attachment", cd)
	}
	if len(w.Body.Bytes()) >= len(content) {
		t.Fatalf("compressed body (%d bytes) not smaller than original (%d bytes)", w.Body.Len(), len(content))
	}

	gr, err := gzip.NewReader(w.Body)
	if err != nil {
		t.Fatal(err)
	}
	got, err := io.ReadAll(gr)
	if err != nil {
		t.Fatal(err)
	}
	if gr.Name != "notes.txt" {
		t.Fatalf("gzip header name = %q, want notes.txt", gr.Name)
	}
	if !bytes.Equal(got, content) {
		t.Fatalf("decompressed content mismatch: %d bytes, want %d", len(got), len(content))
	}
}
