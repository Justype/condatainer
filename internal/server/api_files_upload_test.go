package server

import (
	"bytes"
	"mime/multipart"
	"net/http"
	"net/http/httptest"
	"net/url"
	"os"
	"path/filepath"
	"testing"
)

// TestUploadPartialFileCleanedUp aborts an upload mid-file (truncated
// multipart body, as an XHR abort produces) and checks the half-written
// target is removed instead of left truncated on disk.
func TestUploadPartialFileCleanedUp(t *testing.T) {
	dest := t.TempDir()

	var buf bytes.Buffer
	mw := multipart.NewWriter(&buf)
	fw, err := mw.CreateFormFile("sub/big.bin", "big.bin")
	if err != nil {
		t.Fatal(err)
	}
	if _, err := fw.Write(make([]byte, 4<<20)); err != nil {
		t.Fatal(err)
	}
	if err := mw.Close(); err != nil {
		t.Fatal(err)
	}
	full := buf.Bytes()

	r := httptest.NewRequest(http.MethodPost,
		"/api/fs/upload?dest="+url.QueryEscape(dest)+"&overwrite=true",
		bytes.NewReader(full[:len(full)/2])) // cut off mid-file
	r.Header.Set("Content-Type", mw.FormDataContentType())
	w := httptest.NewRecorder()

	(&srv{}).handleFSUpload(w, r)

	if w.Code == http.StatusOK {
		t.Fatalf("expected an error status for a truncated body, got %d", w.Code)
	}
	if _, err := os.Lstat(filepath.Join(dest, "sub", "big.bin")); !os.IsNotExist(err) {
		t.Fatalf("partial file left behind (stat err: %v)", err)
	}
}
