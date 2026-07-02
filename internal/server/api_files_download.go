package server

import (
	"archive/zip"
	"fmt"
	"io"
	"io/fs"
	"net/http"
	"os"
	"path/filepath"
)

const (
	// maxDownloadEntries caps the number of files a directory download may
	// contain, checked by a pre-walk before any bytes are streamed.
	maxDownloadEntries = 20000
	// maxDownloadBytes caps the total uncompressed size a directory download
	// may contain, checked by the same pre-walk.
	maxDownloadBytes = 5 << 30 // 5 GiB
)

// handleFSDownload serves GET /api/fs/download?path=...&inline=1 — serves a
// single file as-is (following symlinks, matching handleFS's listing
// behavior), or a directory as a streamed zip archive.
//
// By default a file is served with Content-Disposition: attachment, forcing
// a save-as download (the dashboard's explicit Download button always uses
// this). Passing inline=1 serves it with Content-Disposition: inline
// instead, so the browser renders it directly — used when a file row in
// the listing is clicked directly, as opposed to the Download button.
// inline is ignored for directories, which always download as a zip.
//
// Directory downloads do NOT follow symlinks (unlike the file case and
// unlike handleFS's listing): the zip walk uses Lstat-based traversal and
// skips symlink entries entirely, to avoid symlink loops or a link pointing
// outside the intended subtree silently ballooning the archive. Before
// streaming begins, the directory is pre-walked to check it doesn't exceed
// maxDownloadEntries/maxDownloadBytes; oversized directories are rejected
// with 413 rather than truncating an in-progress response.
func (s *srv) handleFSDownload(w http.ResponseWriter, r *http.Request) {
	path := r.URL.Query().Get("path")
	if path == "" {
		http.Error(w, "path is required", http.StatusBadRequest)
		return
	}

	info, err := os.Stat(path)
	if err != nil {
		http.Error(w, err.Error(), http.StatusBadRequest)
		return
	}

	if !info.IsDir() {
		inline := r.URL.Query().Get("inline") == "1"
		downloadFile(w, r, path, info, inline)
		return
	}

	downloadDir(w, path)
}

// downloadFile streams a single file, following symlinks transparently
// (path/info are already resolved by the caller's os.Stat). disposition is
// "inline" when inline is true, "attachment" otherwise; the browser decides
// whether it can actually render an inline response based on the detected
// Content-Type (set by http.ServeContent from the extension or content
// sniffing) — unrenderable types still fall back to a download prompt even
// with inline set.
func downloadFile(w http.ResponseWriter, r *http.Request, path string, info os.FileInfo, inline bool) {
	f, err := os.Open(path)
	if err != nil {
		http.Error(w, err.Error(), http.StatusBadRequest)
		return
	}
	defer f.Close()

	disposition := "attachment"
	if inline {
		disposition = "inline"
	}
	w.Header().Set("Content-Disposition", fmt.Sprintf("%s; filename=%q", disposition, filepath.Base(path)))
	http.ServeContent(w, r, filepath.Base(path), info.ModTime(), f)
}

// downloadDir zips dirPath and streams the archive with a Content-Disposition
// attachment header. Symlinks inside the tree are skipped, not followed.
func downloadDir(w http.ResponseWriter, dirPath string) {
	var entryCount int
	var totalBytes int64
	walkErr := filepath.WalkDir(dirPath, func(p string, d fs.DirEntry, err error) error {
		if err != nil {
			return err
		}
		if d.IsDir() || d.Type()&os.ModeSymlink != 0 {
			return nil
		}
		entryCount++
		if entryCount > maxDownloadEntries {
			return fmt.Errorf("directory contains more than %d entries", maxDownloadEntries)
		}
		info, err := d.Info()
		if err != nil {
			return err
		}
		totalBytes += info.Size()
		if totalBytes > maxDownloadBytes {
			return fmt.Errorf("directory exceeds %d bytes", maxDownloadBytes)
		}
		return nil
	})
	if walkErr != nil {
		http.Error(w, "directory too large to download: "+walkErr.Error(), http.StatusRequestEntityTooLarge)
		return
	}

	base := filepath.Base(filepath.Clean(dirPath))
	w.Header().Set("Content-Disposition", fmt.Sprintf("attachment; filename=%q", base+".zip"))
	w.Header().Set("Content-Type", "application/zip")

	zw := zip.NewWriter(w)
	defer zw.Close()

	filepath.WalkDir(dirPath, func(p string, d fs.DirEntry, err error) error { //nolint:errcheck
		if err != nil || d.IsDir() || d.Type()&os.ModeSymlink != 0 {
			return nil
		}
		rel, err := filepath.Rel(dirPath, p)
		if err != nil {
			return nil
		}
		zipPath := filepath.ToSlash(filepath.Join(base, rel))
		entryWriter, err := zw.Create(zipPath)
		if err != nil {
			return nil
		}
		f, err := os.Open(p)
		if err != nil {
			return nil
		}
		defer f.Close()
		io.Copy(entryWriter, f) //nolint:errcheck
		return nil
	})
}
