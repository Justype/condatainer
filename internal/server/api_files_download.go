package server

import (
	"archive/zip"
	"compress/gzip"
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

// handleFSDownload serves GET /api/fs/download?path=... — a single file
// as-is (following symlinks), or a directory as a streamed zip archive
// (symlinks inside skipped, not followed; oversized directories rejected
// with 413). File-only query params: inline=1 serves with
// Content-Disposition: inline instead of attachment; gz=1 streams the file
// gzip-compressed as <name>.gz and takes precedence over inline. check=1
// only validates the download (200 or 413), streaming nothing.
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
	// FIFOs/sockets/devices are not downloadable — opening one can block forever.
	if !info.IsDir() && !info.Mode().IsRegular() {
		http.Error(w, "not a regular file", http.StatusBadRequest)
		return
	}

	// the dashboard pre-checks with check=1 because its <a download>
	// trigger gets no error feedback from a failed real download
	if r.URL.Query().Get("check") == "1" {
		if info.IsDir() {
			if err := checkDirDownload(path); err != nil {
				http.Error(w, err.Error(), http.StatusRequestEntityTooLarge)
				return
			}
		}
		writeJSON(w, map[string]bool{"ok": true})
		return
	}

	if !info.IsDir() {
		if r.URL.Query().Get("gz") == "1" {
			downloadFileGz(w, path, info)
			return
		}
		inline := r.URL.Query().Get("inline") == "1"
		downloadFile(w, r, path, info, inline)
		return
	}

	downloadDir(w, path)
}

// downloadFile streams a single file, following symlinks. Served inline
// (the browser renders it if it can) when inline is true, as an attachment
// otherwise.
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

// downloadFileGz streams a single file gzip-compressed, as an attachment
// named <base>.gz. The gzip header carries the original name and mtime, so
// gunzip restores them.
func downloadFileGz(w http.ResponseWriter, path string, info os.FileInfo) {
	f, err := os.Open(path)
	if err != nil {
		http.Error(w, err.Error(), http.StatusBadRequest)
		return
	}
	defer f.Close()

	base := filepath.Base(path)
	w.Header().Set("Content-Disposition", fmt.Sprintf("attachment; filename=%q", base+".gz"))
	w.Header().Set("Content-Type", "application/gzip")

	gz := gzip.NewWriter(w)
	gz.Name = base
	gz.ModTime = info.ModTime()
	defer gz.Close()
	io.Copy(gz, f) //nolint:errcheck
}

// checkDirDownload pre-walks dirPath and returns an error if it exceeds
// maxDownloadEntries or maxDownloadBytes (symlinks and other non-regular
// files skipped, matching the zip walk).
func checkDirDownload(dirPath string) error {
	var entryCount int
	var totalBytes int64
	return filepath.WalkDir(dirPath, func(p string, d fs.DirEntry, err error) error {
		if err != nil {
			return err
		}
		if d.IsDir() || !d.Type().IsRegular() {
			return nil
		}
		entryCount++
		if entryCount > maxDownloadEntries {
			return fmt.Errorf("folder is too large to download as zip: more than %d files", maxDownloadEntries)
		}
		info, err := d.Info()
		if err != nil {
			return err
		}
		totalBytes += info.Size()
		if totalBytes > maxDownloadBytes {
			return fmt.Errorf("folder is too large to download as zip: over %d GiB in total", maxDownloadBytes>>30)
		}
		return nil
	})
}

// downloadDir zips dirPath and streams the archive with a Content-Disposition
// attachment header. Symlinks and other non-regular files inside the tree
// are skipped (symlinks are not followed).
func downloadDir(w http.ResponseWriter, dirPath string) {
	if err := checkDirDownload(dirPath); err != nil {
		http.Error(w, err.Error(), http.StatusRequestEntityTooLarge)
		return
	}

	base := filepath.Base(filepath.Clean(dirPath))
	w.Header().Set("Content-Disposition", fmt.Sprintf("attachment; filename=%q", base+".zip"))
	w.Header().Set("Content-Type", "application/zip")

	zw := zip.NewWriter(w)
	defer zw.Close()

	filepath.WalkDir(dirPath, func(p string, d fs.DirEntry, err error) error { //nolint:errcheck
		if err != nil || d.IsDir() || !d.Type().IsRegular() {
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
