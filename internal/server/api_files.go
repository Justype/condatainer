package server

import (
	"errors"
	"net/http"
	"os"
	"path/filepath"
	"sync"
	"syscall"
	"time"
)

// fsErrorResponse writes err for path as a human-readable message with a
// fitting status (403 permission, 404 missing, 400 otherwise).
func fsErrorResponse(w http.ResponseWriter, path string, err error) {
	switch {
	case os.IsPermission(err):
		http.Error(w, "no permission to access "+path, http.StatusForbidden)
	case os.IsNotExist(err):
		http.Error(w, path+" does not exist", http.StatusNotFound)
	case errors.Is(err, syscall.ENOTDIR):
		http.Error(w, path+" is not a directory", http.StatusBadRequest)
	default:
		http.Error(w, err.Error(), http.StatusBadRequest)
	}
}

// maxFSEntries is the maximum number of directory entries returned per request.
// Prevents loading thousands of entries from large shared dirs like /home or /scratch.
const maxFSEntries = 2000

// handleFS serves GET /api/fs?path= — directory listing for the file
// browser — DELETE /api/fs?path= — deleting a file or directory (see
// handleFSDelete in api_files_delete.go) — and PATCH /api/fs?path= —
// renaming or moving a file or directory (see handleFSRename in
// api_files_move.go). Sets X-FS-Truncated: true when the directory
// exceeds maxFSEntries entries.
func (s *srv) handleFS(w http.ResponseWriter, r *http.Request) {
	switch r.Method {
	case http.MethodDelete:
		s.handleFSDelete(w, r)
		return
	case http.MethodPatch:
		s.handleFSRename(w, r)
		return
	}

	path := r.URL.Query().Get("path")
	if path == "" {
		if home, err := os.UserHomeDir(); err == nil {
			path = home
		} else {
			path = "/"
		}
	}

	entries, err := os.ReadDir(path)
	if err != nil {
		fsErrorResponse(w, path, err)
		return
	}

	if len(entries) > maxFSEntries {
		entries = entries[:maxFSEntries]
		w.Header().Set("X-FS-Truncated", "true")
	}

	type fsEntry struct {
		Name       string    `json:"name"`
		IsDir      bool      `json:"is_dir"`
		Path       string    `json:"path"`
		Size       int64     `json:"size,omitempty"`
		ModifiedAt time.Time `json:"modified_at,omitempty"`
	}
	// Each entry needs a stat for size/mtime (and a follow for symlinks) —
	// on a network filesystem those are per-entry round trips, so they run
	// concurrently; done serially, a full listing takes many seconds.
	out := make([]fsEntry, len(entries))
	sem := make(chan struct{}, 64)
	var wg sync.WaitGroup
	for i, e := range entries {
		out[i] = fsEntry{
			Name:  e.Name(),
			IsDir: e.IsDir(),
			Path:  filepath.Join(path, e.Name()),
		}
		wg.Add(1)
		sem <- struct{}{}
		go func(fe *fsEntry, e os.DirEntry) {
			defer wg.Done()
			defer func() { <-sem }()
			// For symlinks, follow the link to determine if the target is a directory.
			if e.Type()&os.ModeSymlink != 0 {
				if info, err := os.Stat(fe.Path); err == nil {
					fe.IsDir = info.IsDir()
					fe.Size = info.Size()
					fe.ModifiedAt = info.ModTime()
				}
			} else if info, err := e.Info(); err == nil {
				fe.Size = info.Size()
				fe.ModifiedAt = info.ModTime()
			}
		}(&out[i], e)
	}
	wg.Wait()
	writeJSON(w, out)
}
