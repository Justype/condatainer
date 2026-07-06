package server

import (
	"net/http"
	"os"
	"path/filepath"
	"time"
)

// maxFSEntries is the maximum number of directory entries returned per request.
// Prevents loading thousands of entries from large shared dirs like /home or /scratch.
const maxFSEntries = 500

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
		http.Error(w, err.Error(), http.StatusBadRequest)
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
	var out []fsEntry
	for _, e := range entries {
		full := filepath.Join(path, e.Name())
		fe := fsEntry{
			Name:  e.Name(),
			IsDir: e.IsDir(),
			Path:  full,
		}
		// For symlinks, follow the link to determine if the target is a directory.
		if e.Type()&os.ModeSymlink != 0 {
			if info, err := os.Stat(full); err == nil {
				fe.IsDir = info.IsDir()
				fe.Size = info.Size()
				fe.ModifiedAt = info.ModTime()
			}
		} else if info, err := e.Info(); err == nil {
			fe.Size = info.Size()
			fe.ModifiedAt = info.ModTime()
		}
		out = append(out, fe)
	}
	writeJSON(w, out)
}
