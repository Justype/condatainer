package server

import (
	"net/http"
	"os"
)

// handleFSDelete serves DELETE /api/fs?path=... — deletes a single file or a
// directory (recursively). Dispatched from handleFS in api_files.go. Uses
// os.Lstat so a symlink entry itself is removed rather than followed and its
// target deleted. No confirmation or undo at this layer — the dashboard's
// own delete button prompts the user before issuing the request.
func (s *srv) handleFSDelete(w http.ResponseWriter, r *http.Request) {
	path := r.URL.Query().Get("path")
	if path == "" {
		http.Error(w, "path is required", http.StatusBadRequest)
		return
	}

	info, err := os.Lstat(path)
	if err != nil {
		http.Error(w, err.Error(), http.StatusBadRequest)
		return
	}

	if info.IsDir() {
		err = os.RemoveAll(path)
	} else {
		err = os.Remove(path)
	}
	if err != nil {
		http.Error(w, err.Error(), http.StatusInternalServerError)
		return
	}

	writeJSON(w, map[string]interface{}{"deleted": path})
}
