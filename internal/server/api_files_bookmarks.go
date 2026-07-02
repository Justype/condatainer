package server

import (
	"encoding/json"
	"net/http"
)

// handleFSBookmarks serves GET/POST/DELETE /api/fs/bookmarks — the
// dashboard's server-persisted file-browser bookmark list (see
// file_bookmarks.go). GET lists all bookmarks. POST {path, label?} stars a
// path. DELETE ?path=... unstars a path. All three return the full updated
// list as JSON, so the client can simply replace its local bookmarks array
// and re-render.
func (s *srv) handleFSBookmarks(w http.ResponseWriter, r *http.Request) {
	switch r.Method {
	case http.MethodGet:
		list, err := listFileBookmarks()
		if err != nil {
			http.Error(w, err.Error(), http.StatusInternalServerError)
			return
		}
		writeJSON(w, list)

	case http.MethodPost:
		var req struct {
			Path  string `json:"path"`
			Label string `json:"label"`
		}
		if err := json.NewDecoder(r.Body).Decode(&req); err != nil {
			http.Error(w, "invalid request body", http.StatusBadRequest)
			return
		}
		if req.Path == "" {
			http.Error(w, "path is required", http.StatusBadRequest)
			return
		}
		list, err := addFileBookmark(req.Path, req.Label)
		if err != nil {
			http.Error(w, err.Error(), http.StatusBadRequest)
			return
		}
		writeJSON(w, list)

	case http.MethodDelete:
		path := r.URL.Query().Get("path")
		if path == "" {
			http.Error(w, "path is required", http.StatusBadRequest)
			return
		}
		list, err := removeFileBookmark(path)
		if err != nil {
			http.Error(w, err.Error(), http.StatusBadRequest)
			return
		}
		writeJSON(w, list)

	default:
		http.Error(w, "method not allowed", http.StatusMethodNotAllowed)
	}
}
