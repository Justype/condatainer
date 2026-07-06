package server

import (
	"encoding/json"
	"net/http"
)

// handleHelperBookmarks serves GET/POST/DELETE /api/helpers/bookmarks — the
// dashboard's server-persisted bookmarked-helper list (see
// helper_bookmarks.go). GET lists all bookmarked helper names. POST {name}
// stars a helper. DELETE ?name=... unstars a helper. All three return the
// full updated list as JSON, so the client can simply replace its local
// bookmark set and re-render.
func (s *srv) handleHelperBookmarks(w http.ResponseWriter, r *http.Request) {
	switch r.Method {
	case http.MethodGet:
		list, err := listHelperBookmarks()
		if err != nil {
			http.Error(w, err.Error(), http.StatusInternalServerError)
			return
		}
		writeJSON(w, list)

	case http.MethodPost:
		var req struct {
			Name string `json:"name"`
		}
		if err := json.NewDecoder(r.Body).Decode(&req); err != nil {
			http.Error(w, "invalid request body", http.StatusBadRequest)
			return
		}
		if req.Name == "" {
			http.Error(w, "name is required", http.StatusBadRequest)
			return
		}
		list, err := addHelperBookmark(req.Name)
		if err != nil {
			http.Error(w, err.Error(), http.StatusBadRequest)
			return
		}
		writeJSON(w, list)

	case http.MethodDelete:
		name := r.URL.Query().Get("name")
		if name == "" {
			http.Error(w, "name is required", http.StatusBadRequest)
			return
		}
		list, err := removeHelperBookmark(name)
		if err != nil {
			http.Error(w, err.Error(), http.StatusBadRequest)
			return
		}
		writeJSON(w, list)

	default:
		http.Error(w, "method not allowed", http.StatusMethodNotAllowed)
	}
}
