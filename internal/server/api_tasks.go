package server

import (
	"fmt"
	"net/http"
	"strings"
)

// handleTaskSub routes GET /api/tasks/{id}/stream and POST /api/tasks/{id}/cancel.
func (s *srv) handleTaskSub(w http.ResponseWriter, r *http.Request) {
	sub := strings.TrimPrefix(r.URL.Path, "/api/tasks/")
	parts := strings.SplitN(sub, "/", 2)
	if len(parts) < 2 {
		http.NotFound(w, r)
		return
	}
	id, action := parts[0], parts[1]
	switch action {
	case "stream":
		s.handleTaskStream(w, r, id)
	case "cancel":
		s.handleTaskCancel(w, r, id)
	default:
		http.NotFound(w, r)
	}
}

// handleTaskStream serves GET /api/tasks/{id}/stream as an SSE stream.
// Clients receive {"t":"out","d":"..."} chunks and a final {"t":"done","ok":bool} event.
func (s *srv) handleTaskStream(w http.ResponseWriter, r *http.Request, id string) {
	v, ok := s.tasks.Load(id)
	if !ok {
		http.NotFound(w, r)
		return
	}
	v.(*taskEntry).broker.ServeSSE(w, r)
}

// handleTaskCancel serves POST /api/tasks/{id}/cancel.
func (s *srv) handleTaskCancel(w http.ResponseWriter, r *http.Request, id string) {
	if r.Method != http.MethodPost {
		http.Error(w, "POST required", http.StatusMethodNotAllowed)
		return
	}
	v, ok := s.tasks.Load(id)
	if !ok {
		http.NotFound(w, r)
		return
	}
	v.(*taskEntry).cancel()
	w.WriteHeader(http.StatusOK)
	fmt.Fprint(w, `{"ok":true}`)
}
