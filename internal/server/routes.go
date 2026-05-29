package server

import (
	"net/http"

	"github.com/Justype/condatainer/web"
)

func registerRoutes(mux *http.ServeMux, s *srv) {
	// Static assets
	mux.Handle("/static/", http.StripPrefix("/static/", http.FileServer(http.FS(web.Files))))

	// Dashboard root
	mux.HandleFunc("/", func(w http.ResponseWriter, r *http.Request) {
		if r.URL.Path != "/" {
			http.NotFound(w, r)
			return
		}
		data, err := web.Files.ReadFile("index.html")
		if err != nil {
			http.Error(w, "index.html not found", http.StatusInternalServerError)
			return
		}
		w.Header().Set("Content-Type", "text/html; charset=utf-8")
		w.Write(data) //nolint:errcheck
	})

	// REST API
	mux.HandleFunc("/api/status", s.handleStatus)
	mux.HandleFunc("/api/env/check", s.handleEnvCheck)
	mux.HandleFunc("/api/helpers/available", s.handleHelpersAvailable)
	mux.HandleFunc("/api/helpers/update", s.handleHelpersUpdate)
	mux.HandleFunc("/api/helpers/", s.handleHelpersSub) // /api/helpers/{id}/...
	mux.HandleFunc("/api/helpers", s.handleHelpersList)
	mux.HandleFunc("/api/overlay/create", s.handleOverlayCreate)
	mux.HandleFunc("/api/overlay/edit", s.handleOverlayEdit)
	mux.HandleFunc("/api/env/info", s.handleEnvInfo)
	mux.HandleFunc("/api/tasks/", s.handleTaskSub)
	mux.HandleFunc("/api/overlays", s.handleOverlaysList)
	mux.HandleFunc("/api/fs", s.handleFS)
	mux.HandleFunc("/api/gpu-options", s.handleGpuOptions)

}
