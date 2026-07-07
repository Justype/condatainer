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
	mux.HandleFunc("/api/helpers/bookmarks", s.handleHelperBookmarks)
	mux.HandleFunc("/api/helpers/", s.handleHelpersSub) // /api/helpers/{id}/...
	mux.HandleFunc("/api/helpers", s.handleHelpersList)
	mux.HandleFunc("/api/overlay/create", s.handleOverlayCreate)
	mux.HandleFunc("/api/overlay/edit", s.handleOverlayEdit)
	mux.HandleFunc("/api/env/info", s.handleEnvInfo)
	mux.HandleFunc("/api/tasks/", s.handleTaskSub)
	mux.HandleFunc("/api/overlays", s.handleOverlaysList)
	mux.HandleFunc("/api/avail/inspect", s.handleAvailInspect)
	mux.HandleFunc("/api/avail", s.handleAvail)
	mux.HandleFunc("/api/search", s.handleSearch)
	mux.HandleFunc("/api/create", s.handleCreate)
	mux.HandleFunc("/api/remove", s.handleRemove)
	mux.HandleFunc("/api/fs", s.handleFS)
	mux.HandleFunc("/api/fs/download", s.handleFSDownload)
	mux.HandleFunc("/api/fs/upload/check", s.handleFSUploadCheck)
	mux.HandleFunc("/api/fs/upload", s.handleFSUpload)
	mux.HandleFunc("/api/fs/bookmarks", s.handleFSBookmarks)
	mux.HandleFunc("/api/fs/copy", s.handleFSCopy)
	mux.HandleFunc("/api/gpu-options", s.handleGpuOptions)

}
