package server

import (
	"net/http"
	"os"
	"path/filepath"
	"sort"
	"strings"

	"github.com/Justype/condatainer/internal/container"
	"github.com/Justype/condatainer/internal/utils"
)

// handleOverlaysList serves GET /api/overlays — lists module (.sqf) overlays.
func (s *srv) handleOverlaysList(w http.ResponseWriter, r *http.Request) {
	installed, err := container.InstalledOverlays()
	if err != nil {
		http.Error(w, err.Error(), http.StatusInternalServerError)
		return
	}
	type overlayEntry struct {
		Name string `json:"name"`
		Path string `json:"path"`
		Size int64  `json:"size"`
		Type string `json:"type"`
	}
	var entries []overlayEntry
	for name, path := range installed {
		if !utils.IsSqf(path) {
			continue
		}
		size := int64(0)
		if info, err := os.Stat(path); err == nil {
			size = info.Size()
		}
		ext := strings.ToLower(filepath.Ext(path))
		entries = append(entries, overlayEntry{
			Name: name,
			Path: path,
			Size: size,
			Type: ext[1:],
		})
	}
	sort.Slice(entries, func(i, j int) bool { return entries[i].Name < entries[j].Name })
	writeJSON(w, entries)
}
