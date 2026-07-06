package server

import (
	"context"
	"encoding/json"
	"fmt"
	"net/http"
	"os"
	"path/filepath"
	"sort"
	"strings"
	"time"

	"log/slog"

	"github.com/Justype/condatainer/internal/container"
	cntexec "github.com/Justype/condatainer/internal/exec"
	"github.com/Justype/condatainer/internal/logging"
	"github.com/Justype/condatainer/internal/overlay"
	"github.com/Justype/condatainer/internal/utils"
	"github.com/Justype/condatainer/internal/weblog"
)

// handleOverlayEdit serves POST /api/overlay/edit — resize, remove, or add packages
// on an existing writable .img overlay. Returns a task ID immediately; progress is
// streamed via GET /api/tasks/{id}/stream.
func (s *srv) handleOverlayEdit(w http.ResponseWriter, r *http.Request) {
	if r.Method != http.MethodPost {
		http.Error(w, "POST required", http.StatusMethodNotAllowed)
		return
	}
	type editReq struct {
		Path           string `json:"path"`
		Size           string `json:"size"`            // new size e.g. "30G"; empty = no resize
		AddPackages    string `json:"add_packages"`    // space-separated
		RemovePackages string `json:"remove_packages"` // space-separated
	}
	var req editReq
	if err := json.NewDecoder(r.Body).Decode(&req); err != nil {
		http.Error(w, "bad request body", http.StatusBadRequest)
		return
	}
	if req.Path == "" {
		http.Error(w, "path required", http.StatusBadRequest)
		return
	}
	if !utils.IsImg(req.Path) {
		http.Error(w, "only writable .img overlays can be edited", http.StatusBadRequest)
		return
	}

	var newSizeMB int
	if req.Size != "" {
		var err error
		newSizeMB, err = utils.ParseSizeToMB(req.Size)
		if err != nil {
			http.Error(w, "invalid size: "+err.Error(), http.StatusBadRequest)
			return
		}
	}

	addPkgs := strings.Fields(req.AddPackages)
	removePkgs := strings.Fields(req.RemovePackages)
	if newSizeMB == 0 && len(addPkgs) == 0 && len(removePkgs) == 0 {
		http.Error(w, "nothing to do: specify size, add_packages, or remove_packages", http.StatusBadRequest)
		return
	}

	taskID := fmt.Sprintf("edit-%d", time.Now().UnixNano())
	broker := newSSEBroker()
	ctx, cancel := context.WithCancel(s.ctx)
	s.tasks.Store(taskID, &taskEntry{broker: broker, cancel: cancel})
	writeJSON(w, map[string]string{"id": taskID})

	bw := &brokerWriter{broker}
	ctx = logging.WithLogger(ctx, slog.New(weblog.New(bw)))
	ctx = logging.WithWriter(ctx, bw)
	io := cntexec.IO{Stdout: bw, Stderr: bw}

	go func() {
		defer cancel()
		defer s.scheduleTaskCleanup(taskID)

		if newSizeMB > 0 {
			fmt.Fprintf(bw, "Resizing %s to %s...\n", req.Path, req.Size)
			if err := overlay.Resize(ctx, req.Path, newSizeMB); err != nil {
				broadcastResult(broker, ctx, fmt.Errorf("resize: %w", err))
				return
			}
		}
		if len(removePkgs) > 0 {
			fmt.Fprintf(bw, "Removing packages: %s\n", strings.Join(removePkgs, " "))
			if err := cntexec.RemovePackages(ctx, req.Path, removePkgs, false, io); err != nil {
				broadcastResult(broker, ctx, fmt.Errorf("remove packages: %w", err))
				return
			}
		}
		if len(addPkgs) > 0 {
			fmt.Fprintf(bw, "Installing packages: %s\n", strings.Join(addPkgs, " "))
			if err := cntexec.InstallPackages(ctx, req.Path, addPkgs, false, io); err != nil {
				broadcastResult(broker, ctx, fmt.Errorf("add packages: %w", err))
				return
			}
		}
		fmt.Fprintf(bw, "Done.\n")
		result, _ := json.Marshal(map[string]interface{}{"t": "done", "ok": true, "path": req.Path})
		broker.publishFinal(result)
	}()
}

// handleOverlaysList serves GET /api/overlays — lists module (.sqf) overlays.
func (s *srv) handleOverlaysList(w http.ResponseWriter, r *http.Request) {
	container.InvalidateInstalledOverlaysCache()
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
