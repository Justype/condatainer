package server

import (
	"context"
	"encoding/json"
	"errors"
	"fmt"
	"net/http"
	"os"
	"path/filepath"
	"strings"
	"syscall"
	"time"
)

// handleFSRename serves PATCH /api/fs?path=... — renames and/or moves a
// file or directory. Body takes:
//   - {"new_name": "..."} — rename within the existing parent directory;
//     new_name must be a bare file name with no path separators.
//   - {"dest_dir": "..."} — move to a different (existing) directory,
//     keeping the same base name.
//   - both — move to dest_dir under the name new_name.
//
// Tries a plain os.Rename first and responds synchronously with
// {"path": newPath}. A move across filesystems (EXDEV) instead runs as a
// background copy-then-delete task (see startCrossDeviceMove) and responds
// with {"id": taskID}. Rejects with 409 if the destination already exists.
func (s *srv) handleFSRename(w http.ResponseWriter, r *http.Request) {
	path := r.URL.Query().Get("path")
	if path == "" {
		http.Error(w, "path is required", http.StatusBadRequest)
		return
	}

	var req struct {
		NewName string `json:"new_name"`
		DestDir string `json:"dest_dir"`
	}
	if err := json.NewDecoder(r.Body).Decode(&req); err != nil {
		http.Error(w, "invalid request body", http.StatusBadRequest)
		return
	}

	srcInfo, err := os.Lstat(path)
	if err != nil {
		http.Error(w, err.Error(), http.StatusBadRequest)
		return
	}

	newName := strings.TrimSpace(req.NewName)
	if newName != "" && (newName != filepath.Base(newName) || newName == "." || newName == "..") {
		http.Error(w, "new_name must be a plain file name, not a path", http.StatusBadRequest)
		return
	}

	var newPath string
	switch destDir := strings.TrimSpace(req.DestDir); {
	case destDir != "":
		info, err := os.Stat(destDir)
		if err != nil {
			http.Error(w, err.Error(), http.StatusBadRequest)
			return
		}
		if !info.IsDir() {
			http.Error(w, "dest_dir is not a directory", http.StatusBadRequest)
			return
		}
		// Reject moving a directory into itself or its own subtree — same
		// guard as copy: a plain rename only gets an opaque EINVAL from the
		// kernel, and the EXDEV copy fallback must never attempt it.
		if srcInfo.IsDir() && destInsideSource(path, destDir) {
			http.Error(w, "cannot move a directory into itself", http.StatusBadRequest)
			return
		}
		base := filepath.Base(path)
		if newName != "" {
			base = newName
		}
		newPath = filepath.Join(destDir, base)

	case newName != "":
		newPath = filepath.Join(filepath.Dir(path), newName)

	default:
		http.Error(w, "new_name or dest_dir is required", http.StatusBadRequest)
		return
	}

	if newPath == path {
		http.Error(w, "source and destination are the same", http.StatusBadRequest)
		return
	}
	if _, err := os.Lstat(newPath); err == nil {
		http.Error(w, fmt.Sprintf("%q already exists", filepath.Base(newPath)), http.StatusConflict)
		return
	}

	if err := os.Rename(path, newPath); err != nil {
		if errors.Is(err, syscall.EXDEV) {
			s.startCrossDeviceMove(w, path, newPath)
			return
		}
		http.Error(w, err.Error(), http.StatusInternalServerError)
		return
	}

	writeJSON(w, map[string]interface{}{"path": newPath})
}

// startCrossDeviceMove starts a background task that copies path to
// newPath, then removes the source only after the copy fully succeeds; on
// failure or cancellation the partial destination is removed and the
// source is left in place. Responds with {"id": taskID} immediately.
func (s *srv) startCrossDeviceMove(w http.ResponseWriter, path, newPath string) {
	taskID := fmt.Sprintf("move-%d", time.Now().UnixNano())
	broker := newSSEBroker()
	ctx, cancel := context.WithCancel(s.ctx)
	s.tasks.Store(taskID, &taskEntry{broker: broker, cancel: cancel})
	writeJSON(w, map[string]string{"id": taskID})

	go func() {
		defer cancel()
		defer s.scheduleTaskCleanup(taskID)

		total, err := countTreeEntries(path)
		if err != nil {
			broadcastResult(broker, ctx, fmt.Errorf("scanning %s: %w", path, err))
			return
		}

		copied := 0
		broadcastProgress(broker, 0, total)
		err = copyTree(ctx, path, newPath, true, func() { // preserve times: this is a move
			copied++
			broadcastProgress(broker, copied, total)
		})
		if err != nil {
			os.RemoveAll(newPath) // clean up partial copy; source is untouched
			broadcastResult(broker, ctx, err)
			return
		}

		if err := os.RemoveAll(path); err != nil {
			// the copy is complete and valid; only the source cleanup failed
			broadcastDone(broker, fmt.Errorf("copied to %s but failed to remove source: %w", newPath, err))
			return
		}
		broadcastDone(broker, nil)
	}()
}
