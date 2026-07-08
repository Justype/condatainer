package server

import (
	"context"
	"encoding/json"
	"fmt"
	"net/http"
	"os"
	"path/filepath"
	"strings"
	"time"
)

// handleFSCopy serves POST /api/fs/copy?path=... — copies a file or
// directory to an existing directory. Body:
// {"dest_dir": "...", "new_name": "..."} — new_name is optional (must be a
// bare file name). Without new_name a conflicting destination auto-renames
// to "name_copy" / "name_copyN"; with new_name a conflict is rejected with
// 409. Runs as a
// background task: returns a task ID immediately and streams
// {"t":"progress","current":n,"total":n} events via
// GET /api/tasks/{id}/stream. On failure or cancellation the partial
// destination is removed; the source is never touched.
func (s *srv) handleFSCopy(w http.ResponseWriter, r *http.Request) {
	if r.Method != http.MethodPost {
		http.Error(w, "POST required", http.StatusMethodNotAllowed)
		return
	}
	path := r.URL.Query().Get("path")
	if path == "" {
		http.Error(w, "path is required", http.StatusBadRequest)
		return
	}
	srcInfo, err := os.Lstat(path)
	if err != nil {
		http.Error(w, err.Error(), http.StatusBadRequest)
		return
	}

	var req struct {
		DestDir string `json:"dest_dir"`
		NewName string `json:"new_name"`
	}
	if err := json.NewDecoder(r.Body).Decode(&req); err != nil {
		http.Error(w, "invalid request body", http.StatusBadRequest)
		return
	}
	destDir := strings.TrimSpace(req.DestDir)
	if destDir == "" {
		http.Error(w, "dest_dir is required", http.StatusBadRequest)
		return
	}
	newName := strings.TrimSpace(req.NewName)
	if newName != "" && (newName != filepath.Base(newName) || newName == "." || newName == "..") {
		http.Error(w, "new_name must be a plain file name, not a path", http.StatusBadRequest)
		return
	}
	info, err := os.Stat(destDir)
	if err != nil {
		http.Error(w, err.Error(), http.StatusBadRequest)
		return
	}
	if !info.IsDir() {
		http.Error(w, "dest_dir is not a directory", http.StatusBadRequest)
		return
	}

	// Reject copying a directory into itself or its own subtree — copyTree
	// would otherwise walk the source while writing into it.
	if srcInfo.IsDir() && destInsideSource(path, destDir) {
		http.Error(w, "cannot copy a directory into itself", http.StatusBadRequest)
		return
	}

	base := filepath.Base(path)
	if newName != "" {
		base = newName
	}
	newPath := filepath.Join(destDir, base)
	if newName == "" {
		// no explicit name: a conflict (including copying into the source's
		// own directory) auto-renames to "name_copy" instead of erroring
		newPath = uniquePath(newPath, srcInfo.IsDir())
	} else {
		if newPath == path {
			http.Error(w, "source and destination are the same", http.StatusBadRequest)
			return
		}
		if _, err := os.Lstat(newPath); err == nil {
			http.Error(w, fmt.Sprintf("%q already exists", filepath.Base(newPath)), http.StatusConflict)
			return
		}
	}

	taskID := fmt.Sprintf("copy-%d", time.Now().UnixNano())
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
		err = copyTree(ctx, path, newPath, false, func() {
			copied++
			broadcastProgress(broker, copied, total)
		})
		if err != nil {
			os.RemoveAll(newPath) // clean up partial copy on failure/cancellation
			broadcastResult(broker, ctx, err)
			return
		}
		broadcastDone(broker, nil)
	}()
}
