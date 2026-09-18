package server

import (
	"context"
	"encoding/json"
	"fmt"
	"net/http"
	"os"
	"time"

	"log/slog"

	"github.com/Justype/condatainer/internal/logging"
	"github.com/Justype/condatainer/internal/logging/weblog"
	"github.com/Justype/condatainer/internal/project"
	"github.com/Justype/condatainer/internal/project/orchestrate"
)

// resolveProjectCwd falls back to the server process's own directory when
// cwd is empty, mirroring cmd/project.go's projectRootOrInit — the server
// has no ambient "current directory" of its own, so every project-scoped
// handler takes cwd explicitly, the same convention handleEnvCheck and
// handleHelperStart already follow.
func resolveProjectCwd(cwd string) string {
	if cwd != "" {
		return cwd
	}
	wd, _ := os.Getwd()
	return wd
}

// handleProjectStatus serves GET /api/project/status?cwd=... — a
// synchronous, read-only read, since project.Status does no scanning-and-
// writing the way project lock does.
func (s *srv) handleProjectStatus(w http.ResponseWriter, r *http.Request) {
	cwd := resolveProjectCwd(r.URL.Query().Get("cwd"))
	status, err := project.Status(cwd)
	if err != nil {
		http.Error(w, err.Error(), http.StatusInternalServerError)
		return
	}
	writeJSON(w, status)
}

// handleProjectLock serves POST /api/project/lock — creates or updates
// cwd's project lock. Returns a task ID immediately; progress is streamed
// via GET /api/tasks/{id}/stream, the same pattern handleOverlayEdit uses.
func (s *srv) handleProjectLock(w http.ResponseWriter, r *http.Request) {
	if r.Method != http.MethodPost {
		http.Error(w, "POST required", http.StatusMethodNotAllowed)
		return
	}
	type lockReq struct {
		CWD string `json:"cwd"`
	}
	var req lockReq
	if r.Body != nil {
		_ = json.NewDecoder(r.Body).Decode(&req)
	}
	root := resolveProjectCwd(req.CWD)

	taskID := fmt.Sprintf("project-lock-%d", time.Now().UnixNano())
	broker := newSSEBroker()
	ctx, cancel := context.WithCancel(s.ctx)
	s.tasks.Store(taskID, &taskEntry{broker: broker, cancel: cancel})
	writeJSON(w, map[string]string{"id": taskID})

	bw := &brokerWriter{broker}
	ctx = logging.WithLogger(ctx, slog.New(weblog.New(bw)))
	ctx = logging.WithWriter(ctx, bw)

	go func() {
		defer cancel()
		defer s.scheduleTaskCleanup(taskID)

		result, err := orchestrate.Lock(ctx, root, orchestrate.LockOptions{})
		if err != nil {
			broadcastResult(broker, ctx, err)
			return
		}
		if len(result.Failed) > 0 {
			for _, e := range result.Failed {
				fmt.Fprintf(bw, "%v\n", e)
			}
		}
		if len(result.UnpublishedFrozenEnv) > 0 {
			fmt.Fprintln(bw, "A frozen environment cannot be rebuilt; run `condatainer project push` so another checkout can restore this project")
		}
		if len(result.UnpinnedHelperOverlays) > 0 {
			fmt.Fprintln(bw, "Overlays used by helpers but not pinned:")
			for _, name := range result.UnpinnedHelperOverlays {
				fmt.Fprintf(bw, "  %s\n", name)
			}
		}
		fmt.Fprintf(bw, "%d pin(s), %d failed.\n", len(result.Current.Pins), len(result.Failed))
		type doneResult struct {
			T                      string              `json:"t"`
			OK                     bool                `json:"ok"`
			Root                   string              `json:"root"`
			PinCount               int                 `json:"pin_count"`
			UnpinnedHelperOverlays []string            `json:"unpinned_helper_overlays,omitempty"`
			ManualPinUsage         map[string][]string `json:"manual_pin_usage,omitempty"`
		}
		data, _ := json.Marshal(doneResult{
			T: "done", OK: len(result.Failed) == 0, Root: result.Root,
			PinCount: len(result.Current.Pins), UnpinnedHelperOverlays: result.UnpinnedHelperOverlays,
			ManualPinUsage: result.ManualPinUsage,
		})
		broker.publishFinal(data)
	}()
}
