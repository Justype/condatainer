package server

import (
	"context"
	"encoding/json"
	"errors"
	"fmt"
	"net/http"
	"os"
	"os/exec"
	"strings"
	"syscall"
	"time"

	"github.com/Justype/condatainer/internal/config"
	"github.com/Justype/condatainer/internal/container"
	"github.com/Justype/condatainer/internal/utils"
)

// createExecutable returns the binary to run for install tasks.
// Overridable in tests.
var createExecutable = os.Executable

// handleCreate serves POST /api/create — installs an overlay by running the
// condatainer CLI itself (`create <name>`) as a subprocess and streaming its
// output. The CLI handles build scripts, conda fallback, dependencies, and
// scheduler submission (exit code config.ExitCodeJobsSubmitted). Returns a
// task ID immediately; progress is streamed via GET /api/tasks/{id}/stream.
//
// Request: {"name":"cellranger/9.0.1","answers":["..."]} — answers are piped
// to the CLI's stdin, one per #INTERACTIVE prompt; without answers the CLI
// runs with --yes (empty responses).
func (s *srv) handleCreate(w http.ResponseWriter, r *http.Request) {
	if r.Method != http.MethodPost {
		http.Error(w, "POST required", http.StatusMethodNotAllowed)
		return
	}
	type createReq struct {
		Name    string   `json:"name"`
		Answers []string `json:"answers"`
	}
	var req createReq
	if err := json.NewDecoder(r.Body).Decode(&req); err != nil {
		http.Error(w, "bad request body", http.StatusBadRequest)
		return
	}
	raw := strings.TrimSpace(req.Name)
	if raw == "" {
		http.Error(w, "name required", http.StatusBadRequest)
		return
	}
	// Reject flag-looking input before normalization ("--x" would normalize to "/x").
	if strings.HasPrefix(raw, "-") {
		http.Error(w, "invalid name", http.StatusBadRequest)
		return
	}
	name := utils.NormalizeNameVersion(raw)

	exe, err := createExecutable()
	if err != nil {
		http.Error(w, "cannot locate condatainer binary: "+err.Error(), http.StatusInternalServerError)
		return
	}

	taskID := fmt.Sprintf("create-%d", time.Now().UnixNano())
	broker := newSSEBroker()
	ctx, cancel := context.WithCancel(s.ctx)
	s.tasks.Store(taskID, &taskEntry{broker: broker, cancel: cancel})
	writeJSON(w, map[string]string{"id": taskID})

	bw := &brokerWriter{broker}

	go func() {
		defer cancel()
		defer s.scheduleTaskCleanup(taskID)

		args := []string{"create", name}
		if len(req.Answers) == 0 {
			args = append(args, "--yes")
		}
		fmt.Fprintf(bw, "Running: condatainer %s\n", strings.Join(args, " "))

		cmd := exec.CommandContext(ctx, exe, args...)
		cmd.Stdout = bw
		cmd.Stderr = bw
		if len(req.Answers) > 0 {
			cmd.Stdin = strings.NewReader(strings.Join(req.Answers, "\n") + "\n")
		}
		// Run in its own process group so cancellation also kills build
		// children (apptainer, mksquashfs); the CLI cleans up on SIGTERM.
		cmd.SysProcAttr = &syscall.SysProcAttr{Setpgid: true}
		cmd.Cancel = func() error {
			return syscall.Kill(-cmd.Process.Pid, syscall.SIGTERM)
		}
		cmd.WaitDelay = 15 * time.Second

		err := cmd.Run()
		if err == nil {
			container.InvalidateInstalledOverlaysCache()
			fmt.Fprintf(bw, "Done.\n")
			result, _ := json.Marshal(map[string]interface{}{"t": "done", "ok": true, "name": name})
			broker.publishFinal(result)
			return
		}
		var exitErr *exec.ExitError
		if errors.As(err, &exitErr) && exitErr.ExitCode() == config.ExitCodeJobsSubmitted {
			fmt.Fprintf(bw, "Build job(s) submitted to the scheduler — the overlay will appear once the job finishes.\n")
			result, _ := json.Marshal(map[string]interface{}{"t": "done", "ok": true, "name": name, "submitted": true})
			broker.publishFinal(result)
			return
		}
		broadcastResult(broker, ctx, fmt.Errorf("create %s: %w", name, err))
	}()
}
