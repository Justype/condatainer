package server

import (
	"context"
	"encoding/json"
	"fmt"
	"io"
	"net/http"
	"os"
	"sort"
	"strconv"
	"strings"
	"sync"
	"time"

	"github.com/Justype/condatainer/internal/helper"
	"github.com/Justype/condatainer/internal/logging"
	"github.com/Justype/condatainer/internal/scheduler"
)

// historyCache reduces repeated LoadHistory() disk reads across the frequent polling endpoints.
// Entries are reused for up to historyCacheTTL; the watcher's own writes bypass this cache
// (they go directly to disk) so stale data is bounded by the TTL.
const historyCacheTTL = 3 * time.Second

var (
	historyCacheMu   sync.Mutex
	historyCacheRuns []*helper.HelperRun
	historyCacheAt   time.Time
)

// cachedLoadHistory returns a cached history slice, re-reading from disk at most once per historyCacheTTL.
func cachedLoadHistory() ([]*helper.HelperRun, error) {
	historyCacheMu.Lock()
	defer historyCacheMu.Unlock()
	if historyCacheRuns != nil && time.Since(historyCacheAt) < historyCacheTTL {
		return historyCacheRuns, nil
	}
	runs, err := helper.LoadHistory()
	if err != nil {
		return nil, err
	}
	historyCacheRuns = runs
	historyCacheAt = time.Now()
	return runs, nil
}

// invalidateHistoryCache clears the cached history so the next API request reads fresh data.
// Called by the watcher after each history write so status transitions appear promptly.
func invalidateHistoryCache() {
	historyCacheMu.Lock()
	historyCacheRuns = nil
	historyCacheMu.Unlock()
}

// handleHelpersList serves GET /api/helpers?q=&limit=&status=
func (s *srv) handleHelpersList(w http.ResponseWriter, r *http.Request) {
	q := r.URL.Query().Get("q")
	status := r.URL.Query().Get("status")

	runs, err := cachedLoadHistory()
	if err != nil {
		http.Error(w, "reading history", http.StatusInternalServerError)
		return
	}

	limit := 100
	if l := r.URL.Query().Get("limit"); l != "" {
		if n, err := strconv.Atoi(l); err == nil && n > 0 {
			limit = n
		}
	}

	var out []*helper.HelperRun
	for i := len(runs) - 1; i >= 0; i-- {
		run := runs[i]
		if status == "active" {
			if !isActiveRunStatus(run.Status) {
				continue
			}
		} else if status != "" && run.Status != status {
			continue
		}
		if q != "" && !strings.Contains(run.Name, q) && !strings.Contains(run.CWD, q) {
			continue
		}
		out = append(out, run)
		if len(out) >= limit {
			break
		}
	}

	var resp []helperListEntry
	for _, run := range out {
		resp = append(resp, s.enrichRun(run))
	}
	writeJSON(w, resp)
}

// handleHelperLogs serves GET /api/helpers/{id}/logs as SSE.
func (s *srv) handleHelperLogs(w http.ResponseWriter, r *http.Request, id string) {
	broker := s.brokers.get(id)

	flusher, ok := w.(http.Flusher)
	if !ok {
		http.Error(w, "streaming not supported", http.StatusInternalServerError)
		return
	}
	w.Header().Set("Content-Type", "text/event-stream")
	w.Header().Set("Cache-Control", "no-cache")
	w.Header().Set("Connection", "keep-alive")

	const maxBackfill = 200
	msgs, _, _ := helper.ReadMessages(id, 0)
	if len(msgs) > maxBackfill {
		msgs = msgs[len(msgs)-maxBackfill:]
	}
	for _, m := range msgs {
		data, _ := json.Marshal(m)
		fmt.Fprintf(w, "data: %s\n\n", data) //nolint:errcheck
	}
	flusher.Flush()

	ch, _, _ := broker.subscribe() // job-log broker has no "finished" concept, unlike task brokers
	defer broker.unsubscribe(ch)
	ticker := time.NewTicker(25 * time.Second)
	defer ticker.Stop()
	for {
		select {
		case <-r.Context().Done():
			return
		case data := <-ch:
			fmt.Fprintf(w, "data: %s\n\n", data) //nolint:errcheck
			flusher.Flush()
		case <-ticker.C:
			fmt.Fprintf(w, ": keepalive\n\n") //nolint:errcheck
			flusher.Flush()
		}
	}
}

// defaultJobLogTail is the default number of lines returned by handleHelperJobLog.
// Prevents loading arbitrarily large log files into the browser.
const defaultJobLogTail = 1000

// handleHelperJobLog serves GET /api/helpers/{id}/joblog — returns the last N lines
// of the job output file (?tail=N, default 1000). Sets X-Log-Truncated: true when
// the file contains more lines than the requested tail.
func (s *srv) handleHelperJobLog(w http.ResponseWriter, r *http.Request, id string) {
	logPath := helper.JobLogFilePath(id)
	tailN := defaultJobLogTail
	if q := r.URL.Query().Get("tail"); q != "" {
		if n, err := strconv.Atoi(q); err == nil && n > 0 {
			tailN = n
		}
	}
	data, truncated, err := tailFile(logPath, tailN)
	if err != nil {
		http.Error(w, "job log not found", http.StatusNotFound)
		return
	}
	w.Header().Set("Content-Type", "text/plain; charset=utf-8")
	if truncated {
		w.Header().Set("X-Log-Truncated", "true")
	}
	w.Write(data) //nolint:errcheck
}

// tailFile returns the last maxLines lines of the file at path, reading backwards
// from EOF in 64 KiB chunks. truncated is true when the file has more lines than maxLines.
func tailFile(path string, maxLines int) (data []byte, truncated bool, err error) {
	f, err := os.Open(path)
	if err != nil {
		return nil, false, err
	}
	defer f.Close()

	size, err := f.Seek(0, io.SeekEnd)
	if err != nil {
		return nil, false, err
	}
	if size == 0 {
		return nil, false, nil
	}

	const chunkSize = 65536
	pos := size
	newlinesSeen := 0

	for pos > 0 {
		n := int64(chunkSize)
		if pos < n {
			n = pos
		}
		pos -= n
		chunk := make([]byte, n)
		if _, err = f.ReadAt(chunk, pos); err != nil {
			return nil, false, err
		}
		for i := int(n) - 1; i >= 0; i-- {
			if chunk[i] == '\n' {
				newlinesSeen++
				if newlinesSeen > maxLines {
					startAt := pos + int64(i) + 1
					out := make([]byte, size-startAt)
					_, err = f.ReadAt(out, startAt)
					return out, true, err
				}
			}
		}
	}

	// File has ≤ maxLines lines — return all of it.
	if _, err = f.Seek(0, io.SeekStart); err != nil {
		return nil, false, err
	}
	data, err = io.ReadAll(f)
	return data, false, err
}

// handleHelperStop serves POST /api/helpers/{id}/stop.
func (s *srv) handleHelperStop(w http.ResponseWriter, r *http.Request, id string) {
	if r.Method != http.MethodPost {
		http.Error(w, "POST required", http.StatusMethodNotAllowed)
		return
	}
	log := logging.FromContext(s.ctx)
	if run := helper.HistoryEntryForID(id); run == nil {
		log.Warn("server: stop — history entry not found", "id", id)
	} else if run.JobID == "" {
		if err := helper.KillHeadlessProcess(id); err != nil {
			log.Debug("server: stop headless helper failed", "id", id, "err", err)
		}
	} else if sched := scheduler.ActiveScheduler(); sched == nil {
		log.Warn("server: stop — no scheduler available", "id", id, "jobID", run.JobID)
	} else {
		cancelCtx, cancelCancel := context.WithTimeout(context.Background(), 30*time.Second)
		if err := sched.CancelJob(cancelCtx, run.JobID); err != nil {
			log.Warn("server: cancel job failed", "id", id, "jobID", run.JobID, "err", err)
		}
		cancelCancel()
	}
	// closeDone cleans up the proxy, updates history, decrements runningCount, and
	// moves the ID to the done-map so the watcher skips it on future ticks.
	// This ensures the count is correct immediately even if the trap never fires
	// (e.g. SIGKILL, node crash, NFS issue).
	s.watcher.closeDone(id, 0, time.Now())
	w.WriteHeader(http.StatusOK)
	fmt.Fprintf(w, `{"ok":true}`)
}

// handleHelperDelete serves DELETE /api/helpers/{id}/delete — removes a finished history entry.
// Active entries (pending/starting/running) are rejected with 409 Conflict.
func (s *srv) handleHelperDelete(w http.ResponseWriter, r *http.Request, id string) {
	if r.Method != http.MethodDelete {
		http.Error(w, "DELETE required", http.StatusMethodNotAllowed)
		return
	}
	if run := helper.HistoryEntryForID(id); run != nil && isActiveRunStatus(run.Status) {
		http.Error(w, "cannot delete an active helper; stop it first", http.StatusConflict)
		return
	}
	if err := helper.DeleteHistoryEntry(id); err != nil {
		http.Error(w, err.Error(), http.StatusInternalServerError)
		return
	}
	invalidateHistoryCache()
	writeJSON(w, map[string]bool{"ok": true})
}

// handleGpuOptions serves GET /api/gpu-options — returns unique GPU type names from the active scheduler.
func (s *srv) handleGpuOptions(w http.ResponseWriter, r *http.Request) {
	w.Header().Set("Content-Type", "application/json")
	sched := scheduler.ActiveScheduler()
	if sched == nil {
		json.NewEncoder(w).Encode([]string{}) //nolint:errcheck
		return
	}
	ctx, cancel := context.WithTimeout(r.Context(), 8*time.Second)
	defer cancel()
	info, err := sched.GetClusterInfo(ctx)
	if err != nil || info == nil {
		json.NewEncoder(w).Encode([]string{}) //nolint:errcheck
		return
	}
	seen := make(map[string]bool)
	for _, g := range info.AvailableGpus {
		if g.Type != "" && g.Type != "gpu" {
			seen[g.Type] = true
		}
	}
	opts := make([]string, 0, len(seen))
	for t := range seen {
		opts = append(opts, t)
	}
	sort.Strings(opts)
	json.NewEncoder(w).Encode(opts) //nolint:errcheck
}
