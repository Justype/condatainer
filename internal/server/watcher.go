package server

import (
	"context"
	"os"
	"path/filepath"
	"sync"
	"sync/atomic"
	"time"

	"github.com/Justype/condatainer/internal/config"
	"github.com/Justype/condatainer/internal/helper"
	"github.com/Justype/condatainer/internal/logging"
	"github.com/Justype/condatainer/internal/scheduler"
)

// autoDoneBuffer is added to the estimated end before synthesizing a done file.
// Covers the SIGTERM→SIGKILL grace period (~30s) and NFS write latency.
const autoDoneBuffer = 3 * time.Minute

// schedPollInterval is the minimum time between scheduler queries for pending/starting helpers.
const schedPollInterval = 30 * time.Second

// runningSchedPollInterval is the minimum time between scheduler queries for running helpers.
// Used as a backstop for SIGKILL/OOM-kill/node-crash; normal exits are caught via the done-event trap.
const runningSchedPollInterval = 5 * time.Minute

// tunnelReopenInterval is the minimum time between tunnel reopen attempts for a helper
// whose SSH connection has dropped. Open() dials in a background goroutine so this
// just prevents hammering the SSH server on persistent auth failures.
const tunnelReopenInterval = 30 * time.Second

// helperInfo caches scheduler job ID, status, timing, and service address for an active helper.
// Avoids re-reading history.jsonl on every tick.
type helperInfo struct {
	jobID     string
	name      string // helper script name (e.g. "code-server")
	status    string
	startedAt time.Time
	walltime  time.Duration
	port      int    // service port; 0 = URL-only (no tunnel needed)
	node      string // compute node hostname
}

// watcher polls NFS helper state directories every 2 seconds.
// It opens SSH tunnels when ready files appear and closes them when done files appear.
// When a helper's walltime expires and no done file appears (e.g. SIGKILL), the watcher
// synthesizes one so history and tunnels are cleaned up automatically.
type watcher struct {
	s               *srv
	offsets         map[string]int64      // per-helper message file byte offset
	done            map[string]struct{}   // IDs that have reached a terminal state
	lastSchedPoll   map[string]time.Time  // last scheduler query time per helper
	lastOpenAttempt map[string]time.Time  // last proxies.Open() attempt time per helper
	info            map[string]helperInfo // cached status+timing for active helpers
	lastMsg         sync.Map              // string → string; last log message per helper (read by HTTP handlers)
	runningCount    atomic.Int32          // count of non-terminal helpers; safe to read from HTTP handlers
}

func newWatcher(s *srv) *watcher {
	return &watcher{
		s:               s,
		offsets:         make(map[string]int64),
		done:            make(map[string]struct{}),
		lastSchedPoll:   make(map[string]time.Time),
		lastOpenAttempt: make(map[string]time.Time),
		info:            make(map[string]helperInfo),
	}
}

// RunningCount returns the number of helpers not yet in a terminal state.
// Safe to call from HTTP handler goroutines.
func (w *watcher) RunningCount() int {
	return int(w.runningCount.Load())
}

// LastMsg returns the cached last log message for a helper, or "".
// Safe to call from HTTP handler goroutines.
func (w *watcher) LastMsg(id string) string {
	if v, ok := w.lastMsg.Load(id); ok {
		return v.(string)
	}
	return ""
}

// prePopulateDone seeds the done-map and info-cache from history on startup.
// Completed helpers go into the done-map so they are skipped on every tick.
// Active helpers go into the info cache so history.jsonl is never re-read in the hot path.
func (w *watcher) prePopulateDone() {
	runs, err := helper.LoadHistory()
	if err != nil {
		return
	}
	for _, r := range runs {
		if r.Status == "done" || r.Status == "failed" || r.Status == "cancelled" {
			w.done[r.ID] = struct{}{}
		} else {
			w.info[r.ID] = helperInfo{
				jobID:     r.JobID,
				name:      r.Name,
				status:    r.Status,
				startedAt: r.StartedAt,
				walltime:  r.Walltime,
				port:      r.Port,
				node:      r.Node,
			}
			w.runningCount.Add(1)
		}
	}
}

// reconcileOnStart does a one-time check for all non-terminal jobs at server startup.
// It catches jobs that died during a server outage: walltime-elapsed is checked first
// (free), then the scheduler is queried once per job while it may still have a record.
func (w *watcher) reconcileOnStart(ctx context.Context) {
	sched := scheduler.ActiveScheduler()

	// Collect IDs first to avoid mutating w.info while iterating.
	ids := make([]string, 0, len(w.info))
	for id := range w.info {
		ids = append(ids, id)
	}

	for _, id := range ids {
		if _, done := w.done[id]; done {
			continue
		}
		inf := w.info[id]

		if w.walltimeExpired(inf) {
			logging.FromContext(w.s.ctx).Debug("server: startup reconcile — walltime expired", "id", id)
			w.synthesizeDone(id)
			w.closeDone(id, 0, inf.startedAt.Add(inf.walltime))
			continue
		}

		if inf.jobID == "" || sched == nil {
			continue
		}
		qCtx, cancel := context.WithTimeout(ctx, 10*time.Second)
		status, err := sched.GetJobStatus(qCtx, inf.jobID)
		cancel()
		if err != nil {
			continue // scheduler unavailable — conservative, don't mark done
		}
		switch status {
		case scheduler.JobStatusDone, scheduler.JobStatusFailed:
			logging.FromContext(w.s.ctx).Debug("server: startup reconcile — scheduler says done", "id", id, "jobID", inf.jobID)
			w.synthesizeDone(id)
			w.closeDone(id, 0, time.Now())
		}
	}
}

// walltimeExpired reports whether the helper's known walltime has elapsed (plus autoDoneBuffer).
func (w *watcher) walltimeExpired(inf helperInfo) bool {
	return inf.walltime > 0 && !inf.startedAt.IsZero() &&
		time.Now().After(inf.startedAt.Add(inf.walltime).Add(autoDoneBuffer))
}

// run is the main watcher loop. Blocks until ctx is cancelled.
func (w *watcher) run(ctx context.Context) {
	w.prePopulateDone()
	w.reconcileOnStart(ctx)
	w.restoreExisting()

	ticker := time.NewTicker(2 * time.Second)
	defer ticker.Stop()
	for {
		select {
		case <-ctx.Done():
			return
		case <-ticker.C:
			w.tick()
		}
	}
}

// tick scans all known helper state directories.
func (w *watcher) tick() {
	baseDir := filepath.Join(config.GetUserStateDir(), "helper")
	entries, err := os.ReadDir(baseDir)
	if err != nil {
		return
	}
	for _, e := range entries {
		if !e.IsDir() {
			continue
		}
		id := e.Name()
		if _, isDone := w.done[id]; isDone {
			continue
		}
		w.pollHelper(id)
	}
}

func (w *watcher) pollHelper(id string) {
	if _, isDone := w.done[id]; isDone {
		return
	}
	// Seed info cache on first encounter (helpers submitted after server start).
	if _, known := w.info[id]; !known {
		if run := helper.HistoryEntryForID(id); run != nil {
			w.info[id] = helperInfo{
				jobID:     run.JobID,
				name:      run.Name,
				status:    run.Status,
				startedAt: run.StartedAt,
				walltime:  run.Walltime,
				port:      run.Port,
				node:      run.Node,
			}
			w.runningCount.Add(1)
		}
	}

	r := helper.PollOnce(id, w.offsets[id])
	w.offsets[id] = r.NewOffset

	for _, m := range r.NewMessages {
		w.s.brokers.publish(id, m.Level, m.Text)
		w.lastMsg.Store(id, m.Text)
	}

	if r.Done != nil {
		w.closeDone(id, r.Done.ExitCode, r.Done.Ts)
		return
	}

	// Walltime-elapsed: synthesize done without querying the scheduler.
	if inf, ok := w.info[id]; ok && w.walltimeExpired(inf) {
		logging.FromContext(w.s.ctx).Debug("server: helper walltime expired, synthesizing done", "id", id)
		w.synthesizeDone(id)
		w.closeDone(id, 0, inf.startedAt.Add(inf.walltime))
		return
	}

	if w.syncSchedulerStatus(id) {
		return
	}

	// Reopen tunnel if it dropped (SSH disconnect or initial dial failure).
	// Open() is called in a goroutine so the watcher is not blocked by the SSH dial.
	// Rate-limited to avoid hammering the SSH server on persistent auth failures.
	if inf, ok := w.info[id]; ok && inf.status == "running" && inf.port > 0 && w.s.proxies.Get(id) == nil {
		if last, seen := w.lastOpenAttempt[id]; !seen || time.Since(last) >= tunnelReopenInterval {
			w.lastOpenAttempt[id] = time.Now()
			go w.s.proxies.Open(id, inf.name, inf.node, inf.port)
		}
	}

	if r.Ready != nil {
		rs := r.Ready
		walltime := time.Duration(rs.WalltimeSec) * time.Second
		if rs.Port > 0 && w.s.proxies.Get(id) == nil {
			w.s.proxies.Open(id, w.info[id].name, rs.Node, rs.Port)
		}
		// Update cached info with actual start time, walltime, and service address from the ready event.
		if inf, ok := w.info[id]; ok {
			w.info[id] = helperInfo{
				jobID:     inf.jobID,
				name:      inf.name,
				status:    "running",
				startedAt: rs.Timestamp,
				walltime:  walltime,
				port:      rs.Port,
				node:      rs.Node,
			}
		}
		_ = helper.UpdateHistoryRun(id, func(run *helper.HelperRun) {
			run.Status = "running"
			run.Node = rs.Node
			run.Port = rs.Port
			run.Label = rs.Label
			run.URLPath = rs.URLPath
			run.ExternalURL = rs.ExternalURL
			run.StartedAt = rs.Timestamp
			run.Walltime = walltime
		})
		invalidateHistoryCache()
	}
}

func (w *watcher) syncSchedulerStatus(id string) bool {
	inf, ok := w.info[id]
	if !ok || inf.jobID == "" {
		return false
	}

	var interval time.Duration
	switch inf.status {
	case "pending", "starting":
		interval = schedPollInterval
	case "running":
		interval = runningSchedPollInterval
	default:
		return false
	}

	if last, ok := w.lastSchedPoll[id]; ok && time.Since(last) < interval {
		return false
	}
	sched := scheduler.ActiveScheduler()
	if sched == nil {
		return false
	}
	w.lastSchedPoll[id] = time.Now()
	status, err := sched.GetJobStatus(context.Background(), inf.jobID)
	if err != nil {
		logging.FromContext(w.s.ctx).Debug("server: GetJobStatus failed", "jobID", inf.jobID, "err", err)
		return false
	}
	switch status {
	case scheduler.JobStatusRunning:
		if inf.status == "pending" {
			w.info[id] = helperInfo{jobID: inf.jobID, name: inf.name, status: "starting", startedAt: inf.startedAt, walltime: inf.walltime}
			_ = helper.UpdateHistoryRun(id, func(run *helper.HelperRun) {
				run.Status = "starting"
			})
			invalidateHistoryCache()
		}
	case scheduler.JobStatusDone, scheduler.JobStatusFailed:
		w.synthesizeDone(id)
		w.closeDone(id, 0, time.Now())
		return true
	}
	return false
}

// closeDone closes the tunnel and marks history done/failed for the given helper.
// endedAt should be the time the job actually ended (done-event Ts, or StartedAt+Walltime
// for synthesised expiry) so the stored duration is accurate even after a server outage.
func (w *watcher) closeDone(id string, exitCode int, endedAt time.Time) {
	if w.s.proxies.Get(id) != nil {
		w.s.proxies.Close(id)
		logging.FromContext(w.s.ctx).Debug("server: helper done, tunnel closed", "id", id, "exit", exitCode)
	}
	status := "done"
	if exitCode != 0 {
		status = "failed"
	}
	_ = helper.UpdateHistoryStatus(id, status, endedAt)
	invalidateHistoryCache()
	if _, wasActive := w.info[id]; wasActive {
		w.runningCount.Add(-1)
	}
	w.done[id] = struct{}{}
	delete(w.offsets, id)
	delete(w.lastSchedPoll, id)
	delete(w.lastOpenAttempt, id)
	delete(w.info, id)
	w.lastMsg.Delete(id)
}

// synthesizeDone appends a done event on behalf of a helper that was killed without
// recording one itself (e.g. SIGKILL at walltime, node crash).
func (w *watcher) synthesizeDone(id string) {
	ec := 0
	if err := helper.AppendEvent(id, helper.StateEvent{Type: "done", ExitCode: &ec}); err != nil {
		logging.FromContext(w.s.ctx).Debug("server: synthesize done event failed", "id", id, "err", err)
	}
}

// restoreExisting opens tunnels for helpers that were already running before this server started.
func (w *watcher) restoreExisting() {
	runs, err := helper.RunningHelpers("")
	if err != nil {
		return
	}
	for _, r := range runs {
		rs, err := helper.ReadReadyEvent(r.ID)
		if err != nil || rs == nil || rs.Port == 0 {
			continue
		}
		w.s.proxies.Open(r.ID, r.Name, rs.Node, rs.Port)
		logging.FromContext(w.s.ctx).Debug("server: restored tunnel", "id", r.ID, "node", rs.Node, "port", rs.Port)
	}
}
