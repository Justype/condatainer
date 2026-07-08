package server

import (
	"encoding/json"
	"fmt"
	"net/http"
	"os"
	"time"

	"github.com/Justype/condatainer/internal/config"
	"github.com/Justype/condatainer/internal/helper"
)

// helperListEntry wraps HelperRun with server-computed fields for the API response.
type helperListEntry struct {
	helper.HelperRun
	AccessURL   string `json:"access_url,omitempty"`
	EndsAtMs    int64  `json:"ends_at_ms,omitempty"`
	LastMessage string `json:"last_message,omitempty"`
	WalltimeStr string `json:"walltime_str,omitempty"`
}

// fmtDuration formats a time.Duration as "HH:MM:SS".
func fmtDuration(d time.Duration) string {
	if d == 0 {
		return ""
	}
	h := int(d.Hours())
	m := int(d.Minutes()) % 60
	s := int(d.Seconds()) % 60
	return fmt.Sprintf("%02d:%02d:%02d", h, m, s)
}

// fmtUptime formats the duration since t as a human-readable string.
func fmtUptime(since time.Time) string {
	d := time.Since(since)
	h := int(d.Hours())
	m := int(d.Minutes()) % 60
	if h == 0 {
		return fmt.Sprintf("%dm", m)
	}
	if h < 24 {
		return fmt.Sprintf("%dh %dm", h, m)
	}
	days := h / 24
	h = h % 24
	return fmt.Sprintf("%dd %dh", days, h)
}

// enrichRun wraps a HelperRun with computed fields using the server's port.
func (s *srv) enrichRun(run *helper.HelperRun) helperListEntry {
	e := helperListEntry{HelperRun: *run}
	e.AccessURL = run.AccessURL(s.port)
	if est := run.EstimatedEnd(); !est.IsZero() {
		e.EndsAtMs = est.UnixMilli()
	}
	e.LastMessage = s.watcher.LastMsg(run.ID)
	e.WalltimeStr = fmtDuration(run.Walltime)
	return e
}

func isActiveRunStatus(status string) bool {
	return status == "pending" || status == "starting" || status == "running"
}

// handleStatus serves GET /api/status — server health and running count.
// Uses the watcher's in-memory counter to avoid re-reading history.jsonl on every tick.
func (s *srv) handleStatus(w http.ResponseWriter, r *http.Request) {
	home, _ := os.UserHomeDir()
	scratch := os.Getenv("SCRATCH")
	hostname, _ := os.Hostname()
	writeJSON(w, map[string]interface{}{
		"port":         s.port,
		"hostname":     hostname,
		"pid":          os.Getpid(),
		"running":      s.watcher.RunningCount(),
		"uptime":       fmtUptime(s.startTime),
		"version":      config.Global.Version,
		"home":         home,
		"scratch":      scratch,
		"notification": config.Global.Notification,
	})
}

func writeJSON(w http.ResponseWriter, v interface{}) {
	w.Header().Set("Content-Type", "application/json")
	json.NewEncoder(w).Encode(v) //nolint:errcheck
}
