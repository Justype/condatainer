package helper

import (
	"bufio"
	"encoding/json"
	"fmt"
	"os"
	"path/filepath"
	"strconv"
	"strings"
	"syscall"
	"time"

	"github.com/Justype/condatainer/internal/config"
	"github.com/Justype/condatainer/internal/utils"
)

// HelperRun records one execution of a helper (stored in JSONL history).
type HelperRun struct {
	ID          string            `json:"id"`            // "{name}-YYYYMMDD-HHMMSS"
	Name        string            `json:"name"`          // helper script name
	JobID       string            `json:"job_id"`        // scheduler job ID or ""
	Node        string            `json:"node"`          // compute node hostname
	Port        int               `json:"port"`          // 0 = URL-only service (vscode-tunnel)
	Label       string            `json:"label"`         // human-readable label
	CWD         string            `json:"cwd"`           // working directory
	Walltime    time.Duration     `json:"walltime_secs"` // from #TIME: / -t flag
	CPUs        int               `json:"cpus,omitempty"`
	Mem         string            `json:"mem,omitempty"`         // e.g. "32GB", "1536MB" — via utils.FormatMemoryMB
	GPU         string            `json:"gpu,omitempty"`         // e.g. "a100:1"
	BaseImage   string            `json:"base_image,omitempty"`  // overridden with -b
	EnvOverlay  string            `json:"env_overlay,omitempty"` // writable .img (-e)
	Overlays    []string          `json:"overlays,omitempty"`    // read-only .sqf (-o)
	Params      map[string]string `json:"params,omitempty"`      // resolved #PARAM: values
	URLPath     string            `json:"url_path"`              // appended to proxy URL (e.g. "?token=abc123")
	ExternalURL string            `json:"external_url"`          // shown as-is, no proxy (e.g. vscode-tunnel URL)
	StartedAt   time.Time         `json:"started_at"`
	EndedAt     *time.Time        `json:"ended_at,omitempty"`
	Status      string            `json:"status"` // "running"|"done"|"failed"
}

// AccessURL returns the link the user should open.
// ExternalURL takes priority; otherwise returns the subdomain proxy URL
// (http://{id}.localhost:{port}/{urlPath}).
func (r *HelperRun) AccessURL(serverPort int) string {
	if r.ExternalURL != "" {
		return r.ExternalURL
	}
	if serverPort > 0 && r.Port > 0 {
		urlPath := strings.TrimPrefix(r.URLPath, "/")
		return fmt.Sprintf("http://%s.localhost:%d/%s", r.ID, serverPort, urlPath)
	}
	return ""
}

// EstimatedEnd returns StartedAt + Walltime; zero value if walltime is unknown.
func (r *HelperRun) EstimatedEnd() time.Time {
	if r.Walltime == 0 {
		return time.Time{}
	}
	return r.StartedAt.Add(r.Walltime)
}

// StateEvent is one entry in the unified events JSONL file per helper run.
// Type is one of "msg", "ready", or "done".
type StateEvent struct {
	Type string    `json:"type"`
	Ts   time.Time `json:"ts"`

	// "msg" fields
	Level string `json:"level,omitempty"` // "info"|"warn"|"error"
	Text  string `json:"text,omitempty"`

	// "ready" fields
	Port        int    `json:"port,omitempty"`
	Node        string `json:"node,omitempty"`
	Label       string `json:"label,omitempty"`
	WalltimeSec int64  `json:"walltime_secs,omitempty"`
	JobID       string `json:"job_id,omitempty"`
	URLPath     string `json:"url_path,omitempty"`
	ExternalURL string `json:"external_url,omitempty"`

	// "done" fields — pointer so exit_code:0 is preserved in JSON (not omitted by omitempty)
	ExitCode *int `json:"exit_code,omitempty"`
}

// ReadyState is derived from a "ready" StateEvent; used by callers that need
// the ready information as a typed struct.
type ReadyState struct {
	Port        int       `json:"port"`
	Node        string    `json:"node"`
	Label       string    `json:"label"`
	Timestamp   time.Time `json:"ts"`
	WalltimeSec int64     `json:"walltime_secs"`
	JobID       string    `json:"job_id"`
	URLPath     string    `json:"url_path"`
	ExternalURL string    `json:"external_url"`
}

// DoneState is derived from a "done" StateEvent.
type DoneState struct {
	ExitCode int       `json:"exit_code"`
	Ts       time.Time `json:"ts"`
}

// MessageLine is derived from a "msg" StateEvent.
type MessageLine struct {
	Level string    `json:"level"`
	Text  string    `json:"text"`
	Ts    time.Time `json:"ts"`
}

// StateDir returns the NFS state directory for a helper run.
//
// When called from inside a helper job (e.g. _server_message, _server_ready,
// _server_done), the wrapper has already set CNT_HELPER_STATE_DIR to the
// exact absolute path. Using it directly avoids any XDG/HOME resolution
// differences inside the Apptainer container — so messages are always written
// to the same directory the server watcher reads, even when the server is down.
func StateDir(id string) string {
	if envID := os.Getenv("CNT_HELPER_ID"); envID != "" && envID == id {
		if envDir := os.Getenv("CNT_HELPER_STATE_DIR"); envDir != "" {
			return envDir
		}
	}
	return config.GetHelperStateDir(id)
}

// EventsFilePath returns the path to the events JSONL file for a helper run.
func EventsFilePath(id string) string {
	return filepath.Join(StateDir(id), "events")
}

// JobLogFilePath returns the path to the job output log.
func JobLogFilePath(id string) string {
	return filepath.Join(StateDir(id), "job.log")
}

// PidFilePath returns the path to the ephemeral PID file for a headless helper run.
func PidFilePath(id string) string {
	return filepath.Join(StateDir(id), "pid")
}

// WriteHelperPid writes the bash wrapper's PID to {stateDir}/pid.
func WriteHelperPid(id string, pid int) error {
	return os.WriteFile(PidFilePath(id), []byte(strconv.Itoa(pid)+"\n"), utils.PermFile)
}

// ReadHelperPid reads the PID from {stateDir}/pid.
// Returns an error if the file does not exist or contains an invalid PID.
func ReadHelperPid(id string) (int, error) {
	data, err := os.ReadFile(PidFilePath(id))
	if err != nil {
		return 0, err
	}
	pid, err := strconv.Atoi(strings.TrimSpace(string(data)))
	if err != nil || pid <= 0 {
		return 0, fmt.Errorf("helper %s: invalid pid file content", id)
	}
	return pid, nil
}

// KillHeadlessProcess sends SIGTERM to the process group of the headless helper.
// PGID == bash wrapper PID because ExecutePlan sets Setpgid: true.
// The wrapper's SIGTERM trap calls condatainer _server_done and exits cleanly.
func KillHeadlessProcess(id string) error {
	pid, err := ReadHelperPid(id)
	if err != nil {
		return fmt.Errorf("reading pid file: %w", err)
	}
	return syscall.Kill(-pid, syscall.SIGTERM)
}

// IsHeadlessProcessAlive reports whether the process in {stateDir}/pid is still
// running. Returns false if the pid file is absent, unreadable, or the process
// no longer exists.
func IsHeadlessProcessAlive(id string) bool {
	pid, err := ReadHelperPid(id)
	if err != nil {
		return false
	}
	p, err := os.FindProcess(pid)
	if err != nil {
		return false
	}
	return p.Signal(syscall.Signal(0)) == nil
}

// AppendEvent appends a StateEvent to the events file for the given helper ID.
// Sets Ts to now if zero. The state directory must already exist (created by newHelperID).
func AppendEvent(id string, ev StateEvent) error {
	if ev.Ts.IsZero() {
		ev.Ts = time.Now()
	}
	path := EventsFilePath(id)
	f, err := os.OpenFile(path, os.O_APPEND|os.O_CREATE|os.O_WRONLY, utils.PermFile)
	if err != nil {
		return err
	}
	defer f.Close()
	data, err := json.Marshal(ev)
	if err != nil {
		return err
	}
	_, err = f.Write(append(data, '\n'))
	return err
}

// ReadEventsFrom reads StateEvents from the events file starting at byte offset.
// Returns the events, the new byte offset, and any error.
// Returns offset unchanged and nil slice if the file does not exist yet.
func ReadEventsFrom(id string, offset int64) ([]StateEvent, int64, error) {
	path := EventsFilePath(id)
	f, err := os.Open(path)
	if err != nil {
		if os.IsNotExist(err) {
			return nil, offset, nil
		}
		return nil, offset, err
	}
	defer f.Close()
	if offset > 0 {
		if _, err := f.Seek(offset, 0); err != nil {
			return nil, offset, err
		}
	}
	var events []StateEvent
	scanner := bufio.NewScanner(f)
	for scanner.Scan() {
		var ev StateEvent
		if err := json.Unmarshal(scanner.Bytes(), &ev); err == nil {
			events = append(events, ev)
		}
	}
	pos, _ := f.Seek(0, 1)
	return events, pos, scanner.Err()
}

// ReadReadyEvent returns the first "ready" event for the given helper, or nil if none yet.
func ReadReadyEvent(id string) (*ReadyState, error) {
	events, _, err := ReadEventsFrom(id, 0)
	if err != nil {
		return nil, err
	}
	for _, ev := range events {
		if ev.Type == "ready" {
			return &ReadyState{
				Port:        ev.Port,
				Node:        ev.Node,
				Label:       ev.Label,
				Timestamp:   ev.Ts,
				WalltimeSec: ev.WalltimeSec,
				JobID:       ev.JobID,
				URLPath:     ev.URLPath,
				ExternalURL: ev.ExternalURL,
			}, nil
		}
	}
	return nil, nil
}

// ReadDoneEvent returns the first "done" event for the given helper, or nil if none yet.
func ReadDoneEvent(id string) (*DoneState, error) {
	events, _, err := ReadEventsFrom(id, 0)
	if err != nil {
		return nil, err
	}
	for _, ev := range events {
		if ev.Type == "done" {
			exitCode := 0
			if ev.ExitCode != nil {
				exitCode = *ev.ExitCode
			}
			return &DoneState{ExitCode: exitCode, Ts: ev.Ts}, nil
		}
	}
	return nil, nil
}

// ReadMessages reads "msg" events from the events file starting at byte offset.
// Returns message lines, the new offset, and any error.
func ReadMessages(id string, offset int64) ([]MessageLine, int64, error) {
	events, newOffset, err := ReadEventsFrom(id, offset)
	var msgs []MessageLine
	for _, ev := range events {
		if ev.Type == "msg" {
			msgs = append(msgs, MessageLine{Level: ev.Level, Text: ev.Text, Ts: ev.Ts})
		}
	}
	return msgs, newOffset, err
}

// HistoryPath returns the JSONL history file path.
func HistoryPath() string {
	return config.GetHelperHistoryPath()
}

// AppendHistory appends a HelperRun record to the JSONL history file.
func AppendHistory(run *HelperRun) error {
	path := HistoryPath()
	if path == "" {
		return fmt.Errorf("cannot determine history path")
	}
	if err := os.MkdirAll(filepath.Dir(path), utils.PermDir); err != nil {
		return err
	}
	f, err := os.OpenFile(path, os.O_APPEND|os.O_CREATE|os.O_WRONLY, utils.PermFile)
	if err != nil {
		return err
	}
	defer f.Close()
	data, err := json.Marshal(run)
	if err != nil {
		return err
	}
	_, err = f.Write(append(data, '\n'))
	return err
}

// LoadHistory reads the JSONL history file and returns all records (oldest first).
func LoadHistory() ([]*HelperRun, error) {
	path := HistoryPath()
	if path == "" {
		return nil, nil
	}
	f, err := os.Open(path)
	if err != nil {
		if os.IsNotExist(err) {
			return nil, nil
		}
		return nil, err
	}
	defer f.Close()

	var runs []*HelperRun
	scanner := bufio.NewScanner(f)
	for scanner.Scan() {
		var run HelperRun
		if err := json.Unmarshal(scanner.Bytes(), &run); err == nil {
			runs = append(runs, &run)
		}
	}
	return runs, scanner.Err()
}

// UpdateHistoryStatus rewrites the JSONL history, setting status and endedAt for the given ID.
// endedAt should reflect when the job actually ended (e.g. the done-event timestamp or
// StartedAt+Walltime for synthesised expiry), not merely when the server noticed.
func UpdateHistoryStatus(id, status string, endedAt time.Time) error {
	return UpdateHistoryRun(id, func(r *HelperRun) {
		r.Status = status
		r.EndedAt = &endedAt
	})
}

// withHistoryLock acquires an exclusive flock on path+".lock" for the duration of fn.
func withHistoryLock(path string, fn func() error) error {
	lf, err := os.OpenFile(path+".lock", os.O_CREATE|os.O_WRONLY, utils.PermFile)
	if err != nil {
		return err
	}
	defer lf.Close()
	if err := syscall.Flock(int(lf.Fd()), syscall.LOCK_EX); err != nil {
		return err
	}
	defer syscall.Flock(int(lf.Fd()), syscall.LOCK_UN) //nolint:errcheck
	return fn()
}

// UpdateHistoryRun rewrites the JSONL history, applying fn to the record for id.
// If no record with that ID exists the call is a no-op.
// Uses an exclusive file lock to prevent concurrent writers from losing updates.
func UpdateHistoryRun(id string, fn func(*HelperRun)) error {
	path := HistoryPath()
	if path == "" {
		return fmt.Errorf("cannot determine history path")
	}
	return withHistoryLock(path, func() error {
		runs, err := LoadHistory()
		if err != nil {
			return err
		}
		updated := false
		for _, r := range runs {
			if r.ID == id {
				fn(r)
				updated = true
			}
		}
		if !updated {
			return nil
		}
		return rewriteHistory(path, runs)
	})
}

// DeleteHistoryEntry removes the entry with the given ID from the JSONL history file
// and deletes the associated state directory (events, job.log, etc.).
// It is a no-op if the ID is not found. Active (in-flight) entries are rejected
// so callers must stop/cancel the job before deleting its history record.
func DeleteHistoryEntry(id string) error {
	path := HistoryPath()
	if path == "" {
		return fmt.Errorf("cannot determine history path")
	}
	if err := withHistoryLock(path, func() error {
		runs, err := LoadHistory()
		if err != nil {
			return err
		}
		filtered := runs[:0]
		for _, r := range runs {
			if r.ID != id {
				filtered = append(filtered, r)
			}
		}
		if len(filtered) == len(runs) {
			return nil // ID not found — no-op
		}
		return rewriteHistory(path, filtered)
	}); err != nil {
		return err
	}
	// Remove the per-run state directory (events, job.log, …); best-effort.
	stateDir := config.GetHelperStateDir(id)
	if stateDir != "" {
		os.RemoveAll(stateDir) //nolint:errcheck
	}
	return nil
}

// rewriteHistory writes all runs back to the JSONL file atomically (tmp then rename).
func rewriteHistory(path string, runs []*HelperRun) error {
	tmp := path + ".tmp"
	f, err := utils.CreateFileWritable(tmp)
	if err != nil {
		return err
	}
	for _, r := range runs {
		data, err := json.Marshal(r)
		if err != nil {
			f.Close()
			os.Remove(tmp)
			return err
		}
		if _, err := f.Write(append(data, '\n')); err != nil {
			f.Close()
			os.Remove(tmp)
			return err
		}
	}
	f.Close()
	return os.Rename(tmp, path)
}

// isActiveStatus reports whether a status is in-flight (not yet finished).
func isActiveStatus(status string) bool {
	return status == "pending" || status == "starting" || status == "running"
}

// RunningHelpers returns all HelperRun records with an active status, verifying
// each is still alive via the events file (no done event yet).
func RunningHelpers(name string) ([]*HelperRun, error) {
	runs, err := LoadHistory()
	if err != nil {
		return nil, err
	}
	var alive []*HelperRun
	for _, r := range runs {
		if !isActiveStatus(r.Status) {
			continue
		}
		if name != "" && r.Name != name {
			continue
		}
		done, _ := ReadDoneEvent(r.ID)
		if done != nil {
			status := "done"
			if done.ExitCode != 0 {
				status = "failed"
			}
			_ = UpdateHistoryStatus(r.ID, status, done.Ts)
			continue
		}
		// Walltime elapsed: treat as done even without a done event (e.g. SIGKILL).
		if r.Walltime > 0 && !r.StartedAt.IsZero() && time.Now().After(r.StartedAt.Add(r.Walltime)) {
			_ = UpdateHistoryStatus(r.ID, "done", r.StartedAt.Add(r.Walltime))
			continue
		}
		// Headless: if the pid file exists but the process is gone, treat as done.
		// Catches SIGKILL'd jobs when the server is not running to detect them.
		if r.JobID == "" {
			if _, pidErr := ReadHelperPid(r.ID); pidErr == nil && !IsHeadlessProcessAlive(r.ID) {
				_ = UpdateHistoryStatus(r.ID, "done", time.Now())
				continue
			}
		}
		alive = append(alive, r)
	}
	return alive, nil
}
