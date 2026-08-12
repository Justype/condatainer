// Package producer coordinates commands that create or replace one image.
//
// Its lock is separate from image.AcquireLock: a producer lock prevents two
// builds or pulls from targeting the same pathname, while the inode flock
// prevents replacing an image that a running container is using.
package producer

import (
	"context"
	"encoding/json"
	"fmt"
	"os"
	"strconv"
	"strings"
	"syscall"
	"time"

	"github.com/Justype/condatainer/internal/scheduler"
	"github.com/Justype/condatainer/internal/utils"
)

// PreparedSuffix marks a producer's in-progress output.
const PreparedSuffix = ".part"

// Info is the JSON metadata stored in a producer lock file.
type Info struct {
	Runner    string `json:"runner"`
	JobID     string `json:"job_id"`
	Node      string `json:"node"`
	PID       int    `json:"pid"`
	CreatedAt string `json:"created_at"`
}

// Guard is a local producer lock held for one target. Release removes only the
// lock this guard created.
type Guard struct {
	path string
	info Info
}

// AcquireLocal claims target for the current process. A definitely stale lock
// and its owner-derived partial output are removed before one bounded retry.
func AcquireLocal(target string) (*Guard, error) {
	info := Info{
		Runner:    "local",
		Node:      ShortHostname(),
		PID:       os.Getpid(),
		CreatedAt: time.Now().Format(time.RFC3339),
	}
	path := Path(target)
	for attempt := 0; attempt < 2; attempt++ {
		if err := Acquire(path, info); err == nil {
			return &Guard{path: path, info: info}, nil
		} else if !os.IsExist(err) {
			return nil, fmt.Errorf("cannot create producer lock %s: %w", path, err)
		}

		existing, err := Read(path)
		if err != nil {
			return nil, fmt.Errorf("cannot read producer lock %s: %w", path, err)
		}
		stale, _, checkErr := IsStale(existing)
		if !stale {
			if checkErr != nil {
				return nil, fmt.Errorf("producer lock found at %s: %w", path, checkErr)
			}
			return nil, fmt.Errorf("another build or pull is producing %s (lock: %s)", target, path)
		}
		if err := os.Remove(path); err != nil && !os.IsNotExist(err) {
			return nil, fmt.Errorf("cannot remove stale producer lock %s: %w", path, err)
		}
		os.Remove(PreparedPath(target, existing)) //nolint:errcheck
	}
	return nil, fmt.Errorf("another build or pull claimed %s", target)
}

// Info returns the metadata recorded by this guard.
func (g *Guard) Info() Info { return g.info }

// Release gives up the producer lock.
func (g *Guard) Release() error {
	if g == nil || g.path == "" {
		return nil
	}
	path := g.path
	g.path = ""
	err := os.Remove(path)
	if os.IsNotExist(err) {
		return nil
	}
	return err
}

// Path returns the producer lock path for target.
func Path(target string) string { return target + ".lock" }

// PreparedPath returns a producer-private output path beside target. Because it
// is derived from Info, stale-lock cleanup can find the abandoned output.
func PreparedPath(target string, info Info) string {
	return target + "." + Tag(info) + PreparedSuffix
}

// Tag renders the filesystem-safe producer identity shared by prepared outputs
// and private build workspaces.
func Tag(info Info) string {
	owner := info.Runner
	if owner == "" {
		owner = "local"
	}
	tag := owner
	if info.JobID != "" {
		tag += "-" + info.JobID
	} else {
		tag += "-" + info.Node + "-" + strconv.Itoa(info.PID)
	}
	return sanitizeTag(tag)
}

func sanitizeTag(tag string) string {
	return strings.Map(func(r rune) rune {
		switch {
		case r >= 'a' && r <= 'z', r >= 'A' && r <= 'Z', r >= '0' && r <= '9', r == '-', r == '_':
			return r
		default:
			return '-'
		}
	}, tag)
}

// LocalInfo identifies the current process for provisional workspace paths.
// Build lock acquisition produces the same tag unless a scheduler lock is
// adopted, in which case the workspace is re-sited to the scheduler job tag.
func LocalInfo() Info {
	return Info{Runner: "local", Node: ShortHostname(), PID: os.Getpid()}
}

// ShortHostname returns the unqualified hostname used in local lock metadata.
func ShortHostname() string {
	h, _ := os.Hostname()
	if idx := strings.Index(h, "."); idx > 0 {
		return h[:idx]
	}
	return h
}

// Acquire atomically creates path and records info. Callers check os.IsExist to
// distinguish an active producer from other filesystem errors.
func Acquire(path string, info Info) error {
	data, err := json.Marshal(info)
	if err != nil {
		return fmt.Errorf("failed to marshal producer lock: %w", err)
	}
	f, err := os.OpenFile(path, os.O_CREATE|os.O_EXCL|os.O_WRONLY, utils.PermFile)
	if err != nil {
		return err
	}
	defer f.Close()
	if _, err = f.Write(data); err != nil {
		return err
	}
	utils.ShareWithParentGroup(path)
	return nil
}

// Overwrite updates a lock already held by the caller.
func Overwrite(path string, info Info) error {
	data, err := json.Marshal(info)
	if err != nil {
		return fmt.Errorf("failed to marshal producer lock: %w", err)
	}
	return os.WriteFile(path, data, utils.PermFile)
}

// Read parses a producer lock. An empty legacy lock returns zero Info and is
// consequently stale.
func Read(path string) (Info, error) {
	data, err := os.ReadFile(path)
	if err != nil {
		return Info{}, err
	}
	if len(data) == 0 {
		return Info{}, nil
	}
	var info Info
	if err := json.Unmarshal(data, &info); err != nil {
		return Info{}, fmt.Errorf("corrupt lock file: %w", err)
	}
	return info, nil
}

// IsStale reports whether the recorded producer is definitely gone. An
// unverifiable remote owner is conservatively treated as alive and returned
// with an explanatory error.
func IsStale(info Info) (stale bool, status scheduler.JobStatus, err error) {
	if info.Runner == "" {
		return true, scheduler.JobStatusUnknown, nil
	}
	if info.Runner != "local" {
		if info.JobID == "" {
			return true, scheduler.JobStatusUnknown, nil
		}
		sched := scheduler.ActiveScheduler()
		if sched == nil {
			return false, scheduler.JobStatusUnknown, fmt.Errorf("scheduler unavailable, cannot verify job %s", info.JobID)
		}
		st, err := sched.GetJobStatus(context.Background(), info.JobID)
		if err != nil {
			return false, scheduler.JobStatusUnknown, fmt.Errorf("cannot check job %s: %w", info.JobID, err)
		}
		if st == scheduler.JobStatusUnknown {
			return false, st, fmt.Errorf("cannot determine status of job %s", info.JobID)
		}
		return !st.IsAlive(), st, nil
	}

	host := ShortHostname()
	if info.Node != host {
		return false, scheduler.JobStatusUnknown, fmt.Errorf("lock held by node %q (current: %q); cannot verify remotely", info.Node, host)
	}
	proc, err := os.FindProcess(info.PID)
	if err != nil {
		return true, scheduler.JobStatusUnknown, nil
	}
	if err := proc.Signal(syscall.Signal(0)); err != nil {
		if err == syscall.EPERM {
			return false, scheduler.JobStatusRunning, nil
		}
		return true, scheduler.JobStatusUnknown, nil
	}
	return false, scheduler.JobStatusRunning, nil
}
