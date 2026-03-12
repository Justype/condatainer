package scheduler

import (
	"context"
	"os"
	"os/exec"
	"path/filepath"
	"time"
)

// siblingBin resolves a companion binary name relative to an already-known
// absolute binary path. When base is absolute, it checks the same directory
// (two fast os.Stat calls instead of a full PATH scan). Falls back to
// exec.LookPath when base is relative or not absolute.
func siblingBin(base, name string) string {
	if filepath.IsAbs(base) {
		p := filepath.Join(filepath.Dir(base), name)
		if _, err := os.Stat(p); err == nil {
			return p
		}
		return ""
	}
	p, _ := exec.LookPath(name)
	return p
}

// DefaultCommandTimeout is the maximum time to wait for a scheduler command to respond.
// Set this before calling any scheduler methods (e.g., from cmd/root.go after loading config).
// A value of 0 disables the timeout (command runs until completion).
var DefaultCommandTimeout = 5 * time.Second

// runCommand executes a scheduler CLI command with DefaultCommandTimeout.
// Returns TimeoutError if the command does not respond within the timeout.
// If DefaultCommandTimeout is 0, the command runs without a time limit.
func runCommand(schedulerName, operation, bin string, args ...string) ([]byte, error) {
	if DefaultCommandTimeout == 0 {
		out, err := exec.Command(bin, args...).CombinedOutput()
		return out, err
	}
	ctx, cancel := context.WithTimeout(context.Background(), DefaultCommandTimeout)
	defer cancel()
	out, err := exec.CommandContext(ctx, bin, args...).CombinedOutput()
	if ctx.Err() == context.DeadlineExceeded {
		return nil, &TimeoutError{Scheduler: schedulerName, Operation: operation, Timeout: DefaultCommandTimeout}
	}
	return out, err
}
