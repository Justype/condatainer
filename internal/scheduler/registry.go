package scheduler

import "sync"

var (
	activeScheduler Scheduler
	schedulerMu     sync.RWMutex

	debugMode bool

	slurmMem = true
)

// SetDebugMode enables or disables debug output for scheduler operations.
// Call from cmd/root.go after loading config.
func SetDebugMode(enabled bool) { debugMode = enabled }

// SetSlurmMem controls whether generated SLURM scripts carry --mem/--mem-per-cpu
// (default true). Call from cmd/root.go after loading config.
func SetSlurmMem(enabled bool) { slurmMem = enabled }

// SetActiveScheduler configures the scheduler instance that the application should use.
// Passing nil clears any previously configured scheduler.
func SetActiveScheduler(s Scheduler) {
	schedulerMu.Lock()
	defer schedulerMu.Unlock()
	activeScheduler = s
}

// ActiveScheduler returns the currently configured scheduler instance (may be nil).
func ActiveScheduler() Scheduler {
	schedulerMu.RLock()
	defer schedulerMu.RUnlock()
	return activeScheduler
}

// ClearActiveScheduler resets the active scheduler reference.
func ClearActiveScheduler() {
	SetActiveScheduler(nil)
}
