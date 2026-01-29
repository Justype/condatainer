package scheduler

import "sync"

var (
	activeScheduler Scheduler
	schedulerMu     sync.RWMutex
)

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
