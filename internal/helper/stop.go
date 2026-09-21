package helper

import (
	"context"
	"errors"
	"fmt"
	"time"

	"github.com/Justype/condatainer/internal/scheduler"
)

// ErrNoScheduler is returned by StopRun for a scheduled run when no scheduler
// is available to cancel it.
var ErrNoScheduler = errors.New("no scheduler available")

// StopRun asks a helper run to stop: a scheduled job is cancelled through the
// scheduler, a headless one gets SIGTERM to its process group. A failure is
// returned, but the caller still closes the run out, since the process may
// already be gone.
func StopRun(ctx context.Context, r *HelperRun) error {
	if r.Headless() {
		return KillHeadlessProcess(r.ID)
	}
	sched := scheduler.ActiveScheduler()
	if sched == nil {
		return fmt.Errorf("%w for job %s", ErrNoScheduler, r.JobID)
	}
	ctx, cancel := context.WithTimeout(ctx, 30*time.Second)
	defer cancel()
	return sched.CancelJob(ctx, r.JobID)
}
