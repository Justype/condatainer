package helper

import (
	"errors"
	"fmt"
	"path/filepath"
	"syscall"
	"time"

	"github.com/Justype/condatainer/internal/image/producer"
	"github.com/Justype/condatainer/internal/utils"
)

// ErrOtherHost marks a headless helper that runs on a different host than the
// caller, so its process cannot be signalled from here.
var ErrOtherHost = errors.New("helper runs on another host")

// Liveness is what a headless helper's lock file says about its wrapper.
type Liveness int

const (
	// LivenessUnknown: no usable record, so nothing is claimed either way.
	LivenessUnknown Liveness = iota
	// LivenessAlive: the wrapper's holder process still has the lock.
	LivenessAlive
	// LivenessGone: the lock was recorded and has since been released.
	LivenessGone
)

// LockFilePath returns the lock file a headless helper's wrapper holds.
func LockFilePath(id string) string {
	return filepath.Join(StateDir(id), "lock")
}

// HoldLock creates path, takes an exclusive lock on it and records info. The
// lock lasts until the returned handle is closed or the process exits.
func HoldLock(path string, info producer.Info) (*utils.FileLock, error) {
	if err := producer.Acquire(path, producer.Info{}); err != nil && !errors.Is(err, syscall.EEXIST) {
		return nil, err
	}
	var lock *utils.FileLock
	var err error
	// A liveness probe holds a shared lock for an instant; wait it out.
	for i := 0; i < 20; i++ {
		if lock, err = utils.AcquireFileLock(path, true); !errors.Is(err, utils.ErrLockConflict) {
			break
		}
		time.Sleep(50 * time.Millisecond)
	}
	if err != nil {
		return nil, err
	}
	if err := producer.Overwrite(path, info); err != nil {
		lock.Close()
		return nil, err
	}
	return lock, nil
}

// HeadlessLiveness reports whether a headless helper's wrapper is still
// running, judged by its lock file, which is visible from every host. The
// recorded owner is returned when there is one.
func HeadlessLiveness(id string) (Liveness, producer.Info) {
	path := LockFilePath(id)
	probe, err := utils.AcquireFileLock(path, false)
	if errors.Is(err, utils.ErrLockConflict) {
		info, _ := producer.Read(path)
		return LivenessAlive, info
	}
	if err != nil {
		return LivenessUnknown, producer.Info{}
	}
	probe.Close()
	// Free and never written: the holder has not started yet.
	if info, err := producer.Read(path); err == nil && info.PID > 0 {
		return LivenessGone, info
	}
	return LivenessUnknown, producer.Info{}
}

// KillHeadlessProcess sends SIGTERM to the process group of the headless
// helper's wrapper, which is its session leader. The wrapper's TERM trap
// records done and exits. It refuses a helper owned by another host.
func KillHeadlessProcess(id string) error {
	info, err := producer.Read(LockFilePath(id))
	if err != nil || info.PID <= 0 {
		return fmt.Errorf("helper %s has no lock record", id)
	}
	if host := producer.ShortHostname(); info.Node != host {
		return fmt.Errorf("%w: %s runs on %s, this is %s; stop it from there", ErrOtherHost, id, info.Node, host)
	}
	return syscall.Kill(-info.PID, syscall.SIGTERM)
}
