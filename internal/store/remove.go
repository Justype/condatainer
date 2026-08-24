package store

import (
	"errors"
	"fmt"
	"os"

	"github.com/Justype/condatainer/internal/artifactcache"
	"github.com/Justype/condatainer/internal/image"
)

// resolveIdentityFor locates the artifact to remove. Injected for tests, which
// have no real image for key regeneration to read.
var resolveIdentityFor = ResolveIdentity

// ErrNotStored reports an artifact that resolved to a flat install rather than
// a store entry.
var ErrNotStored = errors.New("artifact is not a store entry")

// Remove deletes one store entry, lock and protection permitting.
//
// Only the store layout. A flat artifact answers to a bare name and is what
// `remove` deals with; deleting one here would silently take a name out of
// service through a command that reads as identity housekeeping.
//
// The exclusive inode lock is the whole safety story, and its failures stay
// distinct because they send the reader to different places: ErrInUse means a
// container is running, ErrProtected means someone pinned this identity by
// clearing its write bit, and a missing file means a concurrent remover won.
func Remove(name string, q IdentityQuery, dirs []string) (Candidate, error) {
	candidate, _, err := resolveIdentityFor(name, q, dirs)
	if err != nil {
		return Candidate{}, err
	}
	if candidate.Layout != LayoutStored {
		return Candidate{}, fmt.Errorf("%w: %s is installed flat; use `condatainer remove`",
			ErrNotStored, candidate.Path)
	}
	if err := removeLocked(candidate.Path); err != nil {
		return Candidate{}, err
	}
	return candidate, nil
}

// removeLocked unlinks one artifact while holding its exclusive inode lock, and
// forgets it from the artifact cache.
//
// The lock is released before the cache is updated but after the unlink, which
// is the only order that cannot leave a cached record for a live file: a reader
// that repopulates the cache in between is reading a path that no longer exists
// and caches nothing.
func removeLocked(path string) error {
	lock, err := image.AcquireLock(path, true)
	if err != nil {
		return err
	}
	removeErr := os.Remove(path)
	lock.Close() //nolint:errcheck
	if removeErr != nil && !os.IsNotExist(removeErr) {
		return removeErr
	}
	artifactcache.Forget(path)
	return nil
}
