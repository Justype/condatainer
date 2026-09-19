package image

import (
	"errors"
	"fmt"
	"io/fs"

	"github.com/Justype/condatainer/internal/utils"
)

// Why a lock attempt failed, so a caller can react rather than only report.
var (
	// ErrProtected marks an image whose write bit is clear. See AcquireLock.
	ErrProtected = errors.New("image is write-protected")
	// ErrInUse marks a conflicting lock held by another process.
	ErrInUse = errors.New("image is in use")
)

// lockError is a message written for the reader that still matches its sentinel
// with errors.Is. Wrapping the sentinel with %w would print its text after the
// message, repeating it.
type lockError struct {
	msg  string
	kind error
}

func (e *lockError) Error() string { return e.msg }
func (e *lockError) Unwrap() error { return e.kind }

// Lock represents a file lock on an overlay image. It must be closed to
// release the lock.
//
// An alias, not a new type: internal/libexec needs the exact same
// non-blocking file lock for its own toolchain-generation sentinel,
// and — to stay importable by internal/image/squashfs without cycling back
// through this package — returns *utils.FileLock directly rather than
// importing internal/image for it. The alias means both are the same type,
// interchangeable at any call site (internal/runtime/exec.Prepare holds a
// single slice of overlay and toolchain locks together) with no conversion.
type Lock = utils.FileLock

// AcquireLock takes a non-blocking whole-file lock on the overlay image:
// exclusive when write is true, shared otherwise. It is the kind of lock
// Apptainer holds on a mounted ext3 image, so a running container is seen.
//
// A write lock opens O_RDWR, as an fcntl write lock requires, and that is also
// the rule: an image with its write bit clear is protected, and CondaTainer never
// modifies or removes it — not even for its owner, who can unlink it through
// the directory and can restore the bit. Clearing it is how an artifact is
// pinned.
func AcquireLock(path string, write bool) (*Lock, error) {
	lock, err := utils.AcquireFileLock(path, write)
	if err != nil {
		if errors.Is(err, utils.ErrLockConflict) {
			if write {
				return nil, &lockError{utils.StylePath(path) + " is currently in use", ErrInUse}
			}
			return nil, &lockError{utils.StylePath(path) + " is open for writing by another process", ErrInUse}
		}
		return nil, openFailure(path, err, write)
	}
	return lock, nil
}

// openFailure names the actual reason the image could not be opened. Only a
// lock conflict is "in use"; a protected or missing image is neither, and
// reporting one as the other sends the reader looking for a container that is
// not running.
func openFailure(path string, err error, write bool) error {
	styled := utils.StylePath(path)
	switch {
	case errors.Is(err, fs.ErrNotExist):
		return fmt.Errorf("%s does not exist", styled)
	case errors.Is(err, fs.ErrPermission) && write:
		return &lockError{styled + " is write-protected; chmod +w to allow changes", ErrProtected}
	case errors.Is(err, fs.ErrPermission):
		return fmt.Errorf("%s is not readable", styled)
	}
	return fmt.Errorf("can't open %s: %w", styled, err)
}

// CheckAvailable reports whether the image can be locked for reading or writing
// right now. It releases immediately, so it is a pre-flight check: a caller that
// must not race another process holds the lock across its own work instead.
func CheckAvailable(path string, write bool) error {
	lock, err := AcquireLock(path, write)
	if err != nil {
		return err
	}
	return lock.Close()
}
