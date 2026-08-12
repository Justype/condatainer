package image

import (
	"errors"
	"fmt"
	"io/fs"
	"os"
	"syscall"

	"github.com/Justype/condatainer/internal/utils"
)

// Why a lock attempt failed, so a caller can react rather than only report.
var (
	// ErrProtected marks an image whose write bit is clear. See AcquireLock.
	ErrProtected = errors.New("image is write-protected")
	// ErrInUse marks a conflicting flock held by another process.
	ErrInUse = errors.New("image is in use")
)

// Lock represents a file lock on an overlay image.
// It must be closed to release the lock.
type Lock struct {
	file  *os.File
	path  string
	write bool
}

// Close releases the lock by closing the file.
func (l *Lock) Close() error {
	if l.file == nil {
		return nil
	}
	err := l.file.Close()
	l.file = nil
	return err
}

// AcquireLock takes a non-blocking flock on the overlay image: exclusive
// (LOCK_EX) when write is true, shared (LOCK_SH) otherwise.
//
// A write lock opens O_RDWR, which is deliberate and not merely what flock
// needs: an image with its write bit clear is protected, and CondaTainer never
// modifies or removes it — not even for its owner, who can unlink it through
// the directory and can restore the bit. Clearing it is how an artifact is
// pinned.
func AcquireLock(path string, write bool) (*Lock, error) {
	flag, op := os.O_RDONLY, syscall.LOCK_SH
	if write {
		flag, op = os.O_RDWR, syscall.LOCK_EX
	}
	f, err := os.OpenFile(path, flag, 0)
	if err != nil {
		return nil, openFailure(path, err, write)
	}
	if err := syscall.Flock(int(f.Fd()), op|syscall.LOCK_NB); err != nil {
		f.Close()
		if write {
			return nil, fmt.Errorf("%s is currently in use: %w", utils.StylePath(path), ErrInUse)
		}
		return nil, fmt.Errorf("%s is currently being written: %w", utils.StylePath(path), ErrInUse)
	}
	return &Lock{file: f, path: path, write: write}, nil
}

// openFailure names the actual reason the image could not be opened. Only a
// flock conflict is "in use"; a protected or missing image is neither, and
// reporting one as the other sends the reader looking for a container that is
// not running.
func openFailure(path string, err error, write bool) error {
	styled := utils.StylePath(path)
	switch {
	case errors.Is(err, fs.ErrNotExist):
		return fmt.Errorf("%s does not exist", styled)
	case errors.Is(err, fs.ErrPermission) && write:
		return fmt.Errorf("%s is write-protected; chmod +w to allow changes: %w", styled, ErrProtected)
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
