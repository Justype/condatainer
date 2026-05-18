package overlay

import (
	"fmt"
	"os"
	"syscall"

	"github.com/Justype/condatainer/internal/utils"
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

// AcquireLock attempts to acquire a lock on the overlay file.
// If write is true, it attempts an exclusive lock (LOCK_EX).
// If write is false, it attempts a shared lock (LOCK_SH).
// Both are non-blocking (LOCK_NB).
func AcquireLock(path string, write bool) (*Lock, error) {
	if write {
		f, err := os.OpenFile(path, os.O_RDWR, 0)
		if err != nil {
			return nil, fmt.Errorf("Can't open %s for writing, currently in use", utils.StylePath(path))
		}
		err = syscall.Flock(int(f.Fd()), syscall.LOCK_EX|syscall.LOCK_NB)
		if err != nil {
			f.Close()
			return nil, fmt.Errorf("Can't open %s for writing, currently in use", utils.StylePath(path))
		}
		return &Lock{file: f, path: path, write: true}, nil
	}
	f, err := os.Open(path)
	if err != nil {
		return nil, fmt.Errorf("Can't open %s for reading, currently in use for writing", utils.StylePath(path))
	}
	err = syscall.Flock(int(f.Fd()), syscall.LOCK_SH|syscall.LOCK_NB)
	if err != nil {
		f.Close()
		return nil, fmt.Errorf("Can't open %s for reading, currently in use for writing", utils.StylePath(path))
	}
	return &Lock{file: f, path: path, write: false}, nil
}

// CheckAvailable checks if the overlay image is available for reading or writing.
func CheckAvailable(path string, write bool) error {
	lock, err := AcquireLock(path, write)
	if err != nil {
		return err
	}
	return lock.Close()
}
