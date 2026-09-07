package utils

import (
	"errors"
	"fmt"
	"os"
	"syscall"
)

// ErrFlockConflict marks a non-blocking flock call that lost to another
// holder: the file opened fine, but a conflicting lock is already held.
var ErrFlockConflict = errors.New("file is locked by another process")

// FlockHandle is an open file held under a non-blocking flock. It must be
// closed to release the lock.
type FlockHandle struct {
	file *os.File
}

// Close releases the lock by closing the file.
func (h *FlockHandle) Close() error {
	if h.file == nil {
		return nil
	}
	err := h.file.Close()
	h.file = nil
	return err
}

// AcquireFlock opens path and takes a non-blocking flock on it: exclusive
// (LOCK_EX) when write is true, shared (LOCK_SH) otherwise. write also picks
// the open mode (O_RDWR/O_RDONLY) — deliberately, not merely what flock
// needs: some callers use the open mode itself as a permission check, ahead
// of and independent from the flock call, and that has to stay true for
// whatever this is built into.
//
// A failure to open returns the raw *fs.PathError from os.OpenFile
// unwrapped, so a caller can classify it with errors.Is against
// fs.ErrNotExist/fs.ErrPermission. A failure to flock — the file opened, but
// a conflicting lock is already held — wraps ErrFlockConflict, so a caller
// can tell the two failure stages apart.
func AcquireFlock(path string, write bool) (*FlockHandle, error) {
	flag, op := os.O_RDONLY, syscall.LOCK_SH
	if write {
		flag, op = os.O_RDWR, syscall.LOCK_EX
	}
	f, err := os.OpenFile(path, flag, 0)
	if err != nil {
		return nil, err
	}
	if err := syscall.Flock(int(f.Fd()), op|syscall.LOCK_NB); err != nil {
		f.Close()
		return nil, fmt.Errorf("%w: %w", ErrFlockConflict, err)
	}
	return &FlockHandle{file: f}, nil
}
