package server

import (
	"encoding/json"
	"os"
	"path/filepath"
	"syscall"

	"github.com/Justype/condatainer/internal/utils"
)

// jsonStore persists a single JSON value of type T to a file, guarded by an
// exclusive flock on a sibling ".lock" file so concurrent HTTP requests
// can't interleave reads and writes. Writes are atomic (tmp file then
// rename). Used for small pieces of dashboard state — bookmarks, pinned
// helpers, and similar — that live under the user's state dir.
//
// Reads (Load) are not locked: writes are atomic via rename, so a reader
// always sees either the old or the new complete file, never a partial one.
type jsonStore[T any] struct {
	path string
}

// newJSONStore returns a store backed by the file at path.
func newJSONStore[T any](path string) *jsonStore[T] {
	return &jsonStore[T]{path: path}
}

// Load reads the current value. Returns the zero value of T, not an error,
// if the file does not exist yet.
func (s *jsonStore[T]) Load() (T, error) {
	var zero T
	data, err := os.ReadFile(s.path)
	if err != nil {
		if os.IsNotExist(err) {
			return zero, nil
		}
		return zero, err
	}
	var v T
	if err := json.Unmarshal(data, &v); err != nil {
		return zero, err
	}
	return v, nil
}

// Update loads the current value, applies fn (which mutates it in place),
// and atomically writes the result back — all under an exclusive lock so
// concurrent updates can't clobber each other. Returns the updated value.
// If fn returns an error, the file is left untouched.
func (s *jsonStore[T]) Update(fn func(*T) error) (T, error) {
	var result T
	err := s.withLock(func() error {
		v, err := s.Load()
		if err != nil {
			return err
		}
		if err := fn(&v); err != nil {
			return err
		}
		if err := s.save(v); err != nil {
			return err
		}
		result = v
		return nil
	})
	return result, err
}

// save atomically writes v to s.path (tmp file then rename).
func (s *jsonStore[T]) save(v T) error {
	if err := os.MkdirAll(filepath.Dir(s.path), utils.PermDir); err != nil {
		return err
	}
	data, err := json.MarshalIndent(v, "", "  ")
	if err != nil {
		return err
	}
	tmp := s.path + ".tmp"
	f, err := utils.CreateFileWritable(tmp)
	if err != nil {
		return err
	}
	if _, err := f.Write(data); err != nil {
		f.Close()
		os.Remove(tmp) //nolint:errcheck
		return err
	}
	if err := f.Close(); err != nil {
		os.Remove(tmp) //nolint:errcheck
		return err
	}
	return os.Rename(tmp, s.path)
}

// withLock acquires an exclusive flock on s.path+".lock" for the duration
// of fn.
func (s *jsonStore[T]) withLock(fn func() error) error {
	if err := os.MkdirAll(filepath.Dir(s.path), utils.PermDir); err != nil {
		return err
	}
	lf, err := os.OpenFile(s.path+".lock", os.O_CREATE|os.O_WRONLY, utils.PermFile)
	if err != nil {
		return err
	}
	defer lf.Close()
	if err := syscall.Flock(int(lf.Fd()), syscall.LOCK_EX); err != nil {
		return err
	}
	defer syscall.Flock(int(lf.Fd()), syscall.LOCK_UN) //nolint:errcheck
	return fn()
}
