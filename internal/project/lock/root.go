package lock

import (
	"errors"
	"fmt"
	"os"
	"path/filepath"
)

// ErrNoProject reports that no cnt-lock/ was found at or above a directory.
var ErrNoProject = errors.New("no project found")

// RootAt returns dir when it directly contains cnt-lock/, and ErrNoProject
// otherwise.
//
// Nothing searches a parent, so a subdirectory of a project is not a project:
// the root is always a directory the caller is standing in. Only cnt-lock/ marks
// one — a Git root is never a fallback, since a repository may hold several
// projects or none.
func RootAt(dir string) (string, error) {
	absolute, err := filepath.Abs(dir)
	if err != nil {
		return "", err
	}
	if info, err := os.Stat(filepath.Join(absolute, DirName)); err == nil && info.IsDir() {
		return absolute, nil
	}
	return "", fmt.Errorf("%w at %s", ErrNoProject, absolute)
}

// RootFor resolves the project root a command should act on: the explicit
// --project directory when given, otherwise the current one. An explicit
// directory is taken as named, never re-derived.
func RootFor(explicit, start string) (string, error) {
	if explicit == "" {
		return RootAt(start)
	}
	dir, err := filepath.Abs(explicit)
	if err != nil {
		return "", err
	}
	info, err := os.Stat(dir)
	if err != nil {
		return "", err
	}
	if !info.IsDir() {
		return "", fmt.Errorf("project path is not a directory: %s", dir)
	}
	return dir, nil
}

// Dir is the lock directory inside a project root.
func Dir(root string) string { return filepath.Join(root, DirName) }

// FilePath is the lock document inside a project root.
func FilePath(root string) string { return filepath.Join(Dir(root), FileName) }

// ArtifactsPath is the vendored artifact directory inside a project root.
func ArtifactsPath(root string) string { return filepath.Join(Dir(root), ArtifactsDir) }

// Load reads and validates a project's lock. A project with a cnt-lock/ but no
// lock.json yet is an empty lock, not an error: `project lock` creates the
// directory before it has anything to record.
func Load(root string) (*Lock, error) {
	data, err := os.ReadFile(FilePath(root))
	if errors.Is(err, os.ErrNotExist) {
		return New(), nil
	}
	if err != nil {
		return nil, err
	}
	return Unmarshal(data)
}
