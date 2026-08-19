// Package project holds what a project is, above the lock that records it.
package project

import (
	"fmt"
	"path/filepath"
)

// WorkDir is the directory a project's job must run in, given whatever working
// directory its script declared. An empty declared value means none was given.
//
// The answer is always the project root, and a declared directory that is not
// the root is refused rather than honored.
//
// Inside a project a relative path means one thing: relative to the project
// root. That is what the scanner keys a `#DEP:` path on and where restore
// materializes it. The runtime resolves the script's *own* relative paths —
// input data, outputs — against the process working directory, so a job that
// runs anywhere else silently splits the two: CondaTainer's overlays resolve
// against the root while the script's own paths resolve somewhere else. One
// consistent anchor is the only arrangement where both are right.
//
// Refusing rather than overriding: a declared `--chdir` is something the author
// wrote down on purpose, and quietly relocating their job would break the
// relative paths they wrote it for. Two intentions are in conflict and only the
// author can resolve it.
//
// The default matters as much as the check. Nothing sets a working directory
// today, so a job takes the scheduler's — the submission directory under SLURM
// and LSF, but `$HOME` under PBS. A project therefore states the root instead
// of inheriting an answer that varies by scheduler.
func WorkDir(root, declared string) (string, error) {
	root, err := filepath.Abs(root)
	if err != nil {
		return "", err
	}
	if declared == "" {
		return root, nil
	}
	absolute, err := filepath.Abs(declared)
	if err != nil {
		return "", err
	}
	if !samePath(absolute, root) {
		return "", fmt.Errorf("the script runs in %s but this project is rooted at %s; a project's relative paths resolve against the root",
			absolute, root)
	}
	return root, nil
}

// samePath compares two absolute paths, resolving symlinks when both exist.
// A path that cannot be resolved is compared as written, which is the same
// answer for everything except a link, and a link that is not there yet cannot
// be followed anyway.
func samePath(a, b string) bool {
	if a == b {
		return true
	}
	resolvedA, errA := filepath.EvalSymlinks(a)
	resolvedB, errB := filepath.EvalSymlinks(b)
	return errA == nil && errB == nil && resolvedA == resolvedB
}
