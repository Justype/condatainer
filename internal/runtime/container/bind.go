package container

import (
	"os"
	"path/filepath"
	"strings"

	"github.com/Justype/condatainer/internal/config"
	"github.com/Justype/condatainer/internal/scheduler"
)

// BoundExecDir and BoundExecPath are where the running executable is bound
// inside a container, so a nested call has a stable path to reach it by.
//
// Not /usr/bin. A bind there puts a mount boundary inside the directory dpkg
// unpacks into most, and dpkg then reads the copied-up /usr as something other
// than a directory and refuses every package — "trying to overwrite '/usr',
// which is also in package login".
//
// Not under the metadata directory either. A pack runs mksquashfs in a container
// with the staged metadata bound at /.cnt and passed as an archive root, so a
// bind below it is walked and the executable lands in the artifact. Its own
// directory is CondaTainer's, no package writes there, and nothing packs it.
const (
	BoundExecDir  = "/.cnt_bin"
	BoundExecPath = BoundExecDir + "/condatainer"
)

// tmpDirEnvVars are the temp directory environment variables read by name, in
// priority order: CNT_TMPDIR is condatainer's own staging root, the rest are what
// the C++ standard library's temp_directory_path() consults. The scheduler's
// node-local tmp is not here — its variable depends on the detected scheduler,
// so TmpDirBinds resolves it through scheduler.ActiveTmpDir.
var tmpDirEnvVars = []string{"CNT_TMPDIR", "TMPDIR", "TMP", "TEMP", "TEMPDIR"}

// TmpDirBinds returns every host temp directory that exists, so a tool inside the
// container can follow whichever variable it reads. Apptainer binds only /tmp by
// default, and an unbound variable names a host path that is not there — which is
// what micromamba dies on when the site profile exports one of these.
func TmpDirBinds() []string {
	dirs := make([]string, 0, len(tmpDirEnvVars)+1)
	for _, env := range tmpDirEnvVars {
		dirs = append(dirs, os.Getenv(env))
	}
	// Scheduler-assigned node-local tmp (fast SSD scratch): the variable it comes
	// from is whatever the detected scheduler names.
	dirs = append(dirs, scheduler.ActiveTmpDir())

	binds := make([]string, 0, len(dirs))
	for _, dir := range dirs {
		if dir == "" {
			continue
		}
		if _, err := os.Stat(dir); err != nil {
			continue
		}
		binds = append(binds, dir)
	}
	return binds
}

// DeduplicateBindPaths resolves, deduplicates, and filters child paths from bind directories.
// It handles formats: "/path", "/path:/container", "/path:/container:ro"
func DeduplicateBindPaths(paths []string) []string {
	// First pass: resolve all paths and deduplicate by host path
	seen := make(map[string]bool)
	resolved := make([]string, 0, len(paths))
	for _, bind := range paths {
		if bind == "" {
			continue
		}
		// Parse bind format
		parts := strings.SplitN(bind, ":", 3)
		hostPath := parts[0]

		// Resolve host path to absolute and follow symlinks
		absHostPath, err := filepath.Abs(hostPath)
		if err != nil {
			absHostPath = hostPath
		}
		if realPath, err := filepath.EvalSymlinks(absHostPath); err == nil {
			absHostPath = realPath
		}

		// Skip if host path already seen
		if seen[absHostPath] {
			continue
		}
		seen[absHostPath] = true

		// Reconstruct the bind string with resolved host path
		var resolvedBind string
		switch len(parts) {
		case 3:
			resolvedBind = absHostPath + ":" + parts[1] + ":" + parts[2]
		case 2:
			resolvedBind = absHostPath + ":" + parts[1]
		default:
			resolvedBind = absHostPath
		}
		resolved = append(resolved, resolvedBind)
	}

	// Second pass: filter out child paths covered by parent paths.
	// Binds with an explicit container path (e.g. /host:/container) are deliberate
	// remaps and must never be filtered — the parent bind does not provide that mapping.
	filtered := make([]string, 0, len(resolved))
	for _, bind := range resolved {
		parts := strings.SplitN(bind, ":", 3)
		hostPath := parts[0]

		// Explicit remap — always keep
		if len(parts) >= 2 {
			filtered = append(filtered, bind)
			continue
		}

		isChild := false
		for _, otherBind := range resolved {
			otherParts := strings.SplitN(otherBind, ":", 3)
			otherHostPath := otherParts[0]
			if hostPath == otherHostPath {
				continue
			}
			rel, err := filepath.Rel(otherHostPath, hostPath)
			if err == nil && rel != "." && !strings.HasPrefix(rel, "..") {
				isChild = true
				break
			}
		}
		if !isChild {
			filtered = append(filtered, bind)
		}
	}

	return filtered
}

// BindPaths collects all directories suitable for --bind flags.
// Deduplication and child path filtering is handled by DeduplicateBindPaths().
func BindPaths(paths ...string) []string {
	bindPaths := make([]string, 0, len(paths)+10)

	// Add provided paths
	for _, path := range paths {
		if path == "" {
			continue
		}
		if _, err := os.Stat(path); err != nil {
			continue
		}
		bindPaths = append(bindPaths, path)
	}

	// Add current working directory
	if cwd, err := os.Getwd(); err == nil && cwd != "" {
		bindPaths = append(bindPaths, cwd)
	}

	// Add $SCRATCH if set
	if scratch := os.Getenv("SCRATCH"); scratch != "" {
		bindPaths = append(bindPaths, scratch)
	}

	// Add every temp directory the host advertises.
	bindPaths = append(bindPaths, TmpDirBinds()...)

	// Collect all base directories
	baseDirs := []string{}

	// Extra root dir (CNT_EXTRA_ROOT, group/lab layer)
	if dir := config.GetExtraRootDir(); dir != "" {
		baseDirs = append(baseDirs, dir)
	}

	// Root, Scratch, User directories
	if dir := config.GetRootDir(); dir != "" {
		baseDirs = append(baseDirs, dir)
	}
	if dir := config.GetScratchDataDir(); dir != "" {
		baseDirs = append(baseDirs, dir)
	}
	if dir := config.GetUserDataDir(); dir != "" {
		baseDirs = append(baseDirs, dir)
	}

	// Add base directories
	for _, dir := range baseDirs {
		if dir == "" {
			continue
		}
		if _, err := os.Stat(dir); err != nil {
			continue
		}
		bindPaths = append(bindPaths, dir)
	}

	// Bind condatainer executable for nested calls.
	// Always bind the actual exe regardless of install location.
	if execPath, err := os.Executable(); err == nil && execPath != "" {
		bindPaths = append(bindPaths, execPath+":"+BoundExecPath+":ro")
	}

	return bindPaths
}
