package container

import (
	"os"
	"path/filepath"
	"strings"

	"github.com/Justype/condatainer/internal/config"
)

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

	// Second pass: filter out child paths covered by parent paths
	filtered := make([]string, 0, len(resolved))
	for _, bind := range resolved {
		parts := strings.SplitN(bind, ":", 3)
		hostPath := parts[0]

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

	// Helper to check if directory is writable
	isWritable := func(dir string) bool {
		testFile := filepath.Join(dir, ".write-test")
		f, err := os.Create(testFile)
		if err != nil {
			return false
		}
		f.Close()
		os.Remove(testFile)
		return true
	}

	// Collect all base directories
	baseDirs := []string{}

	// Extra base directories
	baseDirs = append(baseDirs, config.GetExtraBaseDirs()...)

	// Portable, Scratch, User directories
	if dir := config.GetPortableDataDir(); dir != "" {
		baseDirs = append(baseDirs, dir)
	}
	if dir := config.GetScratchDataDir(); dir != "" {
		baseDirs = append(baseDirs, dir)
	}
	if dir := config.GetUserDataDir(); dir != "" {
		baseDirs = append(baseDirs, dir)
	}

	// Add base directories with :ro flag if not writable
	for _, dir := range baseDirs {
		if dir == "" {
			continue
		}
		if _, err := os.Stat(dir); err != nil {
			continue
		}
		if !isWritable(dir) {
			bindPaths = append(bindPaths, dir+":"+dir+":ro")
		} else {
			bindPaths = append(bindPaths, dir)
		}
	}

	// Bind condatainer executable for nested calls (non-portable only)
	if !config.IsPortable() {
		if execPath, err := os.Executable(); err == nil && execPath != "" {
			bindPaths = append(bindPaths, execPath+":/usr/bin/condatainer:ro")
		}
	}

	return bindPaths
}
