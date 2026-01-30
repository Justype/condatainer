package exec

import (
	"os"
	"path/filepath"

	"github.com/Justype/condatainer/internal/config"
)

// BindPaths collects all directories suitable for --bind flags.
// Deduplication and child path filtering is handled by apptainer.DeduplicateBindPaths().
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
