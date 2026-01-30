package exec

import (
	"os"
	"path/filepath"
	"strings"

	"github.com/Justype/condatainer/internal/config"
)

// BindPaths returns unique absolute directories suitable for --bind flags while omitting nested duplicates.
func BindPaths(paths ...string) []string {
	absPaths := make([]string, 0, len(paths)+3)
	for _, path := range paths {
		if path == "" {
			continue
		}
		// Check if path exists before adding
		if _, err := os.Stat(path); err != nil {
			continue
		}
		if resolved := resolvePath(path); resolved != "" {
			absPaths = append(absPaths, resolved)
		}
	}

	// Do not autobind the home dir (already done by apptainer; apptainer has --no-home flag)
	if cwd, err := os.Getwd(); err == nil && cwd != "" {
		if resolved := resolvePath(cwd); resolved != "" {
			absPaths = append(absPaths, resolved)
		}
	}
	if scratch := os.Getenv("SCRATCH"); scratch != "" {
		if resolved := resolvePath(scratch); resolved != "" {
			absPaths = append(absPaths, resolved)
		}
	}

	// Auto-bind all possible base directories (from search paths)
	// Check if writable and add :ro flag if read-only
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

	// Get all base directories from search paths
	baseDirs := make([]string, 0)

	// Extra base directories
	for _, dir := range config.GetExtraBaseDirs() {
		if dir != "" {
			baseDirs = append(baseDirs, dir)
		}
	}

	// Portable directory
	if portableDir := config.GetPortableDataDir(); portableDir != "" {
		baseDirs = append(baseDirs, portableDir)
	}

	// Scratch directory
	if scratchDir := config.GetScratchDataDir(); scratchDir != "" {
		baseDirs = append(baseDirs, scratchDir)
	}

	// User directory
	if userDir := config.GetUserDataDir(); userDir != "" {
		baseDirs = append(baseDirs, userDir)
	}

	// Add base directories with appropriate flags
	for _, dir := range baseDirs {
		if resolved := resolvePath(dir); resolved != "" {
			// Check if directory exists and is readable
			if _, err := os.Stat(resolved); err == nil {
				// Add with :ro flag if not writable (format: path:path:ro)
				if !isWritable(resolved) {
					resolved = resolved + ":" + resolved + ":ro"
				}
				absPaths = append(absPaths, resolved)
			}
		}
	}

	unique := make([]string, 0, len(absPaths))
	seen := map[string]struct{}{}
	for _, candidate := range absPaths {
		if _, ok := seen[candidate]; ok {
			continue
		}
		seen[candidate] = struct{}{}
		unique = append(unique, candidate)
	}

	filtered := make([]string, 0, len(unique))
	for _, candidate := range unique {
		skip := false
		// Extract the host path (first component before any colons)
		candidatePath := extractHostPath(candidate)
		for _, parent := range unique {
			if parent == candidate {
				continue
			}
			// Extract the host path from parent too
			parentPath := extractHostPath(parent)
			if isChildPath(candidatePath, parentPath) {
				skip = true
				break
			}
		}
		if !skip {
			filtered = append(filtered, candidate)
		}
	}

	// also bind the condatainer executable to /usr/bin/condatainer to help with nested calls
	// If is portable, the bin dir is already bound via PATH modification
	if !config.IsPortable() {
		execPath, err := os.Executable()
		if err == nil && execPath != "" {
			bindingPath := execPath + ":/usr/bin/condatainer:ro"
			filtered = append(filtered, bindingPath)
		}
	}

	return filtered
}

// extractHostPath extracts the host path from a bind path.
// Handles formats: path, path:ro, path:container_path, path:container_path:ro
func extractHostPath(bindPath string) string {
	// Split by colons
	parts := strings.Split(bindPath, ":")
	if len(parts) == 0 {
		return ""
	}
	// First part is always the host path
	return parts[0]
}

func isChildPath(child, parent string) bool {
	rel, err := filepath.Rel(parent, child)
	if err != nil {
		return false
	}
	return rel != "." && !strings.HasPrefix(rel, "..")
}

func resolvePath(path string) string {
	abs, err := filepath.Abs(path)
	if err != nil {
		return ""
	}
	real, err := filepath.EvalSymlinks(abs)
	if err != nil {
		return abs
	}
	return real
}
