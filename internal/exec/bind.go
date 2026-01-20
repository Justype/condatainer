package exec

import (
	"os"
	"path/filepath"
	"strings"
)

// BindPaths returns unique absolute directories suitable for --bind flags while omitting nested duplicates.
func BindPaths(paths ...string) []string {
	absPaths := make([]string, 0, len(paths)+3)
	for _, path := range paths {
		if path == "" {
			continue
		}
		if resolved := resolvePath(path); resolved != "" {
			absPaths = append(absPaths, resolved)
		}
	}

	if home, err := os.UserHomeDir(); err == nil && home != "" {
		if resolved := resolvePath(home); resolved != "" {
			absPaths = append(absPaths, resolved)
		}
	}
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
		for _, parent := range unique {
			if parent == candidate {
				continue
			}
			if isChildPath(candidate, parent) {
				skip = true
				break
			}
		}
		if !skip {
			filtered = append(filtered, candidate)
		}
	}

	return filtered
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
