package build

import (
	"os"
	"path/filepath"
	"sort"
	"strings"

	"github.com/Justype/condatainer/internal/config"
	"github.com/Justype/condatainer/internal/utils"
)

var (
	defBuiltList map[string]bool
)

// getDefListPathForOverlay derives the .def_list path from an overlay's path
// Returns the .def_list file in the same directory as the overlay
func getDefListPathForOverlay(overlayPath string) string {
	imagesDir := filepath.Dir(overlayPath)
	return filepath.Join(imagesDir, ".def_list")
}

// loadDefListFromPath loads entries from a single .def_list file
// Returns a map of normalized names, or empty map if file doesn't exist
func loadDefListFromPath(listPath string) map[string]bool {
	defList := make(map[string]bool)

	data, err := os.ReadFile(listPath)
	if err != nil {
		// File doesn't exist or can't be read - return empty map
		return defList
	}

	// Parse entries
	lines := strings.Split(string(data), "\n")
	for _, line := range lines {
		line = strings.TrimSpace(line)
		if line != "" && !strings.HasPrefix(line, "#") {
			defList[line] = true
		}
	}

	return defList
}

// saveDefListToPath saves entries to a specific .def_list file
// Creates parent directory if needed
func saveDefListToPath(listPath string, entries []string) error {
	// Ensure directory exists
	dir := filepath.Dir(listPath)
	if err := os.MkdirAll(dir, utils.PermDir); err != nil {
		return err
	}

	// Sort entries for consistent output
	sort.Strings(entries)

	// Write to file
	content := "# def-built overlays list (auto-generated)\n"
	content += strings.Join(entries, "\n")
	if len(entries) > 0 {
		content += "\n"
	}

	return os.WriteFile(listPath, []byte(content), utils.PermFile)
}

// loadAllDefLists loads the def list from all image directories
// Searches all image directories (extra → portable → scratch → user)
// and merges all .def_list files found. Higher-priority directories are processed first.
func loadAllDefLists() map[string]bool {
	defList := make(map[string]bool)

	// Search all image directories for .def_list files
	// Directories are returned in priority order (extra → portable → scratch → user)
	for _, imagesDir := range config.GetImageSearchPaths() {
		listPath := filepath.Join(imagesDir, ".def_list")
		data, err := os.ReadFile(listPath)
		if err != nil {
			// File doesn't exist in this directory, continue to next
			continue
		}

		// Parse and merge entries from this .def_list
		// Since we process higher-priority directories first, their entries
		// are added first (though for boolean map, order doesn't affect the result)
		lines := strings.Split(string(data), "\n")
		for _, line := range lines {
			line = strings.TrimSpace(line)
			if line != "" && !strings.HasPrefix(line, "#") {
				defList[line] = true
			}
		}
	}

	return defList
}

// overlayPathInImageSearchPaths reports whether overlayPath is located under any configured image search path.
func overlayPathInImageSearchPaths(overlayPath string) bool {
	if overlayPath == "" {
		return false
	}

	absOverlay, err := filepath.Abs(overlayPath)
	if err != nil {
		absOverlay = filepath.Clean(overlayPath)
	}

	for _, imagesDir := range config.GetImageSearchPaths() {
		if imagesDir == "" {
			continue
		}
		absDir, err := filepath.Abs(imagesDir)
		if err != nil {
			absDir = imagesDir
		}

		rel, err := filepath.Rel(absDir, absOverlay)
		if err != nil {
			continue
		}
		// rel == "." means same dir; otherwise ensure overlay is a descendant (not starting with "..")
		if rel == "." || !strings.HasPrefix(rel, "..") {
			return true
		}
	}

	return false
}

// UpdateDefBuiltList adds an overlay name to the def list in the same directory as the overlay
// This should be called after building a .def overlay
func UpdateDefBuiltList(nameVersion, overlayPath string) {
	// Only update .def_list when the overlay is inside one of the configured image search paths
	if !overlayPathInImageSearchPaths(overlayPath) {
		utils.PrintDebug(".def_list not updated for %s — path %s is not in image search paths", nameVersion, overlayPath)
		return
	}

	normalized := utils.NormalizeNameVersion(nameVersion)

	// Get the .def_list path in the same directory as the overlay
	listPath := getDefListPathForOverlay(overlayPath)

	// Load existing entries from that specific .def_list file
	entries := loadDefListFromPath(listPath)

	// Add the new entry
	entries[normalized] = true

	// Convert map to sorted slice
	entrySlice := make([]string, 0, len(entries))
	for name := range entries {
		entrySlice = append(entrySlice, name)
	}

	// Save back to the same .def_list file
	if err := saveDefListToPath(listPath, entrySlice); err != nil {
		utils.PrintDebug("Failed to update .def_list at %s: %v", listPath, err)
	}

	// Clear global cache to force reload on next GetDefBuiltList() call
	defBuiltList = nil
}

// GetDefBuiltList returns the def list, loading from all directories if needed
// Reads .def_list from ALL image directories (extra → portable → scratch → user) and merges them.
// Writes only happen to the first writable directory when UpdateDefBuiltList() or RemoveFromDefBuiltList() is called.
func GetDefBuiltList() map[string]bool {
	if defBuiltList != nil {
		return defBuiltList
	}

	// Load from all image directories following priority order (extra → portable → scratch → user)
	defBuiltList = loadAllDefLists()

	// Merge in entries from build script metadata (e.g., .def files in build-scripts/)
	// Note: We don't auto-save here to avoid persisting entries from read-only directories
	// into the writable directory. Entries are only saved to the first writable directory
	// (following the same priority order) when explicitly added via UpdateDefBuiltList()
	// or removed via RemoveFromDefBuiltList()
	scripts, err := GetAllBuildScripts(true)
	if err == nil && len(scripts) > 0 {
		for name, info := range scripts {
			if info.IsContainer {
				defBuiltList[name] = true
			}
		}
	}

	return defBuiltList
}

// RemoveFromDefBuiltList removes an overlay name from the def list in the same directory as the overlay
// This should be called after removing a .def overlay
func RemoveFromDefBuiltList(nameVersion, overlayPath string) {
	// Only modify .def_list when the overlay is inside one of the configured image search paths
	if !overlayPathInImageSearchPaths(overlayPath) {
		utils.PrintDebug(".def_list not modified for %s — path %s is not in image search paths", nameVersion, overlayPath)
		return
	}

	normalized := utils.NormalizeNameVersion(nameVersion)

	// Get the .def_list path in the same directory as the overlay
	listPath := getDefListPathForOverlay(overlayPath)

	// Load existing entries from that specific .def_list file
	entries := loadDefListFromPath(listPath)

	// Remove the entry if it exists
	if _, exists := entries[normalized]; !exists {
		// Entry not in this .def_list file, nothing to do
		return
	}

	delete(entries, normalized)

	// Convert map to sorted slice
	entrySlice := make([]string, 0, len(entries))
	for name := range entries {
		entrySlice = append(entrySlice, name)
	}

	// Save back to the same .def_list file
	if err := saveDefListToPath(listPath, entrySlice); err != nil {
		utils.PrintDebug("Failed to update .def_list at %s: %v", listPath, err)
	}

	// Clear global cache to force reload on next GetDefBuiltList() call
	defBuiltList = nil
}
