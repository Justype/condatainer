package build

import (
	"os"
	"path/filepath"
	"sort"
	"strings"

	"github.com/Justype/condatainer/internal/utils"
)

var (
	defBuiltWhitelist     map[string]bool
	defBuiltWhitelistPath string
)

// getDefBuiltWhitelistPath returns the path to the persistent whitelist file
func getDefBuiltWhitelistPath() string {
	if defBuiltWhitelistPath != "" {
		return defBuiltWhitelistPath
	}

	// Try user config directory first
	userConfigDir, err := os.UserConfigDir()
	if err == nil {
		defBuiltWhitelistPath = filepath.Join(userConfigDir, "condatainer", ".def_built_overlays")
		return defBuiltWhitelistPath
	}

	// Fallback to temp directory
	defBuiltWhitelistPath = filepath.Join(os.TempDir(), ".condatainer_def_built_overlays")
	return defBuiltWhitelistPath
}

// loadDefBuiltWhitelist loads the whitelist from persistent file
func loadDefBuiltWhitelist() map[string]bool {
	whitelist := make(map[string]bool)
	whitelistPath := getDefBuiltWhitelistPath()
	data, err := os.ReadFile(whitelistPath)
	if err != nil {
		return whitelist
	}

	// Each line is a normalized overlay name
	lines := strings.Split(string(data), "\n")
	for _, line := range lines {
		line = strings.TrimSpace(line)
		if line != "" && !strings.HasPrefix(line, "#") {
			whitelist[line] = true
		}
	}

	return whitelist
}

// saveDefBuiltWhitelist saves the whitelist to persistent file
func saveDefBuiltWhitelist(whitelist map[string]bool) {
	whitelistPath := getDefBuiltWhitelistPath()

	// Ensure directory exists
	dir := filepath.Dir(whitelistPath)
	if err := os.MkdirAll(dir, 0755); err != nil {
		return
	}

	// Create sorted list of names
	names := make([]string, 0, len(whitelist))
	for name := range whitelist {
		names = append(names, name)
	}
	sort.Strings(names)

	// Write to file
	content := "# .def-built overlays (auto-generated)\n"
	content += strings.Join(names, "\n") + "\n"
	os.WriteFile(whitelistPath, []byte(content), 0644)
}

// UpdateDefBuiltWhitelist adds an overlay name to the whitelist and saves it
// This should be called after building a .def overlay
func UpdateDefBuiltWhitelist(nameVersion string) {
	normalized := utils.NormalizeNameVersion(nameVersion)

	// Load current whitelist
	if defBuiltWhitelist == nil {
		defBuiltWhitelist = loadDefBuiltWhitelist()
	}

	// Add new entry
	defBuiltWhitelist[normalized] = true

	// Save to file
	saveDefBuiltWhitelist(defBuiltWhitelist)
}

// GetDefBuiltWhitelist returns the whitelist, loading from file if needed
func GetDefBuiltWhitelist() map[string]bool {
	if defBuiltWhitelist != nil {
		return defBuiltWhitelist
	}

	// Always load persistent file first (contains entries from -s/-n builds)
	defBuiltWhitelist = loadDefBuiltWhitelist()

	// Merge in entries from build script metadata
	scripts, err := GetAllBuildScripts(true)
	if err == nil && len(scripts) > 0 {
		for name, info := range scripts {
			if info.IsContainer {
				defBuiltWhitelist[name] = true
			}
		}
		// Save merged result to persistent file
		saveDefBuiltWhitelist(defBuiltWhitelist)
	}

	return defBuiltWhitelist
}
