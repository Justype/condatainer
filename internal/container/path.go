package container

import (
	"fmt"
	"path/filepath"
	"strings"

	"github.com/Justype/condatainer/internal/utils"
)

// FormatOverlayMount formats an overlay path with the appropriate :ro or :rw suffix
func FormatOverlayMount(path string, writable bool) string {
	// Check if path already has :ro or :rw suffix
	if strings.HasSuffix(path, ":ro") || strings.HasSuffix(path, ":rw") {
		return path
	}

	if utils.IsImg(path) {
		// For .img files, add :ro or :rw suffix based on writable flag
		if writable {
			return path + ":rw"
		}
		return path + ":ro"
	} else if utils.IsSqf(path) {
		// For .sqf files, always add :ro suffix (they're always read-only)
		return path + ":ro"
	}
	// For .sif and other files, no suffix needed
	return path
}

// BuildPathEnv constructs the PATH environment variable based on overlays
func BuildPathEnv(overlays []string) string {
	paths := []string{"/usr/sbin", "/usr/bin"}

	for _, overlay := range overlays {
		// Strip :ro or :rw suffix for path checking
		cleanOverlay := strings.TrimSuffix(strings.TrimSuffix(overlay, ":ro"), ":rw")
		name := strings.TrimSuffix(filepath.Base(cleanOverlay), filepath.Ext(cleanOverlay))
		normalized := utils.NormalizeNameVersion(name)
		if normalized == "" {
			continue
		}
		if strings.Count(normalized, "/") > 1 {
			continue
		}
		var relative string
		if utils.IsImg(cleanOverlay) {
			relative = "/ext3/env/bin"
		} else if utils.IsSqf(cleanOverlay) {
			relative = fmt.Sprintf("/cnt/%s/bin", normalized)
		} else {
			utils.PrintWarning("Unknown overlay file extension for %s. Skipping PATH addition.", utils.StylePath(overlay))
			continue
		}
		paths = append([]string{relative}, paths...)
	}

	return strings.Join(paths, ":")
}
