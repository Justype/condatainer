package image

import (
	"fmt"
	"os"
	"path/filepath"
	"strings"

	"github.com/Justype/condatainer/internal/config"
	"github.com/Justype/condatainer/internal/utils"
)

// ScanOptions selects what ScanOverlays walks and how it keys the result.
type ScanOptions struct {
	Dirs []string // search paths; nil means config.GetImageSearchPaths()

	// Aliases also keys a <base>/<name> image by its bare <name>. Off by
	// default: a dependency or a removal names an image exactly.
	Aliases bool
}

// ScanOverlays walks the image search paths and returns every installed overlay
// keyed by normalized name (samtools--1.21.sqf → "samtools/1.21"), with its
// copies in search order. An unreadable directory is reported with a usable map.
func ScanOverlays(opts ScanOptions) (map[string][]string, error) {
	dirs := opts.Dirs
	if dirs == nil {
		dirs = config.GetImageSearchPaths()
	}

	var alias string
	if opts.Aliases {
		if base := config.ResolvedDefaultDistro(); base != "" {
			alias = base + "/"
		}
	}

	overlays := make(map[string][]string)
	var firstErr error
	for _, dir := range dirs {
		entries, err := os.ReadDir(dir)
		if err != nil {
			if !os.IsNotExist(err) && firstErr == nil {
				firstErr = fmt.Errorf("unable to list overlays in %s: %w", dir, err)
			}
			continue
		}
		for _, entry := range entries {
			if entry.IsDir() || !utils.IsOverlay(entry.Name()) {
				continue
			}
			name := strings.ReplaceAll(strings.TrimSuffix(entry.Name(), filepath.Ext(entry.Name())), "--", "/")
			path := filepath.Join(dir, entry.Name())
			overlays[name] = append(overlays[name], path)
			if bare, ok := strings.CutPrefix(name, alias); alias != "" && ok && bare != "" {
				overlays[bare] = append(overlays[bare], path)
			}
		}
	}
	return overlays, firstErr
}

// FirstPaths reduces a scan to the highest-priority copy of each overlay.
func FirstPaths(scan map[string][]string) map[string]string {
	paths := make(map[string]string, len(scan))
	for name, copies := range scan {
		if len(copies) > 0 {
			paths[name] = copies[0]
		}
	}
	return paths
}

// Names reduces a scan to the set of installed overlay names.
func Names(scan map[string][]string) map[string]bool {
	names := make(map[string]bool, len(scan))
	for name := range scan {
		names[name] = true
	}
	return names
}
