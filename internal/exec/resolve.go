package exec

import (
	"fmt"
	"os"
	"path/filepath"
	"strings"

	"github.com/Justype/condatainer/internal/config"
	"github.com/Justype/condatainer/internal/utils"
)

// InstalledOverlays scans all image search paths and returns a mapping
// from normalized overlay names to their absolute file path.
// Searches user → scratch → legacy → system directories.
// First match wins (user overlays shadow system ones).
func InstalledOverlays() (map[string]string, error) {
	overlays := map[string]string{}
	dirs := config.GetImageSearchPaths()
	utils.PrintDebug("[RESOLVE] Scanning for installed overlays in: %v", dirs)
	for _, dir := range dirs {
		if err := populateOverlays(dir, overlays); err != nil {
			return nil, err
		}
	}
	utils.PrintDebug("[RESOLVE] Found %d installed overlays", len(overlays))
	return overlays, nil
}

func populateOverlays(dir string, store map[string]string) error {
	entries, err := os.ReadDir(dir)
	if err != nil {
		if os.IsNotExist(err) {
			return nil
		}
		return fmt.Errorf("unable to list overlays in %s: %w", dir, err)
	}

	for _, entry := range entries {
		if entry.IsDir() {
			continue
		}
		name := entry.Name()
		if !utils.IsOverlay(name) {
			continue
		}
		normalized := strings.TrimSuffix(name, filepath.Ext(name))
		normalized = strings.ReplaceAll(normalized, "--", "/")
		store[normalized] = filepath.Join(dir, name)
	}
	return nil
}

// ResolveOverlayPaths coerces user input into absolute overlay paths that exist on disk.
func ResolveOverlayPaths(inputs []string) ([]string, error) {
	if len(inputs) == 0 {
		return nil, nil
	}

	installed, err := InstalledOverlays()
	if err != nil {
		return nil, err
	}

	resolved := make([]string, 0, len(inputs))
	for _, entry := range inputs {
		if entry == "" {
			continue
		}

		// Check for :ro or :rw suffix and strip it temporarily for path resolution
		var suffix string
		pathToResolve := entry
		if strings.HasSuffix(entry, ":ro") {
			suffix = ":ro"
			pathToResolve = strings.TrimSuffix(entry, ":ro")
		} else if strings.HasSuffix(entry, ":rw") {
			suffix = ":rw"
			pathToResolve = strings.TrimSuffix(entry, ":rw")
		}

		if utils.IsOverlay(pathToResolve) {
			absPath, err := filepath.Abs(pathToResolve)
			if err != nil {
				return nil, fmt.Errorf("failed to normalize overlay path %s: %w", pathToResolve, err)
			}
			if !utils.FileExists(absPath) {
				return nil, fmt.Errorf("overlay file %s not found", absPath)
			}
			resolved = append(resolved, absPath+suffix)
			continue
		}

		normalized := utils.NormalizeNameVersion(pathToResolve)
		if normalized == "" {
			return nil, fmt.Errorf("invalid overlay specification %q", entry)
		}

		utils.PrintDebug("[RESOLVE] Looking up overlay %s in installed map", normalized)
		if mapped, ok := installed[normalized]; ok {
			utils.PrintDebug("[RESOLVE] Found overlay %s -> %s", normalized, mapped)
			resolved = append(resolved, mapped+suffix)
			continue
		}

		utils.PrintDebug("[RESOLVE] Overlay %s not in installed map, trying buildOverlayPathFromSpec", normalized)
		pathFromSpec, err := buildOverlayPathFromSpec(normalized)
		if err != nil {
			return nil, err
		}
		resolved = append(resolved, pathFromSpec+suffix)
	}

	return resolved, nil
}

func buildOverlayPathFromSpec(normalized string) (string, error) {
	// Search all image directories for the overlay
	overlayName := fmt.Sprintf("%s.sqf", strings.ReplaceAll(normalized, "/", "--"))

	for _, dir := range config.GetImageSearchPaths() {
		path := filepath.Join(dir, overlayName)
		if utils.FileExists(path) {
			return path, nil
		}
	}

	return "", fmt.Errorf("overlay %s not found (searched: %v)", normalized, config.GetImageSearchPaths())
}

func ensureSingleImage(paths []string) error {
	count := 0
	for _, path := range paths {
		if utils.IsImg(path) {
			count++
		}
	}
	if count > 1 {
		return fmt.Errorf("only one .img overlay can be used at a time")
	}
	return nil
}

func putImgToLast(paths []string) []string {
	var imgPaths []string
	var other []string

	for _, path := range paths {
		if utils.IsImg(path) {
			imgPaths = append(imgPaths, path)
			continue
		}
		other = append(other, path)
	}

	return append(other, imgPaths...)
}
