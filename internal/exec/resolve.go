package exec

import (
	"fmt"
	"os"
	"path/filepath"
	"strings"

	"github.com/Justype/condatainer/internal/config"
	"github.com/Justype/condatainer/internal/utils"
)

// InstalledOverlays scans the images directory and returns a mapping
// from normalized overlay names to their absolute file path.
// All overlays (including ref overlays) are stored in a single ImagesDir.
func InstalledOverlays() (map[string]string, error) {
	overlays := map[string]string{}
	if err := populateOverlays(config.Global.ImagesDir, overlays); err != nil {
		return nil, err
	}
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
		if utils.IsOverlay(entry) {
			absPath, err := filepath.Abs(entry)
			if err != nil {
				return nil, fmt.Errorf("failed to normalize overlay path %s: %w", entry, err)
			}
			if !utils.FileExists(absPath) {
				return nil, fmt.Errorf("overlay file %s not found", absPath)
			}
			resolved = append(resolved, absPath)
			continue
		}

		normalized := utils.NormalizeNameVersion(entry)
		if normalized == "" {
			return nil, fmt.Errorf("invalid overlay specification %q", entry)
		}
		if mapped, ok := installed[normalized]; ok {
			resolved = append(resolved, mapped)
			continue
		}

		pathFromSpec, err := buildOverlayPathFromSpec(normalized)
		if err != nil {
			return nil, err
		}
		resolved = append(resolved, pathFromSpec)
	}

	return resolved, nil
}

func buildOverlayPathFromSpec(normalized string) (string, error) {
	// All overlays are stored in a single ImagesDir (including ref overlays)
	overlayName := fmt.Sprintf("%s.sqf", strings.ReplaceAll(normalized, "/", "--"))
	path := filepath.Join(config.Global.ImagesDir, overlayName)

	if utils.FileExists(path) {
		return path, nil
	}
	return "", fmt.Errorf("overlay %s not found at %s", normalized, path)
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
