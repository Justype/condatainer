package container

import (
	"fmt"
	"log/slog"
	"path/filepath"
	"strings"

	"github.com/Justype/condatainer/catalog"
	"github.com/Justype/condatainer/internal/config"
	"github.com/Justype/condatainer/internal/image"
	"github.com/Justype/condatainer/internal/utils"
)

// cachedInstalledOverlays caches the result of the last InstalledOverlays scan.
// Nil means the cache is cold or has been invalidated.
var cachedInstalledOverlays map[string]string

// InstalledOverlays maps every installed overlay name to its highest-priority
// copy, plus the bare-name alias of each <base>/<name> image. Cached for the
// process; call InvalidateInstalledOverlaysCache after installing a new image.
func InstalledOverlays() (map[string]string, error) {
	if cachedInstalledOverlays != nil {
		return cachedInstalledOverlays, nil
	}
	scan, err := image.ScanOverlays(image.ScanOptions{Aliases: true})
	if err != nil {
		// Resolution feeds container launch, so an unreadable image directory is
		// reported rather than silently resolving to a shorter list.
		return nil, err
	}
	overlays := image.FirstPaths(scan)
	slog.Default().Debug("found installed overlays", "count", len(overlays))
	cachedInstalledOverlays = overlays
	return overlays, nil
}

// InvalidateInstalledOverlaysCache clears the cached overlay map so the next call
// to InstalledOverlays rescans the image directories.
func InvalidateInstalledOverlaysCache() {
	cachedInstalledOverlays = nil
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

		// If it's an absolute path to an existing file, use it directly (regardless of extension)
		if filepath.IsAbs(pathToResolve) && utils.FileExists(pathToResolve) {
			resolved = append(resolved, pathToResolve+suffix)
			continue
		}

		if utils.IsOverlay(pathToResolve) || utils.IsSif(pathToResolve) {
			absPath, err := filepath.Abs(pathToResolve)
			if err != nil {
				return nil, fmt.Errorf("failed to normalize overlay path %s: %w", pathToResolve, err)
			}
			if !utils.FileExists(absPath) {
				return nil, missingOverlayError(absPath)
			}
			resolved = append(resolved, absPath+suffix)
			continue
		}

		dep, err := catalog.ParseDep(pathToResolve)
		if err != nil {
			return nil, fmt.Errorf("invalid overlay specification %q", entry)
		}
		normalized := dep.NameVersion()

		// With a constraint: pick the highest installed version that satisfies it,
		// up to and including the preferred version (upper bound).
		if dep.Op != "" {
			prefix := dep.Name + "/"
			bestVer := ""
			bestPath := ""
			for key, path := range installed {
				if strings.HasPrefix(key, prefix) {
					ver := strings.TrimPrefix(key, prefix)
					if dep.Satisfies(ver) {
						if bestVer == "" || catalog.CompareVersions(ver, bestVer) > 0 {
							bestVer = ver
							bestPath = path
						}
					}
				}
			}
			if bestPath != "" {
				slog.Default().Debug("constraint satisfied", "spec", dep.String(), "by", dep.Name+"/"+bestVer)
				resolved = append(resolved, bestPath+suffix)
				continue
			}
			// No satisfying version installed — fall through to exact preferred version.
		}

		slog.Default().Debug("looking up overlay", "name", normalized)
		if mapped, ok := installed[normalized]; ok {
			slog.Default().Debug("found overlay", "name", normalized, "path", mapped)
			resolved = append(resolved, mapped+suffix)
			continue
		}

		slog.Default().Debug("overlay not in installed map, trying buildOverlayPathFromSpec", "name", normalized)
		pathFromSpec, err := buildOverlayPathFromSpec(normalized)
		if err != nil {
			return nil, err
		}
		resolved = append(resolved, pathFromSpec+suffix)
	}

	return resolved, nil
}

// missingOverlayError reports an overlay path that does not exist. Naming a
// specific .img by path is a request for a writable mount by that exact name —
// guessing that a read-only .sqf beside it was meant instead would be exactly
// the silent substitution CLAUDE.md's refusal-over-guessing stance rules out.
// When the .img would have paired with a snapshot (LookupSnapshot), the
// message names it and the two ways to proceed, so the fix is obvious instead
// of requiring the reader to already know overlay snapshots exist.
func missingOverlayError(absPath string) error {
	if !utils.IsImg(absPath) {
		return fmt.Errorf("overlay file %s not found", absPath)
	}
	lookup := LookupSnapshot(absPath)
	if lookup.Path == "" {
		return fmt.Errorf("overlay file %s not found", absPath)
	}
	return fmt.Errorf("overlay file %s not found; a snapshot %s exists here — `overlay create %s` to continue from it, or `exec -o %s` to run against it read-only",
		absPath, lookup.Path, filepath.Base(absPath), lookup.Path)
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

	// Fallback: bare name (no "/") → try <base>--<name>.sqf
	if !strings.Contains(normalized, "/") && config.ResolvedDefaultDistro() != "" {
		prefixed := config.ResolvedDefaultDistro() + "--" + normalized + ".sqf"
		for _, dir := range config.GetImageSearchPaths() {
			path := filepath.Join(dir, prefixed)
			if utils.FileExists(path) {
				slog.Default().Debug("found overlay via OS-prefix fallback", "path", path)
				return path, nil
			}
		}
	}

	return "", fmt.Errorf("overlay %s not found (searched: %v)", normalized, config.GetImageSearchPaths())
}
