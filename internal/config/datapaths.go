package config

import (
	"fmt"
	"os"
	"path/filepath"
	"slices"
	"strings"
	"sync"

	"github.com/Justype/condatainer/internal/utils"
	"github.com/spf13/viper"
)

// =============================================================================
// Types & Globals
// =============================================================================

// DataPaths holds the search paths for data directories.
// Search order: extra → portable → scratch → user (first match wins for lookups).
// Write operations go to the first writable path in the same order.
type DataPaths struct {
	ImagesDirs        []string // Search paths for images
	BuildScriptsDirs  []string // Search paths for build scripts
	HelperScriptsDirs []string // Search paths for helper scripts
}

// GlobalDataPaths holds the computed data paths.
var GlobalDataPaths DataPaths

//=============================================================================
// Config Getters (env/viper)
// =============================================================================

// splitPipe splits an env var value on "|".
// Use when ":" is ambiguous — e.g. values contain "://" (URLs) or ":ro"/":rw" markers.
func splitPipe(val string) []string {
	var result []string
	for _, entry := range strings.Split(val, "|") {
		if entry = strings.TrimSpace(entry); entry != "" {
			result = append(result, entry)
		}
	}
	return result
}

// splitPipeOrColon splits an env var value on "|" or ":" ("|" takes precedence).
// Use for plain paths or plain strings where ":" is unambiguous (no markers).
func splitPipeOrColon(val string) []string {
	sep := ":"
	if strings.Contains(val, "|") {
		sep = "|"
	}
	var result []string
	for _, entry := range strings.Split(val, sep) {
		if entry = strings.TrimSpace(entry); entry != "" {
			result = append(result, entry)
		}
	}
	return result
}

// getEnvSlice reads a string-slice config key, checking the corresponding env var first.
// The env var name is derived from the key: "extra_image_dirs" → "CNT_EXTRA_IMAGE_DIRS".
// The split function controls how the env var value is parsed.
func getEnvSlice(key string, split func(string) []string) []string {
	envKey := "CNT_" + strings.ToUpper(key)
	if envVal := os.Getenv(envKey); envVal != "" {
		if vals := split(envVal); len(vals) > 0 {
			return vals
		}
	}
	return viper.GetStringSlice(key)
}

// ParseDirEntry splits a dir entry ("path", "path:ro", or "path:rw") into (path, readOnly).
// :ro = force search-only (never written). :rw = explicit writable annotation (same as no marker).
// Default (no marker) is writable if filesystem perms allow — matches all other dir types.
func ParseDirEntry(entry string) (string, bool) {
	if after, ok := strings.CutSuffix(entry, ":ro"); ok {
		return after, true
	}
	if after, ok := strings.CutSuffix(entry, ":rw"); ok {
		return after, false
	}
	return entry, false
}

// GetChannels returns the conda channels from config or environment.
// CNT_CHANNELS supports "|" and ":" as separators ("|" takes precedence).
func GetChannels() []string { return getEnvSlice("channels", splitPipeOrColon) }

// GetExtraBaseDirs returns extra base directories from config or environment.
// CNT_EXTRA_BASE_DIRS supports "|" and ":" as separators ("|" takes precedence).
func GetExtraBaseDirs() []string { return getEnvSlice("extra_base_dirs", splitPipeOrColon) }

// GetExtraImageDirs returns explicit extra image directories from config or environment.
// CNT_EXTRA_IMAGE_DIRS uses "|" as separator; entries support ":ro"/":rw" markers.
func GetExtraImageDirs() []string { return getEnvSlice("extra_image_dirs", splitPipe) }

// GetExtraBuildDirs returns explicit extra build-scripts directories from config or environment.
// CNT_EXTRA_BUILD_DIRS supports "|" and ":" as separators ("|" takes precedence).
func GetExtraBuildDirs() []string { return getEnvSlice("extra_build_dirs", splitPipeOrColon) }

// GetExtraHelperDirs returns explicit extra helper-scripts directories from config or environment.
// CNT_EXTRA_HELPER_DIRS uses "|" as separator; entries support ":ro"/":rw" markers.
func GetExtraHelperDirs() []string { return getEnvSlice("extra_helper_dirs", splitPipe) }

// GetExtraScriptsLinks returns extra remote build script source URLs.
// CNT_EXTRA_SCRIPTS_LINKS uses "|" as separator (URLs contain "://").
func GetExtraScriptsLinks() []string { return getEnvSlice("extra_scripts_links", splitPipe) }

// =============================================================================
// Dir Getters (single path, no IO)
// =============================================================================

// GetUserDataDir returns the user's data directory following XDG spec.
// Returns $XDG_DATA_HOME/condatainer or ~/.local/share/condatainer
func GetUserDataDir() string {
	if dataHome := os.Getenv("XDG_DATA_HOME"); dataHome != "" {
		return filepath.Join(dataHome, "condatainer")
	}
	if home, err := os.UserHomeDir(); err == nil {
		return filepath.Join(home, ".local", "share", "condatainer")
	}
	return ""
}

// GetUserCacheDir returns the user's cache directory following XDG spec.
// Returns $XDG_CACHE_HOME/condatainer or ~/.cache/condatainer
func GetUserCacheDir() string {
	if cacheHome := os.Getenv("XDG_CACHE_HOME"); cacheHome != "" {
		return filepath.Join(cacheHome, "condatainer")
	}
	if home, err := os.UserHomeDir(); err == nil {
		return filepath.Join(home, ".cache", "condatainer")
	}
	return ""
}

// GetUserStateDir returns the user's state directory following XDG spec.
// Returns $XDG_STATE_HOME/condatainer or ~/.local/state/condatainer
func GetUserStateDir() string {
	if stateHome := os.Getenv("XDG_STATE_HOME"); stateHome != "" {
		return filepath.Join(stateHome, "condatainer")
	}
	if home, err := os.UserHomeDir(); err == nil {
		return filepath.Join(home, ".local", "state", "condatainer")
	}
	return ""
}

// GetScratchDataDir returns the scratch data directory for HPC systems.
// Returns $SCRATCH/condatainer if SCRATCH is set.
func GetScratchDataDir() string {
	if scratch := os.Getenv("SCRATCH"); scratch != "" {
		return filepath.Join(scratch, "condatainer")
	}
	return ""
}

// GetPortableDataDir returns the portable data directory for group/shared use.
// Only detects truly portable installations - directories next to the executable
// that are NOT in standard user locations (~/.config, ~/.local, etc.).
// This is useful for shared group directories or portable installations.
//
// Detection logic:
//  1. If executable is in <dir>/bin/ and <dir> is not a generic parent ($HOME, /usr, etc.), use <dir>
//  2. Otherwise check for <exeDir>/condatainer/ directory
func GetPortableDataDir() string {
	portableDirOnce.Do(func() {
		portableDirCache = detectPortableDataDir()
	})
	return portableDirCache
}

// IsPortable returns true if the current installation is portable (not in standard user locations).
func IsPortable() bool {
	return GetPortableDataDir() != ""
}

// GetPortableBinDir returns the bin directory of a portable installation, or empty string if not portable.
func GetPortableBinDir() string {
	portableBinOnce.Do(func() {
		if portableDir := GetPortableDataDir(); portableDir != "" {
			binDir := filepath.Join(portableDir, "bin")
			if stat, err := os.Stat(binDir); err == nil && stat.IsDir() {
				portableBinCache = binDir
			}
		}
	})
	return portableBinCache
}

// =============================================================================
// Portable Detection (internal)
// =============================================================================

var (
	portableDirOnce  sync.Once
	portableDirCache string
	portableBinOnce  sync.Once
	portableBinCache string
)

func detectPortableDataDir() string {
	exe, err := os.Executable()
	if err != nil {
		return ""
	}

	exeDir := filepath.Dir(exe)
	if realExe, err := filepath.EvalSymlinks(exe); err == nil {
		exeDir = filepath.Dir(realExe)
	}

	// Executable must be in <install_dir>/bin/
	if filepath.Base(exeDir) == "bin" {
		parentDir := filepath.Dir(exeDir)
		if !isNonPortableParent(parentDir) {
			return parentDir
		}
	}

	return ""
}

// isNonPortableParent checks if a directory is a common user or system bin parent
// that should not be treated as a portable installation root.
func isNonPortableParent(dir string) bool {
	home, _ := os.UserHomeDir()
	excludedParents := []string{
		home,
		filepath.Join(home, ".local"),
		"/usr",
		"/usr/local",
		"/opt",
	}
	return slices.Contains(excludedParents, dir)
}

// =============================================================================
// Init
// =============================================================================

// InitDataPaths initializes GlobalDataPaths based on environment and config.
// This should be called after LoadDefaults and LoadFromViper.
func InitDataPaths() {
	GlobalDataPaths = DataPaths{
		ImagesDirs:        buildImageSearchPaths(),
		BuildScriptsDirs:  buildScriptSearchPaths("build-scripts"),
		HelperScriptsDirs: buildScriptSearchPaths("helper-scripts"),
	}
}

// =============================================================================
// Search Path Builders (internal, read)
// =============================================================================

// buildImageSearchPaths builds the search paths for images.
// Priority: extra_image_dirs → extra_base_dirs → portable → scratch → user
func buildImageSearchPaths() []string {
	var paths []string
	seen := make(map[string]bool)

	addPath := func(dir string) {
		if dir == "" {
			return
		}
		dir = os.ExpandEnv(dir)
		absDir, err := filepath.Abs(dir)
		if err != nil {
			absDir = dir
		}
		if !seen[absDir] {
			seen[absDir] = true
			paths = append(paths, absDir)
		}
	}

	for _, entry := range GetExtraImageDirs() {
		path, _ := ParseDirEntry(entry)
		addPath(path)
	}
	for _, baseDir := range GetExtraBaseDirs() {
		addPath(filepath.Join(baseDir, "images"))
	}
	if portableDir := GetPortableDataDir(); portableDir != "" {
		addPath(filepath.Join(portableDir, "images"))
	}
	if scratchDir := GetScratchDataDir(); scratchDir != "" {
		addPath(filepath.Join(scratchDir, "images"))
	}
	if userDir := GetUserDataDir(); userDir != "" {
		addPath(filepath.Join(userDir, "images"))
	}

	return paths
}

// buildScriptSearchPaths builds the search paths for scripts (build-scripts or helper-scripts).
// Priority: extra_build_dirs (build-scripts only) → extra_base_dirs → portable → scratch → user
func buildScriptSearchPaths(subdir string) []string {
	var paths []string
	seen := make(map[string]bool)

	addPath := func(dir string) {
		if dir == "" {
			return
		}
		dir = os.ExpandEnv(dir)
		absDir, err := filepath.Abs(dir)
		if err != nil {
			absDir = dir
		}
		if !seen[absDir] {
			seen[absDir] = true
			paths = append(paths, absDir)
		}
	}

	switch subdir {
	case "build-scripts":
		for _, path := range GetExtraBuildDirs() {
			addPath(path)
		}
	case "helper-scripts":
		for _, entry := range GetExtraHelperDirs() {
			path, _ := ParseDirEntry(entry)
			addPath(path)
		}
	}
	for _, baseDir := range GetExtraBaseDirs() {
		addPath(filepath.Join(baseDir, subdir))
	}
	if portableDir := GetPortableDataDir(); portableDir != "" {
		addPath(filepath.Join(portableDir, subdir))
		// Auto-detect cnt-scripts subdir for build-scripts (split-repo layout)
		if subdir == "build-scripts" && DirExists(filepath.Join(portableDir, "cnt-scripts")) {
			addPath(filepath.Join(portableDir, "cnt-scripts", subdir))
		}
	}
	if scratchDir := GetScratchDataDir(); scratchDir != "" {
		addPath(filepath.Join(scratchDir, subdir))
	}
	if userDir := GetUserDataDir(); userDir != "" {
		addPath(filepath.Join(userDir, subdir))
	}

	return paths
}

// =============================================================================
// Search Path Getters (public, read)
// =============================================================================

// GetImageSearchPaths returns all paths to search for images.
func GetImageSearchPaths() []string {
	if len(GlobalDataPaths.ImagesDirs) == 0 {
		InitDataPaths()
	}
	return GlobalDataPaths.ImagesDirs
}

// GetBuildScriptSearchPaths returns all paths to search for build scripts.
func GetBuildScriptSearchPaths() []string {
	if len(GlobalDataPaths.BuildScriptsDirs) == 0 {
		InitDataPaths()
	}
	return GlobalDataPaths.BuildScriptsDirs
}

// GetHelperScriptSearchPaths returns all paths to search for helper scripts.
func GetHelperScriptSearchPaths() []string {
	if len(GlobalDataPaths.HelperScriptsDirs) == 0 {
		InitDataPaths()
	}
	return GlobalDataPaths.HelperScriptsDirs
}

// GetCacheSearchPaths returns all personal cache directories to search.
// Cache is always personal (scratch → user) to avoid cross-user pollution.
func GetCacheSearchPaths() []string {
	dirs := cacheWriteDirs()
	paths := make([]string, len(dirs))
	for i, d := range dirs {
		paths[i] = d.Path
	}
	return paths
}

// =============================================================================
// Write Dir Helpers (internal)
// =============================================================================

// SearchDir is a candidate directory for write operations.
// Personal dirs (scratch, XDG user) are created on first use via IsWritableDir.
// Shared dirs (portable, extra_base_dirs, extra_*_dirs) are only probed via CanWriteToDir.
type SearchDir struct {
	Path     string
	Personal bool
}

// firstWritableDir returns the first writable path from the slice.
// Personal dirs are created if needed (IsWritableDir); shared dirs are probed only (CanWriteToDir).
func firstWritableDir(dirs []SearchDir) string {
	for _, d := range dirs {
		if d.Personal {
			if utils.IsWritableDir(d.Path) {
				return d.Path
			}
		} else {
			if utils.CanWriteToDir(d.Path) {
				return d.Path
			}
		}
	}
	return ""
}

// imageWriteDirs returns the ordered write candidates for image directories.
// Shared: extra_image_dirs (non-:ro), extra_base_dirs, portable.
// Personal: scratch, user.
func imageWriteDirs() []SearchDir {
	var dirs []SearchDir
	for _, entry := range GetExtraImageDirs() {
		path, readOnly := ParseDirEntry(entry)
		if !readOnly && path != "" {
			dirs = append(dirs, SearchDir{Path: filepath.Clean(os.ExpandEnv(path))})
		}
	}
	for _, base := range GetExtraBaseDirs() {
		dirs = append(dirs, SearchDir{Path: filepath.Join(os.ExpandEnv(base), "images")})
	}
	if p := GetPortableDataDir(); p != "" {
		dirs = append(dirs, SearchDir{Path: filepath.Join(p, "images")})
	}
	if s := GetScratchDataDir(); s != "" {
		dirs = append(dirs, SearchDir{Path: filepath.Join(s, "images"), Personal: true})
	}
	if u := GetUserDataDir(); u != "" {
		dirs = append(dirs, SearchDir{Path: filepath.Join(u, "images"), Personal: true})
	}
	return dirs
}

// helperWriteDirs returns the ordered write candidates for helper-scripts directories.
// Shared: extra_helper_dirs (non-:ro), extra_base_dirs, portable.
// Personal: scratch, user.
func helperWriteDirs() []SearchDir {
	var dirs []SearchDir
	for _, entry := range GetExtraHelperDirs() {
		path, readOnly := ParseDirEntry(entry)
		if !readOnly && path != "" {
			dirs = append(dirs, SearchDir{Path: filepath.Clean(os.ExpandEnv(path))})
		}
	}
	for _, base := range GetExtraBaseDirs() {
		dirs = append(dirs, SearchDir{Path: filepath.Join(os.ExpandEnv(base), "helper-scripts")})
	}
	if p := GetPortableDataDir(); p != "" {
		dirs = append(dirs, SearchDir{Path: filepath.Join(p, "helper-scripts")})
	}
	if s := GetScratchDataDir(); s != "" {
		dirs = append(dirs, SearchDir{Path: filepath.Join(s, "helper-scripts"), Personal: true})
	}
	if u := GetUserDataDir(); u != "" {
		dirs = append(dirs, SearchDir{Path: filepath.Join(u, "helper-scripts"), Personal: true})
	}
	return dirs
}

// cacheWriteDirs returns the ordered write candidates for cache directories.
// Cache is always personal to avoid cross-user pollution.
func cacheWriteDirs() []SearchDir {
	var dirs []SearchDir
	if s := GetScratchDataDir(); s != "" {
		dirs = append(dirs, SearchDir{Path: filepath.Join(s, "cache"), Personal: true})
	}
	if u := GetUserCacheDir(); u != "" {
		dirs = append(dirs, SearchDir{Path: u, Personal: true})
	}
	return dirs
}

// =============================================================================
// Write Dir Resolvers (public)
// =============================================================================

// GetWritableImagesDir returns the first writable images directory.
// Explicit extra_image_dirs entries marked :ro are skipped.
// Shared dirs (extra_image_dirs, extra_base_dirs, portable) are probed only.
// Personal dirs (scratch, user) are created on first use.
func GetWritableImagesDir() (string, error) {
	dirs := imageWriteDirs()
	if dir := firstWritableDir(dirs); dir != "" {
		return dir, nil
	}
	paths := make([]string, len(dirs))
	for i, d := range dirs {
		paths[i] = d.Path
	}
	return "", fmt.Errorf("no writable images directory found (searched: %v)", paths)
}

// GetWritableHelperScriptsDir returns the first writable helper scripts directory.
// Explicit extra_helper_dirs entries marked :ro are skipped.
// Shared dirs (extra_helper_dirs, extra_base_dirs, portable) are probed only.
// Personal dirs (scratch, user) are created on first use.
func GetWritableHelperScriptsDir() (string, error) {
	dirs := helperWriteDirs()
	if dir := firstWritableDir(dirs); dir != "" {
		return dir, nil
	}
	paths := make([]string, len(dirs))
	for i, d := range dirs {
		paths[i] = d.Path
	}
	return "", fmt.Errorf("no writable helper scripts directory found (searched: %v)", paths)
}

// GetWritableCacheDir returns the first writable personal cache directory.
// Always personal (scratch → user cache) — never writes to shared dirs.
func GetWritableCacheDir() (string, error) {
	dirs := cacheWriteDirs()
	if dir := firstWritableDir(dirs); dir != "" {
		return dir, nil
	}
	return "", fmt.Errorf("no writable cache directory found")
}

// GetBaseImageWritePath returns the path where a new base image should be written.
func GetBaseImageWritePath() (string, error) {
	dir, err := GetWritableImagesDir()
	if err != nil {
		return "", err
	}
	return filepath.Join(dir, BaseImageSifName()), nil
}

// =============================================================================
// Finders
// =============================================================================

// FindImage searches all image paths for an image by name.
// Returns the full path to the image if found.
func FindImage(name string) (string, error) {
	baseName := strings.TrimSuffix(strings.TrimSuffix(name, ".sif"), ".sqf")

	for _, dir := range GetImageSearchPaths() {
		candidate := filepath.Join(dir, name)
		if _, err := os.Stat(candidate); err == nil {
			return candidate, nil
		}
		candidate = filepath.Join(dir, baseName+".sif")
		if _, err := os.Stat(candidate); err == nil {
			return candidate, nil
		}
		candidate = filepath.Join(dir, baseName+".sqf")
		if _, err := os.Stat(candidate); err == nil {
			return candidate, nil
		}
	}

	return "", fmt.Errorf("image not found: %s (searched: %v)", name, GetImageSearchPaths())
}

// FindHelperScript searches all helper script paths for a script by name.
// Returns the full path to the script if found.
func FindHelperScript(name string) (string, error) {
	for _, dir := range GetHelperScriptSearchPaths() {
		candidate := filepath.Join(dir, name)
		if _, err := os.Stat(candidate); err == nil {
			return candidate, nil
		}
	}
	return "", fmt.Errorf("helper script not found: %s", name)
}

// FindBaseImage searches all image paths for the base image.
// Returns the full path if found, empty string otherwise.
func FindBaseImage() string {
	sifName := BaseImageSifName()
	for _, dir := range GetImageSearchPaths() {
		if candidate := filepath.Join(dir, sifName); fileExists(candidate) {
			return candidate
		}
	}
	return ""
}

// =============================================================================
// Utilities
// =============================================================================

// DirExists checks if a directory exists and is a directory.
func DirExists(path string) bool {
	info, err := os.Stat(path)
	if err != nil {
		return false
	}
	return info.IsDir()
}
