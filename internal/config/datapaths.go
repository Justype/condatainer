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

// DataPaths holds the search paths for data directories.
// Search order: extra → portable → scratch → user (first match wins for lookups).
// Write operations go to the first writable path in the same order.
type DataPaths struct {
	ImagesDirs        []string // Search paths for images
	BuildScriptsDirs  []string // Search paths for build scripts
	HelperScriptsDirs []string // Search paths for helper scripts
}

// GlobalDataPaths holds the computed data paths
var GlobalDataPaths DataPaths

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

// GetUserDataDir returns the user's data directory following XDG spec.
// Returns $XDG_DATA_HOME/condatainer or ~/.local/share/condatainer
func GetUserDataDir() string {
	// XDG_DATA_HOME takes priority
	if dataHome := os.Getenv("XDG_DATA_HOME"); dataHome != "" {
		return filepath.Join(dataHome, "condatainer")
	}

	// Default to ~/.local/share/condatainer
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
	// XDG_STATE_HOME takes priority
	if stateHome := os.Getenv("XDG_STATE_HOME"); stateHome != "" {
		return filepath.Join(stateHome, "condatainer")
	}

	// Default to ~/.local/state/condatainer
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
	// Resolve symlinks to get real path
	if realExe, err := filepath.EvalSymlinks(exe); err == nil {
		exeDir = filepath.Dir(realExe)
	}

	// Skip if executable is in standard user locations
	if isStandardUserPath(exeDir) {
		return ""
	}

	// Executable is in <install_dir>/bin/, data is in <install_dir>/
	// e.g., /project/group/condatainer/bin/condatainer → /project/group/condatainer/
	// e.g., $HOME/condatainer/bin/condatainer → $HOME/condatainer/
	// Skip generic bin dirs whose parent is $HOME, $HOME/.local, /usr, /usr/local, /opt
	if filepath.Base(exeDir) == "bin" {
		parentDir := filepath.Dir(exeDir)
		if !isNonPortableParent(parentDir) {
			return parentDir
		}
	}

	return ""
}

// IsPortable returns true if the current installation is portable (not in standard user locations).
// This is useful for determining if nested PATH updates should include the installation's bin directory.
func IsPortable() bool {
	return GetPortableDataDir() != ""
}

// GetPortableBinDir returns the bin directory of a portable installation, or empty string if not portable.
// This is typically used to add the installation's bin directory to PATH for nested usage.
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

// isStandardUserPath checks if a path is under standard user directories
// (home dir config/data locations that should not be treated as portable).
func isStandardUserPath(path string) bool {
	home, err := os.UserHomeDir()
	if err != nil {
		return false
	}

	// Standard XDG and home subdirectories
	standardPrefixes := []string{
		filepath.Join(home, ".config"),
		filepath.Join(home, ".local"),
		filepath.Join(home, ".cache"),
	}

	// Also check XDG environment variables
	if xdgConfig := os.Getenv("XDG_CONFIG_HOME"); xdgConfig != "" {
		standardPrefixes = append(standardPrefixes, xdgConfig)
	}
	if xdgData := os.Getenv("XDG_DATA_HOME"); xdgData != "" {
		standardPrefixes = append(standardPrefixes, xdgData)
	}

	for _, prefix := range standardPrefixes {
		if strings.HasPrefix(path, prefix) {
			return true
		}
	}
	return false
}

// isNonPortableParent checks if a directory is a common parent directory
// that should not be treated as a portable installation root.
// This includes $HOME, $HOME/.local, /usr, /usr/local, and /opt.
func isNonPortableParent(dir string) bool {
	home, _ := os.UserHomeDir()
	excludedParents := []string{
		home,                          // $HOME/bin
		filepath.Join(home, ".local"), // $HOME/.local/bin
		"/usr",                        // /usr/bin
		"/usr/local",                  // /usr/local/bin
		"/opt",                        // /opt/bin
	}

	return slices.Contains(excludedParents, dir)
}

// InitDataPaths initializes GlobalDataPaths based on environment and config.
// This should be called after LoadDefaults and LoadFromViper.
func InitDataPaths() {
	GlobalDataPaths = DataPaths{
		ImagesDirs:        buildImageSearchPaths(),
		BuildScriptsDirs:  buildScriptSearchPaths("build-scripts"),
		HelperScriptsDirs: buildScriptSearchPaths("helper-scripts"),
	}
}

// buildImageSearchPaths builds the search paths for images.
// Priority: extra_image_dirs → extra_base_dirs → portable → scratch → user
func buildImageSearchPaths() []string {
	var paths []string
	seen := make(map[string]bool)

	addPath := func(dir string) {
		if dir == "" {
			return
		}
		// Expand environment variables
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

	// 0a. Explicit extra image directories (highest priority, strip :ro/:rw marker)
	for _, entry := range GetExtraImageDirs() {
		path, _ := ParseDirEntry(entry)
		addPath(path)
	}

	// 0b. Extra base directories - auto-append /images
	for _, baseDir := range GetExtraBaseDirs() {
		addPath(filepath.Join(baseDir, "images"))
	}

	// 1. Portable directory (for group/shared use - preferred for both read and write)
	if portableDir := GetPortableDataDir(); portableDir != "" {
		addPath(filepath.Join(portableDir, "images"))
	}

	// 2. Scratch directory (HPC systems, large storage)
	if scratchDir := GetScratchDataDir(); scratchDir != "" {
		addPath(filepath.Join(scratchDir, "images"))
	}

	// 3. User directory (personal storage - fallback)
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
		// Expand environment variables
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

	// 0a. Explicit extra script directories (highest priority)
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

	// 0b. Extra base directories - auto-append subdir
	for _, baseDir := range GetExtraBaseDirs() {
		addPath(filepath.Join(baseDir, subdir))
	}

	// 1. Portable directory (for group/shared use - preferred)
	if portableDir := GetPortableDataDir(); portableDir != "" {
		addPath(filepath.Join(portableDir, subdir))

		// Auto-detect cnt-scripts subdir for build-scripts only (split-repo layout):
		// clone cnt-scripts next to the binary and it is found automatically
		if subdir == "build-scripts" && DirExists(filepath.Join(portableDir, "cnt-scripts")) {
			addPath(filepath.Join(portableDir, "cnt-scripts", subdir))
		}
	}

	// 2. Scratch directory (HPC systems)
	if scratchDir := GetScratchDataDir(); scratchDir != "" {
		addPath(filepath.Join(scratchDir, subdir))
	}

	// 3. User directory (personal storage - fallback)
	if userDir := GetUserDataDir(); userDir != "" {
		addPath(filepath.Join(userDir, subdir))
	}

	return paths
}

// buildCacheSearchPaths builds search paths for cache directories.
// Priority: extra_base_dirs → portable → scratch
// Note: Per-user cache is intentionally excluded for shared CondaTainer caches.
func buildCacheSearchPaths() []string {
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

	for _, baseDir := range GetExtraBaseDirs() {
		addPath(filepath.Join(baseDir, "cache"))
	}

	if portableDir := GetPortableDataDir(); portableDir != "" {
		addPath(filepath.Join(portableDir, "cache"))
	}

	if scratchDir := GetScratchDataDir(); scratchDir != "" {
		addPath(filepath.Join(scratchDir, "cache"))
	}

	return paths
}

// GetCacheSearchPaths returns all paths to search for cache directories.
// Includes the user cache dir as the final fallback.
func GetCacheSearchPaths() []string {
	paths := buildCacheSearchPaths()
	if userCacheDir := GetUserCacheDir(); userCacheDir != "" {
		seen := make(map[string]bool)
		for _, p := range paths {
			seen[p] = true
		}
		if !seen[userCacheDir] {
			paths = append(paths, userCacheDir)
		}
	}
	return paths
}

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

// FindImage searches all image paths for an image by name.
// Returns the full path to the image if found.
func FindImage(name string) (string, error) {
	// Normalize the name (remove .sif or .sqf extension for searching)
	baseName := strings.TrimSuffix(strings.TrimSuffix(name, ".sif"), ".sqf")

	for _, dir := range GetImageSearchPaths() {
		// Try exact name first
		candidate := filepath.Join(dir, name)
		if _, err := os.Stat(candidate); err == nil {
			return candidate, nil
		}

		// Try with .sif extension
		candidate = filepath.Join(dir, baseName+".sif")
		if _, err := os.Stat(candidate); err == nil {
			return candidate, nil
		}

		// Try with .sqf extension (overlay)
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

// GetWritableImagesDir returns the first writable images directory.
// Explicit extra_image_dirs entries marked :ro are skipped.
// Probes existing directories first (no side effects); creates only the last
// (user-owned) directory if none are currently writable.
func GetWritableImagesDir() (string, error) {
	// Phase 0: explicit extra_image_dirs — skip :ro entries
	for _, entry := range GetExtraImageDirs() {
		path, readOnly := ParseDirEntry(entry)
		if !readOnly && utils.CanWriteToDir(path) {
			return path, nil
		}
	}
	// Phase 1: remaining search paths (extra_base_dirs → portable → scratch → user)
	// Build a deduplicated list excluding explicit extra_image_dirs paths already checked
	extraSet := make(map[string]bool)
	for _, entry := range GetExtraImageDirs() {
		path, _ := ParseDirEntry(entry)
		if p, err := filepath.Abs(os.ExpandEnv(path)); err == nil {
			extraSet[p] = true
		}
	}
	paths := GetImageSearchPaths()
	for _, dir := range paths {
		if extraSet[dir] {
			continue // already checked in phase 0
		}
		if utils.CanWriteToDir(dir) {
			return dir, nil
		}
	}
	if len(paths) > 0 {
		last := paths[len(paths)-1]
		if !extraSet[last] && utils.IsWritableDir(last) {
			return last, nil
		}
	}
	return "", fmt.Errorf("no writable images directory found (searched: %v)", paths)
}

// GetWritableCacheDir returns the first writable cache directory.
// Probes existing directories first (no side effects). Falls back to the XDG
// user cache directory (~/.cache/condatainer) if no shared cache is writable,
// creating it on first use. This ensures cache always works even when the
// portable install dir is read-only and $SCRATCH is unset.
func GetWritableCacheDir() (string, error) {
	for _, dir := range buildCacheSearchPaths() {
		if utils.CanWriteToDir(dir) {
			return dir, nil
		}
	}
	if userCacheDir := GetUserCacheDir(); userCacheDir != "" {
		if utils.IsWritableDir(userCacheDir) {
			return userCacheDir, nil
		}
	}
	return "", fmt.Errorf("no writable shared cache directory found")
}

// GetWritableHelperScriptsDir returns the first writable helper scripts directory.
// Explicit extra_helper_dirs entries marked :ro are skipped.
// Probes existing directories first; creates only the last (user-owned) directory
// if none are currently writable.
func GetWritableHelperScriptsDir() (string, error) {
	// Phase 0: explicit extra_helper_dirs — skip :ro entries
	for _, entry := range GetExtraHelperDirs() {
		path, readOnly := ParseDirEntry(entry)
		if !readOnly && utils.CanWriteToDir(path) {
			return path, nil
		}
	}
	// Phase 1: remaining search paths, excluding extra_helper_dirs already checked
	extraSet := make(map[string]bool)
	for _, entry := range GetExtraHelperDirs() {
		path, _ := ParseDirEntry(entry)
		if p, err := filepath.Abs(os.ExpandEnv(path)); err == nil {
			extraSet[p] = true
		}
	}
	paths := GetHelperScriptSearchPaths()
	for _, dir := range paths {
		if extraSet[dir] {
			continue
		}
		if utils.CanWriteToDir(dir) {
			return dir, nil
		}
	}
	if len(paths) > 0 {
		last := paths[len(paths)-1]
		if !extraSet[last] && utils.IsWritableDir(last) {
			return last, nil
		}
	}
	return "", fmt.Errorf("no writable helper scripts directory found")
}

// ImageInfo holds information about an image found during enumeration.
type ImageInfo struct {
	Name   string // Filename (e.g., "python.sqf")
	Path   string // Full absolute path
	Source string // Source directory type: "extra", "portable", "scratch", "user", or "legacy"
}

// ScriptFileInfo holds information about a script found during enumeration.
type ScriptFileInfo struct {
	Name   string // Filename
	Path   string // Full absolute path
	Source string // Source directory type
}

// // categorizeSource determines the source type of a directory based on its index and path.
// func categorizeSource(dir string, index int) string {
// 	portableDir := GetPortableDataDir()
// 	scratchDir := GetScratchDataDir()
// 	userDir := GetUserDataDir()

// 	if portableDir != "" && strings.HasPrefix(dir, portableDir) {
// 		return "portable"
// 	}
// 	if scratchDir != "" && strings.HasPrefix(dir, scratchDir) {
// 		return "scratch"
// 	}
// 	if userDir != "" && strings.HasPrefix(dir, userDir) {
// 		return "user"
// 	}

// 	// Fallback based on index (matching the search order)
// 	switch index {
// 	case 0:
// 		return "extra"
// 	case 1:
// 		return "portable"
// 	case 2:
// 		return "scratch"
// 	case 3:
// 		return "user"
// 	default:
// 		return "legacy"
// 	}
// }

// DirExists checks if a directory exists and is a directory.
func DirExists(path string) bool {
	info, err := os.Stat(path)
	if err != nil {
		return false
	}
	return info.IsDir()
}

// BaseImageDefName is the filename for the base image definition file.
const BaseImageDefName = "base_image.def"

// FindBaseImage searches all image paths for the base image.
// Returns the full path if found, empty string otherwise.
// slug-named .sif (e.g. "ubuntu24--base_image.sif")
func FindBaseImage() string {
	sifName := BaseImageSifName()
	for _, dir := range GetImageSearchPaths() {
		if candidate := filepath.Join(dir, sifName); fileExists(candidate) {
			return candidate
		}
	}
	return ""
}

// GetBaseImageWritePath returns the path where a new base image should be written.
// Returns the first writable images directory + slug-based .sif name (e.g. ubuntu24--base_image.sif).
func GetBaseImageWritePath() (string, error) {
	dir, err := GetWritableImagesDir()
	if err != nil {
		return "", err
	}
	return filepath.Join(dir, BaseImageSifName()), nil
}
