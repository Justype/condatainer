package config

import (
	"fmt"
	"os"
	"path/filepath"
	"strings"

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

// GetExtraBaseDirs returns extra base directories from config or environment.
// Environment variable CONDATAINER_EXTRA_BASE_DIRS uses colon-separated paths (Unix convention).
// Config file uses YAML array format.
func GetExtraBaseDirs() []string {
	// Check environment variable first (colon-separated, like PATH)
	if envDirs := os.Getenv("CONDATAINER_EXTRA_BASE_DIRS"); envDirs != "" {
		var dirs []string
		for _, dir := range strings.Split(envDirs, ":") {
			dir = strings.TrimSpace(dir)
			if dir != "" {
				dirs = append(dirs, dir)
			}
		}
		if len(dirs) > 0 {
			return dirs
		}
	}

	// Fall back to viper config
	return viper.GetStringSlice("extra_base_dirs")
}

// GetSystemDataDir returns the system-wide data directory.
// Checks XDG_DATA_DIRS first, then standard paths.
func GetSystemDataDir() string {
	// Check XDG_DATA_DIRS first (colon-separated list)
	if dataDirs := os.Getenv("XDG_DATA_DIRS"); dataDirs != "" {
		for _, dir := range strings.Split(dataDirs, ":") {
			if dir == "" {
				continue
			}
			candidate := filepath.Join(dir, "condatainer")
			if stat, err := os.Stat(candidate); err == nil && stat.IsDir() {
				return candidate
			}
		}
	}

	// Standard fallbacks
	standardPaths := []string{
		"/usr/local/share/condatainer",
		"/usr/share/condatainer",
		"/opt/condatainer",
	}
	for _, dir := range standardPaths {
		if stat, err := os.Stat(dir); err == nil && stat.IsDir() {
			return dir
		}
	}

	return ""
}

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
//  1. If executable is in <dir>/bin/, check if <dir> has images/ subdirectory
//  2. Otherwise check for <exeDir>/condatainer/ directory
func GetPortableDataDir() string {
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

	// Case 1: Executable is in <install_dir>/bin/, data is in <install_dir>/
	// e.g., /project/group/condatainer/bin/condatainer → /project/group/condatainer/
	if filepath.Base(exeDir) == "bin" {
		parentDir := filepath.Dir(exeDir)
		// Check if parent has images/ subdirectory (indicator of portable install)
		if stat, err := os.Stat(filepath.Join(parentDir, "images")); err == nil && stat.IsDir() {
			return parentDir
		}
	}

	// Case 2: Check for sibling condatainer/ directory
	// e.g., /opt/bin/condatainer → /opt/bin/condatainer/
	candidate := filepath.Join(exeDir, "condatainer")
	if stat, err := os.Stat(candidate); err == nil && stat.IsDir() {
		return candidate
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
	if portableDir := GetPortableDataDir(); portableDir != "" {
		binDir := filepath.Join(portableDir, "bin")
		if stat, err := os.Stat(binDir); err == nil && stat.IsDir() {
			return binDir
		}
	}
	return ""
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

	for _, excluded := range excludedParents {
		if dir == excluded {
			return true
		}
	}
	return false
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
// Priority: extra_base_dirs → portable → scratch → user
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

	// 0. Extra base directories (highest priority) - auto-append /images
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
// Priority: extra_base_dirs → portable → scratch → user
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

	// 0. Extra base directories (highest priority) - auto-append subdir
	for _, baseDir := range GetExtraBaseDirs() {
		addPath(filepath.Join(baseDir, subdir))
	}

	// 1. Portable directory (for group/shared use - preferred)
	if portableDir := GetPortableDataDir(); portableDir != "" {
		addPath(filepath.Join(portableDir, subdir))
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

// FindBuildScriptFile searches all build script paths for a script by name.
// Returns the full path to the script if found.
func FindBuildScriptFile(name string) (string, error) {
	for _, dir := range GetBuildScriptSearchPaths() {
		// Try exact name
		candidate := filepath.Join(dir, name)
		if _, err := os.Stat(candidate); err == nil {
			return candidate, nil
		}

		// Try with .def extension
		candidate = filepath.Join(dir, name+".def")
		if _, err := os.Stat(candidate); err == nil {
			return candidate, nil
		}
	}

	return "", fmt.Errorf("build script not found: %s", name)
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
// Creates the directory if it doesn't exist.
func GetWritableImagesDir() (string, error) {
	for _, dir := range GetImageSearchPaths() {
		// Try to create directory if it doesn't exist
		if err := os.MkdirAll(dir, utils.PermDir); err != nil {
			continue
		}

		// Test write permission
		testFile := filepath.Join(dir, ".write-test")
		f, err := os.Create(testFile)
		if err != nil {
			continue
		}
		f.Close()
		os.Remove(testFile)
		return dir, nil
	}

	return "", fmt.Errorf("no writable images directory found (searched: %v)", GetImageSearchPaths())
}

// GetWritableBuildScriptsDir returns the first writable build scripts directory.
func GetWritableBuildScriptsDir() (string, error) {
	for _, dir := range GetBuildScriptSearchPaths() {
		if err := os.MkdirAll(dir, utils.PermDir); err != nil {
			continue
		}

		testFile := filepath.Join(dir, ".write-test")
		f, err := os.Create(testFile)
		if err != nil {
			continue
		}
		f.Close()
		os.Remove(testFile)
		return dir, nil
	}

	return "", fmt.Errorf("no writable build scripts directory found")
}

// GetWritableHelperScriptsDir returns the first writable helper scripts directory.
func GetWritableHelperScriptsDir() (string, error) {
	for _, dir := range GetHelperScriptSearchPaths() {
		if err := os.MkdirAll(dir, utils.PermDir); err != nil {
			continue
		}

		testFile := filepath.Join(dir, ".write-test")
		f, err := os.Create(testFile)
		if err != nil {
			continue
		}
		f.Close()
		os.Remove(testFile)
		return dir, nil
	}

	return "", fmt.Errorf("no writable helper scripts directory found")
}

// ImageInfo holds information about an image found during enumeration.
type ImageInfo struct {
	Name   string // Filename (e.g., "python.sqf")
	Path   string // Full absolute path
	Source string // Source directory type: "extra", "portable", "scratch", "user", or "legacy"
}

// ListAllImages returns all images from all search paths.
// Images in higher-priority directories shadow those in lower-priority ones.
func ListAllImages() ([]ImageInfo, error) {
	seen := make(map[string]bool)
	var results []ImageInfo

	for i, dir := range GetImageSearchPaths() {
		entries, err := os.ReadDir(dir)
		if err != nil {
			if os.IsNotExist(err) {
				continue
			}
			// Log but continue with other directories
			continue
		}

		source := categorizeSource(dir, i)

		for _, e := range entries {
			if e.IsDir() {
				continue
			}

			name := e.Name()
			// Skip non-image files
			if !strings.HasSuffix(name, ".sif") && !strings.HasSuffix(name, ".sqf") {
				continue
			}

			if seen[name] {
				continue // Already found in higher-priority dir
			}
			seen[name] = true

			results = append(results, ImageInfo{
				Name:   name,
				Path:   filepath.Join(dir, name),
				Source: source,
			})
		}
	}

	return results, nil
}

// ScriptFileInfo holds information about a script found during enumeration.
type ScriptFileInfo struct {
	Name   string // Filename
	Path   string // Full absolute path
	Source string // Source directory type
}

// ListAllBuildScripts returns all build scripts from all search paths.
func ListAllBuildScripts() ([]ScriptFileInfo, error) {
	seen := make(map[string]bool)
	var results []ScriptFileInfo

	for i, dir := range GetBuildScriptSearchPaths() {
		source := categorizeSource(dir, i)

		err := filepath.Walk(dir, func(path string, info os.FileInfo, err error) error {
			if err != nil {
				if os.IsNotExist(err) {
					return nil
				}
				return nil // Continue with other files
			}

			if info.IsDir() {
				return nil
			}

			name := info.Name()

			// Skip helper scripts (.py and .sh files) and templates
			if strings.HasSuffix(name, ".py") || strings.HasSuffix(name, ".sh") {
				return nil
			}
			if strings.Contains(path, "template") {
				return nil
			}

			// Use relative path from dir as the key
			relPath, _ := filepath.Rel(dir, path)
			if relPath == "" {
				relPath = name
			}

			if seen[relPath] {
				return nil // Already found in higher-priority dir
			}
			seen[relPath] = true

			results = append(results, ScriptFileInfo{
				Name:   relPath,
				Path:   path,
				Source: source,
			})

			return nil
		})

		if err != nil && !os.IsNotExist(err) {
			continue
		}
	}

	return results, nil
}

// ListAllHelperScripts returns all helper scripts from all search paths.
func ListAllHelperScripts() ([]ScriptFileInfo, error) {
	seen := make(map[string]bool)
	var results []ScriptFileInfo

	for i, dir := range GetHelperScriptSearchPaths() {
		entries, err := os.ReadDir(dir)
		if err != nil {
			if os.IsNotExist(err) {
				continue
			}
			continue
		}

		source := categorizeSource(dir, i)

		for _, e := range entries {
			if e.IsDir() {
				continue
			}

			name := e.Name()
			// Skip hidden files
			if strings.HasPrefix(name, ".") {
				continue
			}

			if seen[name] {
				continue
			}
			seen[name] = true

			results = append(results, ScriptFileInfo{
				Name:   name,
				Path:   filepath.Join(dir, name),
				Source: source,
			})
		}
	}

	return results, nil
}

// categorizeSource determines the source type of a directory based on its index and path.
func categorizeSource(dir string, index int) string {
	portableDir := GetPortableDataDir()
	scratchDir := GetScratchDataDir()
	userDir := GetUserDataDir()

	if portableDir != "" && strings.HasPrefix(dir, portableDir) {
		return "portable"
	}
	if scratchDir != "" && strings.HasPrefix(dir, scratchDir) {
		return "scratch"
	}
	if userDir != "" && strings.HasPrefix(dir, userDir) {
		return "user"
	}

	// Fallback based on index (matching the search order)
	switch index {
	case 0:
		return "extra"
	case 1:
		return "portable"
	case 2:
		return "scratch"
	case 3:
		return "user"
	default:
		return "legacy"
	}
}

// DirExists checks if a directory exists and is a directory.
func DirExists(path string) bool {
	info, err := os.Stat(path)
	if err != nil {
		return false
	}
	return info.IsDir()
}

// GetExistingImagePaths returns only the image search paths that actually exist.
func GetExistingImagePaths() []string {
	var existing []string
	for _, dir := range GetImageSearchPaths() {
		if DirExists(dir) {
			existing = append(existing, dir)
		}
	}
	return existing
}

// GetExistingBuildScriptPaths returns only the build script search paths that actually exist.
func GetExistingBuildScriptPaths() []string {
	var existing []string
	for _, dir := range GetBuildScriptSearchPaths() {
		if DirExists(dir) {
			existing = append(existing, dir)
		}
	}
	return existing
}

// GetExistingHelperScriptPaths returns only the helper script search paths that actually exist.
func GetExistingHelperScriptPaths() []string {
	var existing []string
	for _, dir := range GetHelperScriptSearchPaths() {
		if DirExists(dir) {
			existing = append(existing, dir)
		}
	}
	return existing
}

// BaseImageDefName is the filename for the base image definition file.
const BaseImageDefName = "base_image.def"

// FindBaseImage searches all image paths for the base image.
// Returns the full path if found, empty string otherwise.
// slug-named .sqf (e.g. "ubuntu24--base_image.sqf")
func FindBaseImage() string {
	sqfName := BaseImageSqfName()
	for _, dir := range GetImageSearchPaths() {
		if candidate := filepath.Join(dir, sqfName); fileExists(candidate) {
			return candidate
		}
	}
	return ""
}

// GetBaseImageWritePath returns the path where a new base image should be written.
// Returns the first writable images directory + slug-based .sqf name (e.g. ubuntu24--base_image.sqf).
func GetBaseImageWritePath() (string, error) {
	dir, err := GetWritableImagesDir()
	if err != nil {
		return "", err
	}
	return filepath.Join(dir, BaseImageSqfName()), nil
}

// FindBaseImageDef searches all build-scripts paths for the base image definition file.
// Prefers <slug>/base_image.def (e.g. ubuntu24/base_image.def), falls back to base_image.def.
// Returns the full path if found, empty string otherwise.
func FindBaseImageDef() string {
	slug := Global.DefaultDistro
	if slug == "" {
		return ""
	}
	for _, dir := range GetBuildScriptSearchPaths() {
		// slug subdirectory: ubuntu24/base_image.def
		if candidate := filepath.Join(dir, slug, BaseImageDefName); fileExists(candidate) {
			return candidate
		}
	}
	return ""
}

// GetBaseImageDefWritePath returns the path where a new base image def should be written.
// This is the first writable build-scripts directory + base_image.def
func GetBaseImageDefWritePath() (string, error) {
	dir, err := GetWritableBuildScriptsDir()
	if err != nil {
		return "", err
	}
	return filepath.Join(dir, BaseImageDefName), nil
}
