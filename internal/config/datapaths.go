package config

import (
	"encoding/json"
	"fmt"
	"os"
	"path/filepath"
	"slices"
	"strconv"
	"strings"
	"sync"

	"github.com/Justype/condatainer/internal/utils"
	"github.com/spf13/viper"
)

// =============================================================================
// Types & Globals
// =============================================================================

// DataPaths holds the search paths for data directories.
// Read order: scratch → user → extra-root → root (first match wins for lookups).
// Write order is the opposite — extra-root → root → scratch → user, first writable —
// so reads resolve nearest and writes land furthest out. See searchPaths.
type DataPaths struct {
	ImagesDirs        []string // Search paths for images
	HelperScriptsDirs []string // Search paths for helper scripts
}

// GlobalDataPaths holds the computed data paths.
var GlobalDataPaths DataPaths

//=============================================================================
// Config Getters (env/viper)
// =============================================================================

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
// The env var name is derived from the key: "channels" → "CNT_CHANNELS". The split
// function controls how the env var value is parsed. An env var replaces the config
// value outright.
func getEnvSlice(key string, split func(string) []string) []string {
	envKey := "CNT_" + strings.ToUpper(key)
	if envVal := os.Getenv(envKey); envVal != "" {
		if vals := split(envVal); len(vals) > 0 {
			return vals
		}
	}
	return viper.GetStringSlice(key)
}

// GetChannels returns the conda channels from config or environment.
// CNT_CHANNELS supports "|" and ":" as separators ("|" takes precedence).
// Overwrite semantics: if set in config, highest-priority layer wins (no merge).
func GetChannels() []string { return getEnvSlice("channels", splitPipeOrColon) }

// GetExtraRootDir returns the extra root directory from CNT_EXTRA_ROOT env var.
// Single value, env only (no config key). Used for group/lab shared installations.
// Returns empty string if not set. Cached via sync.Once.
func GetExtraRootDir() string {
	extraRootOnce.Do(func() {
		root := os.Getenv("CNT_EXTRA_ROOT")
		if root == "" {
			return
		}
		root = os.ExpandEnv(root)
		if abs, err := filepath.Abs(root); err == nil {
			extraRootCache = abs
		} else {
			extraRootCache = root
		}
	})
	return extraRootCache
}

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

// GetUserConfigDir returns the user's config directory following XDG spec.
// Returns $XDG_CONFIG_HOME/condatainer or ~/.config/condatainer
func GetUserConfigDir() string {
	if configHome := os.Getenv("XDG_CONFIG_HOME"); configHome != "" {
		return filepath.Join(configHome, "condatainer")
	}
	if home, err := os.UserHomeDir(); err == nil {
		return filepath.Join(home, ".config", "condatainer")
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

// GetRootDir returns the installation root directory.
// CNT_ROOT env var takes priority; otherwise auto-detected from the executable location.
//
// Detection logic:
//  1. CNT_ROOT env var — explicit root (set by sysadmin in module file)
//  2. Executable in <dir>/bin/ and <dir> is not a standard user/system parent → use <dir>
//
// Cached via sync.Once.
func GetRootDir() string {
	rootDirOnce.Do(func() {
		rootDirCache = detectRootDir()
	})
	return rootDirCache
}

// =============================================================================
// Root Detection (internal)
// =============================================================================

var (
	rootDirOnce    sync.Once
	rootDirCache   string
	extraRootOnce  sync.Once
	extraRootCache string
)

func detectRootDir() string {
	// Explicit override: CNT_ROOT sets the root dir directly,
	// bypassing the executable-location heuristic below.
	if root := os.Getenv("CNT_ROOT"); root != "" {
		root = os.ExpandEnv(root)
		if abs, err := filepath.Abs(root); err == nil {
			return abs
		}
		return root
	}

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
		if !isNonRootParent(parentDir) {
			return parentDir
		}
	}

	return ""
}

// isNonRootParent returns true for standard user/system directories that should
// not be treated as an installation root.
func isNonRootParent(dir string) bool {
	home, _ := os.UserHomeDir()
	excludedParents := []string{
		home,
		filepath.Join(home, ".local"),
		"/usr",
		"/usr/local",
		"/opt",
		"/",
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
		HelperScriptsDirs: helperScriptSearchPaths(),
	}
}

// =============================================================================
// Search Path Builders (internal, read)
// =============================================================================

// searchPaths builds the read search order for subdir across every data tier:
// personal (scratch → user) before shared (CNT_EXTRA_ROOT → root).
//
// Reads go nearest-first and writes go furthest-first (see imageWriteDirs), so a
// build lands as far out as permissions allow — one copy serving the whole group —
// while the person who wants their own version of a name just builds it into their
// own directory and has it win for them. Order *within* a tier is the same in both
// directions: scratch before user, and a group's own extra-root before a site root.
func searchPaths(subdir string) []string {
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

	if scratchDir := GetScratchDataDir(); scratchDir != "" {
		addPath(filepath.Join(scratchDir, subdir))
	}
	if userDir := GetUserDataDir(); userDir != "" {
		addPath(filepath.Join(userDir, subdir))
	}
	if extraRoot := GetExtraRootDir(); extraRoot != "" {
		addPath(filepath.Join(extraRoot, subdir))
	}
	if rootDir := GetRootDir(); rootDir != "" {
		addPath(filepath.Join(rootDir, subdir))
	}

	return paths
}

// buildImageSearchPaths builds the search paths for images.
// Priority: scratch → user → CNT_EXTRA_ROOT → root
func buildImageSearchPaths() []string { return searchPaths("images") }

// helperScriptSearchPaths builds the search paths for helper scripts.
// Priority: scratch → user → CNT_EXTRA_ROOT → root
//
// Recipes have no equivalent: they come from the catalog's `sources`, not from a
// data directory.
func helperScriptSearchPaths() []string { return searchPaths("helper-scripts") }

// =============================================================================
// Data Layers
// =============================================================================

// DataLayer identifies which tier a data directory belongs to. The short forms
// (u/r/e) match the -l/--layer vocabulary used by the config commands, so a
// single mental model covers both config files and data directories.
type DataLayer string

const (
	LayerExtraRoot DataLayer = "extra-root" // $CNT_EXTRA_ROOT           (e)
	LayerAppRoot   DataLayer = "app-root"   // $CNT_ROOT or <install-dir> (r)
	LayerUser      DataLayer = "user"       // $SCRATCH, else XDG data dir (u)
	LayerUnknown   DataLayer = "unknown"    // not under any tier root
)

// ParseDataLayer resolves a user-supplied layer name or its short form.
func ParseDataLayer(s string) (DataLayer, error) {
	switch strings.ToLower(strings.TrimSpace(s)) {
	case "u", "user":
		return LayerUser, nil
	case "r", "app-root", "approot", "root":
		return LayerAppRoot, nil
	case "e", "extra-root", "extraroot", "extra":
		return LayerExtraRoot, nil
	}
	return "", fmt.Errorf("unknown layer %q (want: u/user, r/app-root, e/extra-root)", s)
}

// ClassifyDataDir reports which layer a data directory belongs to, by comparing
// its parent against each tier root. A directory shared by two tiers (e.g.
// CNT_ROOT == $SCRATCH/condatainer) is labelled with the shared one, which is what
// decides whether writes may land there — not with whichever tier reads it first.
// Every search path belongs to a tier; LayerUnknown is for a path from elsewhere.
func ClassifyDataDir(path string) DataLayer {
	parent := filepath.Clean(filepath.Dir(filepath.Clean(path)))
	match := func(root string) bool {
		return root != "" && filepath.Clean(root) == parent
	}
	switch {
	case match(GetExtraRootDir()):
		return LayerExtraRoot
	case match(GetRootDir()):
		return LayerAppRoot
	case match(GetScratchDataDir()), match(GetUserDataDir()):
		return LayerUser
	}
	return LayerUnknown
}

// FilterDirsByLayer returns the subset of dirs belonging to layer, preserving order.
func FilterDirsByLayer(dirs []string, layer DataLayer) []string {
	var out []string
	for _, d := range dirs {
		if ClassifyDataDir(d) == layer {
			out = append(out, d)
		}
	}
	return out
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

// GetHelperScriptSearchPaths returns all paths to search for helper scripts.
func GetHelperScriptSearchPaths() []string {
	if len(GlobalDataPaths.HelperScriptsDirs) == 0 {
		InitDataPaths()
	}
	return GlobalDataPaths.HelperScriptsDirs
}

// GetCacheSearchPaths returns all personal cache directories to search.
// Cache is always the personal XDG cache dir to avoid cross-user pollution.
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
// Personal dirs (scratch, XDG user) are always created on first use via EnsureWritableDir.
// Shared dirs (CNT_EXTRA_ROOT, root, extra_*_dirs): subdirs are created if the parent
// exists, but the parent itself is never auto-created.
type SearchDir struct {
	Path     string
	Personal bool
}

// deduplicateWriteDirs removes duplicate paths, upgrading Personal=false → true when
// the same path appears in both a shared tier (root) and a personal tier (scratch).
// This handles the case where CNT_ROOT or root dir equals $SCRATCH/condatainer.
func deduplicateWriteDirs(dirs []SearchDir) []SearchDir {
	seen := make(map[string]int) // path → index in result
	result := make([]SearchDir, 0, len(dirs))
	for _, d := range dirs {
		if idx, ok := seen[d.Path]; ok {
			if d.Personal && !result[idx].Personal {
				result[idx].Personal = true
			}
		} else {
			seen[d.Path] = len(result)
			result = append(result, d)
		}
	}
	return result
}

// firstWritableDir returns the first writable path from the slice.
// Personal dirs are always created on first use (EnsureWritableDir).
// Shared dirs are created only if the parent directory already exists — the parent
// itself is never auto-created, but subdirs (images/, helper-scripts/) are.
func firstWritableDir(dirs []SearchDir) string {
	for _, d := range dirs {
		if d.Personal {
			if utils.EnsureWritableDir(d.Path) {
				return d.Path
			}
		} else {
			if utils.DirExists(filepath.Dir(d.Path)) && utils.EnsureWritableDir(d.Path) {
				return d.Path
			}
		}
	}
	return ""
}

// imageWriteDirs returns the ordered write candidates for image directories.
// Shared: CNT_EXTRA_ROOT, root. Personal: scratch, user.
func imageWriteDirs() []SearchDir {
	var dirs []SearchDir
	if extraRoot := GetExtraRootDir(); extraRoot != "" {
		dirs = append(dirs, SearchDir{Path: filepath.Join(extraRoot, "images")})
	}
	if p := GetRootDir(); p != "" {
		dirs = append(dirs, SearchDir{Path: filepath.Join(p, "images")})
	}
	if s := GetScratchDataDir(); s != "" {
		dirs = append(dirs, SearchDir{Path: filepath.Join(s, "images"), Personal: true})
	}
	if u := GetUserDataDir(); u != "" {
		dirs = append(dirs, SearchDir{Path: filepath.Join(u, "images"), Personal: true})
	}
	return deduplicateWriteDirs(dirs)
}

// helperWriteDirs returns the ordered write candidates for helper-scripts directories.
// Shared: CNT_EXTRA_ROOT, root. Personal: scratch, user.
func helperWriteDirs() []SearchDir {
	var dirs []SearchDir
	if extraRoot := GetExtraRootDir(); extraRoot != "" {
		dirs = append(dirs, SearchDir{Path: filepath.Join(extraRoot, "helper-scripts")})
	}
	if p := GetRootDir(); p != "" {
		dirs = append(dirs, SearchDir{Path: filepath.Join(p, "helper-scripts")})
	}
	if s := GetScratchDataDir(); s != "" {
		dirs = append(dirs, SearchDir{Path: filepath.Join(s, "helper-scripts"), Personal: true})
	}
	if u := GetUserDataDir(); u != "" {
		dirs = append(dirs, SearchDir{Path: filepath.Join(u, "helper-scripts"), Personal: true})
	}
	return deduplicateWriteDirs(dirs)
}

// cacheWriteDirs returns the write candidates for cache directories.
// Cache always lives in the personal XDG cache dir ($XDG_CACHE_HOME/condatainer
// or ~/.cache/condatainer) — never in data directories or shared roots.
func cacheWriteDirs() []SearchDir {
	if u := GetUserCacheDir(); u != "" {
		return []SearchDir{{Path: u, Personal: true}}
	}
	return nil
}

// =============================================================================
// Write Dir Resolvers (public)
// =============================================================================

// GetWritableImagesDir returns the first writable images directory.
// Shared dirs (CNT_EXTRA_ROOT, root) are probed only.
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
// Shared dirs (CNT_EXTRA_ROOT, root) are probed only.
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

// GetWritableCacheDir returns the writable personal cache directory
// ($XDG_CACHE_HOME/condatainer or ~/.cache/condatainer) — never a shared dir.
func GetWritableCacheDir() (string, error) {
	dirs := cacheWriteDirs()
	if dir := firstWritableDir(dirs); dir != "" {
		return dir, nil
	}
	return "", fmt.Errorf("no writable cache directory found")
}

// GetBaseImageWritePath returns the path where a new base image should be written.
func GetBaseImageWritePath() (string, error) {
	sifName := BaseImageSifName()
	if sifName == "" {
		return "", fmt.Errorf("no base configured: set `base` in config")
	}
	dir, err := GetWritableImagesDir()
	if err != nil {
		return "", err
	}
	return filepath.Join(dir, sifName), nil
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
	if sifName == "" {
		return ""
	}
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

// GetHelperStateDir returns the per-ID helper state directory on NFS.
// Each helper run writes its ready/messages/done files under this path.
// Path: ~/.local/state/condatainer/helper/{id}/
func GetHelperStateDir(id string) string {
	stateDir := GetUserStateDir()
	if stateDir == "" {
		return ""
	}
	return filepath.Join(stateDir, "helper", id)
}

// GetHelperHistoryPath returns the path to the JSONL helper run history file.
// Path: ~/.local/state/condatainer/helper-history.jsonl
func GetHelperHistoryPath() string {
	stateDir := GetUserStateDir()
	if stateDir == "" {
		return ""
	}
	return filepath.Join(stateDir, "helper-history.jsonl")
}

// GetFileBookmarksPath returns the path to the JSON file storing the dashboard's
// file-browser bookmarks.
// Path: ~/.local/state/condatainer/file-bookmarks.json
func GetFileBookmarksPath() string {
	stateDir := GetUserStateDir()
	if stateDir == "" {
		return ""
	}
	return filepath.Join(stateDir, "file-bookmarks.json")
}

// GetHelperBookmarksPath returns the path to the JSON file storing the
// dashboard's bookmarked (starred) helper names.
// Path: ~/.local/state/condatainer/helper-bookmarks.json
func GetHelperBookmarksPath() string {
	stateDir := GetUserStateDir()
	if stateDir == "" {
		return ""
	}
	return filepath.Join(stateDir, "helper-bookmarks.json")
}

// ServerState is the JSON written to the server PID file (NFS-visible).
type ServerState struct {
	Host   string `json:"host"`
	Port   int    `json:"port"`
	Token  string `json:"token"`
	PID    int    `json:"pid"`
	Daemon bool   `json:"daemon"` // Saved for 'server restart' using same daemon mode.
}

// GetServerPidFilePath returns the path to the server PID/state file for the current host.
// Path: ~/.local/state/condatainer/server-{hostname}.pid
// Per-host files allow NFS home dirs shared across login nodes without conflicts.
func GetServerPidFilePath() string {
	stateDir := GetUserStateDir()
	if stateDir == "" {
		return ""
	}
	host, _ := os.Hostname()
	if host == "" {
		host = "local"
	}
	return filepath.Join(stateDir, "server-"+host+".pid")
}

// GetServerLogPath returns the path to the server log file for the current host.
// Path: ~/.local/state/condatainer/server-{hostname}.log
func GetServerLogPath() string {
	host, _ := os.Hostname()
	if host == "" {
		host = "local"
	}
	return GetServerLogPathForHost(host)
}

// GetServerLogPathForHost returns the path to the server log file for the given host.
// Path: ~/.local/state/condatainer/server-{host}.log
// Since the state dir is on NFS, this works for any login node.
func GetServerLogPathForHost(host string) string {
	stateDir := GetUserStateDir()
	if stateDir == "" {
		return ""
	}
	return filepath.Join(stateDir, "server-"+host+".log")
}

// GetRunningServerPort returns the port of the currently running server on this host,
// reading it directly from the PID file. Returns 0 if the server is not running.
func GetRunningServerPort() int {
	pidFile := GetServerPidFilePath()
	if pidFile == "" {
		return 0
	}
	data, err := os.ReadFile(pidFile)
	if err != nil {
		return 0
	}
	var ss ServerState
	if err := json.Unmarshal(data, &ss); err != nil {
		return 0
	}
	return ss.Port
}

// GetServerSavedPortPath returns the path to the shared saved-port file.
// Path: ~/.local/state/condatainer/server.port
// Written by "server start --port" so the same port is reused on subsequent starts
// across all login nodes.
func GetServerSavedPortPath() string {
	stateDir := GetUserStateDir()
	if stateDir == "" {
		return ""
	}
	return filepath.Join(stateDir, "server.port")
}

// ReadSavedServerPort reads the port saved by a previous "server start --port" call.
// Returns 0 if no port has been saved or the file cannot be read.
func ReadSavedServerPort() int {
	path := GetServerSavedPortPath()
	if path == "" {
		return 0
	}
	data, err := os.ReadFile(path)
	if err != nil {
		return 0
	}
	n, err := strconv.Atoi(strings.TrimSpace(string(data)))
	if err != nil || n <= 0 {
		return 0
	}
	return n
}

// SaveServerPort writes port to the shared saved-port file so it is reused on
// subsequent "server start" calls without --port.
func SaveServerPort(port int) error {
	path := GetServerSavedPortPath()
	if path == "" {
		return fmt.Errorf("cannot determine state directory")
	}
	dir := filepath.Dir(path)
	if err := utils.MkdirAllShared(dir); err != nil {
		return err
	}
	if err := os.WriteFile(path, []byte(strconv.Itoa(port)+"\n"), utils.PermFile); err != nil {
		return err
	}
	utils.ShareWithParentGroup(path)
	return nil
}

// ListServerPidFiles returns paths to all server-*.pid files in the user state directory.
// Each file corresponds to a login node that may have (or have had) a server running.
// On NFS home dirs shared across login nodes, this reveals servers on other nodes.
func ListServerPidFiles() []string {
	stateDir := GetUserStateDir()
	if stateDir == "" {
		return nil
	}
	entries, err := os.ReadDir(stateDir)
	if err != nil {
		return nil
	}
	var files []string
	for _, e := range entries {
		if !e.IsDir() && strings.HasPrefix(e.Name(), "server-") && strings.HasSuffix(e.Name(), ".pid") {
			files = append(files, filepath.Join(stateDir, e.Name()))
		}
	}
	return files
}
