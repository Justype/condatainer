package build

import (
	"compress/gzip"
	"crypto/sha256"
	"encoding/json"
	"fmt"
	"io"
	"net/http"
	"os"
	"path/filepath"
	"reflect"
	"strings"
	"time"

	"github.com/Justype/condatainer/internal/config"
	"github.com/Justype/condatainer/internal/utils"
)

// ForceRefresh skips the on-disk metadata cache and fetches fresh data from remote.
// Set by CLI commands (e.g. `update` command).
var ForceRefresh bool

// RemoteMetadataCache is the on-disk envelope for the cached remote build script metadata.
type RemoteMetadataCache struct {
	FetchedAt    time.Time                    `json:"fetched_at"`
	SourceURL    string                       `json:"source_url"`
	PrebuiltLink string                       `json:"prebuilt_link,omitempty"`
	Metadata     map[string]RemoteScriptEntry `json:"metadata"`
}

type localScriptsDirFingerprint struct {
	Exists        bool           `json:"exists"`
	ModTimeUnix   int64          `json:"mod_time_unix"`
	ImmediateDirs map[string]int `json:"immediate_dirs,omitempty"`
}

type LocalScriptsCache struct {
	FetchedAt    time.Time                             `json:"fetched_at"`
	SearchPaths  []string                              `json:"search_paths"`
	Fingerprints map[string]localScriptsDirFingerprint `json:"fingerprints"`
	Scripts      map[string]ScriptInfo                 `json:"scripts"`
}

// remoteMetadataCachePathForURL returns the absolute path for the metadata cache file for a given base URL.
// Uses a 12-char SHA-256 prefix so each URL gets a stable, unique filename.
func remoteMetadataCachePathForURL(baseURL string) (string, error) {
	cacheDir, err := config.GetWritableCacheDir()
	if err != nil {
		return "", err
	}
	h := sha256.Sum256([]byte(baseURL))
	name := fmt.Sprintf("remote-scripts-%x.json.gz", h[:6]) // 12 hex chars
	return filepath.Join(cacheDir, name), nil
}

func localScriptsCachePath() (string, error) {
	cacheDir, err := config.GetWritableCacheDir()
	if err != nil {
		return "", err
	}
	return filepath.Join(cacheDir, "local-scripts-cache.json.gz"), nil
}

func collectLocalScriptsFingerprints(paths []string) map[string]localScriptsDirFingerprint {
	out := make(map[string]localScriptsDirFingerprint, len(paths))
	for _, dir := range paths {
		fp := localScriptsDirFingerprint{Exists: false}
		info, err := os.Stat(dir)
		if err != nil || !info.IsDir() {
			out[dir] = fp
			continue
		}

		fp.Exists = true
		fp.ModTimeUnix = info.ModTime().UnixNano()
		fp.ImmediateDirs = map[string]int{}

		entries, err := os.ReadDir(dir)
		if err == nil {
			for _, entry := range entries {
				if !entry.IsDir() {
					continue
				}
				subPath := filepath.Join(dir, entry.Name())
				subInfo, statErr := os.Stat(subPath)
				if statErr != nil {
					continue
				}
				fp.ImmediateDirs[entry.Name()] = int(subInfo.ModTime().Unix())
			}
		}

		out[dir] = fp
	}
	return out
}

func loadLocalScriptsCache(path string) (*LocalScriptsCache, error) {
	var cache LocalScriptsCache
	if err := utils.ReadGzipJSONFile(path, &cache); err != nil {
		return nil, err
	}
	return &cache, nil
}

func saveLocalScriptsCache(path string, cache *LocalScriptsCache) error {
	return utils.WriteGzipJSONFileAtomic(path, cache)
}

func loadRemoteMetadataCacheAny(path, sourceURL string, ttl time.Duration, checkTTL bool) (*RemoteMetadataCache, error) {
	var cache RemoteMetadataCache
	if err := utils.ReadGzipJSONFile(path, &cache); err != nil {
		return nil, fmt.Errorf("corrupt cache: %w", err)
	}
	if cache.SourceURL != sourceURL {
		return nil, fmt.Errorf("source URL changed")
	}
	if checkTTL && time.Since(cache.FetchedAt) > ttl {
		return nil, fmt.Errorf("cache expired")
	}
	return &cache, nil
}

// loadRemoteMetadataCache reads and validates the on-disk cache.
// Returns an error if the file is missing, corrupt, the URL has changed, or the TTL has expired.
func loadRemoteMetadataCache(path, sourceURL string, ttl time.Duration) (*RemoteMetadataCache, error) {
	return loadRemoteMetadataCacheAny(path, sourceURL, ttl, true)
}

// loadRemoteMetadataCacheIgnoreExpiry is like loadRemoteMetadataCache but skips the TTL check.
// Used as a stale fallback when the network is unavailable.
func loadRemoteMetadataCacheIgnoreExpiry(path, sourceURL string) (*RemoteMetadataCache, error) {
	return loadRemoteMetadataCacheAny(path, sourceURL, 0, false)
}

// saveRemoteMetadataCache writes the cache envelope to disk atomically.
func saveRemoteMetadataCache(path string, cache *RemoteMetadataCache) error {
	return utils.WriteGzipJSONFileAtomic(path, cache)
}

// fetchRemoteMetadata performs the HTTP+gzip+JSON fetch from the given metadata URL.
func fetchRemoteMetadata(url string) (map[string]RemoteScriptEntry, error) {
	client := &http.Client{Timeout: 30 * time.Second}
	resp, err := client.Get(url)
	if err != nil {
		return nil, fmt.Errorf("failed to fetch remote metadata: %w", err)
	}
	defer resp.Body.Close()

	if resp.StatusCode != http.StatusOK {
		return nil, fmt.Errorf("failed to fetch remote metadata: HTTP %d", resp.StatusCode)
	}

	gzReader, err := gzip.NewReader(resp.Body)
	if err != nil {
		return nil, fmt.Errorf("failed to decompress metadata: %w", err)
	}
	defer gzReader.Close()

	data, err := io.ReadAll(gzReader)
	if err != nil {
		return nil, fmt.Errorf("failed to read metadata: %w", err)
	}

	var metadata map[string]RemoteScriptEntry
	if err := json.Unmarshal(data, &metadata); err != nil {
		return nil, fmt.Errorf("failed to parse metadata: %w", err)
	}
	return metadata, nil
}

// convertToScriptInfoMap converts a raw metadata map to ScriptInfo map, tagging each entry with sourceURL and prebuiltLink.
func convertToScriptInfoMap(metadata map[string]RemoteScriptEntry, sourceURL, prebuiltLink string) map[string]ScriptInfo {
	scripts := make(map[string]ScriptInfo, len(metadata))
	for name, entry := range metadata {
		scripts[name] = ScriptInfo{
			Name:         name,
			Path:         entry.RelativePath,
			IsContainer:  strings.HasSuffix(entry.RelativePath, ".def"),
			IsRemote:     true,
			SourceURL:    sourceURL,
			PrebuiltLink: prebuiltLink,
		}
	}
	return scripts
}

// fetchPrebuiltLink retrieves the prebuilt base URL for a remote scripts source.
// It fetches {baseURL}/metadata/prebuilt_link (plain text, one line).
// Returns empty string if the file is absent or unreachable — meaning no prebuilt for that source.
func fetchPrebuiltLink(baseURL string) string {
	client := &http.Client{Timeout: 10 * time.Second}
	resp, err := client.Get(baseURL + "/metadata/prebuilt_link")
	if err != nil || resp.StatusCode != http.StatusOK {
		return ""
	}
	defer resp.Body.Close()
	data, err := io.ReadAll(io.LimitReader(resp.Body, 512))
	if err != nil {
		return ""
	}
	return strings.TrimSpace(string(data))
}

// PruneOrphanedCaches removes cache files (build and helper) for URLs no longer in ScriptsLinks.
// Errors are non-fatal.
func PruneOrphanedCaches() {
	cacheDir, err := config.GetWritableCacheDir()
	if err != nil {
		utils.PrintDebug("PruneOrphanedCaches: cannot get cache dir: %v", err)
		return
	}

	active := make(map[string]bool, len(config.Global.ScriptsLinks))
	for _, u := range config.Global.ScriptsLinks {
		active[u] = true
	}

	entries, err := os.ReadDir(cacheDir)
	if err != nil {
		utils.PrintDebug("PruneOrphanedCaches: cannot read cache dir: %v", err)
		return
	}

	for _, e := range entries {
		if e.IsDir() {
			continue
		}
		name := e.Name()
		if !strings.HasPrefix(name, "remote-scripts-") && !strings.HasPrefix(name, "helper-scripts-") {
			continue
		}
		if !strings.HasSuffix(name, ".json.gz") {
			continue
		}
		path := filepath.Join(cacheDir, name)
		var envelope struct {
			SourceURL string `json:"source_url"`
		}
		if err := utils.ReadGzipJSONFile(path, &envelope); err != nil {
			continue
		}
		if !active[envelope.SourceURL] {
			utils.PrintDebug("Removing orphaned cache file: %s (URL: %s)", name, envelope.SourceURL)
			if err := os.Remove(path); err != nil {
				utils.PrintDebug("Failed to remove orphaned cache: %v", err)
			}
		}
	}
}

// PreferRemote controls the precedence of build script resolution.
// When true, remote scripts take precedence over local scripts.
// When false (default), local scripts take precedence over remote scripts.
// Set by CLI commands (--remote flag) or config (prefer_remote: true).
var PreferRemote bool

// in-process caches — valid for the lifetime of one CLI invocation
var (
	cachedLocalScripts       map[string]ScriptInfo
	cachedRemoteScriptsByURL = map[string]map[string]ScriptInfo{} // per-URL
	cachedMergedRemote       map[string]ScriptInfo                // merged across all URLs
)

// ScriptInfo holds information about a build script
type ScriptInfo struct {
	Name         string // name/version format (e.g., "samtools/1.21")
	Path         string // full path to build script (local) or relative path (remote)
	IsContainer  bool   // true if .def file
	IsRemote     bool   // true if from remote repository
	SourceURL    string // base URL of the remote this script came from (empty for local)
	PrebuiltLink string // base URL for prebuilt downloads (from metadata/prebuilt_link); empty = no prebuilt
}

// RemoteScriptEntry represents a single script entry in the metadata
type RemoteScriptEntry struct {
	RelativePath string   `json:"relative_path"`
	Deps         []string `json:"deps"`
	Whatis       string   `json:"whatis"`
}

// GetLocalBuildScripts scans all build-scripts directories and returns a map of name -> ScriptInfo.
// Searches user → scratch → legacy → system directories.
// First match wins (user scripts shadow system ones).
func GetLocalBuildScripts() (map[string]ScriptInfo, error) {
	if ForceRefresh {
		cachedLocalScripts = nil
	}

	if !ForceRefresh && cachedLocalScripts != nil {
		return cachedLocalScripts, nil
	}

	searchPaths := config.GetBuildScriptSearchPaths()
	fingerprints := collectLocalScriptsFingerprints(searchPaths)
	if !ForceRefresh {
		if cachePath, err := localScriptsCachePath(); err == nil {
			if cache, cacheErr := loadLocalScriptsCache(cachePath); cacheErr == nil {
				if reflect.DeepEqual(cache.SearchPaths, searchPaths) && reflect.DeepEqual(cache.Fingerprints, fingerprints) {
					utils.PrintDebug("Using cached local build scripts (fetched %s ago)", time.Since(cache.FetchedAt).Round(time.Minute))
					cachedLocalScripts = cache.Scripts
					if cachedLocalScripts == nil {
						cachedLocalScripts = map[string]ScriptInfo{}
					}
					return cachedLocalScripts, nil
				}
			}
		}
	}

	scripts := make(map[string]ScriptInfo)

	// Scan all build script directories in priority order
	for _, buildScriptsDir := range searchPaths {
		if !utils.DirExists(buildScriptsDir) {
			continue
		}

		err := filepath.WalkDir(buildScriptsDir, func(path string, d os.DirEntry, err error) error {
			if err != nil {
				return nil // Skip errors, continue with other files
			}

			// Skip directories
			if d.IsDir() {
				return nil
			}

			// Skip .py and .sh files (helper scripts)
			if strings.HasSuffix(d.Name(), ".py") || strings.HasSuffix(d.Name(), ".sh") {
				return nil
			}

			// Skip template files
			if strings.Contains(path, "template") {
				return nil
			}

			// Generate key (relative path)
			relPath, err := filepath.Rel(buildScriptsDir, path)
			if err != nil {
				return nil
			}

			// Handle .def files
			isContainer := strings.HasSuffix(relPath, ".def")
			if isContainer {
				relPath = strings.TrimSuffix(relPath, ".def")
			}

			// First match wins - don't overwrite if already found in higher-priority dir
			if _, exists := scripts[relPath]; !exists {
				scripts[relPath] = ScriptInfo{
					Name:        relPath,
					Path:        path,
					IsContainer: isContainer,
					IsRemote:    false,
				}
			}
			return nil
		})

		if err != nil {
			// Log but continue with other directories
			utils.PrintDebug("Error scanning %s: %v", buildScriptsDir, err)
		}
	}

	cachedLocalScripts = scripts
	if cachePath, err := localScriptsCachePath(); err == nil {
		envelope := &LocalScriptsCache{
			FetchedAt:    time.Now(),
			SearchPaths:  searchPaths,
			Fingerprints: fingerprints,
			Scripts:      scripts,
		}
		if writeErr := saveLocalScriptsCache(cachePath, envelope); writeErr != nil {
			utils.PrintDebug("Failed to write local scripts cache: %v", writeErr)
		}
	}
	return scripts, nil
}

// getRemoteBuildScriptsForURL fetches (or loads from per-URL disk cache) the build scripts for one base URL.
func getRemoteBuildScriptsForURL(baseURL string) (map[string]ScriptInfo, error) {
	// In-process cache
	if scripts, ok := cachedRemoteScriptsByURL[baseURL]; ok {
		return scripts, nil
	}

	metaURL := baseURL + "/metadata/build-scripts.json.gz"
	cachePath, cacheErr := remoteMetadataCachePathForURL(baseURL)

	ttl := config.Global.MetadataCacheTTL
	if !ForceRefresh && cacheErr == nil && ttl > 0 {
		if cache, err := loadRemoteMetadataCache(cachePath, baseURL, ttl); err == nil {
			utils.PrintDebug("Using cached remote metadata for %s (fetched %s ago)", baseURL, time.Since(cache.FetchedAt).Round(time.Minute))
			scripts := convertToScriptInfoMap(cache.Metadata, baseURL, cache.PrebuiltLink)
			cachedRemoteScriptsByURL[baseURL] = scripts
			return scripts, nil
		}
	}

	// Network fetch
	utils.PrintDebug("Fetching remote build metadata from %s...", baseURL)
	metadata, err := fetchRemoteMetadata(metaURL)
	if err != nil {
		// Stale fallback
		if cacheErr == nil {
			if cache, staleErr := loadRemoteMetadataCacheIgnoreExpiry(cachePath, baseURL); staleErr == nil {
				utils.PrintWarning("Network unavailable for %s; using cached metadata (fetched %s ago)",
					baseURL, time.Since(cache.FetchedAt).Round(time.Minute))
				scripts := convertToScriptInfoMap(cache.Metadata, baseURL, cache.PrebuiltLink)
				cachedRemoteScriptsByURL[baseURL] = scripts
				return scripts, nil
			}
		}
		return nil, err
	}
	utils.PrintDebug("Fetched remote build metadata from %s (%d entries)", baseURL, len(metadata))

	// Fetch per-source prebuilt base URL (non-fatal; empty = no prebuilt for this source)
	prebuiltLink := fetchPrebuiltLink(baseURL)

	// Persist to disk (non-fatal)
	if cacheErr == nil && ttl > 0 {
		envelope := &RemoteMetadataCache{
			FetchedAt:    time.Now(),
			SourceURL:    baseURL,
			PrebuiltLink: prebuiltLink,
			Metadata:     metadata,
		}
		if writeErr := saveRemoteMetadataCache(cachePath, envelope); writeErr != nil {
			utils.PrintDebug("Failed to write metadata cache for %s: %v", baseURL, writeErr)
		}
	}

	scripts := convertToScriptInfoMap(metadata, baseURL, prebuiltLink)
	cachedRemoteScriptsByURL[baseURL] = scripts
	return scripts, nil
}

// GetRemoteBuildScripts fetches and merges build script metadata from all configured remote URLs.
// Earlier URLs in ScriptsLinks take precedence (first match wins).
func GetRemoteBuildScripts() (map[string]ScriptInfo, error) {
	// In-process merged cache
	if cachedMergedRemote != nil {
		return cachedMergedRemote, nil
	}

	if ForceRefresh {
		cachedRemoteScriptsByURL = map[string]map[string]ScriptInfo{}
	}

	merged := map[string]ScriptInfo{}
	var firstErr error

	for _, baseURL := range config.Global.ScriptsLinks {
		scripts, err := getRemoteBuildScriptsForURL(baseURL)
		if err != nil {
			utils.PrintDebug("Failed to fetch remote scripts from %s: %v", baseURL, err)
			if firstErr == nil {
				firstErr = err
			}
			continue
		}
		for name, info := range scripts {
			if _, exists := merged[name]; !exists {
				merged[name] = info // earlier URL wins
			}
		}
	}

	if len(merged) == 0 && firstErr != nil {
		return nil, firstErr
	}

	cachedMergedRemote = merged
	return merged, nil
}

// FindBuildScript looks for a build script by name/version.
// By default, local scripts take precedence over remote.
// When PreferRemote is true, remote scripts take precedence over local.
// Returns the ScriptInfo and a boolean indicating if it was found.
func FindBuildScript(nameVersion string) (ScriptInfo, bool) {
	normalized := utils.NormalizeNameVersion(nameVersion)

	if PreferRemote {
		// Remote first: remote scripts take precedence over local
		remoteScripts, err := GetRemoteBuildScripts()
		if err == nil {
			if info, found := remoteScripts[normalized]; found {
				return info, true
			}
		}

		localScripts, err := GetLocalBuildScripts()
		if err == nil {
			if info, found := localScripts[normalized]; found {
				return info, true
			}
		}
	} else {
		// Default: local scripts take precedence over remote
		localScripts, err := GetLocalBuildScripts()
		if err == nil {
			if info, found := localScripts[normalized]; found {
				return info, true
			}
		}

		remoteScripts, err := GetRemoteBuildScripts()
		if err == nil {
			if info, found := remoteScripts[normalized]; found {
				return info, true
			}
		}
	}

	return ScriptInfo{}, false
}

// DownloadRemoteScript downloads a remote build script to the local tmp directory.
// Uses info.SourceURL (set during metadata fetch) as the base for the download URL.
func DownloadRemoteScript(info ScriptInfo, tmpDir string) (string, error) {
	if !info.IsRemote {
		return info.Path, nil
	}

	// Use SourceURL if available; fall back to primary scripts_link for backward compat
	baseURL := info.SourceURL
	if baseURL == "" {
		baseURL = config.Global.ScriptsLink
	}

	rawURL := fmt.Sprintf("%s/build-scripts/%s", baseURL, info.Path)

	// Create HTTP client with timeout
	client := &http.Client{
		Timeout: 60 * time.Second,
	}

	// Fetch the script
	resp, err := client.Get(rawURL)
	if err != nil {
		return "", fmt.Errorf("failed to download script from %s: %w", rawURL, err)
	}
	defer resp.Body.Close()

	if resp.StatusCode != http.StatusOK {
		return "", fmt.Errorf("failed to download script: HTTP %d from %s", resp.StatusCode, rawURL)
	}

	// Create local path with normalized name (replace / with --)
	// Format: tmp/remote--name--version.def or tmp/remote--name--version.sh
	normalizedName := strings.ReplaceAll(info.Name, "/", "--")
	localPath := filepath.Join(tmpDir, "remote--"+normalizedName)
	if info.IsContainer {
		localPath += ".def"
	} else {
		localPath += ".sh"
	}

	// Ensure tmp directory exists
	if err := os.MkdirAll(tmpDir, utils.PermDir); err != nil {
		return "", fmt.Errorf("failed to create tmp directory: %w", err)
	}

	// Write to file
	file, err := utils.CreateFileWritable(localPath)
	if err != nil {
		return "", fmt.Errorf("failed to create file: %w", err)
	}
	defer file.Close()

	if _, err := io.Copy(file, resp.Body); err != nil {
		return "", fmt.Errorf("failed to write script: %w", err)
	}

	// Set permissions for downloaded script (664 for group-writable)
	if err := os.Chmod(localPath, utils.PermFile); err != nil {
		utils.PrintDebug("Failed to set permissions on %s: %v", localPath, err)
	}

	return localPath, nil
}
