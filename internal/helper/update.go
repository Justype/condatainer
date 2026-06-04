package helper

import (
	"compress/gzip"
	"context"
	"crypto/sha256"
	"encoding/json"
	"fmt"
	"io"
	"net/http"
	"os"
	"path/filepath"
	"sync"
	"time"

	"github.com/Justype/condatainer/internal/config"
	"github.com/Justype/condatainer/internal/utils"
)

// RemoteScriptEntry represents a helper script entry from remote metadata.
type RemoteScriptEntry struct {
	Path      string `json:"path"`
	SourceURL string `json:"-"`
}

type remoteHelperEntry struct {
	Path string `json:"path"`
}

// RemoteMetadataCache is the on-disk envelope for cached helper script metadata.
type RemoteMetadataCache struct {
	FetchedAt time.Time                    `json:"fetched_at"`
	SourceURL string                       `json:"source_url"`
	Metadata  map[string]remoteHelperEntry `json:"metadata"`
}

var (
	remoteMetadataMu       sync.Mutex
	remoteMetadataBySource = map[string]map[string]RemoteScriptEntry{}
)

// RefreshRemoteMetadata refreshes or loads remote helper script metadata.
// Earlier URLs in ScriptsLinks take precedence.
func RefreshRemoteMetadata(ctx context.Context, force bool, w io.Writer) (map[string]RemoteScriptEntry, error) {
	if force {
		remoteMetadataMu.Lock()
		remoteMetadataBySource = map[string]map[string]RemoteScriptEntry{}
		remoteMetadataMu.Unlock()
	}

	merged := map[string]RemoteScriptEntry{}
	var firstErr error

	for _, baseURL := range config.Global.ScriptsLinks {
		meta, err := fetchRemoteMetadataForSource(ctx, baseURL, force, w)
		if err != nil {
			if w != nil {
				fmt.Fprintf(w, "WARN: failed to fetch helper metadata from %s: %v\n", baseURL, err)
			}
			if firstErr == nil {
				firstErr = err
			}
			continue
		}
		for name, entry := range meta {
			if _, exists := merged[name]; !exists {
				merged[name] = entry
			}
		}
	}

	if len(merged) == 0 && firstErr != nil {
		return nil, firstErr
	}
	return merged, nil
}

// UpdateRemoteScripts syncs helper scripts from the configured remote script URLs.
// When name is non-empty, only that helper is updated.
func UpdateRemoteScripts(ctx context.Context, name string, forceMetadata bool, w io.Writer) error {
	entries, err := RefreshRemoteMetadata(ctx, forceMetadata, w)
	if err != nil {
		return fmt.Errorf("failed to fetch remote helper metadata: %w", err)
	}
	if name != "" {
		entry, ok := entries[name]
		if !ok {
			return fmt.Errorf("helper script %q not found in remote metadata", name)
		}
		entries = map[string]RemoteScriptEntry{name: entry}
	}

	helperScriptsDir, err := config.GetWritableHelperScriptsDir()
	if err != nil {
		return err
	}
	if err := os.MkdirAll(helperScriptsDir, utils.PermDir); err != nil {
		return fmt.Errorf("failed to create helper scripts directory: %w", err)
	}

	if name == "" && w != nil {
		fmt.Fprintln(w, "Updating all helper scripts...")
	}
	for scriptName, entry := range entries {
		if entry.Path == "" {
			continue
		}
		sourceURL := entry.SourceURL
		if sourceURL == "" {
			sourceURL = config.Global.ScriptsLink
		}
		url := fmt.Sprintf("%s/%s", sourceURL, entry.Path)
		dest := filepath.Join(helperScriptsDir, filepath.Base(entry.Path))
		if w != nil {
			fmt.Fprintf(w, "Updating %s\n", scriptName)
		}
		if err := downloadRemoteExecutable(ctx, url, dest); err != nil {
			if w != nil {
				fmt.Fprintf(w, "WARN: failed to update %s: %v\n", scriptName, err)
			}
			continue
		}
	}
	if w != nil {
		fmt.Fprintln(w, "Helper update finished.")
	}
	return nil
}

func fetchRemoteMetadataForSource(ctx context.Context, baseURL string, force bool, w io.Writer) (map[string]RemoteScriptEntry, error) {
	remoteMetadataMu.Lock()
	if cached, ok := remoteMetadataBySource[baseURL]; ok {
		remoteMetadataMu.Unlock()
		return cached, nil
	}
	remoteMetadataMu.Unlock()

	metaURL := baseURL + "/metadata/helper-scripts.json.gz"
	cachePath, cacheErr := remoteMetadataCachePathForSource(baseURL)
	ttl := config.Global.MetadataCacheTTL

	if !force && cacheErr == nil && ttl > 0 {
		if cache, err := loadRemoteMetadataCache(cachePath, baseURL, ttl, true); err == nil {
			if w != nil {
				fmt.Fprintf(w, "Using cached helper metadata for %s (fetched %s ago)\n", baseURL, time.Since(cache.FetchedAt).Round(time.Minute))
			}
			result := injectRemoteSource(cache.Metadata, baseURL)
			remoteMetadataMu.Lock()
			remoteMetadataBySource[baseURL] = result
			remoteMetadataMu.Unlock()
			return result, nil
		}
	}

	if w != nil {
		fmt.Fprintf(w, "Fetching helper metadata from %s ...\n", metaURL)
	}
	rawMeta, err := fetchRawRemoteMetadata(ctx, metaURL)
	if err != nil {
		if cacheErr == nil {
			if cache, staleErr := loadRemoteMetadataCache(cachePath, baseURL, 0, false); staleErr == nil {
				if w != nil {
					fmt.Fprintf(w, "WARN: network unavailable for %s; using cached helper metadata (fetched %s ago)\n",
						metaURL, time.Since(cache.FetchedAt).Round(time.Minute))
				}
				result := injectRemoteSource(cache.Metadata, baseURL)
				remoteMetadataMu.Lock()
				remoteMetadataBySource[baseURL] = result
				remoteMetadataMu.Unlock()
				return result, nil
			}
		}
		return nil, err
	}

	if cacheErr == nil && ttl > 0 {
		envelope := RemoteMetadataCache{
			FetchedAt: time.Now(),
			SourceURL: baseURL,
			Metadata:  rawMeta,
		}
		if err := utils.WriteGzipJSONFileAtomic(cachePath, envelope); err != nil && w != nil {
			fmt.Fprintf(w, "WARN: failed to write helper metadata cache for %s: %v\n", baseURL, err)
		}
	}

	result := injectRemoteSource(rawMeta, baseURL)
	remoteMetadataMu.Lock()
	remoteMetadataBySource[baseURL] = result
	remoteMetadataMu.Unlock()
	return result, nil
}

func remoteMetadataCachePathForSource(baseURL string) (string, error) {
	cacheDir, err := config.GetWritableCacheDir()
	if err != nil {
		return "", err
	}
	h := sha256.Sum256([]byte(baseURL))
	name := fmt.Sprintf("helper-scripts-%x.json.gz", h[:6])
	return filepath.Join(cacheDir, name), nil
}

func loadRemoteMetadataCache(path, baseURL string, ttl time.Duration, checkTTL bool) (*RemoteMetadataCache, error) {
	var cache RemoteMetadataCache
	if err := utils.ReadGzipJSONFile(path, &cache); err != nil {
		return nil, err
	}
	if cache.SourceURL != baseURL {
		return nil, fmt.Errorf("source URL changed")
	}
	if checkTTL && time.Since(cache.FetchedAt) > ttl {
		return nil, fmt.Errorf("cache expired")
	}
	return &cache, nil
}

func injectRemoteSource(raw map[string]remoteHelperEntry, baseURL string) map[string]RemoteScriptEntry {
	out := make(map[string]RemoteScriptEntry, len(raw))
	for name, entry := range raw {
		out[name] = RemoteScriptEntry{Path: entry.Path, SourceURL: baseURL}
	}
	return out
}

func fetchRawRemoteMetadata(ctx context.Context, url string) (map[string]remoteHelperEntry, error) {
	req, err := http.NewRequestWithContext(ctx, http.MethodGet, url, nil)
	if err != nil {
		return nil, err
	}
	resp, err := http.DefaultClient.Do(req)
	if err != nil {
		return nil, err
	}
	defer resp.Body.Close()

	if resp.StatusCode != http.StatusOK {
		return nil, fmt.Errorf("failed to fetch metadata: HTTP %d", resp.StatusCode)
	}
	gzReader, err := gzip.NewReader(resp.Body)
	if err != nil {
		return nil, fmt.Errorf("failed to decompress metadata: %w", err)
	}
	defer gzReader.Close()

	var metadata map[string]remoteHelperEntry
	if err := json.NewDecoder(gzReader).Decode(&metadata); err != nil {
		return nil, err
	}
	return metadata, nil
}

func downloadRemoteExecutable(ctx context.Context, url, destPath string) error {
	req, err := http.NewRequestWithContext(ctx, http.MethodGet, url, nil)
	if err != nil {
		return err
	}
	resp, err := http.DefaultClient.Do(req)
	if err != nil {
		return err
	}
	defer resp.Body.Close()

	if resp.StatusCode != http.StatusOK {
		return fmt.Errorf("HTTP %d", resp.StatusCode)
	}
	if err := os.MkdirAll(filepath.Dir(destPath), utils.PermDir); err != nil {
		return err
	}
	out, err := utils.CreateFileWritable(destPath)
	if err != nil {
		return err
	}
	defer out.Close()

	if _, err := io.Copy(out, resp.Body); err != nil {
		return err
	}
	return os.Chmod(destPath, utils.PermExec)
}
