package helper

import (
	"bytes"
	"compress/gzip"
	"context"
	"encoding/json"
	"fmt"
	"io"
	"net/http"
	"os"
	"path/filepath"
	"strings"

	"github.com/Justype/condatainer/catalog"
	"github.com/Justype/condatainer/internal/config"
	"github.com/Justype/condatainer/internal/utils"
)

// helperIndexPath is where a collection publishes its helper listing.
const helperIndexPath = "index/helpers.json"

// RemoteScriptEntry is one helper a source offers.
type RemoteScriptEntry struct {
	Path   string `json:"path"`
	Source string `json:"-"` // the source base it came from
}

// RefreshRemoteMetadata reads the helper index from every configured source,
// merged with earlier sources winning.
//
// Fetching and caching belong to the catalog, which already bounds staleness and
// serves a cached copy when a source is unreachable. A source with no helper
// index simply contributes nothing.
func RefreshRemoteMetadata(ctx context.Context, force bool, w io.Writer) (map[string]RemoteScriptEntry, error) {
	if force {
		if err := config.RefreshCatalogCache(); err != nil && w != nil {
			fmt.Fprintf(w, "WARN: failed to clear the source cache: %v\n", err)
		}
	}

	cat, err := config.OpenCatalog(ctx)
	if err != nil {
		return nil, err
	}

	merged := map[string]RemoteScriptEntry{}
	var firstErr error
	for _, src := range cat {
		meta, err := readHelperIndex(ctx, cat, src, w)
		if err != nil {
			if helperIndexNotFound(err) {
				continue
			}
			if w != nil {
				fmt.Fprintf(w, "WARN: no helper index in %s: %v\n", src.Name, err)
			}
			if firstErr == nil {
				firstErr = err
			}
			continue
		}
		if src.Stale && w != nil {
			fmt.Fprintf(w, "WARN: %s was unreachable; using its cached helper index\n", src.Name)
		}
		for name, entry := range meta {
			if _, exists := merged[name]; !exists {
				entry.Source = src.Base
				merged[name] = entry
			}
		}
	}

	if len(merged) == 0 && firstErr != nil {
		return nil, firstErr
	}
	return merged, nil
}

// readHelperIndex reads index/helpers.json.gz, falling back to the plain file.

// helperIndexNotFound identifies an optional index that a source does not
// publish. Catalog HTTP errors include the response status in their stable
// message; local sources return an os.ErrNotExist wrapper.
func helperIndexNotFound(err error) bool {
	return os.IsNotExist(err) || strings.HasSuffix(err.Error(), ": 404 Not Found")
}
func readHelperIndex(ctx context.Context, cat catalog.Catalog, src *catalog.Source, w io.Writer) (map[string]RemoteScriptEntry, error) {
	if data, err := cat.ReadPath(ctx, src, helperIndexPath+".gz"); err == nil {
		if plain, err := gunzip(data); err == nil {
			return decodeHelperIndex(plain)
		}
	}
	data, err := cat.ReadPath(ctx, src, helperIndexPath)
	if err != nil {
		return nil, err
	}
	return decodeHelperIndex(data)
}

func decodeHelperIndex(data []byte) (map[string]RemoteScriptEntry, error) {
	var meta map[string]RemoteScriptEntry
	if err := json.Unmarshal(data, &meta); err != nil {
		return nil, fmt.Errorf("failed to parse helper index: %w", err)
	}
	return meta, nil
}

func gunzip(data []byte) ([]byte, error) {
	zr, err := gzip.NewReader(bytes.NewReader(data))
	if err != nil {
		return nil, err
	}
	defer zr.Close()
	return io.ReadAll(zr)
}

// UpdateRemoteScripts syncs helper scripts from the configured sources.
// When name is non-empty, only that helper is updated.
func UpdateRemoteScripts(ctx context.Context, name string, forceMetadata bool, w io.Writer) error {
	entries, err := RefreshRemoteMetadata(ctx, forceMetadata, w)
	if err != nil {
		return fmt.Errorf("failed to read helper indexes: %w", err)
	}
	if name != "" {
		entry, ok := entries[name]
		if !ok {
			return fmt.Errorf("helper script %q not found in any source", name)
		}
		entries = map[string]RemoteScriptEntry{name: entry}
	}

	helperScriptsDir, err := config.GetWritableHelperScriptsDir()
	if err != nil {
		return err
	}
	if err := utils.MkdirAllShared(helperScriptsDir); err != nil {
		return fmt.Errorf("failed to create helper scripts directory: %w", err)
	}

	if name == "" && w != nil {
		fmt.Fprintln(w, "Updating all helper scripts...")
	}
	for scriptName, entry := range entries {
		if entry.Path == "" || entry.Source == "" {
			continue
		}
		url := fmt.Sprintf("%s/%s", entry.Source, entry.Path)
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
	destDir := filepath.Dir(destPath)
	if err := utils.MkdirAllShared(destDir); err != nil {
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
	if err := utils.MakeExecutable(destPath); err != nil {
		return err
	}
	return nil
}
