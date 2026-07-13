package overlay

import (
	"os"
	"path/filepath"
	"sync"

	"github.com/Justype/condatainer/internal/config"
	"github.com/Justype/condatainer/internal/utils"
)

// osTypeCacheName is the cache file name inside the condatainer cache directory.
const osTypeCacheName = "sqf-ostype.json.gz"

// osTypeEntry records the IsOSType verdict for a .sqf file at a given size/mtime.
type osTypeEntry struct {
	Size  int64 `json:"size"`
	MTime int64 `json:"mtime"` // ModTime in unix nanoseconds
	OS    bool  `json:"os"`
}

// osTypeCache persists IsOSType verdicts across processes so repeated checks
// (e.g. shell completion, one process per keystroke) skip the unsquashfs probe.
// Entries are keyed by absolute path and validated against current size/mtime,
// so rebuilt overlays (create, build --update) are re-probed automatically.
type osTypeCache struct {
	mu       sync.Mutex
	pathFn   func() string // resolves the cache file location; "" disables persistence
	filePath string        // resolved by load()
	entries  map[string]osTypeEntry
	loaded   bool
}

var globalOSTypeCache = &osTypeCache{pathFn: defaultOSTypeCachePath}

// defaultOSTypeCachePath returns the cache file path, or "" if no cache dir is writable.
func defaultOSTypeCachePath() string {
	dir, err := config.GetWritableCacheDir()
	if err != nil {
		return ""
	}
	return filepath.Join(dir, osTypeCacheName)
}

// load reads the cache file once per process. A missing or unreadable file
// leaves the cache empty. Callers must hold c.mu.
func (c *osTypeCache) load() {
	if c.loaded {
		return
	}
	c.loaded = true
	c.entries = map[string]osTypeEntry{}
	if c.pathFn != nil {
		c.filePath = c.pathFn()
	}
	if c.filePath == "" {
		return
	}
	if err := utils.ReadGzipJSONFile(c.filePath, &c.entries); err != nil || c.entries == nil {
		c.entries = map[string]osTypeEntry{}
	}
}

// save rewrites the cache file. Callers must hold c.mu.
func (c *osTypeCache) save() {
	if c.filePath == "" {
		return
	}
	_ = utils.WriteGzipJSONFileAtomic(c.filePath, c.entries)
}

// lookup returns the cached verdict for path if an entry matches fi's size and mtime.
func (c *osTypeCache) lookup(path string, fi os.FileInfo) (verdict bool, ok bool) {
	c.mu.Lock()
	defer c.mu.Unlock()
	c.load()
	e, ok := c.entries[path]
	if !ok || e.Size != fi.Size() || e.MTime != fi.ModTime().UnixNano() {
		return false, false
	}
	return e.OS, true
}

// store records a verdict for path at fi's size/mtime and persists the cache.
func (c *osTypeCache) store(path string, fi os.FileInfo, isOS bool) {
	c.mu.Lock()
	defer c.mu.Unlock()
	c.load()
	c.entries[path] = osTypeEntry{Size: fi.Size(), MTime: fi.ModTime().UnixNano(), OS: isOS}
	c.save()
}

// forget drops the entry for path, persisting the cache only if an entry existed.
func (c *osTypeCache) forget(path string) {
	c.mu.Lock()
	defer c.mu.Unlock()
	c.load()
	if _, ok := c.entries[path]; !ok {
		return
	}
	delete(c.entries, path)
	c.save()
}

// ForgetOSType drops the cached IsOSType verdict for an overlay path.
// Call after deleting an overlay file so the cache does not keep stale entries.
func ForgetOSType(path string) {
	if abs, err := filepath.Abs(path); err == nil {
		path = abs
	}
	globalOSTypeCache.forget(path)
}
