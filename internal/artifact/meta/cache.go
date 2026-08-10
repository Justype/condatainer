package meta

import (
	"os"
	"path/filepath"
	"sync"

	"github.com/Justype/condatainer/internal/config"
	"github.com/Justype/condatainer/internal/utils"
)

// cacheName is the cache file inside the CondaTainer cache directory.
const cacheName = "cnt-runtime.json.gz"

// cacheEntry records what one image's runtime read produced, at a given
// size/mtime. Runtime is nil when the image has none — negative verdicts are
// cached too.
type cacheEntry struct {
	Size    int64    `json:"size"`
	MTime   int64    `json:"mtime"` // ModTime in unix nanoseconds
	Runtime *Runtime `json:"runtime,omitempty"`
}

// runtimeCache persists runtime reads across processes. Entries are keyed by
// absolute path and validated against the file's size and mtime, so a rebuilt
// image is re-read automatically. See the README's Manifests.
type runtimeCache struct {
	mu       sync.Mutex
	pathFn   func() string // resolves the cache file location; "" disables persistence
	filePath string        // resolved by load()
	entries  map[string]cacheEntry
	loaded   bool
}

var globalCache = &runtimeCache{pathFn: defaultCachePath}

// defaultCachePath returns the cache file path, or "" if no cache dir is writable.
func defaultCachePath() string {
	dir, err := config.GetWritableCacheDir()
	if err != nil {
		return ""
	}
	return filepath.Join(dir, cacheName)
}

// load reads the cache file once per process. A missing or unreadable file
// leaves the cache empty. Callers must hold c.mu.
func (c *runtimeCache) load() {
	if c.loaded {
		return
	}
	c.loaded = true
	c.entries = map[string]cacheEntry{}
	if c.pathFn != nil {
		c.filePath = c.pathFn()
	}
	if c.filePath == "" {
		return
	}
	if err := utils.ReadGzipJSONFile(c.filePath, &c.entries); err != nil || c.entries == nil {
		c.entries = map[string]cacheEntry{}
	}
}

// save rewrites the cache file. Callers must hold c.mu.
func (c *runtimeCache) save() {
	if c.filePath == "" {
		return
	}
	_ = utils.WriteGzipJSONFileAtomic(c.filePath, c.entries)
}

// lookup returns the cached result for path when an entry matches fi's size and
// mtime. found reports whether the cache knew; has reports whether the image
// had a runtime document.
func (c *runtimeCache) lookup(path string, fi os.FileInfo) (rt Runtime, has bool, found bool) {
	c.mu.Lock()
	defer c.mu.Unlock()
	c.load()
	e, ok := c.entries[path]
	if !ok || e.Size != fi.Size() || e.MTime != fi.ModTime().UnixNano() {
		return Runtime{}, false, false
	}
	if e.Runtime == nil {
		return Runtime{}, false, true
	}
	return *e.Runtime, true, true
}

// store records a result for path at fi's size/mtime and persists the cache.
// A nil runtime records that the image has none.
func (c *runtimeCache) store(path string, fi os.FileInfo, rt *Runtime) {
	c.mu.Lock()
	defer c.mu.Unlock()
	c.load()
	c.entries[path] = cacheEntry{Size: fi.Size(), MTime: fi.ModTime().UnixNano(), Runtime: rt}
	c.save()
}

// forget drops the cached result for an image path. remove calls this so a
// deleted image does not leave a stale entry behind for a later file that
// happens to land on the same path.
func (c *runtimeCache) forget(path string) {
	c.mu.Lock()
	defer c.mu.Unlock()
	c.load()
	if _, ok := c.entries[path]; !ok {
		return
	}
	delete(c.entries, path)
	c.save()
}

// Forget drops the cached runtime document for an image path.
func Forget(imagePath string) {
	abs, err := filepath.Abs(imagePath)
	if err != nil {
		abs = imagePath
	}
	globalCache.forget(abs)
}
