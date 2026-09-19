package toolpath

import (
	"encoding/json"
	"os"
	"path/filepath"

	"github.com/Justype/condatainer/internal/config"
)

// A small per-user file remembers what Resolve worked out about host binaries,
// so a new process does not run them again: where an apptainer's bundled tools
// live, and whether a tool cleared its version floor. Each entry records the
// size and modification time of the binary it describes and is ignored when
// they differ, so an unloaded module, an upgrade or another node's binary
// simply misses. The file is never read or written for correctness; any
// failure means recomputing.

const cacheFileName = "toolpath.json"

// insideContainer is config.IsInsideContainer, replaceable by tests.
var insideContainer = config.IsInsideContainer

// cachePath returns the cache file, or "" when there is none to use. Replaceable
// so tests never touch a real cache.
var cachePath = defaultCachePath

// defaultCachePath is the personal cache file. Inside a container there is
// none: the same path there can name a different binary than on the host, so
// the two would overwrite each other's entries.
func defaultCachePath() string {
	if insideContainer() {
		return ""
	}
	if dir := config.GetUserCacheDir(); dir != "" {
		return filepath.Join(dir, cacheFileName)
	}
	return ""
}

type cacheEntry struct {
	Size  int64  `json:"size"`
	ModNS int64  `json:"mod_ns"`
	Value string `json:"value"`
}

type cacheFile struct {
	Entries map[string]cacheEntry `json:"entries"`
}

func identity(path string) (size, modNS int64, ok bool) {
	info, err := os.Stat(path)
	if err != nil {
		return 0, 0, false
	}
	return info.Size(), info.ModTime().UnixNano(), true
}

func loadCache() cacheFile {
	c := cacheFile{Entries: map[string]cacheEntry{}}
	if p := cachePath(); p != "" {
		if data, err := os.ReadFile(p); err == nil {
			if json.Unmarshal(data, &c) != nil || c.Entries == nil {
				c = cacheFile{Entries: map[string]cacheEntry{}}
			}
		}
	}
	return c
}

// cacheGet returns the value stored under kind for the binary at path, when
// that binary is unchanged since it was stored.
func cacheGet(kind, path string) (string, bool) {
	size, mod, ok := identity(path)
	if !ok {
		return "", false
	}
	e, found := loadCache().Entries[kind+"|"+path]
	if !found || e.Size != size || e.ModNS != mod {
		return "", false
	}
	return e.Value, true
}

// cachePut stores value under kind for the binary at path. Best effort.
func cachePut(kind, path, value string) {
	p := cachePath()
	size, mod, ok := identity(path)
	if p == "" || !ok {
		return
	}
	if os.MkdirAll(filepath.Dir(p), 0o755) != nil {
		return
	}
	c := loadCache()
	c.Entries[kind+"|"+path] = cacheEntry{Size: size, ModNS: mod, Value: value}
	data, err := json.Marshal(c)
	if err != nil {
		return
	}
	tmp, err := os.CreateTemp(filepath.Dir(p), cacheFileName+".*")
	if err != nil {
		return
	}
	_, werr := tmp.Write(data)
	if cerr := tmp.Close(); werr != nil || cerr != nil || os.Rename(tmp.Name(), p) != nil {
		os.Remove(tmp.Name())
	}
}
