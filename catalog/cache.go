package catalog

import (
	"crypto/sha256"
	"encoding/hex"
	"os"
	"path/filepath"
	"time"
)

// Cache is a directory and a staleness bound supplied by the caller. The layout
// beneath Dir belongs to this package, keyed per source base.
//
// It bounds fetched bytes only. A filesystem source is authoritative and free
// to stat, so nothing local is cached and nothing local can go stale.
type Cache struct {
	Dir string
	TTL time.Duration
}

func (c Cache) enabled() bool { return c.Dir != "" }

// get returns the cached bytes for a source path, and whether they are within
// the TTL. Bytes past the TTL are still returned: a fetch that fails can serve
// them, which is the normal state of a compute node with no route out.
func (c Cache) get(base, path string) (data []byte, fresh bool) {
	if !c.enabled() {
		return nil, false
	}
	name := c.file(base, path)
	info, err := os.Stat(name)
	if err != nil {
		return nil, false
	}
	data, err = os.ReadFile(name)
	if err != nil {
		return nil, false
	}
	return data, time.Since(info.ModTime()) < c.TTL
}

// put writes bytes for a source path, replacing any previous copy in one step
// so a concurrent reader never sees a half-written file.
func (c Cache) put(base, path string, data []byte) error {
	if !c.enabled() {
		return nil
	}
	name := c.file(base, path)
	if err := os.MkdirAll(filepath.Dir(name), 0o755); err != nil {
		return err
	}
	tmp, err := os.CreateTemp(filepath.Dir(name), ".tmp-*")
	if err != nil {
		return err
	}
	defer os.Remove(tmp.Name())
	if _, err := tmp.Write(data); err != nil {
		tmp.Close()
		return err
	}
	if err := tmp.Close(); err != nil {
		return err
	}
	return os.Rename(tmp.Name(), name)
}

// file is where a source path lands: one directory per source base, keyed by a
// digest so two bases cannot collide, with the path kept readable below it.
func (c Cache) file(base, path string) string {
	sum := sha256.Sum256([]byte(base))
	return filepath.Join(c.Dir, hex.EncodeToString(sum[:8]), filepath.FromSlash(path))
}
