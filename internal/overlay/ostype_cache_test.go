package overlay

import (
	"os"
	"path/filepath"
	"testing"
	"time"
)

// newTestCache returns an osTypeCache persisting to a file in a temp dir.
func newTestCache(t *testing.T) (*osTypeCache, string) {
	t.Helper()
	cachePath := filepath.Join(t.TempDir(), osTypeCacheName)
	return &osTypeCache{pathFn: func() string { return cachePath }}, cachePath
}

// writeTempSqf creates a dummy .sqf file and returns its path and FileInfo.
func writeTempSqf(t *testing.T, dir, name, content string) (string, os.FileInfo) {
	t.Helper()
	path := filepath.Join(dir, name)
	if err := os.WriteFile(path, []byte(content), 0o644); err != nil {
		t.Fatal(err)
	}
	fi, err := os.Stat(path)
	if err != nil {
		t.Fatal(err)
	}
	return path, fi
}

func TestOSTypeCacheStoreLookup(t *testing.T) {
	cache, cachePath := newTestCache(t)
	sqf, fi := writeTempSqf(t, t.TempDir(), "os.sqf", "content")

	if _, ok := cache.lookup(sqf, fi); ok {
		t.Fatal("expected miss on empty cache")
	}

	cache.store(sqf, fi, true)
	verdict, ok := cache.lookup(sqf, fi)
	if !ok || !verdict {
		t.Fatalf("expected cached true verdict, got verdict=%v ok=%v", verdict, ok)
	}

	// A fresh cache instance must read the verdict back from disk.
	reloaded := &osTypeCache{pathFn: func() string { return cachePath }}
	verdict, ok = reloaded.lookup(sqf, fi)
	if !ok || !verdict {
		t.Fatalf("expected persisted true verdict, got verdict=%v ok=%v", verdict, ok)
	}
}

func TestOSTypeCacheInvalidatesOnChange(t *testing.T) {
	cache, _ := newTestCache(t)
	dir := t.TempDir()
	sqf, fi := writeTempSqf(t, dir, "app.sqf", "v1")
	cache.store(sqf, fi, false)

	// Same size, different mtime → miss.
	if err := os.Chtimes(sqf, time.Now(), fi.ModTime().Add(time.Second)); err != nil {
		t.Fatal(err)
	}
	fi2, err := os.Stat(sqf)
	if err != nil {
		t.Fatal(err)
	}
	if _, ok := cache.lookup(sqf, fi2); ok {
		t.Fatal("expected miss after mtime change")
	}

	// Different size → miss.
	if err := os.WriteFile(sqf, []byte("longer content"), 0o644); err != nil {
		t.Fatal(err)
	}
	fi3, err := os.Stat(sqf)
	if err != nil {
		t.Fatal(err)
	}
	if _, ok := cache.lookup(sqf, fi3); ok {
		t.Fatal("expected miss after size change")
	}
}

func TestOSTypeCacheForget(t *testing.T) {
	cache, cachePath := newTestCache(t)
	sqf, fi := writeTempSqf(t, t.TempDir(), "os.sqf", "content")
	cache.store(sqf, fi, true)

	cache.forget(sqf)
	if _, ok := cache.lookup(sqf, fi); ok {
		t.Fatal("expected miss after forget")
	}

	reloaded := &osTypeCache{pathFn: func() string { return cachePath }}
	if _, ok := reloaded.lookup(sqf, fi); ok {
		t.Fatal("expected forget to be persisted")
	}
}

func TestOSTypeCacheNoPersistPath(t *testing.T) {
	// pathFn returning "" disables persistence but keeps in-process caching.
	cache := &osTypeCache{pathFn: func() string { return "" }}
	sqf, fi := writeTempSqf(t, t.TempDir(), "os.sqf", "content")

	cache.store(sqf, fi, true)
	verdict, ok := cache.lookup(sqf, fi)
	if !ok || !verdict {
		t.Fatalf("expected in-process verdict, got verdict=%v ok=%v", verdict, ok)
	}
}

func TestOSTypeCacheCorruptFile(t *testing.T) {
	cache, cachePath := newTestCache(t)
	if err := os.WriteFile(cachePath, []byte("not gzip"), 0o644); err != nil {
		t.Fatal(err)
	}
	sqf, fi := writeTempSqf(t, t.TempDir(), "os.sqf", "content")

	// Corrupt cache file is treated as empty; store must recover it.
	if _, ok := cache.lookup(sqf, fi); ok {
		t.Fatal("expected miss with corrupt cache file")
	}
	cache.store(sqf, fi, true)

	reloaded := &osTypeCache{pathFn: func() string { return cachePath }}
	verdict, ok := reloaded.lookup(sqf, fi)
	if !ok || !verdict {
		t.Fatalf("expected recovered cache to persist verdict, got verdict=%v ok=%v", verdict, ok)
	}
}
