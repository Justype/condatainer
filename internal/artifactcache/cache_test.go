package artifactcache

import (
	"encoding/json"
	"os"
	"path/filepath"
	"testing"
)

func TestCacheMergesFieldsAcrossProcesses(t *testing.T) {
	dir := t.TempDir()
	cachePath := filepath.Join(dir, FileName)
	image := filepath.Join(dir, "x.sqf")
	if err := os.WriteFile(image, []byte("image"), 0o644); err != nil {
		t.Fatal(err)
	}
	info, err := os.Lstat(image)
	if err != nil {
		t.Fatal(err)
	}
	first := New(func() string { return cachePath })
	second := New(func() string { return cachePath })
	first.Merge(image, info, func(record *Record) {
		record.RuntimeKnown = true
		record.Runtime = json.RawMessage(`{"name":"x"}`)
	})
	second.Merge(image, info, func(record *Record) {
		record.Name = "x/1"
		record.Identity = Key{Scheme: "identity-v1", SHA256: "abc"}
	})

	reloaded := New(func() string { return cachePath })
	record, ok := reloaded.Lookup(image, info)
	if !ok {
		t.Fatal("merged record was not found")
	}
	if !record.RuntimeKnown || len(record.Runtime) == 0 || record.Name != "x/1" || record.Identity.Scheme != "identity-v1" {
		t.Fatalf("fields were not merged: %#v", record)
	}
}

func TestCacheRejectsChangedFingerprintAndSymlink(t *testing.T) {
	dir := t.TempDir()
	image := filepath.Join(dir, "x.sqf")
	if err := os.WriteFile(image, []byte("one"), 0o644); err != nil {
		t.Fatal(err)
	}
	cache := New(func() string { return "" })
	info, _ := os.Lstat(image)
	cache.Merge(image, info, func(record *Record) { record.Name = "x" })
	if err := os.WriteFile(image, []byte("different"), 0o644); err != nil {
		t.Fatal(err)
	}
	changed, _ := os.Lstat(image)
	if _, ok := cache.Lookup(image, changed); ok {
		t.Fatal("changed file produced a cache hit")
	}
	link := filepath.Join(dir, "link.sqf")
	if err := os.Symlink(image, link); err != nil {
		t.Fatal(err)
	}
	linkInfo, _ := os.Lstat(link)
	cache.Merge(link, linkInfo, func(record *Record) { record.Name = "link" })
	if _, ok := cache.Lookup(link, linkInfo); ok {
		t.Fatal("symlink produced a cache hit")
	}
}

func TestBatchDefersAndMergesUpdates(t *testing.T) {
	dir := t.TempDir()
	cachePath := filepath.Join(dir, FileName)
	firstPath := filepath.Join(dir, "first.sqf")
	secondPath := filepath.Join(dir, "second.sqf")
	for _, path := range []string{firstPath, secondPath} {
		if err := os.WriteFile(path, []byte(path), 0o644); err != nil {
			t.Fatal(err)
		}
	}
	firstInfo, _ := os.Lstat(firstPath)
	secondInfo, _ := os.Lstat(secondPath)
	cache := New(func() string { return cachePath })
	batch := cache.NewBatch()
	batch.Merge(firstPath, firstInfo, func(record *Record) { record.Name = "first" })
	batch.Merge(secondPath, secondInfo, func(record *Record) { record.Name = "second" })
	if _, err := os.Stat(cachePath); !os.IsNotExist(err) {
		t.Fatalf("batch wrote before Flush: %v", err)
	}
	batch.Flush()

	reloaded := New(func() string { return cachePath })
	if record, ok := reloaded.Lookup(firstPath, firstInfo); !ok || record.Name != "first" {
		t.Fatalf("first batch record = %#v, %v", record, ok)
	}
	if record, ok := reloaded.Lookup(secondPath, secondInfo); !ok || record.Name != "second" {
		t.Fatalf("second batch record = %#v, %v", record, ok)
	}
}

func TestPruneMissingDropsStaleRecords(t *testing.T) {
	dir := t.TempDir()
	cache := New(func() string { return "" })
	path := filepath.Join(dir, "gone.sqf")
	if err := os.WriteFile(path, []byte("image"), 0o644); err != nil {
		t.Fatal(err)
	}
	info, _ := os.Lstat(path)
	cache.Merge(path, info, func(record *Record) { record.Name = "gone" })
	if err := os.Remove(path); err != nil {
		t.Fatal(err)
	}
	if got := cache.PruneMissing(); got != 1 {
		t.Fatalf("PruneMissing() = %d, want 1", got)
	}
	if _, ok := cache.Lookup(path, info); ok {
		t.Fatal("pruned record remains cached")
	}
}
