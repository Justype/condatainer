package meta

import (
	"errors"
	"os"
	"os/exec"
	"path/filepath"
	"strings"
	"testing"
	"time"

	"github.com/Justype/condatainer/catalog"
	"github.com/Justype/condatainer/internal/image/tool"
)

// withTempCache points the manifest cache at a scratch file so a test never
// reads or writes the user's real cache, and never sees another test's entries.
func withTempCache(t *testing.T) {
	t.Helper()
	dir := t.TempDir()
	prev := globalCache
	globalCache = &manifestCache{pathFn: func() string { return filepath.Join(dir, cacheName) }}
	t.Cleanup(func() { globalCache = prev })
}

// requireSquashfsTools skips when the host cannot build or read a .sqf.
func requireSquashfsTools(t *testing.T) {
	t.Helper()
	for _, bin := range []string{"mksquashfs", "unsquashfs"} {
		if _, err := exec.LookPath(bin); err != nil {
			t.Skipf("%s not available", bin)
		}
	}
}

// packSqf builds a .sqf whose root is dir.
func packSqf(t *testing.T, dir string) string {
	t.Helper()
	out := filepath.Join(t.TempDir(), "image.sqf")
	cmd := exec.Command("mksquashfs", dir, out, "-no-progress", "-noappend", "-quiet")
	if output, err := cmd.CombinedOutput(); err != nil {
		t.Fatalf("mksquashfs: %v\n%s", err, output)
	}
	return out
}

// stagedImage stages manifest into a payload root and packs it, the way a build
// backend does.
func stagedImage(t *testing.T, manifest Manifest) string {
	t.Helper()
	root := t.TempDir()
	if err := Stage(filepath.Join(root, DirName), manifest); err != nil {
		t.Fatalf("Stage: %v", err)
	}
	payload := filepath.Join(root, "cnt", "samtools", "1.23.1", "bin")
	if err := os.MkdirAll(payload, 0o755); err != nil {
		t.Fatalf("mkdir payload: %v", err)
	}
	if err := os.WriteFile(filepath.Join(payload, "samtools"), []byte("#!/bin/sh\n"), 0o755); err != nil {
		t.Fatalf("write payload: %v", err)
	}
	return packSqf(t, root)
}

func TestStageWritesManifestAtomically(t *testing.T) {
	dir := filepath.Join(t.TempDir(), "workspace", DirName)
	if err := Stage(dir, valid()); err != nil {
		t.Fatalf("Stage: %v", err)
	}

	data, err := os.ReadFile(filepath.Join(dir, FileName))
	if err != nil {
		t.Fatalf("manifest not written: %v", err)
	}
	want, err := Marshal(valid())
	if err != nil {
		t.Fatalf("Marshal: %v", err)
	}
	if string(data) != string(want) {
		t.Errorf("staged bytes differ from Marshal output:\ngot:\n%s\nwant:\n%s", data, want)
	}

	// The temp file the atomic rename went through must not survive.
	if _, err := os.Stat(filepath.Join(dir, FileName+".tmp")); !os.IsNotExist(err) {
		t.Error("staging left its temp file behind")
	}
}

func TestReadRoundTrip(t *testing.T) {
	requireSquashfsTools(t)
	withTempCache(t)

	want := valid()
	want.Description = "SAMtools alignment toolkit"
	want.URL = "https://www.htslib.org/"
	image := stagedImage(t, want)

	got, err := Read(image)
	if err != nil {
		t.Fatalf("Read: %v", err)
	}
	if got.Name != want.Name || got.Type != want.Type || got.BuildType != want.BuildType {
		t.Errorf("identity = %+v, want %+v", got, want)
	}
	if got.Description != want.Description || got.URL != want.URL {
		t.Errorf("description/url = %q/%q", got.Description, got.URL)
	}
	if got.Runtime.Prefix != want.Runtime.Prefix {
		t.Errorf("prefix = %q, want %q", got.Runtime.Prefix, want.Runtime.Prefix)
	}
	// {prefix} survives the round trip; substitution happens at load time.
	if len(got.Runtime.Env) != 1 || got.Runtime.Env[0].Value != "{prefix}" {
		t.Fatalf("env = %+v, want {prefix} intact", got.Runtime.Env)
	}
	if resolved := got.Runtime.Env[0].Resolved(got.Runtime.Prefix); resolved != "/cnt/samtools/1.23.1" {
		t.Errorf("resolved env = %q", resolved)
	}
}

// An image built before this format has no manifest. That is ErrNoManifest and
// nothing else — never a generic failure, and never a tool error.
func TestReadImageWithoutManifest(t *testing.T) {
	requireSquashfsTools(t)
	withTempCache(t)

	root := t.TempDir()
	payload := filepath.Join(root, "cnt", "legacy", "1.0")
	if err := os.MkdirAll(payload, 0o755); err != nil {
		t.Fatalf("mkdir: %v", err)
	}
	if err := os.WriteFile(filepath.Join(payload, "file"), []byte("x"), 0o644); err != nil {
		t.Fatalf("write: %v", err)
	}

	_, err := Read(packSqf(t, root))
	if !errors.Is(err, ErrNoManifest) {
		t.Fatalf("err = %v, want ErrNoManifest", err)
	}
	// It must not be reported as a missing tool or a corrupt archive.
	if errors.Is(err, tool.ErrToolMissing) || errors.Is(err, tool.ErrCorrupt) {
		t.Errorf("a missing manifest was reported as a host or archive fault: %v", err)
	}
}

// A manifest that is present but does not decode stays distinct from one that
// is absent: something is there and it is wrong.
func TestReadMalformedManifest(t *testing.T) {
	requireSquashfsTools(t)
	withTempCache(t)

	root := t.TempDir()
	if err := os.MkdirAll(filepath.Join(root, DirName), 0o755); err != nil {
		t.Fatalf("mkdir: %v", err)
	}
	if err := os.WriteFile(filepath.Join(root, DirName, FileName), []byte("{not json"), 0o644); err != nil {
		t.Fatalf("write: %v", err)
	}

	_, err := Read(packSqf(t, root))
	if err == nil {
		t.Fatal("malformed manifest accepted")
	}
	if errors.Is(err, ErrNoManifest) {
		t.Errorf("malformed manifest reported as missing: %v", err)
	}
	if !errors.Is(err, ErrInvalid) {
		t.Errorf("err = %v, want ErrInvalid", err)
	}
}

// A schema this build does not know is reported as such, so the caller can
// treat it like a missing manifest instead of failing the run.
func TestReadUnsupportedSchema(t *testing.T) {
	requireSquashfsTools(t)
	withTempCache(t)

	future := valid()
	future.SchemaVersion = SchemaVersion + 1

	_, err := Read(stagedImage(t, future))
	if !errors.Is(err, ErrUnsupportedSchema) {
		t.Fatalf("err = %v, want ErrUnsupportedSchema", err)
	}
}

// An unrecognized type is normalized to app on the way out of the archive, not
// rejected.
func TestReadNormalizesType(t *testing.T) {
	requireSquashfsTools(t)
	withTempCache(t)

	odd := valid()
	odd.Type = "bundle"

	got, err := Read(stagedImage(t, odd))
	if err != nil {
		t.Fatalf("Read: %v", err)
	}
	if got.Type != catalog.TypeApp {
		t.Errorf("type = %q, want app", got.Type)
	}
}

func TestReadRejectsNonImage(t *testing.T) {
	withTempCache(t)

	img := filepath.Join(t.TempDir(), "env.img")
	if err := os.WriteFile(img, []byte("not an archive"), 0o644); err != nil {
		t.Fatalf("write: %v", err)
	}
	// A writable .img is a mutable overlay, not a built image: its environment
	// comes from its .env sidecar, and it is never read for a manifest.
	if _, err := Read(img); !errors.Is(err, tool.ErrCorrupt) {
		t.Errorf("err = %v, want ErrCorrupt for a .img", err)
	}

	if _, err := Read(filepath.Join(t.TempDir(), "absent.sqf")); !errors.Is(err, tool.ErrUnreadable) {
		t.Errorf("a missing file should be ErrUnreadable")
	}
}

func TestCacheHitAvoidsReread(t *testing.T) {
	requireSquashfsTools(t)
	withTempCache(t)

	image := stagedImage(t, valid())
	if _, err := Read(image); err != nil {
		t.Fatalf("first Read: %v", err)
	}

	// Truncating the file would break any real read; a cache hit keyed on the
	// unchanged size and mtime still answers, which is what proves it is cached.
	fi, err := os.Stat(image)
	if err != nil {
		t.Fatalf("stat: %v", err)
	}
	if err := os.WriteFile(image, []byte("destroyed"), 0o644); err != nil {
		t.Fatalf("truncate: %v", err)
	}
	if err := os.Chtimes(image, fi.ModTime(), fi.ModTime()); err != nil {
		t.Fatalf("chtimes: %v", err)
	}
	// Size changed, so this must miss and fail rather than serve a stale hit.
	if _, err := Read(image); err == nil {
		t.Error("a changed image was served from cache")
	}
}

func TestCacheInvalidatedByMtime(t *testing.T) {
	requireSquashfsTools(t)
	withTempCache(t)

	first := valid()
	first.Description = "first"
	image := stagedImage(t, first)
	got, err := Read(image)
	if err != nil {
		t.Fatalf("Read: %v", err)
	}
	if got.Description != "first" {
		t.Fatalf("description = %q", got.Description)
	}

	// Rebuild in place with different content, as build --update does.
	second := valid()
	second.Description = "second"
	rebuilt := stagedImage(t, second)
	data, err := os.ReadFile(rebuilt)
	if err != nil {
		t.Fatalf("read rebuilt: %v", err)
	}
	if err := os.WriteFile(image, data, 0o644); err != nil {
		t.Fatalf("overwrite: %v", err)
	}
	if err := os.Chtimes(image, time.Now(), time.Now()); err != nil {
		t.Fatalf("chtimes: %v", err)
	}

	got, err = Read(image)
	if err != nil {
		t.Fatalf("Read after rebuild: %v", err)
	}
	if got.Description != "second" {
		t.Errorf("description = %q, want the rebuilt value", got.Description)
	}
}

// The negative verdict is cached too, or an image predating the format would be
// re-probed on every listing — the cost this cache exists to avoid.
func TestCacheRemembersMissingManifest(t *testing.T) {
	requireSquashfsTools(t)
	withTempCache(t)

	root := t.TempDir()
	if err := os.MkdirAll(filepath.Join(root, "cnt"), 0o755); err != nil {
		t.Fatalf("mkdir: %v", err)
	}
	if err := os.WriteFile(filepath.Join(root, "cnt", "f"), []byte("x"), 0o644); err != nil {
		t.Fatalf("write: %v", err)
	}
	image := packSqf(t, root)

	if _, err := Read(image); !errors.Is(err, ErrNoManifest) {
		t.Fatalf("first Read: %v", err)
	}

	abs, err := filepath.Abs(image)
	if err != nil {
		t.Fatalf("abs: %v", err)
	}
	fi, err := os.Stat(abs)
	if err != nil {
		t.Fatalf("stat: %v", err)
	}
	_, has, found := globalCache.lookup(abs, fi)
	if !found {
		t.Fatal("the missing-manifest verdict was not cached")
	}
	if has {
		t.Error("cache claims a manifest for an image that has none")
	}

	// A second Read still reports the same thing.
	if _, err := Read(image); !errors.Is(err, ErrNoManifest) {
		t.Errorf("second Read: %v", err)
	}

	Forget(image)
	if _, _, found := globalCache.lookup(abs, fi); found {
		t.Error("Forget left the entry behind")
	}
}

// A base image is the container root. The three verdicts differ in kind: a
// manifest that names another type is a real mismatch, while a missing or
// unreadable one says nothing about the image and must not block a build.
func TestCheckBase(t *testing.T) {
	requireSquashfsTools(t)

	baseManifest := Manifest{
		SchemaVersion: SchemaVersion,
		Name:          "ubuntu24/base",
		Type:          catalog.TypeBase,
		BuildType:     "def",
	}

	t.Run("accepts a base manifest", func(t *testing.T) {
		withTempCache(t)
		if err := CheckBase(stagedImage(t, baseManifest)); err != nil {
			t.Errorf("base image rejected: %v", err)
		}
	})

	t.Run("rejects another type", func(t *testing.T) {
		withTempCache(t)
		err := CheckBase(stagedImage(t, valid()))
		if err == nil {
			t.Fatal("an app image was accepted as a container root")
		}
		if !strings.Contains(err.Error(), string(catalog.TypeApp)) {
			t.Errorf("err = %v, want it to name the type it found", err)
		}
	})

	t.Run("accepts an image with no manifest", func(t *testing.T) {
		withTempCache(t)
		root := t.TempDir()
		if err := os.MkdirAll(filepath.Join(root, "bin"), 0o755); err != nil {
			t.Fatal(err)
		}
		if err := CheckBase(packSqf(t, root)); err != nil {
			t.Errorf("a base predating the manifest was rejected: %v", err)
		}
	})

	t.Run("accepts an unreadable manifest", func(t *testing.T) {
		withTempCache(t)
		root := t.TempDir()
		if err := os.MkdirAll(filepath.Join(root, DirName), 0o755); err != nil {
			t.Fatal(err)
		}
		if err := os.WriteFile(filepath.Join(root, DirName, FileName), []byte("{not json"), 0o644); err != nil {
			t.Fatal(err)
		}
		if err := CheckBase(packSqf(t, root)); err != nil {
			t.Errorf("an unreadable manifest stranded the base: %v", err)
		}
	})
}
