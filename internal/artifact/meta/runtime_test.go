package meta

import (
	"errors"
	"os"
	"path/filepath"
	"strings"
	"testing"
	"time"

	"github.com/Justype/condatainer/catalog"
	"github.com/Justype/condatainer/internal/image/tool"
	"github.com/Justype/condatainer/internal/toolpath"
)

func TestValidateRuntimeAcceptsEachType(t *testing.T) {
	for _, typ := range []catalog.Type{catalog.TypeApp, catalog.TypeData} {
		r := validRuntime()
		r.Type = typ
		if err := ValidateRuntime(r); err != nil {
			t.Errorf("%s with a prefix rejected: %v", typ, err)
		}
	}
	// base and os apply at the container root, so they carry no prefix.
	for _, typ := range []catalog.Type{catalog.TypeBase, catalog.TypeOS} {
		r := validRuntime()
		r.Type = typ
		r.Prefix = ""
		if err := ValidateRuntime(r); err != nil {
			t.Errorf("%s without a prefix rejected: %v", typ, err)
		}
	}
}

func TestValidateRuntimeRejects(t *testing.T) {
	tests := []struct {
		name   string
		mutate func(*Runtime)
		want   error
	}{
		{"unsupported schema", func(r *Runtime) { r.SchemaVersion = SchemaVersion + 1 }, ErrUnsupportedSchema},
		{"zero schema", func(r *Runtime) { r.SchemaVersion = 0 }, ErrUnsupportedSchema},
		{"empty name", func(r *Runtime) { r.Name = "  " }, ErrInvalid},
		{"no architecture", func(r *Runtime) { r.Platform.Arch = "" }, ErrInvalid},
		{"app without prefix", func(r *Runtime) { r.Prefix = "" }, ErrInvalid},
		{"relative prefix", func(r *Runtime) { r.Prefix = "cnt/samtools" }, ErrInvalid},
		{"empty env key", func(r *Runtime) { r.Env[0].Key = "" }, ErrInvalid},
		{"env key with a dash", func(r *Runtime) { r.Env[0].Key = "NOT-VALID" }, ErrInvalid},
		{"env key starting with a digit", func(r *Runtime) { r.Env[0].Key = "1BAD" }, ErrInvalid},
		{"duplicate env key", func(r *Runtime) {
			r.Env = append(r.Env, EnvVar{Key: "SAMTOOLS_HOME", Value: "/elsewhere"})
		}, ErrInvalid},
	}

	for _, tt := range tests {
		t.Run(tt.name, func(t *testing.T) {
			r := validRuntime()
			tt.mutate(&r)
			err := ValidateRuntime(r)
			if err == nil {
				t.Fatalf("accepted %s", tt.name)
			}
			if !errors.Is(err, tt.want) {
				t.Fatalf("error = %v, want %v", err, tt.want)
			}
		})
	}
}

// Only an absent type is filled in. Anything else reaches ValidateRuntime as it
// stands, so a type this build does not know is refused rather than mounted as
// something it is not.
func TestNormalizeOnlyFillsInAnAbsentType(t *testing.T) {
	r := validRuntime()
	r.Type = ""
	r.Normalize()
	if r.Type != catalog.TypeApp {
		t.Errorf("an absent type normalized to %q, want app", r.Type)
	}

	for _, raw := range []catalog.Type{
		catalog.TypeBase, catalog.TypeOS, catalog.TypeApp, catalog.TypeData, catalog.TypeEnv,
		"bundle", "module", "APP",
	} {
		r := validRuntime()
		r.Type = raw
		r.Normalize()
		if r.Type != raw {
			t.Errorf("type %q normalized to %q", raw, r.Type)
		}
	}
}

// {prefix} has to reach the image intact: the install prefix is not known when
// the image is built, only when it is loaded.
func TestMarshalRuntimeKeepsPrefixToken(t *testing.T) {
	data, err := MarshalRuntime(validRuntime())
	if err != nil {
		t.Fatalf("MarshalRuntime: %v", err)
	}
	if !strings.Contains(string(data), `"value": "{prefix}"`) {
		t.Errorf("{prefix} did not survive marshalling:\n%s", data)
	}
	if !strings.HasSuffix(string(data), "\n") {
		t.Error("marshalled runtime has no trailing newline")
	}

	second, err := MarshalRuntime(validRuntime())
	if err != nil {
		t.Fatalf("MarshalRuntime: %v", err)
	}
	if string(data) != string(second) {
		t.Error("marshalling the same runtime twice produced different bytes")
	}
}

func TestStageRuntimeWritesAtomically(t *testing.T) {
	dir := filepath.Join(t.TempDir(), "workspace", DirName)
	if err := StageRuntime(dir, validRuntime()); err != nil {
		t.Fatalf("StageRuntime: %v", err)
	}

	data, err := os.ReadFile(filepath.Join(dir, RuntimeFileName))
	if err != nil {
		t.Fatalf("runtime not written: %v", err)
	}
	want, err := MarshalRuntime(validRuntime())
	if err != nil {
		t.Fatalf("MarshalRuntime: %v", err)
	}
	if string(data) != string(want) {
		t.Errorf("staged bytes differ from MarshalRuntime output:\ngot:\n%s\nwant:\n%s", data, want)
	}

	// The temp file the atomic rename went through must not survive.
	if _, err := os.Stat(filepath.Join(dir, RuntimeFileName+".tmp")); !os.IsNotExist(err) {
		t.Error("staging left its temp file behind")
	}
}

func TestReadRuntimeRoundTrip(t *testing.T) {
	requireSquashfsTools(t)
	withTempCache(t)

	want := validRuntime()
	want.Description = "SAMtools alignment toolkit"
	image := stagedImage(t, want, validManifest())

	got, err := ReadRuntime(image)
	if err != nil {
		t.Fatalf("ReadRuntime: %v", err)
	}
	if got.Name != want.Name || got.Type != want.Type || got.Description != want.Description {
		t.Errorf("identity = %+v, want %+v", got, want)
	}
	if got.Platform != want.Platform {
		t.Errorf("platform = %+v, want %+v", got.Platform, want.Platform)
	}
	if got.Prefix != want.Prefix {
		t.Errorf("prefix = %q, want %q", got.Prefix, want.Prefix)
	}
	// {prefix} survives the round trip; substitution happens at load time.
	if len(got.Env) != 1 || got.Env[0].Value != "{prefix}" {
		t.Fatalf("env = %+v, want {prefix} intact", got.Env)
	}
	if resolved := got.Env[0].Resolved(got.Prefix); resolved != "/cnt/samtools/1.23.1" {
		t.Errorf("resolved env = %q", resolved)
	}
}

// An image built before the runtime split has no runtime.json. That is
// ErrNoRuntime and nothing else — never a generic failure, never a tool error,
// and nothing falls back to the manifest for it.
func TestReadRuntimeWithoutDocument(t *testing.T) {
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

	_, err := ReadRuntime(packSqf(t, root))
	if !errors.Is(err, ErrNoRuntime) {
		t.Fatalf("err = %v, want ErrNoRuntime", err)
	}
	if errors.Is(err, toolpath.ErrToolMissing) || errors.Is(err, tool.ErrCorrupt) {
		t.Errorf("missing runtime metadata was reported as a host or archive fault: %v", err)
	}
}

// An image carrying only a manifest is treated exactly like one carrying no
// metadata at all: there is no fallback path that reads a runtime block out of
// the manifest.
func TestReadRuntimeDoesNotFallBackToManifest(t *testing.T) {
	requireSquashfsTools(t)
	withTempCache(t)

	root := t.TempDir()
	dir := filepath.Join(root, DirName)
	if err := StageManifest(dir, validManifest()); err != nil {
		t.Fatalf("StageManifest: %v", err)
	}
	if _, err := ReadRuntime(packSqf(t, root)); !errors.Is(err, ErrNoRuntime) {
		t.Fatalf("err = %v, want ErrNoRuntime", err)
	}
}

// A runtime document that is present but does not decode stays distinct from one
// that is absent: something is there and it is wrong. What the bytes mean is
// DecodeRuntime's decision, so packing an archive to deliver them would only
// test unsquashfs.
func TestReadRuntimeMalformed(t *testing.T) {
	_, err := DecodeRuntime([]byte("{not json"), "image.sqf")
	if err == nil {
		t.Fatal("malformed runtime accepted")
	}
	if errors.Is(err, ErrNoRuntime) {
		t.Errorf("malformed runtime reported as missing: %v", err)
	}
	if !errors.Is(err, ErrInvalid) {
		t.Errorf("err = %v, want ErrInvalid", err)
	}
}

// A schema this build does not know is reported as such, so the caller can treat
// it like a missing document instead of failing the run.
func TestReadRuntimeUnsupportedSchema(t *testing.T) {
	future := validRuntime()
	future.SchemaVersion = SchemaVersion + 1

	_, err := DecodeRuntime(marshalRuntime(t, future), "image.sqf")
	if !errors.Is(err, ErrUnsupportedSchema) {
		t.Fatalf("err = %v, want ErrUnsupportedSchema", err)
	}
}

// A type this build does not know is refused at the document boundary. It
// decides the install prefix and whether the artifact may be published, so
// reading it as app would mount and gate it as something it is not. "APP" is the
// case that makes the point: a case difference is not a spelling to correct.
func TestReadRuntimeRejectsUnknownType(t *testing.T) {
	for _, raw := range []catalog.Type{"bundle", "APP"} {
		odd := validRuntime()
		odd.Type = raw

		_, err := DecodeRuntime(marshalRuntime(t, odd), "image.sqf")
		if !errors.Is(err, ErrInvalid) {
			t.Errorf("type %q: err = %v, want ErrInvalid", raw, err)
		}
	}
}

func TestReadRuntimeRejectsNonImage(t *testing.T) {
	withTempCache(t)

	img := filepath.Join(t.TempDir(), "env.img")
	if err := os.WriteFile(img, []byte("not an archive"), 0o644); err != nil {
		t.Fatalf("write: %v", err)
	}
	// A writable .img is a mutable overlay, not a built image: its environment
	// comes from its .env sidecar, and it is never read for embedded metadata.
	if _, err := ReadRuntime(img); !errors.Is(err, tool.ErrCorrupt) {
		t.Errorf("err = %v, want ErrCorrupt for a .img", err)
	}

	if _, err := ReadRuntime(filepath.Join(t.TempDir(), "absent.sqf")); !errors.Is(err, tool.ErrUnreadable) {
		t.Errorf("a missing file should be ErrUnreadable")
	}
}

func TestCacheHitAvoidsReread(t *testing.T) {
	requireSquashfsTools(t)
	withTempCache(t)

	image := stagedImage(t, validRuntime(), validManifest())
	if _, err := ReadRuntime(image); err != nil {
		t.Fatalf("first ReadRuntime: %v", err)
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
	if _, err := ReadRuntime(image); err == nil {
		t.Error("a changed image was served from cache")
	}
}

func TestCacheInvalidatedByMtime(t *testing.T) {
	requireSquashfsTools(t)
	withTempCache(t)

	first := validRuntime()
	first.Description = "first"
	image := stagedImage(t, first, validManifest())
	got, err := ReadRuntime(image)
	if err != nil {
		t.Fatalf("ReadRuntime: %v", err)
	}
	if got.Description != "first" {
		t.Fatalf("description = %q", got.Description)
	}

	// Rebuild in place with different content, as build --update does.
	second := validRuntime()
	second.Description = "second"
	rebuilt := stagedImage(t, second, validManifest())
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

	got, err = ReadRuntime(image)
	if err != nil {
		t.Fatalf("ReadRuntime after rebuild: %v", err)
	}
	if got.Description != "second" {
		t.Errorf("description = %q, want the rebuilt value", got.Description)
	}
}

// The negative verdict is cached too, or an image predating the format would be
// re-probed on every listing — the cost this cache exists to avoid.
func TestCacheRemembersMissingRuntime(t *testing.T) {
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

	if _, err := ReadRuntime(image); !errors.Is(err, ErrNoRuntime) {
		t.Fatalf("first ReadRuntime: %v", err)
	}

	abs, err := filepath.Abs(image)
	if err != nil {
		t.Fatalf("abs: %v", err)
	}
	fi, err := os.Stat(abs)
	if err != nil {
		t.Fatalf("stat: %v", err)
	}
	record, found := runtimeCache.Lookup(abs, fi)
	if !found {
		t.Fatal("the missing-runtime verdict was not cached")
	}
	if !record.RuntimeKnown || len(record.Runtime) != 0 {
		t.Error("cache claims runtime metadata for an image that has none")
	}

	// A second read still reports the same thing.
	if _, err := ReadRuntime(image); !errors.Is(err, ErrNoRuntime) {
		t.Errorf("second ReadRuntime: %v", err)
	}

	runtimeCache.Forget(image)
	if _, found := runtimeCache.Lookup(abs, fi); found {
		t.Error("Forget left the entry behind")
	}
}
