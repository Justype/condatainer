package key

import (
	"errors"
	"os"
	"path/filepath"
	"strings"
	"testing"

	"github.com/Justype/condatainer/catalog"
	"github.com/Justype/condatainer/internal/artifact/meta"
)

func snapshotManifest() meta.Manifest {
	return meta.Manifest{
		SchemaVersion: meta.SchemaVersion,
		Name:          meta.EnvName,
		Type:          catalog.TypeEnv,
		BuildType:     meta.BuildTypeSnapshot,
		Platform:      meta.Platform{Arch: "amd64"},
	}
}

func writeFile(t *testing.T, name string, data []byte) string {
	t.Helper()
	path := filepath.Join(t.TempDir(), name)
	if err := os.WriteFile(path, data, 0o644); err != nil {
		t.Fatal(err)
	}
	return path
}

// The digest is over the file, in the OCI form a lock, the store and a registry
// descriptor all share.
func TestArtifactDigest(t *testing.T) {
	path := writeFile(t, "a.sqf", []byte("payload"))
	got, err := ArtifactDigest(path)
	if err != nil {
		t.Fatal(err)
	}
	// sha256("payload")
	want := "sha256:239f59ed55e737c77147cf55ad0c1b030b6d7ee748a7426952f9b852d5a935e5"
	if got != want {
		t.Errorf("digest = %s, want %s", got, want)
	}
}

// One byte anywhere changes the identity — the whole point of addressing an
// artifact nothing can regenerate by its content.
func TestArtifactDigestChangesWithContent(t *testing.T) {
	a, err := ArtifactDigest(writeFile(t, "a.sqf", []byte("payload")))
	if err != nil {
		t.Fatal(err)
	}
	b, err := ArtifactDigest(writeFile(t, "b.sqf", []byte("payloae")))
	if err != nil {
		t.Fatal(err)
	}
	if a == b {
		t.Fatal("different bytes produced the same digest")
	}
}

func rec(typ byte, mode uint32, p, id string) TreeRecord {
	return TreeRecord{Type: typ, Mode: mode, Path: p, ID: id}
}

func TestPayloadKeyNamesTheScheme(t *testing.T) {
	ref, err := PayloadKey([]TreeRecord{rec('f', 0o644, "/cnt_env/a", "ab")})
	if err != nil {
		t.Fatal(err)
	}
	if ref.Scheme != string(PayloadTreeV1) {
		t.Errorf("scheme = %q, want %q", ref.Scheme, PayloadTreeV1)
	}
	if len(ref.SHA256) != 64 {
		t.Errorf("sha256 = %q", ref.SHA256)
	}
}

// The ordering is the scheme's, not the walker's: an archive read in a different
// order is the same archive.
func TestPayloadKeyIgnoresRecordOrder(t *testing.T) {
	a := []TreeRecord{
		rec('d', 0o755, "/cnt_env", ""),
		rec('f', 0o644, "/cnt_env/a", "aa"),
		rec('f', 0o755, "/cnt_env/b", "bb"),
	}
	shuffled := []TreeRecord{a[2], a[0], a[1]}
	first, err := PayloadKey(a)
	if err != nil {
		t.Fatal(err)
	}
	second, err := PayloadKey(shuffled)
	if err != nil {
		t.Fatal(err)
	}
	if first != second {
		t.Errorf("order changed the identity: %v vs %v", first, second)
	}
}

// Each field is in the preimage because it changes what the environment does.
// Mode is the one worth a test: identical bytes that cannot be executed are a
// different environment.
func TestPayloadKeyDistinguishesEveryField(t *testing.T) {
	base := []TreeRecord{rec('f', 0o755, "/cnt_env/bin/tool", "aa")}
	want, err := PayloadKey(base)
	if err != nil {
		t.Fatal(err)
	}
	for name, changed := range map[string][]TreeRecord{
		"content": {rec('f', 0o755, "/cnt_env/bin/tool", "bb")},
		"mode":    {rec('f', 0o644, "/cnt_env/bin/tool", "aa")},
		"path":    {rec('f', 0o755, "/cnt_env/bin/other", "aa")},
		"type":    {rec('l', 0o755, "/cnt_env/bin/tool", "aa")},
	} {
		t.Run(name, func(t *testing.T) {
			got, err := PayloadKey(changed)
			if err != nil {
				t.Fatal(err)
			}
			if got == want {
				t.Errorf("changing %s did not change the identity", name)
			}
		})
	}
}

// A path holding the field separator must not be able to shift the framing into
// looking like a different tree.
func TestPayloadKeyFramesHostilePaths(t *testing.T) {
	a, err := PayloadKey([]TreeRecord{rec('f', 0o644, "/cnt_env/a b", "cc")})
	if err != nil {
		t.Fatal(err)
	}
	b, err := PayloadKey([]TreeRecord{rec('f', 0o644, "/cnt_env/a", "b\x00cc")})
	if err != nil {
		t.Fatal(err)
	}
	if a == b {
		t.Error("a path containing a separator collided with another tree")
	}
}

func TestPayloadKeyRefusesAnEmptyTree(t *testing.T) {
	if _, err := PayloadKey(nil); err == nil {
		t.Error("an empty archive was given an identity")
	}
}

// A snapshot has no derivable keys, and the entry points say so with a sentinel
// a caller can branch on rather than a message about an unknown build type.
func TestSnapshotIsUnkeyed(t *testing.T) {
	m := snapshotManifest()
	if !IsSnapshot(m) {
		t.Fatal("IsSnapshot did not recognise the build type")
	}
	if _, err := ReadSources(t.TempDir(), m); !errors.Is(err, ErrUnkeyed) {
		t.Errorf("ReadSources error = %v, want ErrUnkeyed", err)
	}
	if _, err := Regenerate(m, nil); !errors.Is(err, ErrUnkeyed) {
		t.Errorf("Regenerate error = %v, want ErrUnkeyed", err)
	}
	if _, _, err := latest(meta.BuildTypeSnapshot); !errors.Is(err, ErrUnkeyed) {
		t.Errorf("latest error = %v, want ErrUnkeyed", err)
	}
}

// A snapshot manifest carries no keys at all, which the manifest rules already
// permit: identity and equivalence must be both present or both absent.
func TestSnapshotManifestValidatesWithoutKeys(t *testing.T) {
	if err := meta.ValidateManifest(snapshotManifest()); err != nil {
		t.Fatalf("a keyless snapshot manifest was rejected: %v", err)
	}
}

// A frozen environment is always called env. It is addressed by path and never
// enters a namespace where a name would distinguish it, so a manifest claiming
// another one is refused rather than carried.
func TestSnapshotManifestRefusesAnotherName(t *testing.T) {
	m := snapshotManifest()
	m.Name = "rnaseq/1.0"
	if err := meta.ValidateManifest(m); err == nil {
		t.Fatal("a frozen environment claiming its own name was accepted")
	}
}

// The metadata directory is left out, so the manifest that carries the key can
// live in the archive it describes; the key does not depend on how many
// goroutines hashed.
func TestTreeOfNamesEntriesAndSkipsTheMetadataDirectory(t *testing.T) {
	root := t.TempDir()
	for name, data := range map[string]string{"cnt/app/tool": "x", "cnt/app/data": "y", ".cnt/manifest.json": "{}"} {
		p := filepath.Join(root, name)
		if err := os.MkdirAll(filepath.Dir(p), 0o755); err != nil {
			t.Fatal(err)
		}
		if err := os.WriteFile(p, []byte(data), 0o644); err != nil {
			t.Fatal(err)
		}
	}

	one, err := TreeOf(t.Context(), root, "", 1)
	if err != nil {
		t.Fatal(err)
	}
	many, err := TreeOf(t.Context(), root, "", 8)
	if err != nil {
		t.Fatal(err)
	}
	if PayloadPreimage(one) != PayloadPreimage(many) {
		t.Error("the key depends on the worker count")
	}
	for _, r := range one {
		if r.Path == "/.cnt" || strings.HasPrefix(r.Path, "/.cnt/") {
			t.Errorf("the metadata directory was keyed: %s", r.Path)
		}
	}
	if len(one) != 4 { // /cnt, /cnt/app, and its two files
		t.Errorf("got %d records, want 4: %v", len(one), one)
	}
}
