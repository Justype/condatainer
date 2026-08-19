package lock

import (
	"errors"
	"strings"
	"testing"
)

const digestA = "sha256:0123456789abcdef0123456789abcdef0123456789abcdef0123456789abcdef"

func TestMarshalIsDeterministicAndSorted(t *testing.T) {
	build := func() *Lock {
		l := New()
		l.Selections["star/2.7.11b"] = Selection{Artifact: "artifacts/star--2.7.11b@a31f902c12ab"}
		l.Selections["cutadapt/5.0"] = Selection{Artifact: "artifacts/cutadapt--5.0@41b5204f99c1"}
		l.Selections["grch38/genome/gencode49"] = Selection{Artifact: "artifacts/grch38--genome--gencode49@8ce02116f302"}
		return l
	}

	first, err := build().Marshal()
	if err != nil {
		t.Fatal(err)
	}
	second, err := build().Marshal()
	if err != nil {
		t.Fatal(err)
	}
	if string(first) != string(second) {
		t.Fatalf("two runs produced different bytes:\n%s\n---\n%s", first, second)
	}
	if !strings.HasSuffix(string(first), "}\n") {
		t.Errorf("output does not end in exactly one newline: %q", string(first[len(first)-3:]))
	}
	if !strings.Contains(string(first), "\n  \"selections\": {") {
		t.Errorf("output is not two-space indented:\n%s", first)
	}

	keys := []string{"cutadapt/5.0", "grch38/genome/gencode49", "star/2.7.11b"}
	at := -1
	for _, key := range keys {
		i := strings.Index(string(first), key)
		if i < at {
			t.Fatalf("selection keys are not sorted:\n%s", first)
		}
		at = i
	}
	if strings.Contains(string(first), "origins") {
		t.Errorf("an empty origin map was serialized:\n%s", first)
	}
}

func TestRoundTripPreservesOrigins(t *testing.T) {
	l := New()
	l.Selections["star/2.7.11b"] = Selection{Artifact: "artifacts/star--2.7.11b@a31f902c12ab"}
	if err := l.AddOrigin("artifacts/star--2.7.11b@a31f902c12ab", Origin{Repository: "ghcr.io/example/star", ManifestDigest: digestA}); err != nil {
		t.Fatal(err)
	}
	data, err := l.Marshal()
	if err != nil {
		t.Fatal(err)
	}
	back, err := Unmarshal(data)
	if err != nil {
		t.Fatal(err)
	}
	got := back.Origins["artifacts/star--2.7.11b@a31f902c12ab"]
	if len(got) != 1 || got[0].Repository != "ghcr.io/example/star" || got[0].ManifestDigest != digestA {
		t.Fatalf("origins did not round-trip: %#v", got)
	}
}

// Order is retry priority, so a repeat must not reorder and a duplicate must
// not accumulate.
func TestAddOriginDedupesAndKeepsOrder(t *testing.T) {
	l := New()
	artifact := "artifacts/star--2.7.11b@a31f902c12ab"
	first := Origin{Repository: "ghcr.io/example/star", ManifestDigest: digestA}
	second := Origin{Repository: "ghcr.io/mirror/star", ManifestDigest: digestA}
	for _, origin := range []Origin{first, second, first} {
		if err := l.AddOrigin(artifact, origin); err != nil {
			t.Fatal(err)
		}
	}
	got := l.Origins[artifact]
	if len(got) != 2 || got[0] != first || got[1] != second {
		t.Fatalf("origins = %#v", got)
	}
}

func TestAddOriginRejectsBadLocations(t *testing.T) {
	artifact := "artifacts/star--2.7.11b@a31f902c12ab"
	tests := []struct {
		why    string
		origin Origin
	}{
		{"a scheme is not part of a repository", Origin{Repository: "oci://ghcr.io/example/star", ManifestDigest: digestA}},
		{"a tag is mutable and is not an origin", Origin{Repository: "ghcr.io/example/star:latest", ManifestDigest: digestA}},
		{"an empty repository addresses nothing", Origin{Repository: "", ManifestDigest: digestA}},
		{"a tag is not a digest", Origin{Repository: "ghcr.io/example/star", ManifestDigest: "latest"}},
		{"a short digest is not sha256", Origin{Repository: "ghcr.io/example/star", ManifestDigest: "sha256:abcd"}},
		{"a digest is hex", Origin{Repository: "ghcr.io/example/star", ManifestDigest: "sha256:" + strings.Repeat("z", 64)}},
	}
	for _, tt := range tests {
		t.Run(tt.why, func(t *testing.T) {
			if err := New().AddOrigin(artifact, tt.origin); err == nil {
				t.Fatalf("accepted %#v", tt.origin)
			}
		})
	}
}

// A registry may carry a port, which is a colon that does not introduce a tag.
func TestAddOriginAcceptsARegistryPort(t *testing.T) {
	if err := New().AddOrigin("artifacts/star--2.7.11b@a31f902c12ab",
		Origin{Repository: "localhost:5000/example/star", ManifestDigest: digestA}); err != nil {
		t.Fatalf("rejected a ported registry: %v", err)
	}
}

// Every lock path is hostile input: it names a directory a later step reads.
func TestArtifactPathsAreConstrained(t *testing.T) {
	tests := []struct {
		why  string
		path string
	}{
		{"absolute", "/etc/passwd"},
		{"traversal", "artifacts/../../etc"},
		{"outside artifacts", "elsewhere/star@abc"},
		{"unclean", "artifacts/./star@abc"},
		{"nested", "artifacts/deep/star@abc"},
		{"empty", ""},
		{"backslash", `artifacts\star@abc`},
	}
	for _, tt := range tests {
		t.Run(tt.why, func(t *testing.T) {
			l := New()
			l.Selections["star/2.7.11b"] = Selection{Artifact: tt.path}
			if _, err := l.Marshal(); !errors.Is(err, ErrInvalid) {
				t.Fatalf("accepted %q: %v", tt.path, err)
			}
		})
	}
}

func TestUnmarshalRejectsUnknownFieldsAndSchemas(t *testing.T) {
	tests := []struct {
		why  string
		data string
		want error
	}{
		{"unknown top-level field", `{"schema_version":1,"selections":{},"surprise":true}`, ErrInvalid},
		{"unknown selection field", `{"schema_version":1,"selections":{"a":{"artifact":"artifacts/a@b","pinned":true}}}`, ErrInvalid},
		{"future schema", `{"schema_version":99,"selections":{}}`, ErrSchema},
		{"trailing content", `{"schema_version":1,"selections":{}} {}`, ErrInvalid},
		{"origin outside artifacts", `{"schema_version":1,"selections":{},"origins":{"/abs":[]}}`, ErrInvalid},
	}
	for _, tt := range tests {
		t.Run(tt.why, func(t *testing.T) {
			if _, err := Unmarshal([]byte(tt.data)); !errors.Is(err, tt.want) {
				t.Fatalf("error = %v, want %v", err, tt.want)
			}
		})
	}
}

// A duplicated origin in a hand-edited lock is a conflict, not something to
// silently collapse on the next write.
func TestUnmarshalRejectsDuplicateOrigins(t *testing.T) {
	data := `{"schema_version":1,"selections":{},"origins":{"artifacts/a@b":[` +
		`{"repository":"ghcr.io/x/y","manifest_digest":"` + digestA + `"},` +
		`{"repository":"ghcr.io/x/y","manifest_digest":"` + digestA + `"}]}}`
	if _, err := Unmarshal([]byte(data)); !errors.Is(err, ErrInvalid) {
		t.Fatalf("error = %v, want ErrInvalid", err)
	}
}

// Closure directories are reached through manifest edges, so the selection set
// is only the root set. Pruning against it alone would delete the closure.
func TestSelectedArtifactsIsRootsOnly(t *testing.T) {
	l := New()
	l.Selections["star/2.7.11b"] = Selection{Artifact: "artifacts/star--2.7.11b@a31f902c12ab"}
	l.Selections["also/star"] = Selection{Artifact: "artifacts/star--2.7.11b@a31f902c12ab"}
	if got := l.SelectedArtifacts(); len(got) != 1 || !got["artifacts/star--2.7.11b@a31f902c12ab"] {
		t.Fatalf("SelectedArtifacts = %#v", got)
	}
}
