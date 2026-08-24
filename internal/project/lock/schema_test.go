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
		l.Selections["star/2.7.11b"] = Selection{Artifact: "provenance/star--2.7.11b@a31f902c12ab"}
		l.Selections["cutadapt/5.0"] = Selection{Artifact: "provenance/cutadapt--5.0@41b5204f99c1"}
		l.Selections["grch38/genome/gencode49"] = Selection{Artifact: "provenance/grch38--genome--gencode49@8ce02116f302"}
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
	if strings.Contains(string(first), "remotes") {
		t.Errorf("an empty remote map was serialized:\n%s", first)
	}
}

func TestRoundTripPreservesRemotes(t *testing.T) {
	l := New()
	l.Selections["star/2.7.11b"] = Selection{Artifact: "provenance/star--2.7.11b@a31f902c12ab"}
	if err := l.AddRemote("provenance/star--2.7.11b@a31f902c12ab", Remote{Repository: "ghcr.io/example/star", ManifestDigest: digestA}); err != nil {
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
	got := back.Remotes["provenance/star--2.7.11b@a31f902c12ab"]
	if len(got) != 1 || got[0].Repository != "ghcr.io/example/star" || got[0].ManifestDigest != digestA {
		t.Fatalf("remotes did not round-trip: %#v", got)
	}
}

// Order is retry priority, so a repeat must not reorder and a duplicate must
// not accumulate.
func TestAddRemoteDedupesAndKeepsOrder(t *testing.T) {
	l := New()
	artifact := "provenance/star--2.7.11b@a31f902c12ab"
	first := Remote{Repository: "ghcr.io/example/star", ManifestDigest: digestA}
	second := Remote{Repository: "ghcr.io/mirror/star", ManifestDigest: digestA}
	for _, remote := range []Remote{first, second, first} {
		if err := l.AddRemote(artifact, remote); err != nil {
			t.Fatal(err)
		}
	}
	got := l.Remotes[artifact]
	if len(got) != 2 || got[0] != first || got[1] != second {
		t.Fatalf("remotes = %#v", got)
	}
}

func TestAddRemoteRejectsBadLocations(t *testing.T) {
	artifact := "provenance/star--2.7.11b@a31f902c12ab"
	tests := []struct {
		why    string
		remote Remote
	}{
		{"a scheme is not part of a repository", Remote{Repository: "oci://ghcr.io/example/star", ManifestDigest: digestA}},
		{"a tag is mutable and is not an remote", Remote{Repository: "ghcr.io/example/star:latest", ManifestDigest: digestA}},
		{"an empty repository addresses nothing", Remote{Repository: "", ManifestDigest: digestA}},
		{"a tag is not a digest", Remote{Repository: "ghcr.io/example/star", ManifestDigest: "latest"}},
		{"a short digest is not sha256", Remote{Repository: "ghcr.io/example/star", ManifestDigest: "sha256:abcd"}},
		{"a digest is hex", Remote{Repository: "ghcr.io/example/star", ManifestDigest: "sha256:" + strings.Repeat("z", 64)}},
	}
	for _, tt := range tests {
		t.Run(tt.why, func(t *testing.T) {
			if err := New().AddRemote(artifact, tt.remote); err == nil {
				t.Fatalf("accepted %#v", tt.remote)
			}
		})
	}
}

// A registry may carry a port, which is a colon that does not introduce a tag.
func TestAddRemoteAcceptsARegistryPort(t *testing.T) {
	if err := New().AddRemote("provenance/star--2.7.11b@a31f902c12ab",
		Remote{Repository: "localhost:5000/example/star", ManifestDigest: digestA}); err != nil {
		t.Fatalf("rejected a ported registry: %v", err)
	}
}

// Every lock path is hostile input: it names a directory a later step reads.
func TestEntryPathsAreConstrained(t *testing.T) {
	tests := []struct {
		why  string
		path string
	}{
		{"absolute", "/etc/passwd"},
		{"traversal", "provenance/../../etc"},
		{"outside artifacts", "elsewhere/star@abc"},
		{"unclean", "provenance/./star@abc"},
		{"nested", "provenance/deep/star@abc"},
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
		{"unknown selection field", `{"schema_version":1,"selections":{"a":{"artifact":"provenance/a@b","pinned":true}}}`, ErrInvalid},
		{"future schema", `{"schema_version":99,"selections":{}}`, ErrSchema},
		{"trailing content", `{"schema_version":1,"selections":{}} {}`, ErrInvalid},
		{"remote outside artifacts", `{"schema_version":1,"selections":{},"remotes":{"/abs":[]}}`, ErrInvalid},
	}
	for _, tt := range tests {
		t.Run(tt.why, func(t *testing.T) {
			if _, err := Unmarshal([]byte(tt.data)); !errors.Is(err, tt.want) {
				t.Fatalf("error = %v, want %v", err, tt.want)
			}
		})
	}
}

// A duplicated remote in a hand-edited lock is a conflict, not something to
// silently collapse on the next write.
func TestUnmarshalRejectsDuplicateRemotes(t *testing.T) {
	data := `{"schema_version":1,"selections":{},"remotes":{"provenance/a@b":[` +
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
	l.Selections["star/2.7.11b"] = Selection{Artifact: "provenance/star--2.7.11b@a31f902c12ab"}
	l.Selections["also/star"] = Selection{Artifact: "provenance/star--2.7.11b@a31f902c12ab"}
	if got := l.SelectedArtifacts(); len(got) != 1 || !got["provenance/star--2.7.11b@a31f902c12ab"] {
		t.Fatalf("SelectedArtifacts = %#v", got)
	}
}
