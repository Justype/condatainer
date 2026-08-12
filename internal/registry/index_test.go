package registry

import (
	"encoding/json"
	"strings"
	"testing"

	"github.com/opencontainers/go-digest"
	ocispec "github.com/opencontainers/image-spec/specs-go/v1"
)

func platformDescriptor(os, arch, hex string) ocispec.Descriptor {
	return ocispec.Descriptor{
		MediaType: ocispec.MediaTypeImageManifest,
		Digest:    digest.Digest(digestPrefix + strings.Repeat(hex, 64)),
		Platform:  &ocispec.Platform{OS: os, Architecture: arch},
	}
}

func TestPlatformKey(t *testing.T) {
	tests := []struct {
		in   *ocispec.Platform
		want string
	}{
		{nil, ""},
		{&ocispec.Platform{OS: "linux", Architecture: "amd64"}, "linux/amd64"},
		{&ocispec.Platform{OS: "linux", Architecture: "arm64", Variant: "v8"}, "linux/arm64/v8"},
	}
	for _, tt := range tests {
		if got := platformKey(tt.in); got != tt.want {
			t.Errorf("platformKey(%+v) = %q, want %q", tt.in, got, tt.want)
		}
	}
}

// The index platform is Go's architecture spelling, not the uname spelling
// meta.NativeArch reports. Two vocabularies, and the wire format picks this one.
func TestPlatformUsesTheOCISpelling(t *testing.T) {
	plat, ok := platform()
	if !ok {
		t.Skip("not a distribution architecture")
	}
	if plat.OS != "linux" {
		t.Errorf("OS = %q, want linux", plat.OS)
	}
	if plat.Architecture == "x86_64" || plat.Architecture == "aarch64" {
		t.Errorf("architecture = %q, which is the uname spelling", plat.Architecture)
	}
}

func TestMergeAndContainsIndexEntries(t *testing.T) {
	amd := platformDescriptor("linux", "amd64", "a")
	arm := platformDescriptor("linux", "arm64", "b")
	entries := map[string]ocispec.Descriptor{"linux/amd64": amd}

	if containsIndexEntries(entries, map[string]ocispec.Descriptor{"linux/arm64": arm}) {
		t.Fatal("a missing arm64 entry was reported present")
	}
	mergeIndexEntries(entries, map[string]ocispec.Descriptor{"linux/arm64": arm})
	if !containsIndexEntries(entries, map[string]ocispec.Descriptor{
		"linux/amd64": amd,
		"linux/arm64": arm,
	}) {
		t.Fatal("merged platform entries not found")
	}
}

// A platform present at a different digest is somebody else's push landing on
// top of this one, not a match.
func TestContainsIndexEntriesRejectsAReplacedPlatform(t *testing.T) {
	old := platformDescriptor("linux", "amd64", "a")
	replacement := platformDescriptor("linux", "amd64", "c")
	if containsIndexEntries(
		map[string]ocispec.Descriptor{"linux/amd64": replacement},
		map[string]ocispec.Descriptor{"linux/amd64": old},
	) {
		t.Fatal("a different digest for the same platform was reported present")
	}
}

// The merge direction is what makes the local push win: the pusher inserts its
// own child first, then folds in what the registry already had.
func TestMergeIndexEntriesKeepsTheDestination(t *testing.T) {
	mine := platformDescriptor("linux", "amd64", "a")
	published := platformDescriptor("linux", "amd64", "c")

	entries := map[string]ocispec.Descriptor{"linux/amd64": mine}
	mergeIndexEntries(entries, map[string]ocispec.Descriptor{"linux/amd64": published})
	if entries["linux/amd64"].Digest != mine.Digest {
		t.Error("the published child overwrote the one being pushed")
	}
}

// The same children must always encode to the same index digest, or a re-push of
// unchanged content looks like a new artifact.
func TestBuildIndexIsDeterministic(t *testing.T) {
	amd := platformDescriptor("linux", "amd64", "a")
	arm := platformDescriptor("linux", "arm64", "b")
	ann := map[string]string{AnnTitle: "ubuntu24/base"}

	first, err := json.Marshal(buildIndex(map[string]ocispec.Descriptor{
		"linux/amd64": amd, "linux/arm64": arm,
	}, ann, ArtifactTypeBase))
	if err != nil {
		t.Fatal(err)
	}
	second, err := json.Marshal(buildIndex(map[string]ocispec.Descriptor{
		"linux/arm64": arm, "linux/amd64": amd,
	}, ann, ArtifactTypeBase))
	if err != nil {
		t.Fatal(err)
	}
	if digest.FromBytes(first) != digest.FromBytes(second) {
		t.Errorf("insertion order changed the index digest:\n%s\n%s", first, second)
	}

	var idx ocispec.Index
	if err := json.Unmarshal(first, &idx); err != nil {
		t.Fatal(err)
	}
	if idx.MediaType != ocispec.MediaTypeImageIndex || idx.SchemaVersion != 2 {
		t.Errorf("index = %+v, want an OCI image index at schema 2", idx)
	}
	if idx.ArtifactType != ArtifactTypeBase {
		t.Errorf("artifact type = %q, want %q", idx.ArtifactType, ArtifactTypeBase)
	}
	if got := []string{
		platformKey(idx.Manifests[0].Platform), platformKey(idx.Manifests[1].Platform),
	}; got[0] != "linux/amd64" || got[1] != "linux/arm64" {
		t.Errorf("children ordered %v, want sorted by platform key", got)
	}
}

// The index spans architectures, so it may only carry facts that hold for all of
// them. The compressor is a property of one build's bytes; everything else is a
// property of the artifact.
func TestIndexAnnotationsDropOnlyTheCompressor(t *testing.T) {
	child := Annotations(fullManifest(), "zstd")
	index := indexAnnotations(child)

	if _, ok := index[AnnCompression]; ok {
		t.Error("the index claims a compressor its children may disagree on")
	}
	for key, value := range child {
		if key == AnnCompression {
			continue
		}
		if index[key] != value {
			t.Errorf("%s = %q at the index, %q at the child", key, index[key], value)
		}
	}
	if len(index) != len(child)-1 {
		t.Errorf("index annotations = %v, want the child's minus the compressor", index)
	}
}

// The child's annotations must not be edited in place: push reuses the same map
// for the per-arch manifest it already uploaded.
func TestIndexAnnotationsDoNotMutateTheInput(t *testing.T) {
	child := Annotations(fullManifest(), "zstd")
	indexAnnotations(child)
	if child[AnnCompression] != "zstd" {
		t.Error("indexAnnotations stripped the compressor from the child's own annotations")
	}
}

func TestIndexPlatforms(t *testing.T) {
	idx := ocispec.Index{Manifests: []ocispec.Descriptor{
		platformDescriptor("linux", "arm64", "b"),
		platformDescriptor("linux", "amd64", "a"),
		{Digest: digest.Digest(digestPrefix + strings.Repeat("d", 64))}, // no platform
	}}
	got := indexPlatforms(idx)
	if len(got) != 2 || got[0] != "linux/amd64" || got[1] != "linux/arm64" {
		t.Errorf("indexPlatforms = %v, want the two keyed platforms, sorted", got)
	}
}
