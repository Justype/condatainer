package registry

import (
	"bytes"
	"context"
	"fmt"
	"io"
	"os"
	"path/filepath"
	"testing"

	"github.com/opencontainers/go-digest"
	ocispec "github.com/opencontainers/image-spec/specs-go/v1"
)

// memoryBlobs stands in for a repository's blob store. Its Push verifies the
// descriptor against the bytes it actually receives, so every test here asserts
// for free that a chunk's declared digest and size match what was streamed.
type memoryBlobs struct {
	blobs  map[digest.Digest][]byte
	pushed []ocispec.Descriptor
}

func newMemoryBlobs() *memoryBlobs {
	return &memoryBlobs{blobs: map[digest.Digest][]byte{}}
}

func (m *memoryBlobs) Exists(_ context.Context, target ocispec.Descriptor) (bool, error) {
	_, ok := m.blobs[target.Digest]
	return ok, nil
}

func (m *memoryBlobs) Push(_ context.Context, expected ocispec.Descriptor, content io.Reader) error {
	data, err := io.ReadAll(content)
	if err != nil {
		return err
	}
	if int64(len(data)) != expected.Size {
		return fmt.Errorf("declared size %d, streamed %d bytes", expected.Size, len(data))
	}
	if got := digest.FromBytes(data); got != expected.Digest {
		return fmt.Errorf("declared digest %s, streamed %s", expected.Digest, got)
	}
	m.blobs[expected.Digest] = data
	m.pushed = append(m.pushed, expected)
	return nil
}

// withChunkSize runs the body with a chunk size small enough to test cheaply.
func withChunkSize(t *testing.T, size int64) {
	t.Helper()
	original := artifactChunkSize
	artifactChunkSize = size
	t.Cleanup(func() { artifactChunkSize = original })
}

func writeArtifact(t *testing.T, dir, name, content string) string {
	t.Helper()
	path := filepath.Join(dir, name)
	if err := os.WriteFile(path, []byte(content), 0o600); err != nil {
		t.Fatal(err)
	}
	return path
}

func titleOf(desc ocispec.Descriptor) string { return desc.Annotations[ocispec.AnnotationTitle] }

func TestPushArtifactLayersChunks(t *testing.T) {
	const content = "abcdefghij"
	path := writeArtifact(t, t.TempDir(), "sample.sqf", content)
	withChunkSize(t, 4)

	blobs := newMemoryBlobs()
	layers, err := pushArtifactLayers(context.Background(), blobs, path, MediaTypeOverlayBlob)
	if err != nil {
		t.Fatalf("pushArtifactLayers: %v", err)
	}
	if len(layers) != 3 {
		t.Fatalf("layer count = %d, want 3", len(layers))
	}

	// Offset order is reassembly order, and the names must say so.
	var assembled bytes.Buffer
	for i, want := range []int64{4, 4, 2} {
		if layers[i].Size != want {
			t.Errorf("layer %d size = %d, want %d", i, layers[i].Size, want)
		}
		if got, want := titleOf(layers[i]), fmt.Sprintf("sample.sqf.part%06d", i); got != want {
			t.Errorf("layer %d title = %q, want %q", i, got, want)
		}
		if layers[i].MediaType != MediaTypeOverlayBlob {
			t.Errorf("layer %d media type = %q", i, layers[i].MediaType)
		}
		assembled.Write(blobs.blobs[layers[i].Digest])
	}
	if assembled.String() != content {
		t.Errorf("layers concatenate to %q, want %q", assembled.String(), content)
	}
}

// A chunk's digest must be the digest of that byte range of the source, or a
// pull reassembles bytes the registry verified against the wrong claim.
func TestPushArtifactLayersDigestsTheRange(t *testing.T) {
	const content = "abcdefghij"
	path := writeArtifact(t, t.TempDir(), "sample.sqf", content)
	withChunkSize(t, 4)

	layers, err := pushArtifactLayers(context.Background(), newMemoryBlobs(), path, MediaTypeOverlayBlob)
	if err != nil {
		t.Fatal(err)
	}
	for i, want := range []string{"abcd", "efgh", "ij"} {
		if got := digest.FromString(want); layers[i].Digest != got {
			t.Errorf("layer %d digest = %s, want %s (the digest of %q)", i, layers[i].Digest, got, want)
		}
	}
}

// The whole reason for the section-reader push: a 60 GB artifact must not first
// write 60 GB of chunks to a scratch disk the user never named.
func TestPushArtifactLayersStagesNothing(t *testing.T) {
	// The artifact's own directory is claimed before TMPDIR is redirected, so
	// anything left in staging came from the push.
	path := writeArtifact(t, t.TempDir(), "sample.sqf", "abcdefghijklmnop")
	staging := t.TempDir()
	t.Setenv("TMPDIR", staging)
	withChunkSize(t, 4)

	if _, err := pushArtifactLayers(context.Background(), newMemoryBlobs(), path, MediaTypeOverlayBlob); err != nil {
		t.Fatal(err)
	}
	left, err := os.ReadDir(staging)
	if err != nil {
		t.Fatal(err)
	}
	if len(left) != 0 {
		t.Errorf("push left %d entries in TMPDIR; a chunk is a byte range, not a copy", len(left))
	}
}

// An artifact at or under the chunk size stays one layer named for the file, so
// a pull needs no reassembly at all.
func TestPushArtifactLayersSingleLayer(t *testing.T) {
	path := writeArtifact(t, t.TempDir(), "sample.sqf", "abcd")
	withChunkSize(t, 4)

	layers, err := pushArtifactLayers(context.Background(), newMemoryBlobs(), path, MediaTypeOverlayBlob)
	if err != nil {
		t.Fatal(err)
	}
	if len(layers) != 1 {
		t.Fatalf("layer count = %d, want 1", len(layers))
	}
	if got := titleOf(layers[0]); got != "sample.sqf" {
		t.Errorf("title = %q, want the filename with no part suffix", got)
	}
}

// A chunk the registry already holds is not sent again, which is what makes a
// retried push cheap instead of a second full upload.
func TestPushArtifactLayersSkipsBlobsAlreadyPresent(t *testing.T) {
	path := writeArtifact(t, t.TempDir(), "sample.sqf", "abcdefghij")
	withChunkSize(t, 4)
	ctx := context.Background()

	blobs := newMemoryBlobs()
	first, err := pushArtifactLayers(ctx, blobs, path, MediaTypeOverlayBlob)
	if err != nil {
		t.Fatal(err)
	}
	sent := len(blobs.pushed)

	second, err := pushArtifactLayers(ctx, blobs, path, MediaTypeOverlayBlob)
	if err != nil {
		t.Fatal(err)
	}
	if len(blobs.pushed) != sent {
		t.Errorf("re-push sent %d more blobs, want 0", len(blobs.pushed)-sent)
	}
	// Skipping an upload must still produce the descriptor the manifest needs.
	for i := range first {
		if first[i].Digest != second[i].Digest || titleOf(first[i]) != titleOf(second[i]) {
			t.Errorf("layer %d differs between pushes: %+v vs %+v", i, first[i], second[i])
		}
	}
}

func TestPushArtifactLayersEmptyFile(t *testing.T) {
	path := writeArtifact(t, t.TempDir(), "empty.sqf", "")
	withChunkSize(t, 4)

	layers, err := pushArtifactLayers(context.Background(), newMemoryBlobs(), path, MediaTypeOverlayBlob)
	if err != nil {
		t.Fatalf("pushArtifactLayers: %v", err)
	}
	if len(layers) != 1 || layers[0].Size != 0 {
		t.Fatalf("layers = %+v, want one empty layer", layers)
	}
}

func TestAssemblePulledArtifact(t *testing.T) {
	dir := t.TempDir()
	for name, data := range map[string]string{
		"sample.sqf.part000000": "abc",
		"sample.sqf.part000001": "def",
		"sample.sqf.part000002": "ghi",
	} {
		writeArtifact(t, dir, name, data)
	}

	path, err := assemblePulledArtifact(dir, []string{
		"sample.sqf.part000000",
		"sample.sqf.part000001",
		"sample.sqf.part000002",
	})
	if err != nil {
		t.Fatalf("assemblePulledArtifact: %v", err)
	}
	got, err := os.ReadFile(path)
	if err != nil {
		t.Fatal(err)
	}
	if string(got) != "abcdefghi" {
		t.Fatalf("assembled payload = %q", got)
	}
}

// A single layer is the artifact already; there is nothing to concatenate and no
// second copy to make.
func TestAssemblePulledArtifactSingleLayerIsReturnedAsIs(t *testing.T) {
	dir := t.TempDir()
	writeArtifact(t, dir, "sample.sqf", "payload")

	path, err := assemblePulledArtifact(dir, []string{"sample.sqf"})
	if err != nil {
		t.Fatal(err)
	}
	if path != filepath.Join(dir, "sample.sqf") {
		t.Errorf("path = %q, want the downloaded layer itself", path)
	}
}

// A gap in the sequence would otherwise assemble into a file that looks whole.
func TestAssemblePulledArtifactRejectsAGapInTheSequence(t *testing.T) {
	dir := t.TempDir()
	for _, name := range []string{"sample.part000000", "sample.part000002"} {
		writeArtifact(t, dir, name, name)
	}
	if _, err := assemblePulledArtifact(dir, []string{"sample.part000000", "sample.part000002"}); err == nil {
		t.Fatal("a chunk sequence missing part000001 was assembled anyway")
	}
}

// The manifest's layer count is the contract: a directory holding more or fewer
// files than it declares is not the artifact that was published.
func TestAssemblePulledArtifactRejectsAFileCountMismatch(t *testing.T) {
	dir := t.TempDir()
	writeArtifact(t, dir, "sample.sqf.part000000", "abc")
	writeArtifact(t, dir, "sample.sqf.part000001", "def")
	writeArtifact(t, dir, "stowaway", "x")

	if _, err := assemblePulledArtifact(dir, []string{"sample.sqf.part000000", "sample.sqf.part000001"}); err == nil {
		t.Fatal("an extra file in the staging directory was ignored")
	}
}

// Push and pull must agree on the naming, or a chunked artifact pushed by this
// build cannot be reassembled by it.
func TestChunkNamesRoundTrip(t *testing.T) {
	const content = "abcdefghij"
	src := t.TempDir()
	path := writeArtifact(t, src, "sample.sqf", content)
	withChunkSize(t, 4)

	blobs := newMemoryBlobs()
	layers, err := pushArtifactLayers(context.Background(), blobs, path, MediaTypeOverlayBlob)
	if err != nil {
		t.Fatal(err)
	}

	// Stand in for the pull's file store, which writes each layer under its title.
	staging := t.TempDir()
	names := make([]string, 0, len(layers))
	for _, layer := range layers {
		name := titleOf(layer)
		writeArtifact(t, staging, name, string(blobs.blobs[layer.Digest]))
		names = append(names, name)
	}

	assembled, err := assemblePulledArtifact(staging, names)
	if err != nil {
		t.Fatalf("assemblePulledArtifact: %v", err)
	}
	got, err := os.ReadFile(assembled)
	if err != nil {
		t.Fatal(err)
	}
	if string(got) != content {
		t.Errorf("round trip produced %q, want %q", got, content)
	}
}
