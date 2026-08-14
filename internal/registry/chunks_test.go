package registry

import (
	"bytes"
	"context"
	"errors"
	"fmt"
	"io"
	"os"
	"path/filepath"
	"strings"
	"testing"

	"github.com/opencontainers/go-digest"
	ocispec "github.com/opencontainers/image-spec/specs-go/v1"
	"oras.land/oras-go/v2/errdef"
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

	blobs := newMemoryBlobs()
	layers, err := pushArtifactLayers(context.Background(), blobs, path, MediaTypeOverlayBlob, 4)
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

	layers, err := pushArtifactLayers(context.Background(), newMemoryBlobs(), path, MediaTypeOverlayBlob, 4)
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

	if _, err := pushArtifactLayers(context.Background(), newMemoryBlobs(), path, MediaTypeOverlayBlob, 4); err != nil {
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

	layers, err := pushArtifactLayers(context.Background(), newMemoryBlobs(), path, MediaTypeOverlayBlob, 4)
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
	ctx := context.Background()

	blobs := newMemoryBlobs()
	first, err := pushArtifactLayers(ctx, blobs, path, MediaTypeOverlayBlob, 4)
	if err != nil {
		t.Fatal(err)
	}
	sent := len(blobs.pushed)

	second, err := pushArtifactLayers(ctx, blobs, path, MediaTypeOverlayBlob, 4)
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

	layers, err := pushArtifactLayers(context.Background(), newMemoryBlobs(), path, MediaTypeOverlayBlob, 4)
	if err != nil {
		t.Fatalf("pushArtifactLayers: %v", err)
	}
	if len(layers) != 1 || layers[0].Size != 0 {
		t.Fatalf("layers = %+v, want one empty layer", layers)
	}
}

// flakyBlobs refuses the first refusals pushes with a rate limit, then behaves.
// It wraps memoryBlobs, whose Push verifies the streamed bytes against the
// descriptor, so a retry that resent the wrong bytes fails there rather than
// here.
type flakyBlobs struct {
	*memoryBlobs
	refusals int
	attempts int
	exists   int
}

func (f *flakyBlobs) Exists(ctx context.Context, target ocispec.Descriptor) (bool, error) {
	f.exists++
	return f.memoryBlobs.Exists(ctx, target)
}

func (f *flakyBlobs) Push(ctx context.Context, expected ocispec.Descriptor, content io.Reader) error {
	f.attempts++
	if f.attempts <= f.refusals {
		// Read the body first, the way a registry that refuses at finalization
		// does: the bytes are gone by the time the error arrives, so the next
		// attempt has to produce them again from somewhere.
		if _, err := io.Copy(io.Discard, content); err != nil {
			return err
		}
		return codedErr(403, "DENIED", ghcrSecondaryLimitMessage)
	}
	return f.memoryBlobs.Push(ctx, expected, content)
}

// The reason retry lives in pushRange rather than in the HTTP client: an
// io.SectionReader has no GetBody, so every attempt must open a new one. An
// attempt that reused the reader would stream nothing and the digest check would
// catch it.
func TestPushArtifactLayersRetriesWithAFreshReader(t *testing.T) {
	const content = "abcdefghij"
	path := writeArtifact(t, t.TempDir(), "sample.sqf", content)
	blobs := &flakyBlobs{memoryBlobs: newMemoryBlobs(), refusals: 2}

	ctx := withRetryPolicy(context.Background(), newTestPolicy().retryPolicy)
	layers, err := pushArtifactLayers(ctx, blobs, path, MediaTypeOverlayBlob, 4)
	if err != nil {
		t.Fatalf("pushArtifactLayers: %v", err)
	}

	// Three layers, the first refused twice: 2 wasted attempts plus 3 real ones.
	if blobs.attempts != 5 {
		t.Errorf("push attempts = %d, want 5", blobs.attempts)
	}
	var assembled bytes.Buffer
	for i, want := range []string{"abcd", "efgh", "ij"} {
		if got := digest.FromString(want); layers[i].Digest != got {
			t.Errorf("layer %d digest = %s, want the digest of %q", i, layers[i].Digest, want)
		}
		assembled.Write(blobs.blobs[layers[i].Digest])
	}
	if assembled.String() != content {
		t.Errorf("a retried push assembled to %q, want %q", assembled.String(), content)
	}
}

// A pause at chunk 34 resumes at chunk 34: the layers already committed stay in
// the descriptor slice and are never revisited, which is what makes an
// in-process retry cost one request rather than one per completed layer.
func TestRetryDoesNotRescanCommittedLayers(t *testing.T) {
	path := writeArtifact(t, t.TempDir(), "sample.sqf", "abcdefghij")
	blobs := &flakyBlobs{memoryBlobs: newMemoryBlobs(), refusals: 1}

	ctx := withRetryPolicy(context.Background(), newTestPolicy().retryPolicy)
	if _, err := pushArtifactLayers(ctx, blobs, path, MediaTypeOverlayBlob, 4); err != nil {
		t.Fatal(err)
	}
	// One probe on the first layer, plus the one recheck after the refusal.
	// Anything more means a retry went back over work already done.
	if want := 1 + 1; blobs.exists != want {
		t.Errorf("Exists calls = %d, want %d", blobs.exists, want)
	}
}

// A fresh push asks once. Layers upload in order from the first, so a first
// layer the registry does not have means none of this artifact is committed, and
// every further probe is a request spent to be told 404.
func TestFreshPushProbesOnce(t *testing.T) {
	path := writeArtifact(t, t.TempDir(), "sample.sqf", "abcdefghijklmnop")
	blobs := &flakyBlobs{memoryBlobs: newMemoryBlobs()}

	layers, err := pushArtifactLayers(context.Background(), blobs, path, MediaTypeOverlayBlob, 4)
	if err != nil {
		t.Fatal(err)
	}
	if len(layers) != 4 {
		t.Fatalf("layers = %d, want 4", len(layers))
	}
	if blobs.exists != 1 {
		t.Errorf("Exists calls = %d, want exactly one", blobs.exists)
	}
	if blobs.attempts != 4 {
		t.Errorf("uploads = %d, want every layer sent", blobs.attempts)
	}
}

// The case that pays for probing: a push resumed after a failure. The first
// layer is present, so every layer is checked and nothing is re-uploaded.
func TestResumedPushProbesEveryLayer(t *testing.T) {
	path := writeArtifact(t, t.TempDir(), "sample.sqf", "abcdefghijklmnop")
	blobs := &flakyBlobs{memoryBlobs: newMemoryBlobs()}
	ctx := context.Background()

	if _, err := pushArtifactLayers(ctx, blobs, path, MediaTypeOverlayBlob, 4); err != nil {
		t.Fatal(err)
	}
	sent, probed := blobs.attempts, blobs.exists

	if _, err := pushArtifactLayers(ctx, blobs, path, MediaTypeOverlayBlob, 4); err != nil {
		t.Fatal(err)
	}
	if blobs.attempts != sent {
		t.Errorf("a resumed push re-uploaded %d layers", blobs.attempts-sent)
	}
	if got := blobs.exists - probed; got != 4 {
		t.Errorf("probes on the resumed push = %d, want one per layer", got)
	}
}

// Skipping a probe is never a correctness risk: a registry accepts a blob it
// already holds and deduplicates by digest. The cost of guessing wrong is
// bandwidth, in a case the heuristic has established is unlikely.
func TestBlindUploadOfAPresentBlobSucceeds(t *testing.T) {
	path := writeArtifact(t, t.TempDir(), "sample.sqf", "abcdefghijklmnop")
	blobs := &alreadyHasEverything{memoryBlobs: newMemoryBlobs()}

	layers, err := pushArtifactLayers(context.Background(), blobs, path, MediaTypeOverlayBlob, 4)
	if err != nil {
		t.Fatalf("a blind upload of an existing blob failed: %v", err)
	}
	if len(layers) != 4 {
		t.Errorf("layers = %d, want 4", len(layers))
	}
}

// alreadyHasEverything reports nothing present — so probing is elided — and then
// rejects every upload the way a registry does when it already holds the blob.
type alreadyHasEverything struct{ *memoryBlobs }

func (a *alreadyHasEverything) Exists(context.Context, ocispec.Descriptor) (bool, error) {
	return false, nil
}

func (a *alreadyHasEverything) Push(_ context.Context, _ ocispec.Descriptor, content io.Reader) error {
	if _, err := io.Copy(io.Discard, content); err != nil {
		return err
	}
	return errdef.ErrAlreadyExists
}

// A registry can commit the blob and lose the response — GHCR's refusal arrives
// from an intermediary in front of storage that already accepted it.
func TestPushArtifactLayersStopsWhenTheBlobLandedAnyway(t *testing.T) {
	path := writeArtifact(t, t.TempDir(), "sample.sqf", "abcd")
	blobs := &committingBlobs{memoryBlobs: newMemoryBlobs()}

	ctx := withRetryPolicy(context.Background(), newTestPolicy().retryPolicy)
	layers, err := pushArtifactLayers(ctx, blobs, path, MediaTypeOverlayBlob, 4)
	if err != nil {
		t.Fatalf("pushArtifactLayers: %v", err)
	}
	if blobs.attempts != 1 {
		t.Errorf("push attempts = %d, want the committed blob to end it after one", blobs.attempts)
	}
	if len(layers) != 1 || layers[0].Digest != digest.FromString("abcd") {
		t.Errorf("layers = %+v, want the descriptor produced anyway", layers)
	}
}

// committingBlobs stores the blob and then reports a rate limit, which is the
// lost-response case.
type committingBlobs struct {
	*memoryBlobs
	attempts int
}

func (c *committingBlobs) Push(ctx context.Context, expected ocispec.Descriptor, content io.Reader) error {
	c.attempts++
	if err := c.memoryBlobs.Push(ctx, expected, content); err != nil {
		return err
	}
	return codedErr(403, "DENIED", ghcrSecondaryLimitMessage)
}

// Push says it too, and for the same reason. It has nothing local to clean up —
// a chunk is a range of a file that already exists — but the layers the registry
// already accepted are what make running the command again cheap, so the count
// is worth stating.
func TestPushReportsCancellation(t *testing.T) {
	ctx, logs := captureLogs(t)
	ctx, cancel := context.WithCancel(ctx)
	path := writeArtifact(t, t.TempDir(), "sample.sqf", "abcdefghijklmnop")

	blobs := &cancellingBlobs{memoryBlobs: newMemoryBlobs(), cancel: cancel, after: 2}
	_, err := pushArtifactLayers(ctx, blobs, path, MediaTypeOverlayBlob, 4)
	if !errors.Is(err, context.Canceled) {
		t.Fatalf("err = %v, want the cancellation", err)
	}

	out := logs.String()
	if !strings.Contains(out, verbUpload+" cancelled") {
		t.Errorf("a cancelled push said nothing:\n%s", out)
	}
	if !strings.Contains(out, "accepted=2/4") {
		t.Errorf("the cancellation does not say what the registry kept:\n%s", out)
	}
}

// cancellingBlobs interrupts the command once the registry has accepted after
// layers, the way Ctrl-C lands partway through a transfer.
type cancellingBlobs struct {
	*memoryBlobs
	cancel context.CancelFunc
	after  int
	pushed int
}

func (c *cancellingBlobs) Push(ctx context.Context, expected ocispec.Descriptor, content io.Reader) error {
	if c.pushed >= c.after {
		c.cancel()
		return context.Canceled
	}
	c.pushed++
	return c.memoryBlobs.Push(ctx, expected, content)
}
