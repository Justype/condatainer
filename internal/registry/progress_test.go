package registry

import (
	"bytes"
	"context"
	"io"
	"log/slog"
	"strings"
	"testing"

	ocispec "github.com/opencontainers/image-spec/specs-go/v1"
	"oras.land/oras-go/v2"

	"github.com/Justype/condatainer/internal/logging"
)

// captureLogs returns a context carrying a logger and the buffer it writes to.
func captureLogs(t *testing.T) (context.Context, *bytes.Buffer) {
	t.Helper()
	var buf bytes.Buffer
	logger := slog.New(slog.NewTextHandler(&buf, &slog.HandlerOptions{Level: slog.LevelDebug}))
	return logging.WithLogger(context.Background(), logger), &buf
}

func TestProgressReaderPassesBytesThroughAndReportsOnce(t *testing.T) {
	ctx, logs := captureLogs(t)
	payload := bytes.Repeat([]byte("x"), 1024)

	reader := newProgressReader(ctx, bytes.NewReader(payload), int64(len(payload)), verbUpload)
	got, err := io.ReadAll(reader)
	if err != nil {
		t.Fatal(err)
	}
	if !bytes.Equal(got, payload) {
		t.Fatal("the progress reader changed the payload")
	}

	// The caller reports success independently of the reader hitting EOF; only
	// one of the two is guaranteed to happen, and both must not log.
	reader.finish()
	if n := strings.Count(logs.String(), verbUpload+" progress"); n != 1 {
		t.Errorf("logged completion %d times, want 1:\n%s", n, logs)
	}
}

// A transfer that never reaches EOF still reports what it moved, so a caller can
// see how far a failed push got.
func TestProgressReaderReportsWhatItSaw(t *testing.T) {
	ctx, logs := captureLogs(t)
	reader := newProgressReader(ctx, bytes.NewReader(bytes.Repeat([]byte("x"), 100)), 1000, verbDownload)

	if _, err := io.CopyN(io.Discard, reader, 100); err != nil {
		t.Fatal(err)
	}
	reader.finish()
	if !strings.Contains(logs.String(), verbDownload+" progress") {
		t.Errorf("no completion report:\n%s", logs)
	}
}

// discardTarget accepts any push and drops the bytes. Only Push is implemented;
// nothing here calls the rest of oras.Target.
type discardTarget struct{ oras.Target }

func (discardTarget) Push(_ context.Context, _ ocispec.Descriptor, content io.Reader) error {
	_, err := io.Copy(io.Discard, content)
	return err
}

// A small blob is gone before a report would help, so narrating it is noise.
func TestWithTransferProgressStaysQuietForSmallBlobs(t *testing.T) {
	ctx, logs := captureLogs(t)
	target := withTransferProgress(discardTarget{}, verbUpload)

	desc := ocispec.Descriptor{Size: progressMinSize - 1}
	if err := target.Push(ctx, desc, strings.NewReader("small")); err != nil {
		t.Fatal(err)
	}
	if logs.Len() != 0 {
		t.Errorf("a blob under the threshold was narrated:\n%s", logs)
	}
}

func TestWithTransferProgressNarratesLargeBlobs(t *testing.T) {
	ctx, logs := captureLogs(t)
	target := withTransferProgress(discardTarget{}, verbDownload)

	desc := ocispec.Descriptor{Size: progressMinSize}
	body := io.LimitReader(zeroes{}, progressMinSize)
	if err := target.Push(ctx, desc, body); err != nil {
		t.Fatal(err)
	}
	if !strings.Contains(logs.String(), verbDownload+" progress") {
		t.Errorf("a blob at the threshold was not narrated:\n%s", logs)
	}
}

// zeroes is an endless source, used with io.LimitReader so a multi-megabyte
// transfer costs no allocation.
type zeroes struct{}

func (zeroes) Read(p []byte) (int, error) { return len(p), nil }

func TestLayerProgressReaderReportsOverallLayer(t *testing.T) {
	ctx, logs := captureLogs(t)
	reader := newLayerProgressReader(ctx, strings.NewReader("chunk"), 5, verbUpload, 2, 3)
	if _, err := io.Copy(io.Discard, reader); err != nil {
		t.Fatal(err)
	}

	got := logs.String()
	if !strings.Contains(got, "layer=2/3") {
		t.Errorf("progress does not report the overall layer:\n%s", got)
	}
	if !strings.Contains(got, "final=true") || !strings.Contains(got, "last=false") {
		t.Errorf("intermediate-layer completion is not distinguishable:\n%s", got)
	}
}
