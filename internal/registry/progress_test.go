package registry

import (
	"bytes"
	"context"
	"io"
	"log/slog"
	"strings"
	"testing"
	"time"

	ocispec "github.com/opencontainers/image-spec/specs-go/v1"
	"oras.land/oras-go/v2"

	"github.com/Justype/condatainer/internal/logging"
	"github.com/Justype/condatainer/internal/utils"
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

// The 20 GiB fixture ends in a 4 KiB chunk: 21474840576 bytes in 2 GiB layers is
// ten full ones and a remainder far below the reporting threshold. Without the
// last layer reporting anyway, the display finished at "layer=10/11" and the next
// message printed onto that same open line.
func TestLastLayerClosesTheProgressLineHoweverSmall(t *testing.T) {
	ctx, logs := captureLogs(t)
	dir := t.TempDir()
	// Two layers: one worth narrating, then a tiny remainder.
	body := strings.Repeat("x", int(progressMinSize)+16)
	path := writeArtifact(t, dir, "sample.sqf", body)

	if _, err := pushArtifactLayers(ctx, newMemoryBlobs(), path, MediaTypeOverlayBlob, progressMinSize); err != nil {
		t.Fatal(err)
	}

	out := logs.String()
	if !strings.Contains(out, "layer=2/2") {
		t.Errorf("the final layer never reported:\n%s", out)
	}

	// last=true is what tells the terminal handler to end the line. Exactly one
	// record may carry it, and it has to be the final layer's — a line closed by
	// layer 1 leaves layer 2 drawing onto a line that is already finished.
	var closers []string
	for _, line := range strings.Split(strings.TrimSpace(out), "\n") {
		if strings.Contains(line, "last=true") {
			closers = append(closers, line)
		}
	}
	if len(closers) != 1 {
		t.Fatalf("%d records claimed to close the line, want 1:\n%s", len(closers), out)
	}
	if !strings.Contains(closers[0], "layer=2/2") || !strings.Contains(closers[0], "final=true") {
		t.Errorf("the line was closed by the wrong record: %s", closers[0])
	}
}

// The counter reports the artifact, not a layer: with layers arriving at once
// there is no current one to name, and a completed count read 0/11 for the first
// several minutes while three 2 GiB layers were each half done.
func TestDownloadProgressCountsTheArtifact(t *testing.T) {
	ctx, logs := captureLogs(t)
	p := newDownloadProgress(ctx, 4*progressStep)

	for range 2 {
		if _, err := io.Copy(io.Discard, p.attempt(io.LimitReader(zeroes{}, progressStep))); err != nil {
			t.Fatal(err)
		}
	}
	p.finish()

	out := logs.String()
	if want := `total="` + utils.FormatSize(4*progressStep) + `"`; !strings.Contains(out, want) {
		t.Errorf("no line reports the artifact total (%s):\n%s", want, out)
	}
	if strings.Contains(out, "layer") {
		t.Errorf("the download still names a layer:\n%s", out)
	}
	if got := strings.Count(out, "final=true"); got != 1 {
		t.Errorf("%d lines claimed completion, want 1:\n%s", got, out)
	}
}

// A retried layer rewrites its own range from the start, so what an abandoned
// attempt transferred is not progress — leaving it counted would report an
// artifact as more complete than it is, and could pass its own total.
func TestDownloadProgressWithdrawsAnAbandonedAttempt(t *testing.T) {
	ctx, logs := captureLogs(t)
	p := newDownloadProgress(ctx, 2*progressStep)

	failing := io.MultiReader(io.LimitReader(zeroes{}, progressStep), errReader{})
	if _, err := io.Copy(io.Discard, p.attempt(failing)); err == nil {
		t.Fatal("the failing attempt reported success")
	}
	if p.done != 0 {
		t.Errorf("done = %d after an abandoned attempt, want 0", p.done)
	}

	// The retry, which succeeds, must land the artifact on exactly its total.
	if _, err := io.Copy(io.Discard, p.attempt(io.LimitReader(zeroes{}, 2*progressStep))); err != nil {
		t.Fatal(err)
	}
	p.finish()
	if want := `done="` + utils.FormatSize(2*progressStep) + `"`; !strings.Contains(logs.String(), want) {
		t.Errorf("the retried download did not finish at its total (%s):\n%s", want, logs)
	}
}

type errReader struct{}

func (errReader) Read([]byte) (int, error) { return 0, io.ErrUnexpectedEOF }

// finish reports once however many times it is called: the caller reports
// completion and nothing else is guaranteed to.
func TestDownloadProgressFinishesOnce(t *testing.T) {
	ctx, logs := captureLogs(t)
	p := newDownloadProgress(ctx, progressStep)
	p.finish()
	p.finish()
	if got := strings.Count(logs.String(), "final=true"); got != 1 {
		t.Errorf("finish reported %d times, want 1", got)
	}
}

// An interrupt otherwise leaves the last drawn progress line, with the
// terminal's "^C" on the end of it, as the whole account of what happened.
func TestDownloadProgressReportsCancellation(t *testing.T) {
	ctx, logs := captureLogs(t)
	p := newDownloadProgress(ctx, 4*progressStep)
	if _, err := io.Copy(io.Discard, p.attempt(io.LimitReader(zeroes{}, progressStep))); err != nil {
		t.Fatal(err)
	}
	p.interrupted()

	out := logs.String()
	if !strings.Contains(out, verbDownload+" cancelled") {
		t.Errorf("a cancelled download said nothing:\n%s", out)
	}
	// How far it got is the useful part: it decides whether running it again is
	// worth it, and the figure is gone from the screen the moment it scrolls.
	if want := `done="` + utils.FormatSize(progressStep) + `"`; !strings.Contains(out, want) {
		t.Errorf("the cancellation does not say how far it got (%s):\n%s", want, out)
	}
	if !strings.Contains(out, `total="`+utils.FormatSize(4*progressStep)+`"`) {
		t.Errorf("the cancellation does not say how far there was to go:\n%s", out)
	}
}

// A cancelled transfer never completed, so it must not also report completion.
func TestDownloadProgressCancellationSuppressesCompletion(t *testing.T) {
	ctx, logs := captureLogs(t)
	p := newDownloadProgress(ctx, progressStep)
	p.interrupted()
	p.finish()

	if strings.Contains(logs.String(), "final=true") {
		t.Errorf("a cancelled download also claimed to finish:\n%s", logs)
	}
}

// A running report quotes the speed since the last report, and the final one the
// average over the whole transfer.
func TestDownloadProgressReportsSpeed(t *testing.T) {
	ctx, logs := captureLogs(t)
	p := newDownloadProgress(ctx, 200<<20)

	p.advance(1) // starts the clock
	p.start, p.sampleAt = time.Now().Add(-2*time.Second), time.Now().Add(-time.Second)
	p.advance(100 << 20)
	if !strings.Contains(logs.String(), "speed=") || !strings.Contains(logs.String(), "/s") {
		t.Fatalf("no speed on a running report:\n%s", logs)
	}

	logs.Reset()
	p.finish()
	if !strings.Contains(logs.String(), "speed=") || !strings.Contains(logs.String(), "final=true") {
		t.Fatalf("no speed on the final report:\n%s", logs)
	}
}

// An upload's blob line carries a speed too: the running one from the last
// report, the average on the last.
func TestProgressReaderReportsSpeed(t *testing.T) {
	ctx, logs := captureLogs(t)
	reader := newProgressReader(ctx, bytes.NewReader(bytes.Repeat([]byte("x"), 100<<20)), 100<<20, verbUpload)

	if _, err := io.CopyN(io.Discard, reader, 1); err != nil { // starts the clock
		t.Fatal(err)
	}
	reader.start, reader.sampleAt = time.Now().Add(-2*time.Second), time.Now().Add(-time.Second)
	if _, err := io.Copy(io.Discard, reader); err != nil {
		t.Fatal(err)
	}
	if got := logs.String(); strings.Count(got, "speed=") < 1 || !strings.Contains(got, "final=true") {
		t.Fatalf("no speed on the upload report:\n%s", got)
	}
}
