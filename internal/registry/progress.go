package registry

import (
	"context"
	"fmt"
	"io"
	"log/slog"
	"sync"
	"time"

	ocispec "github.com/opencontainers/image-spec/specs-go/v1"
	"oras.land/oras-go/v2"

	"github.com/Justype/condatainer/internal/logging"
	"github.com/Justype/condatainer/internal/utils"
)

// When a transfer is worth narrating. Below progressMinSize a blob is gone
// before a first report would help; above it, report every progressStep bytes or
// every progressInterval, whichever comes first — the timer is what distinguishes
// a slow link from a stalled one.
const (
	progressMinSize  = int64(4 << 20)
	progressStep     = int64(64 << 20)
	progressInterval = 2 * time.Second
)

// Transfer directions, used as the log message stem. verbPublish covers the
// small writes that finish a push — manifest, index, tags — which move no
// payload and so are never narrated by byte count, only by what happens to them.
const (
	verbUpload   = "upload"
	verbDownload = "download"
	verbPublish  = "publish"
)

// progressTarget narrates large blobs written into an ORAS target, in either
// direction: ORAS models a download as a push too, into a local file store.
//
// Reports go to the context logger, never to stdout — printing would bypass the
// dashboard's SSE log stream, leaving a pull started from the web UI silent for
// the length of a multi-gigabyte download.
type progressTarget struct {
	oras.Target
	verb string
}

// withTransferProgress wraps target so blobs past progressMinSize report as they
// move. verb is [verbUpload] or [verbDownload].
func withTransferProgress(target oras.Target, verb string) oras.Target {
	return &progressTarget{Target: target, verb: verb}
}

func (t *progressTarget) Push(ctx context.Context, desc ocispec.Descriptor, content io.Reader) error {
	if desc.Size < progressMinSize {
		return t.Target.Push(ctx, desc, content)
	}
	reader := newProgressReader(ctx, content, desc.Size, t.verb)
	err := t.Target.Push(ctx, desc, reader)
	if err == nil {
		// A target that stops at the declared size never reads the EOF that would
		// have reported completion.
		reader.finish()
	}
	return err
}

// downloadProgress narrates a whole artifact rather than one layer at a time.
//
// Layers arrive concurrently into their own ranges of one file, so there is no
// "current layer" to name: a completed count reads zero for the first several
// minutes while three 2 GiB layers are each half done. Push reports per layer
// for the opposite reason — one upload is in flight at a time, and a retry names
// it.
type downloadProgress struct {
	log   *slog.Logger
	total int64

	mu         sync.Mutex
	done       int64
	nextReport int64
	lastReport time.Time
	finished   bool
}

func newDownloadProgress(ctx context.Context, total int64) *downloadProgress {
	return &downloadProgress{log: logging.FromContext(ctx), total: total}
}

// attempt returns a reader that counts one transfer's bytes into the artifact's
// total. One per attempt: a layer that is retried rewrites its own range from
// the start, so what an abandoned attempt transferred is not progress.
func (d *downloadProgress) attempt(r io.Reader) io.Reader {
	return &countingReader{parent: d, reader: r}
}

// advance adds n bytes and reports on the usual cadence.
func (d *downloadProgress) advance(n int) {
	d.mu.Lock()
	d.done += int64(n)
	now := time.Now()
	if d.nextReport == 0 {
		d.nextReport, d.lastReport = progressStep, now
	}
	due := d.done >= d.nextReport || now.Sub(d.lastReport) >= progressInterval
	if due {
		d.nextReport, d.lastReport = d.done+progressStep, now
	}
	d.mu.Unlock()

	if due {
		d.report(false)
	}
}

// withdraw removes the bytes an abandoned attempt contributed, so a retried
// layer does not leave the artifact reading as more complete than it is.
func (d *downloadProgress) withdraw(n int64) {
	d.mu.Lock()
	d.done -= n
	if d.done < 0 {
		d.done = 0
	}
	d.mu.Unlock()
}

// interrupted reports how far a cancelled transfer got, and ends the progress
// line. Without it the last drawn line, with the terminal's "^C" on the end, is
// the whole account of what happened.
func (d *downloadProgress) interrupted() {
	d.mu.Lock()
	already, done := d.finished, d.done
	d.finished = true
	d.mu.Unlock()
	if already || d.log == nil {
		return
	}
	d.log.Warn(verbDownload+" cancelled",
		"done", utils.FormatSize(done), "total", utils.FormatSize(d.total))
}

// finish reports completion once, whatever the last layer's size was: the final
// record is what ends the interactive progress line.
func (d *downloadProgress) finish() {
	d.mu.Lock()
	already := d.finished
	d.finished = true
	d.mu.Unlock()
	if !already {
		d.report(true)
	}
}

func (d *downloadProgress) report(final bool) {
	if d.log == nil {
		return
	}
	d.mu.Lock()
	done := d.done
	d.mu.Unlock()
	d.log.Info(verbDownload+" progress", "kind", "progress",
		"done", utils.FormatSize(done), "total", utils.FormatSize(d.total),
		"final", final, "last", final)
}

// countingReader counts one attempt's bytes into the artifact's total, and gives
// them back if the attempt does not finish.
type countingReader struct {
	parent *downloadProgress
	reader io.Reader
	read   int64
}

func (r *countingReader) Read(p []byte) (int, error) {
	n, err := r.reader.Read(p)
	if n > 0 {
		r.read += int64(n)
		r.parent.advance(n)
	}
	if err != nil && err != io.EOF {
		r.parent.withdraw(r.read)
		r.read = 0
	}
	return n, err
}

// progressReader counts bytes on their way through and reports periodically. It
// is read by one goroutine — whichever is draining it — so it holds no lock;
// concurrent transfers stay separate because each has its own reader and slog
// handlers are safe to share.
type progressReader struct {
	reader     io.Reader
	log        *slog.Logger
	verb       string
	total      int64
	done       int64
	nextReport int64
	lastReport time.Time
	finished   bool
	layer      int
	layers     int
}

func newProgressReader(ctx context.Context, r io.Reader, total int64, verb string) *progressReader {
	return &progressReader{reader: r, log: logging.FromContext(ctx), verb: verb, total: total}
}

func newLayerProgressReader(ctx context.Context, r io.Reader, total int64, verb string, layer, layers int) *progressReader {
	reader := newProgressReader(ctx, r, total, verb)
	reader.layer, reader.layers = layer, layers
	return reader
}

func (r *progressReader) Read(p []byte) (int, error) {
	n, err := r.reader.Read(p)
	r.done += int64(n)

	now := time.Now()
	if r.nextReport == 0 {
		r.nextReport = progressStep
		r.lastReport = now
	}
	if r.done >= r.nextReport || now.Sub(r.lastReport) >= progressInterval {
		r.report(false)
		r.nextReport = r.done + progressStep
		r.lastReport = now
	}
	if err == io.EOF {
		r.finish()
	}
	return n, err
}

// finish reports completion once, however many times it is called: the reader
// hits EOF and the caller reports success independently, and only one of the two
// is guaranteed to happen.
func (r *progressReader) finish() {
	if r.finished {
		return
	}
	r.finished = true
	r.report(true)
}

func (r *progressReader) report(final bool) {
	if r.log == nil {
		return
	}
	attrs := []any{"kind", "progress", "done", utils.FormatSize(r.done), "total", utils.FormatSize(r.total)}
	if r.layers > 0 {
		attrs = append(attrs, "layer", fmt.Sprintf("%d/%d", r.layer, r.layers))
	}
	attrs = append(attrs, "final", final, "last", r.layers == 0 || r.layer == r.layers)
	r.log.Info(r.verb+" progress", attrs...)
}
