package registry

import (
	"context"
	"fmt"
	"io"
	"log/slog"
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

// Transfer directions, used as the log message stem.
const (
	verbUpload   = "upload"
	verbDownload = "download"
)

// progressTarget narrates large blobs written into an ORAS target. It covers
// both directions, because ORAS models a download as a push too: a remote
// repository is the target of an upload, a local file store the target of a
// download.
//
// Reports go to the context logger, never to stdout. A transport package that
// prints bypasses the dashboard's SSE log stream, so a pull started from the web
// UI would sit silent for the length of a multi-gigabyte download.
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
	attrs := []any{"kind", "progress", "done", utils.FormatBytes(r.done), "total", utils.FormatBytes(r.total)}
	if r.layers > 0 {
		attrs = append(attrs, "layer", fmt.Sprintf("%d/%d", r.layer, r.layers))
	}
	attrs = append(attrs, "final", final, "last", r.layers == 0 || r.layer == r.layers)
	r.log.Info(r.verb+" progress", attrs...)
}
