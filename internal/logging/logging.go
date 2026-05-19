package logging

import (
	"context"
	"io"
	"log/slog"
)

type loggerKey struct{}
type writerKey struct{}

var discardLogger = slog.New(slog.NewTextHandler(io.Discard, nil))

// WithLogger returns a context carrying logger. A nil logger leaves ctx unchanged.
func WithLogger(ctx context.Context, logger *slog.Logger) context.Context {
	if logger == nil {
		return ctx
	}
	return context.WithValue(ctx, loggerKey{}, logger)
}

// FromContext returns the logger attached to ctx, or a discard logger.
func FromContext(ctx context.Context) *slog.Logger {
	if logger, ok := ctx.Value(loggerKey{}).(*slog.Logger); ok && logger != nil {
		return logger
	}
	return discardLogger
}

// WithWriter attaches w to ctx for raw subprocess output streaming.
// Internal functions that run subprocesses (e.g. resize2fs, e2fsck) write
// their stdout/stderr to this writer when present, enabling web terminal streaming.
func WithWriter(ctx context.Context, w io.Writer) context.Context {
	return context.WithValue(ctx, writerKey{}, w)
}

// WriterFromCtx returns the io.Writer attached to ctx, or nil if none was set.
func WriterFromCtx(ctx context.Context) io.Writer {
	w, _ := ctx.Value(writerKey{}).(io.Writer)
	return w
}
