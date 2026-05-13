package logging

import (
	"context"
	"io"
	"log/slog"
)

type contextKey struct{}

var discardLogger = slog.New(slog.NewTextHandler(io.Discard, nil))

// WithLogger returns a context carrying logger. A nil logger leaves ctx unchanged.
func WithLogger(ctx context.Context, logger *slog.Logger) context.Context {
	if logger == nil {
		return ctx
	}
	return context.WithValue(ctx, contextKey{}, logger)
}

// FromContext returns the logger attached to ctx, or a discard logger.
func FromContext(ctx context.Context) *slog.Logger {
	if logger, ok := ctx.Value(contextKey{}).(*slog.Logger); ok && logger != nil {
		return logger
	}
	return discardLogger
}
