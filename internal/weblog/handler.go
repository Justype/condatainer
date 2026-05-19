// Package weblog provides a slog.Handler that routes structured log records
// to an io.Writer (typically the SSE broker writer), so internal packages using
// logging.FromContext(ctx) produce visible output in the web terminal.
//
// Level → output:
//
//	Debug → "[CNT…] msg"  (only if utils.DebugMode)
//	Info  → "[CNT] msg" or "[CNT✓] msg" / "[CNT▶] msg" / "[CNT◇] msg" via kind attr
//	        (suppressed when utils.QuietMode)
//	Warn  → "[CNT!] msg"
//	Error → "[CNT✗] msg"
package weblog

import (
	"context"
	"fmt"
	"io"
	"log/slog"
	"strings"

	"github.com/Justype/condatainer/internal/utils"
)

// Handler is a slog.Handler that writes plain-text records to an io.Writer.
type Handler struct {
	w     io.Writer
	attrs []slog.Attr
}

// New returns a Handler that writes to w.
func New(w io.Writer) *Handler { return &Handler{w: w} }

// Enabled gates Debug on utils.DebugMode; all other levels always pass.
func (h *Handler) Enabled(_ context.Context, level slog.Level) bool {
	if level == slog.LevelDebug {
		return utils.DebugMode
	}
	return true
}

// Handle formats the record and writes it to h.w.
func (h *Handler) Handle(_ context.Context, r slog.Record) error {
	var kind string
	var extras []string
	visit := func(a slog.Attr) {
		if a.Key == "kind" {
			kind = a.Value.String()
			return
		}
		extras = append(extras, fmt.Sprintf("%s=%v", a.Key, a.Value.Any()))
	}
	for _, a := range h.attrs {
		visit(a)
	}
	r.Attrs(func(a slog.Attr) bool {
		visit(a)
		return true
	})

	msg := r.Message
	if len(extras) > 0 {
		msg = msg + " " + strings.Join(extras, " ")
	}

	var line string
	switch r.Level {
	case slog.LevelDebug:
		line = "[CNT…] " + msg
	case slog.LevelInfo:
		if utils.QuietMode {
			return nil
		}
		switch kind {
		case "success":
			line = "[CNT✓] " + msg
		case "hint":
			line = "[CNT▶] " + msg
		case "note":
			line = "[CNT◇] " + msg
		default:
			line = "[CNT] " + msg
		}
	case slog.LevelWarn:
		line = "[CNT!] " + msg
	case slog.LevelError:
		line = "[CNT✗] " + msg
	default:
		line = "[CNT] " + msg
	}

	fmt.Fprintln(h.w, line)
	return nil
}

// WithAttrs returns a new handler whose records inherit the given attributes.
func (h *Handler) WithAttrs(attrs []slog.Attr) slog.Handler {
	nh := *h
	nh.attrs = append(append([]slog.Attr(nil), h.attrs...), attrs...)
	return &nh
}

// WithGroup is a no-op — group names would clutter the terminal output.
func (h *Handler) WithGroup(_ string) slog.Handler { return h }
