// Package clilog provides a slog.Handler that routes structured log records
// through the CLI's utils.Print* functions, so any internal package using
// logging.FromContext(ctx) produces the familiar [CNT][TAG] terminal output.
//
// Level → printer:
//
//	Debug → utils.PrintDebug
//	Info  → utils.PrintMessage (or PrintSuccess/PrintHint/PrintNote when
//	        an attribute kind=success|hint|note is attached)
//	Warn  → utils.PrintWarning
//	Error → utils.PrintError
package clilog

import (
	"context"
	"fmt"
	"log/slog"
	"strings"

	"github.com/Justype/condatainer/internal/utils"
)

// Handler is a slog.Handler that dispatches records to utils.Print*.
type Handler struct {
	attrs []slog.Attr
}

// New returns a fresh Handler.
func New() *Handler { return &Handler{} }

// Enabled gates Debug records on utils.DebugMode; all other levels pass through
// (Info/Warn already respect QuietMode inside the Print* functions themselves).
func (h *Handler) Enabled(_ context.Context, level slog.Level) bool {
	if level == slog.LevelDebug {
		return utils.DebugMode
	}
	return true
}

// Handle formats the record and dispatches to the matching utils.Print* call.
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

	switch r.Level {
	case slog.LevelDebug:
		utils.PrintDebug("%s", msg)
	case slog.LevelInfo:
		switch kind {
		case "success":
			utils.PrintSuccess("%s", msg)
		case "hint":
			utils.PrintHint("%s", msg)
		case "note":
			utils.PrintNote("%s", msg)
		default:
			utils.PrintMessage("%s", msg)
		}
	case slog.LevelWarn:
		utils.PrintWarning("%s", msg)
	case slog.LevelError:
		utils.PrintError("%s", msg)
	}
	return nil
}

// WithAttrs returns a new handler whose records inherit the given attributes.
func (h *Handler) WithAttrs(attrs []slog.Attr) slog.Handler {
	nh := *h
	nh.attrs = append(append([]slog.Attr(nil), h.attrs...), attrs...)
	return &nh
}

// WithGroup is a no-op for this handler — group names would just clutter the
// terminal output, and attrs are still rendered.
func (h *Handler) WithGroup(_ string) slog.Handler { return h }
