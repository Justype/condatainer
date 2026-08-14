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
	"io"
	"log/slog"
	"os"
	"strings"
	"sync"
	"sync/atomic"

	"golang.org/x/term"

	"github.com/Justype/condatainer/internal/utils"
)

// Handler is a slog.Handler that dispatches records to utils.Print*.
type Handler struct {
	attrs []slog.Attr
	state *progressState
	// out is where the in-place progress line is drawn, and tty whether that
	// destination can render one. Fields rather than os.Stderr and a live
	// terminal check so the line's lifecycle can be tested: a progress line left
	// open across a multi-minute rate-limit wait is a real failure mode and one
	// no test could otherwise reach.
	out io.Writer
	tty bool
}

type progressState struct {
	sync.Mutex
	active bool
}

// current is the handler an interrupt should close the progress line on.
//
// A package-level slot because the signal handler runs far from where the
// handler is built, and there is one of these per process by construction —
// New is called once per command, and the last one is the one drawing.
var current atomic.Pointer[Handler]

// New returns a fresh Handler writing to stderr.
func New() *Handler {
	h := &Handler{
		state: &progressState{},
		out:   os.Stderr,
		tty:   term.IsTerminal(int(os.Stderr.Fd())),
	}
	current.Store(h)
	return h
}

// EndProgressLine closes an in-place progress line if one is open, and does
// nothing otherwise.
//
// For an interrupt, which arrives with no log record to close the line the way
// an ordinary message does. Writing a newline unconditionally instead — as the
// signal handler used to — leaves a blank line whenever a record follows,
// because that record closes the line again.
func EndProgressLine() {
	if h := current.Load(); h != nil {
		h.endLine()
	}
}

// endLine ends an open progress line. The caller must not hold state's lock.
func (h *Handler) endLine() {
	h.state.Lock()
	defer h.state.Unlock()
	h.endLineLocked()
}

func (h *Handler) endLineLocked() {
	if h.state.active {
		fmt.Fprintln(h.out)
		h.state.active = false
	}
}

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
	var final, last bool
	var extras []string
	visit := func(a slog.Attr) {
		if a.Key == "kind" {
			kind = a.Value.String()
			return
		}
		if a.Key == "final" {
			final = a.Value.Bool()
			return
		}
		if a.Key == "last" {
			last = a.Value.Bool()
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

	h.state.Lock()
	defer h.state.Unlock()
	if kind == "progress" && h.tty {
		fmt.Fprintf(h.out, "\r\x1b[2K[CNT] %s", msg)
		h.state.active = true
		if final && last {
			fmt.Fprintln(h.out)
			h.state.active = false
		}
		return nil
	}
	// Any other record ends an open progress line first. This is what keeps a
	// rate-limit pause from leaving a stalled-looking byte count on screen for
	// the length of the wait.
	h.endLineLocked()

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
