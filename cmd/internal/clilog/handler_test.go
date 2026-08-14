package clilog

import (
	"bytes"
	"context"
	"log/slog"
	"strings"
	"sync"
	"testing"
	"time"
)

// newTestHandler returns a handler drawing into buf as though it were a
// terminal, which is the only mode where a progress line exists to leave open.
func newTestHandler() (*Handler, *bytes.Buffer) {
	var buf bytes.Buffer
	return &Handler{state: &progressState{}, out: &buf, tty: true}, &buf
}

func progressRecord(msg string, final, last bool) slog.Record {
	r := slog.NewRecord(time.Now(), slog.LevelInfo, msg, 0)
	r.AddAttrs(
		slog.String("kind", "progress"),
		slog.Bool("final", final),
		slog.Bool("last", last),
	)
	return r
}

// A progress line is drawn in place, so it is still open when the next record
// arrives. A rate-limit pause logs at Warn precisely so that record closes it —
// otherwise a stalled byte count sits on screen for the length of the wait.
func TestProgressLineIsClosedBeforeTheNextMessage(t *testing.T) {
	h, buf := newTestHandler()
	ctx := context.Background()

	if err := h.Handle(ctx, progressRecord("upload progress done=1.25 GB total=2.00 GB layer=4/11", false, false)); err != nil {
		t.Fatal(err)
	}
	if !h.state.active {
		t.Fatal("a progress record left no line open")
	}
	if strings.HasSuffix(buf.String(), "\n") {
		t.Fatal("the progress line was terminated; it must stay open to be redrawn")
	}

	warn := slog.NewRecord(time.Now(), slog.LevelWarn, "upload rate limited", 0)
	warn.AddAttrs(slog.String("layer", "4/11"))
	if err := h.Handle(ctx, warn); err != nil {
		t.Fatal(err)
	}
	if h.state.active {
		t.Error("the progress line is still open across the wait that follows")
	}
	if !strings.HasSuffix(buf.String(), "\n") {
		t.Error("no newline was written to close the progress line")
	}
}

// The final report of the last layer closes its own line, and a second close on
// the next record would print a blank one.
func TestFinalProgressClosesItsOwnLineOnce(t *testing.T) {
	h, buf := newTestHandler()
	ctx := context.Background()

	if err := h.Handle(ctx, progressRecord("upload progress done=2.00 GB total=2.00 GB layer=11/11", true, true)); err != nil {
		t.Fatal(err)
	}
	if h.state.active {
		t.Error("a final report on the last layer left the line open")
	}
	closed := buf.String()

	if err := h.Handle(ctx, slog.NewRecord(time.Now(), slog.LevelInfo, "installed", 0)); err != nil {
		t.Fatal(err)
	}
	if buf.String() != closed {
		t.Errorf("a blank line was printed after an already-closed progress line: %q", buf.String()[len(closed):])
	}
}

// Every message a command prints after a transfer has to come through here, or
// it lands on the end of the open progress line: a finished push once read
// "layer=10/11published /path/to/image.sqf".
func TestASuccessMessageEndsTheProgressLine(t *testing.T) {
	h, buf := newTestHandler()
	ctx := context.Background()

	if err := h.Handle(ctx, progressRecord("upload progress done=2.00 GB total=2.00 GB layer=10/11", false, false)); err != nil {
		t.Fatal(err)
	}
	done := slog.NewRecord(time.Now(), slog.LevelInfo, "published /images/testdata--layers--20g.sqf", 0)
	done.AddAttrs(slog.String("kind", "success"))
	if err := h.Handle(ctx, done); err != nil {
		t.Fatal(err)
	}
	if h.state.active {
		t.Error("the progress line is still open after the success message")
	}
	if !strings.HasSuffix(buf.String(), "\n") {
		t.Error("the success message was appended to the progress line")
	}
}

// A final report on a layer that is not the last keeps the line open, because
// the next layer redraws it in place.
func TestFinalProgressOnAnEarlyLayerKeepsTheLine(t *testing.T) {
	h, _ := newTestHandler()
	if err := h.Handle(context.Background(),
		progressRecord("upload progress done=2.00 GB total=2.00 GB layer=4/11", true, false)); err != nil {
		t.Fatal(err)
	}
	if !h.state.active {
		t.Error("layer 4 of 11 closed the progress line; layer 5 has nowhere to draw")
	}
}

// WithAttrs copies the handler, and the progress state must stay shared or two
// loggers derived from one handler would each think the line was theirs.
func TestWithAttrsSharesProgressState(t *testing.T) {
	h, _ := newTestHandler()
	derived, ok := h.WithAttrs([]slog.Attr{slog.String("artifact", "x")}).(*Handler)
	if !ok {
		t.Fatal("WithAttrs did not return a *Handler")
	}
	if derived.state != h.state {
		t.Error("the derived handler tracks its own progress line")
	}
	if derived.out != h.out || derived.tty != h.tty {
		t.Error("the derived handler lost its destination")
	}
}

// Handlers are shared across goroutines; the state is what guards the line.
func TestHandleIsConcurrencySafe(t *testing.T) {
	h, _ := newTestHandler()
	ctx := context.Background()
	var wg sync.WaitGroup
	for i := range 50 {
		wg.Add(1)
		go func() {
			defer wg.Done()
			if i%2 == 0 {
				_ = h.Handle(ctx, progressRecord("upload progress", false, false))
				return
			}
			_ = h.Handle(ctx, slog.NewRecord(time.Now(), slog.LevelWarn, "upload rate limited", 0))
		}()
	}
	wg.Wait()
}

// An interrupt arrives with no log record to close the progress line the way an
// ordinary message does, so the signal handler asks for it directly.
func TestEndProgressLineClosesAnOpenLine(t *testing.T) {
	h, buf := newTestHandler()
	current.Store(h)
	t.Cleanup(func() { current.Store(nil) })

	if err := h.Handle(context.Background(),
		progressRecord("download progress done=4.38 GB total=20.00 GB", false, false)); err != nil {
		t.Fatal(err)
	}
	EndProgressLine()

	if h.state.active {
		t.Error("the progress line survived the interrupt")
	}
	if !strings.HasSuffix(buf.String(), "\n") {
		t.Error("no newline was written to close the line")
	}

	// The old signal handler wrote a newline unconditionally, so whatever logged
	// next closed the line a second time and left a blank one behind.
	closed := buf.String()
	EndProgressLine()
	if buf.String() != closed {
		t.Errorf("a second call printed a blank line: %q", buf.String()[len(closed):])
	}
}

// Nothing drawing means nothing to close, including before any command runs.
func TestEndProgressLineIsSafeWhenNothingIsDrawing(t *testing.T) {
	current.Store(nil)
	EndProgressLine()

	h, buf := newTestHandler()
	current.Store(h)
	t.Cleanup(func() { current.Store(nil) })
	EndProgressLine()
	if buf.Len() != 0 {
		t.Errorf("closed a line that was never open: %q", buf)
	}
}
