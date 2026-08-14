package registry

import (
	"context"
	"errors"
	"fmt"
	"log/slog"
	"strings"
	"testing"
	"time"

	"github.com/Justype/condatainer/internal/logging"
)

// testPolicy is the real policy with the waiting removed: it records what would
// have been slept so the schedule can be asserted, and jitter pinned so it can.
type testPolicy struct {
	retryPolicy
	waits *[]time.Duration
}

func newTestPolicy() testPolicy {
	waits := &[]time.Duration{}
	p := defaultRetryPolicy()
	p.jitter = func() float64 { return 1 }
	p.sleep = func(_ context.Context, d time.Duration) error {
		*waits = append(*waits, d)
		return nil
	}
	return testPolicy{retryPolicy: p, waits: waits}
}

// rateLimited is the error a registry gives when it wants the client to pause,
// in the shape GHCR actually produces: a plain 403 DENIED whose message is the
// only thing separating it from a permission failure.
func rateLimited() error { return codedErr(403, "DENIED", ghcrSecondaryLimitMessage) }

func TestRetryWaitsOutARateLimit(t *testing.T) {
	p := newTestPolicy()
	attempts := 0

	err := p.run(context.Background(), mutation{verb: verbUpload, attrs: []any{"layer", "34/41"}, do: func(context.Context) error {
		attempts++
		if attempts < 3 {
			return rateLimited()
		}
		return nil
	}})
	if err != nil {
		t.Fatalf("run: %v", err)
	}
	if attempts != 3 {
		t.Errorf("attempts = %d, want 3", attempts)
	}
	// Doubling from one minute, which is the floor because GitHub's secondary
	// limits are measured in minutes.
	if want := []time.Duration{time.Minute, 2 * time.Minute}; !equalDurations(*p.waits, want) {
		t.Errorf("waits = %v, want %v", *p.waits, want)
	}
}

// The failure mode that matters most: an ordinary denial must not spend fifteen
// minutes of backoff to arrive at "your token is wrong".
func TestRetryDoesNotRetryARealDenial(t *testing.T) {
	p := newTestPolicy()
	attempts := 0

	err := p.run(context.Background(), mutation{verb: verbUpload, attrs: []any{"layer", "1/11"}, do: func(context.Context) error {
		attempts++
		return codedErr(403, "DENIED", "permission_denied: write_package")
	}})
	if !errors.Is(err, ErrUnauthorized) {
		t.Errorf("err = %v, want ErrUnauthorized", err)
	}
	if attempts != 1 {
		t.Errorf("attempts = %d, want the denial reported at once", attempts)
	}
	if len(*p.waits) != 0 {
		t.Errorf("waited %v before reporting a permission failure", *p.waits)
	}
}

func TestRetryIsBounded(t *testing.T) {
	p := newTestPolicy()
	attempts := 0

	err := p.run(context.Background(), mutation{verb: verbUpload, attrs: []any{"layer", "34/41"}, do: func(context.Context) error {
		attempts++
		return rateLimited()
	}})
	if !errors.Is(err, ErrRateLimited) {
		t.Errorf("err = %v, want ErrRateLimited", err)
	}
	if want := maxRateLimitRetries + 1; attempts != want {
		t.Errorf("attempts = %d, want %d", attempts, want)
	}
	if want := []time.Duration{1 * time.Minute, 2 * time.Minute, 4 * time.Minute, 8 * time.Minute}; !equalDurations(*p.waits, want) {
		t.Errorf("waits = %v, want %v", *p.waits, want)
	}
}

// A registry can commit a blob and then lose the response — GHCR's refusal comes
// from an intermediary in front of the storage that already accepted it. Without
// the recheck the next attempt retransmits gigabytes that are already there.
func TestRetryStopsWhenTheWriteLandedAnyway(t *testing.T) {
	p := newTestPolicy()
	attempts := 0

	err := p.run(context.Background(), mutation{
		verb:  verbUpload,
		attrs: []any{"layer", "34/41"},
		do: func(context.Context) error {
			attempts++
			return rateLimited()
		},
		committed: func(context.Context) (bool, error) { return true, nil },
	})
	if err != nil {
		t.Fatalf("run: %v", err)
	}
	if attempts != 1 {
		t.Errorf("attempts = %d, want the committed write to end it after one", attempts)
	}
}

// ORAS fixes the Authorization header when it opens an upload session and reuses
// it to finalize the blob, so a token that ages out mid-layer is unrecoverable
// in place — but a fresh attempt mints a new one.
func TestRetryReauthenticatesOnceOnAnExpiredToken(t *testing.T) {
	p := newTestPolicy()
	attempts := 0

	err := p.run(context.Background(), mutation{verb: verbUpload, attrs: []any{"layer", "5/11"}, do: func(context.Context) error {
		attempts++
		if attempts == 1 {
			return codedErr(401, "UNAUTHORIZED", "authentication required")
		}
		return nil
	}})
	if err != nil {
		t.Fatalf("run: %v", err)
	}
	if attempts != 2 {
		t.Errorf("attempts = %d, want one re-authentication", attempts)
	}
	if len(*p.waits) != 0 {
		t.Errorf("waited %v to re-authenticate; nothing is throttling us", *p.waits)
	}
}

// A credential that is wrong stays wrong, so the second 401 is reported.
func TestRetryReportsAPersistentAuthFailure(t *testing.T) {
	p := newTestPolicy()
	attempts := 0

	err := p.run(context.Background(), mutation{verb: verbUpload, attrs: []any{"layer", "5/11"}, do: func(context.Context) error {
		attempts++
		return codedErr(401, "UNAUTHORIZED", "authentication required")
	}})
	if !errors.Is(err, ErrUnauthorized) {
		t.Errorf("err = %v, want ErrUnauthorized", err)
	}
	if attempts != 2 {
		t.Errorf("attempts = %d, want exactly one retry", attempts)
	}
}

// A re-authentication must not consume the rate-limit budget: the two failures
// are unrelated, and an expired token early on should not shorten the wait a
// genuine limit later gets.
func TestReauthenticationDoesNotSpendTheRetryBudget(t *testing.T) {
	p := newTestPolicy()
	attempts := 0

	err := p.run(context.Background(), mutation{verb: verbUpload, attrs: []any{"layer", "5/11"}, do: func(context.Context) error {
		attempts++
		if attempts == 1 {
			return codedErr(401, "UNAUTHORIZED", "authentication required")
		}
		return rateLimited()
	}})
	if !errors.Is(err, ErrRateLimited) {
		t.Fatalf("err = %v, want ErrRateLimited", err)
	}
	if len(*p.waits) != maxRateLimitRetries {
		t.Errorf("got %d waits, want the full budget of %d", len(*p.waits), maxRateLimitRetries)
	}
}

// Ctrl-C during a four-minute backoff must exit in the moment it is pressed.
func TestRetryStopsOnCancellation(t *testing.T) {
	p := defaultRetryPolicy()
	p.sleep = sleepContext // the real one, so the context does the work
	ctx, cancel := context.WithCancel(context.Background())

	attempts := 0
	done := make(chan error, 1)
	go func() {
		done <- p.run(ctx, mutation{verb: verbUpload, attrs: []any{"layer", "34/41"}, do: func(context.Context) error {
			attempts++
			cancel() // fail, then interrupt the wait that follows
			return rateLimited()
		}})
	}()

	select {
	case err := <-done:
		if !errors.Is(err, context.Canceled) {
			t.Errorf("err = %v, want the cancellation", err)
		}
		if attempts != 1 {
			t.Errorf("attempts = %d, want the wait interrupted before a second", attempts)
		}
	case <-time.After(5 * time.Second):
		t.Fatal("run did not return; the backoff ignored cancellation")
	}
}

// directedErr carries a wait the registry asked for, which is what the transport
// will supply once it preserves Retry-After.
type directedErr struct {
	error
	after time.Duration
}

func (d directedErr) RetryAfter() (time.Duration, bool) { return d.after, true }
func (d directedErr) Unwrap() error                     { return d.error }

// What the server asked for beats the local schedule, in both directions: it can
// ask for less than a minute as well as more.
func TestServerDirectedDelayWins(t *testing.T) {
	p := newTestPolicy()
	attempts := 0

	err := p.run(context.Background(), mutation{verb: verbUpload, attrs: []any{"layer", "2/11"}, do: func(context.Context) error {
		attempts++
		if attempts == 1 {
			return directedErr{error: rateLimited(), after: 5 * time.Second}
		}
		return nil
	}})
	if err != nil {
		t.Fatalf("run: %v", err)
	}
	if want := []time.Duration{5 * time.Second}; !equalDurations(*p.waits, want) {
		t.Errorf("waits = %v, want %v", *p.waits, want)
	}
}

// Nothing implements serverDirectedDelay until a transport preserves the header,
// so the local backoff has to stand on its own.
func TestDelayFallsBackWhenNoHeaderSurvives(t *testing.T) {
	p := defaultRetryPolicy()
	p.jitter = func() float64 { return 1 }
	for attempt, want := range []time.Duration{1, 2, 4, 8} {
		if got := p.delay(attempt, rateLimited()); got != want*time.Minute {
			t.Errorf("delay(%d) = %v, want %v", attempt, got, want*time.Minute)
		}
	}
	// Bounded, however many attempts a future policy allows.
	if got := p.delay(20, rateLimited()); got != maxRetryDelay {
		t.Errorf("delay(20) = %v, want it capped at %v", got, maxRetryDelay)
	}
}

// Machines rate limited together must not return at the same instant and trip
// the same limit again.
func TestJitterSpreadsTheWaits(t *testing.T) {
	p := defaultRetryPolicy()
	seen := map[time.Duration]bool{}
	for range 50 {
		d := p.delay(0, rateLimited())
		if d < time.Duration(float64(baseRetryDelay)*(1-retryJitter)) ||
			d > time.Duration(float64(baseRetryDelay)*(1+retryJitter)) {
			t.Fatalf("delay %v is outside the jitter band around %v", d, baseRetryDelay)
		}
		seen[d] = true
	}
	if len(seen) < 2 {
		t.Error("every wait was identical; jitter is not spreading anything")
	}
}

// recordingHandler captures what a pause would print, in the order it prints.
type recordingHandler struct {
	lines *[]string
	attrs []slog.Attr
}

func (h recordingHandler) Enabled(context.Context, slog.Level) bool { return true }
func (h recordingHandler) WithGroup(string) slog.Handler            { return h }

func (h recordingHandler) WithAttrs(attrs []slog.Attr) slog.Handler {
	h.attrs = append(append([]slog.Attr(nil), h.attrs...), attrs...)
	return h
}

func (h recordingHandler) Handle(_ context.Context, r slog.Record) error {
	parts := []string{r.Message}
	r.Attrs(func(a slog.Attr) bool {
		parts = append(parts, fmt.Sprintf("%s=%v", a.Key, a.Value.Any()))
		return true
	})
	*h.lines = append(*h.lines, strings.Join(parts, " "))
	return nil
}

// A paused push has to say so, and say the same thing on the way back, or an
// operator watching a three-hour transfer cannot tell a wait from a hang.
//
// The pause is logged at Warn deliberately: the terminal handler closes an open
// progress line on the first record that is not progress, so a wait cannot leave
// a stalled-looking byte count on screen for four minutes.
func TestRateLimitStatesAreReported(t *testing.T) {
	var lines []string
	ctx := logging.WithLogger(context.Background(),
		slog.New(recordingHandler{lines: &lines}))
	p := newTestPolicy()

	attempts := 0
	if err := p.run(ctx, mutation{verb: verbUpload, attrs: []any{"layer", "4/11"}, do: func(context.Context) error {
		attempts++
		if attempts == 1 {
			return rateLimited()
		}
		return nil
	}}); err != nil {
		t.Fatalf("run: %v", err)
	}

	want := []string{
		"upload rate limited layer=4/11 retry=1/4 wait=1m0s",
		"upload resumed layer=4/11",
	}
	if len(lines) != len(want) {
		t.Fatalf("logged %v, want %v", lines, want)
	}
	for i := range want {
		if lines[i] != want[i] {
			t.Errorf("line %d = %q, want %q", i, lines[i], want[i])
		}
	}
}

// The mutation's attributes are repeated on every line about it, so appending
// per-attempt detail must not write through the shared backing array.
func TestRetryAttrsAreNotSharedAcrossAttempts(t *testing.T) {
	var lines []string
	ctx := logging.WithLogger(context.Background(),
		slog.New(recordingHandler{lines: &lines}))
	p := newTestPolicy()

	_ = p.run(ctx, mutation{verb: verbUpload, attrs: []any{"layer", "34/41"}, do: func(context.Context) error {
		return rateLimited()
	}})
	for _, line := range lines {
		if strings.Count(line, "retry=") > 1 || strings.Count(line, "layer=") > 1 {
			t.Errorf("attributes accumulated across attempts: %q", line)
		}
	}
}

func equalDurations(got, want []time.Duration) bool {
	if len(got) != len(want) {
		return false
	}
	for i := range got {
		if got[i] != want[i] {
			return false
		}
	}
	return true
}
