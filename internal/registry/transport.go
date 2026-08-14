package registry

import (
	"bytes"
	"context"
	"io"
	"net/http"
	"strconv"
	"strings"
	"sync"
	"time"

	"github.com/Justype/condatainer/internal/config"
)

// What the transport preserves, adds, and how much of a body it keeps.
const (
	// maxCapturedBody bounds what is read from a failing response: enough for an
	// OCI error document and the plain text some intermediaries send instead.
	maxCapturedBody = 8 << 10

	// expectContinueMinSize is the body past which a request asks permission
	// before sending. Below it the extra round trip costs more than it saves;
	// above it, a rejection discovered after the body is a layer sent for
	// nothing.
	expectContinueMinSize = 4 << 20
)

// transferFailure is what a registry said when it refused, captured before ORAS
// throws it away: `errcode.ErrorResponse` keeps method, URL, status, and the
// parsed error document, and no headers at all — so Retry-After, the rate-limit
// counters, and the upload session reach no caller without this.
type transferFailure struct {
	StatusCode int
	Host       string
	RetryAfter string
	Reset      string
	Remaining  string
	// cancel abandons the upload session this request belonged to, or is nil.
	// A closure rather than a URL because cancelling needs the credential the
	// failed request carried, and this is the last place that has it.
	cancel func(context.Context)
}

// failureSlot holds the most recent failure of one attempt. One attempt may make
// several requests — a token exchange, then the write — and the last failure is
// the one that ended it.
type failureSlot struct {
	mu   sync.Mutex
	last *transferFailure
}

func (s *failureSlot) set(f *transferFailure) {
	s.mu.Lock()
	defer s.mu.Unlock()
	s.last = f
}

func (s *failureSlot) get() *transferFailure {
	if s == nil {
		return nil
	}
	s.mu.Lock()
	defer s.mu.Unlock()
	return s.last
}

type failureSlotKey struct{}

func withFailureSlot(ctx context.Context, s *failureSlot) context.Context {
	return context.WithValue(ctx, failureSlotKey{}, s)
}

func failureSlotFrom(ctx context.Context) *failureSlot {
	s, _ := ctx.Value(failureSlotKey{}).(*failureSlot)
	return s
}

// delay reports the wait the registry asked for.
//
// Retry-After first, as seconds or an HTTP date. Then a rate-limit reset, but
// only when the matching remaining count is zero: a reset timestamp on its own
// says when a window rolls over, not that anything is being refused now.
func (f *transferFailure) delay(now time.Time) (time.Duration, bool) {
	if v := strings.TrimSpace(f.RetryAfter); v != "" {
		if secs, err := strconv.Atoi(v); err == nil && secs >= 0 {
			return time.Duration(secs) * time.Second, true
		}
		if when, err := http.ParseTime(v); err == nil {
			if d := when.Sub(now); d > 0 {
				return d, true
			}
			return 0, true
		}
	}
	if remaining, err := strconv.Atoi(strings.TrimSpace(f.Remaining)); err != nil || remaining > 0 {
		return 0, false
	}
	if unix, err := strconv.ParseInt(strings.TrimSpace(f.Reset), 10, 64); err == nil {
		if d := time.Unix(unix, 0).Sub(now); d > 0 {
			return d, true
		}
	}
	return 0, false
}

// directedError attaches a server-directed wait to an error, satisfying
// [serverDirectedDelay] so the retry policy prefers it over local backoff.
type directedError struct {
	error
	wait time.Duration
}

func (e directedError) RetryAfter() (time.Duration, bool) { return e.wait, true }
func (e directedError) Unwrap() error                     { return e.error }

// annotate returns err carrying whatever wait the registry asked for, or err
// unchanged. A failure with no usable header leaves the local schedule in charge.
func (f *transferFailure) annotate(err error) error {
	if f == nil || err == nil {
		return err
	}
	if wait, ok := f.delay(time.Now()); ok {
		return directedError{error: err, wait: wait}
	}
	return err
}

// inspectTransport preserves what ORAS discards and adds what it never sends.
// It makes no policy decisions: it records a failure into the slot the retry
// loop put in the context, restores every body it reads, and otherwise leaves
// the exchange alone.
type inspectTransport struct{ next http.RoundTripper }

func (t *inspectTransport) RoundTrip(req *http.Request) (*http.Response, error) {
	// Ask permission before sending a large body. ORAS streams straight into the
	// PUT, so without this an auth failure, quota rejection, proxy body limit, or
	// rate limit is discovered only after a whole layer has crossed the wire. Go
	// falls back to sending after ExpectContinueTimeout, so a registry that does
	// not implement it costs a second.
	if req.Body != nil && req.ContentLength >= expectContinueMinSize &&
		(req.Method == http.MethodPut || req.Method == http.MethodPatch) {
		req.Header.Set("Expect", "100-continue")
	}

	resp, err := t.next.RoundTrip(req)
	if err != nil || resp.StatusCode < 400 {
		return resp, err
	}

	failure := &transferFailure{
		StatusCode: resp.StatusCode,
		Host:       req.URL.Host,
		RetryAfter: resp.Header.Get("Retry-After"),
		Reset:      firstHeader(resp.Header, "X-RateLimit-Reset", "RateLimit-Reset"),
		Remaining:  firstHeader(resp.Header, "X-RateLimit-Remaining", "RateLimit-Remaining"),
	}
	if session := uploadSessionURL(req); session != "" {
		auth := req.Header.Get("Authorization")
		failure.cancel = func(ctx context.Context) { t.cancelUpload(ctx, session, auth) }
	}
	if slot := failureSlotFrom(req.Context()); slot != nil {
		slot.set(failure)
	}

	// Read a bounded prefix and put it back. ORAS decodes the body itself, and a
	// consumed one would turn a specific registry error into "Forbidden".
	if resp.Body != nil {
		prefix, readErr := io.ReadAll(io.LimitReader(resp.Body, maxCapturedBody))
		if readErr == nil {
			rest := resp.Body
			resp.Body = struct {
				io.Reader
				io.Closer
			}{Reader: io.MultiReader(bytes.NewReader(prefix), rest), Closer: rest}
		}
	}
	return resp, nil
}

// cancelUpload abandons a blob upload session, best effort. A failed finalize
// leaves the session open and ORAS never resumes one, so four retries would
// leave four behind. Failures are ignored: this tidies after an error already
// being handled and must never become the error the caller sees.
func (t *inspectTransport) cancelUpload(ctx context.Context, session, authorization string) {
	ctx, stop := context.WithTimeout(ctx, 10*time.Second)
	defer stop()

	req, err := http.NewRequestWithContext(ctx, http.MethodDelete, session, nil)
	if err != nil {
		return
	}
	if authorization != "" {
		req.Header.Set("Authorization", authorization)
	}
	resp, err := t.next.RoundTrip(req)
	if err != nil {
		return
	}
	io.Copy(io.Discard, io.LimitReader(resp.Body, maxCapturedBody)) //nolint:errcheck
	resp.Body.Close()                                               //nolint:errcheck
}

// uploadSessionURL returns the session a request was finalizing, or "". Only a
// PUT has one to abandon; the digest query is dropped because the session is
// addressed without it.
func uploadSessionURL(req *http.Request) string {
	if req.Method != http.MethodPut || !strings.Contains(req.URL.Path, "/blobs/upload") {
		return ""
	}
	session := *req.URL
	session.RawQuery = ""
	return session.String()
}

func firstHeader(h http.Header, names ...string) string {
	for _, name := range names {
		if v := h.Get(name); v != "" {
			return v
		}
	}
	return ""
}

// userAgent identifies CondaTainer to a registry. Providers key abuse
// heuristics, exemptions, and support triage on it, and the "oras-go" default is
// no answer when reporting a secondary rate limit to GitHub.
func userAgent() string {
	return "condatainer/" + config.VERSION + " (oras-go/v2)"
}
