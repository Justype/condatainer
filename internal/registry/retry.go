package registry

import (
	"context"
	"errors"
	"fmt"
	"math/rand/v2"
	"net/http"
	"time"

	"oras.land/oras-go/v2/registry/remote/errcode"

	"github.com/Justype/condatainer/internal/logging"
)

// What a rate-limited push waits, when the registry gives no figure of its own.
//
// A minute is the floor because GitHub's secondary limits are measured in
// minutes ("please wait a few minutes before you try again"), and a client
// returning in one second feeds the problem it is recovering from. Four retries
// then span about fifteen minutes — long enough to outlast a limit, short enough
// that a stuck push is not mistaken for a slow one.
const (
	maxRateLimitRetries = 4
	baseRetryDelay      = time.Minute
	maxRetryDelay       = 15 * time.Minute
	// retryJitter spreads the waits of several machines rate limited together, so
	// they do not all return at once and trip the same limit again.
	retryJitter = 0.2
)

// serverDirectedDelay is implemented by an error carrying a wait the registry
// asked for — a Retry-After or rate-limit reset header. An interface because
// those headers reach no caller on their own: `errcode.ErrorResponse` keeps only
// method, URL, status, and the parsed error document, so capturing them needs
// [inspectTransport].
type serverDirectedDelay interface {
	RetryAfter() (time.Duration, bool)
}

// mutation is one registry write that may be refused for arriving too often.
type mutation struct {
	// verb is the log message stem, matching the transfer it belongs to:
	// verbUpload for a blob, verbPublish for a manifest, tag, or index.
	verb string
	// attrs say which write, and are repeated on every line about it so a pause
	// and the resume that follows can be read as one event:
	//
	//	[CNT] upload rate limited layer=4/11 retry=1/4 wait=1m0s
	//	[CNT] upload resumed layer=4/11
	attrs []any
	// do performs the write, once per attempt. Anything that cannot survive a
	// second attempt — above all a request body — must be created inside it: a
	// blob upload streams from an io.SectionReader with no GetBody, so a reader
	// reused across attempts sends nothing the second time.
	do func(context.Context) error
	// committed reports whether the write already landed, asked after a wait
	// because a failure may have lost the response rather than the write: GHCR
	// can commit a blob and then have an intermediary answer 403, and the next
	// attempt would retransmit gigabytes that are already there. Optional.
	committed func(context.Context) (bool, error)
}

// retryPolicy bounds how long a mutation is retried. The zero value is not
// usable; see defaultRetryPolicy.
type retryPolicy struct {
	maxRetries int
	baseDelay  time.Duration
	maxDelay   time.Duration
	// sleep waits, or reports why it stopped. Injected so tests do not.
	sleep func(context.Context, time.Duration) error
	// jitter returns a factor near 1. Injected so tests are deterministic.
	jitter func() float64
}

// retryPolicyKey carries a policy through a call chain that would otherwise
// thread one past every layer of the push, mirroring how the logger travels. It
// also lets a test remove the waiting without a package-level variable whose
// value depends on which test ran last.
type retryPolicyKey struct{}

// withRetryPolicy returns ctx carrying p for every mutation beneath it.
func withRetryPolicy(ctx context.Context, p retryPolicy) context.Context {
	return context.WithValue(ctx, retryPolicyKey{}, p)
}

// retryPolicyFrom returns the policy ctx carries, or the default.
func retryPolicyFrom(ctx context.Context) retryPolicy {
	if p, ok := ctx.Value(retryPolicyKey{}).(retryPolicy); ok {
		return p
	}
	return defaultRetryPolicy()
}

func defaultRetryPolicy() retryPolicy {
	return retryPolicy{
		maxRetries: maxRateLimitRetries,
		baseDelay:  baseRetryDelay,
		maxDelay:   maxRetryDelay,
		sleep:      sleepContext,
		jitter:     func() float64 { return 1 + retryJitter*(2*rand.Float64()-1) },
	}
}

// sleepContext waits for d, or returns as soon as ctx ends. A push pausing four
// minutes must still answer Ctrl-C in the moment it is pressed.
func sleepContext(ctx context.Context, d time.Duration) error {
	timer := time.NewTimer(d)
	defer timer.Stop()
	select {
	case <-ctx.Done():
		return ctx.Err()
	case <-timer.C:
		return nil
	}
}

// run performs m, waiting and retrying while the registry says it is rate
// limiting this client, and returns the classified error otherwise.
//
// Only a rate limit is retried, and only on positive evidence: spending 1, 2, 4,
// and 8 minutes to arrive at "your token is wrong" would be a worse bug than the
// one this fixes.
//
// A 401 is the single exception, retried once without a wait. ORAS fixes the
// Authorization header when it opens an upload session and reuses it to
// finalize, so a token expiring mid-layer cannot be refreshed in place — but a
// fresh attempt re-runs the session and mints a new one. A second consecutive
// 401 is a credential problem and is reported as one.
func (p retryPolicy) run(ctx context.Context, m mutation) error {
	log := logging.FromContext(ctx)
	authRetried := false

	for attempt := 0; ; attempt++ {
		// Spaced here rather than after a write, so an Exists hit costs nothing:
		// pushRange returns before reaching this.
		if err := pacerFrom(ctx).wait(ctx); err != nil {
			return err
		}
		// A slot per attempt, so what the transport saw is this attempt's.
		slot := &failureSlot{}
		err := classify(m.do(withFailureSlot(ctx, slot)))
		if err == nil {
			return nil
		}
		failure := slot.get()
		err = failure.annotate(err)

		switch {
		case errors.Is(err, ErrUnauthorized) && !authRetried && isExpiredCredential(err):
			// No wait: nothing is throttling us, the credential simply aged out
			// during an upload long enough for that to happen.
			authRetried = true
			log.Debug(m.verb+" re-authenticating, the credential expired mid-transfer", m.attrs...)
			attempt-- // Not a rate-limit attempt; it must not consume the budget.
			continue
		case !errors.Is(err, ErrRateLimited):
			return err
		case attempt >= p.maxRetries:
			return err
		}

		// Abandon the upload session before waiting. ORAS never resumes one, so
		// four retries would otherwise leave four half-open sessions behind.
		if failure != nil && failure.cancel != nil {
			failure.cancel(ctx)
		}

		// Warn, not Info, because it also closes the interactive progress line:
		// one left open across a four-minute sleep leaves a stalled-looking byte
		// count on screen.
		wait := p.delay(attempt, err)
		log.Warn(m.verb+" rate limited", m.logAttrs(
			"retry", fmt.Sprintf("%d/%d", attempt+1, p.maxRetries),
			"wait", wait.Round(time.Second),
		)...)
		if sleepErr := p.sleep(ctx, wait); sleepErr != nil {
			return sleepErr
		}

		if m.committed != nil {
			if done, checkErr := m.committed(ctx); checkErr == nil && done {
				log.Info(m.verb+" already accepted, the error lost the response", m.attrs...)
				return nil
			}
		}
		log.Info(m.verb+" resumed", m.attrs...)
	}
}

// logAttrs returns the mutation's attributes followed by extra, without writing
// through the shared backing array — every attempt appends to the same attrs.
func (m mutation) logAttrs(extra ...any) []any {
	out := make([]any, 0, len(m.attrs)+len(extra))
	return append(append(out, m.attrs...), extra...)
}

// delay reports how long to wait before attempt+1. What the registry asked for
// wins; otherwise the wait doubles from baseDelay with jitter, bounded.
func (p retryPolicy) delay(attempt int, err error) time.Duration {
	var directed serverDirectedDelay
	if errors.As(err, &directed) {
		if d, ok := directed.RetryAfter(); ok && d > 0 {
			return min(d, p.maxDelay)
		}
	}

	wait := p.baseDelay << attempt
	if wait > p.maxDelay || wait <= 0 {
		wait = p.maxDelay
	}
	if p.jitter != nil {
		wait = time.Duration(float64(wait) * p.jitter())
	}
	return wait
}

// isExpiredCredential reports whether a registry answered 401, as opposed to the
// 403 that means the credential is present and not enough. Only the first is
// worth a second attempt, because only the first is fixed by a fresh token.
func isExpiredCredential(err error) bool {
	var resp *errcode.ErrorResponse
	return errors.As(err, &resp) && resp.StatusCode == http.StatusUnauthorized
}
