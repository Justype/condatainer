package registry

import (
	"context"
	"strings"
	"sync"
	"time"
)

// transferProfile is what a registry is known to enforce. Guardrails, not
// protocol branching: it changes only how much is asked for at once.
//
// A zero field means "make no claim", not "no limit applies" — an undocumented
// guess written down here becomes stale policy nobody can correct.
type transferProfile struct {
	// MaxLayerSize is the hard per-layer ceiling, which the layer plan clamps to.
	MaxLayerSize int64
	// UploadTimeout is how long one blob upload may take.
	UploadTimeout time.Duration
	// TokenLifetime bounds a layer as surely as UploadTimeout does: ORAS fixes
	// the credential when it opens an upload session and reuses it to finalize,
	// so a layer that outlives its token cannot be saved in place.
	TokenLifetime time.Duration
	// MinMutationGap is the pause between completed writes.
	MinMutationGap time.Duration
}

// profileFor returns what is known about a registry host, or an empty profile.
// Deliberately short: a profile is a claim someone must be able to point at
// documentation for.
func profileFor(host string) transferProfile {
	switch strings.ToLower(strings.TrimSpace(host)) {
	case ghcrHost:
		return transferProfile{
			// Documented as 10 GB. A clamp rounds down to whole GiB, so this
			// lands at 9 GiB — below the limit rather than at it.
			MaxLayerSize:  10_000_000_000,
			UploadTimeout: 10 * time.Minute,
			// The observed refusal was a request counter, and spacing is the
			// cheapest way to spend it more slowly: eleven seconds across a
			// 20 GiB push.
			MinMutationGap: time.Second,
		}
	default:
		return transferProfile{}
	}
}

// pacer spaces out completed writes. The zero gap makes every method a no-op, so
// an unknown registry pays nothing.
type pacer struct {
	gap   time.Duration
	sleep func(context.Context, time.Duration) error

	mu   sync.Mutex
	last time.Time
}

func newPacer(gap time.Duration) *pacer {
	return &pacer{gap: gap, sleep: sleepContext}
}

// wait blocks until at least gap has passed since the previous write, then
// records this one as the most recent.
//
// The clock starts when the last write *finished*, so a layer that took four
// minutes pays nothing: the gap stops a burst of small requests, it does not add
// latency to work that was already slow.
func (p *pacer) wait(ctx context.Context) error {
	if p == nil || p.gap <= 0 {
		return nil
	}
	p.mu.Lock()
	remaining := time.Duration(0)
	if !p.last.IsZero() {
		remaining = p.gap - time.Since(p.last)
	}
	p.mu.Unlock()

	if remaining > 0 {
		if err := p.sleep(ctx, remaining); err != nil {
			return err
		}
	}
	p.mu.Lock()
	p.last = time.Now()
	p.mu.Unlock()
	return nil
}

type pacerKey struct{}

func withPacer(ctx context.Context, p *pacer) context.Context {
	return context.WithValue(ctx, pacerKey{}, p)
}

// pacerFrom returns the pacer ctx carries, or nil — whose methods are no-ops.
func pacerFrom(ctx context.Context) *pacer {
	p, _ := ctx.Value(pacerKey{}).(*pacer)
	return p
}
