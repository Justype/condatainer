package registry

import (
	"context"
	"fmt"
	"sync"
	"time"

	"github.com/Justype/condatainer/internal/utils"
)

// How the guard decides a layer will not finish in time.
const (
	// throughputSamples is how many recent layers the estimate is drawn from.
	// Recent because a link degrades; few because the slowest of them is used.
	throughputSamples = 3
	// uploadBudgetMargin is the fraction of the limit a layer is expected to fit
	// inside. The rest is room for variance: a projection landing at 99% of the
	// budget is a layer that fails.
	uploadBudgetMargin = 0.8
)

// throughputGuard refuses to start a layer that measurement says cannot finish
// inside the destination's per-upload budget. The alternative is finding out ten
// minutes later: an upload timeout is not a rate limit, so nothing retries it.
//
// A zero budget disables it, which is every registry documenting no limit.
type throughputGuard struct {
	budget time.Duration

	mu    sync.Mutex
	rates []float64 // bytes per second, oldest first
}

// newThroughputGuard builds the guard for a push, or nil when the destination
// claims no per-upload limit. Both of the profile's limits bound a layer: an
// upload that outlives its token cannot be saved in place either, because ORAS
// fixes the credential when it opens the session.
func newThroughputGuard(profile transferProfile) *throughputGuard {
	budget := profile.UploadTimeout
	if profile.TokenLifetime > 0 && (budget == 0 || profile.TokenLifetime < budget) {
		budget = profile.TokenLifetime
	}
	if budget <= 0 {
		return nil
	}
	return &throughputGuard{budget: time.Duration(float64(budget) * uploadBudgetMargin)}
}

// observe records how long a completed layer took.
func (g *throughputGuard) observe(bytes int64, elapsed time.Duration) {
	if g == nil || bytes <= 0 || elapsed <= 0 {
		return
	}
	g.mu.Lock()
	defer g.mu.Unlock()
	g.rates = append(g.rates, float64(bytes)/elapsed.Seconds())
	if len(g.rates) > throughputSamples {
		g.rates = g.rates[len(g.rates)-throughputSamples:]
	}
}

// check reports why a layer of bytes should not be started, or nil. The first
// layer is always allowed: the only honest estimate is a measured one.
func (g *throughputGuard) check(bytes int64) error {
	if g == nil {
		return nil
	}
	g.mu.Lock()
	rate := slowestRate(g.rates)
	g.mu.Unlock()
	if rate <= 0 {
		return nil
	}

	projected := time.Duration(float64(bytes) / rate * float64(time.Second))
	if projected <= g.budget {
		return nil
	}
	needed := float64(bytes) / g.budget.Seconds()
	return fmt.Errorf(
		"this link cannot carry a %s layer inside the registry's upload budget: measured %s/s, need %s/s to finish in %s (projected %s)",
		utils.FormatSize(bytes), utils.FormatSize(int64(rate)), utils.FormatSize(int64(needed)),
		g.budget.Round(time.Second), projected.Round(time.Second))
}

// slowestRate is the conservative estimate: the worst of the recent samples,
// not an average a single good sample could hold up.
func slowestRate(rates []float64) float64 {
	slowest := 0.0
	for _, rate := range rates {
		if slowest == 0 || rate < slowest {
			slowest = rate
		}
	}
	return slowest
}

type throughputGuardKey struct{}

func withThroughputGuard(ctx context.Context, g *throughputGuard) context.Context {
	return context.WithValue(ctx, throughputGuardKey{}, g)
}

// throughputGuardFrom returns the guard ctx carries, or nil — whose methods are
// no-ops.
func throughputGuardFrom(ctx context.Context) *throughputGuard {
	g, _ := ctx.Value(throughputGuardKey{}).(*throughputGuard)
	return g
}
