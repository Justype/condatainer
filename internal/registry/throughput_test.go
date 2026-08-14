package registry

import (
	"context"
	"strings"
	"testing"
	"time"
)

// A registry that documents no per-upload limit gets no guard: guessing a budget
// would be worse than having none.
func TestThroughputGuardOnlyExistsWhereALimitIsDocumented(t *testing.T) {
	if g := newThroughputGuard(profileFor("registry.example.test")); g != nil {
		t.Errorf("an undocumented registry got a guard with a %v budget", g.budget)
	}
	g := newThroughputGuard(profileFor("ghcr.io"))
	if g == nil {
		t.Fatal("GHCR documents a ten-minute upload budget and got no guard")
	}
	// 80% of ten minutes, leaving room for variance: a projection landing at 99%
	// of the budget is a layer that fails.
	if want := 8 * time.Minute; g.budget != want {
		t.Errorf("budget = %v, want %v", g.budget, want)
	}
}

// Both profile limits bound a layer for the same reason — ORAS fixes the
// credential when it opens the session — so the tighter one wins.
func TestThroughputGuardTakesTheTighterLimit(t *testing.T) {
	tests := []struct {
		why     string
		profile transferProfile
		want    time.Duration
	}{
		{"only an upload timeout", transferProfile{UploadTimeout: 10 * time.Minute}, 8 * time.Minute},
		{"only a token lifetime", transferProfile{TokenLifetime: 5 * time.Minute}, 4 * time.Minute},
		{"the token expires first", transferProfile{UploadTimeout: 10 * time.Minute, TokenLifetime: 5 * time.Minute}, 4 * time.Minute},
		{"the upload budget expires first", transferProfile{UploadTimeout: 5 * time.Minute, TokenLifetime: time.Hour}, 4 * time.Minute},
	}
	for _, tt := range tests {
		g := newThroughputGuard(tt.profile)
		if g == nil || g.budget != tt.want {
			t.Errorf("%s: budget = %v, want %v", tt.why, g, tt.want)
		}
	}
}

// The only honest estimate is a measured one, so the first layer always goes.
func TestThroughputGuardNeverRejectsTheFirstLayer(t *testing.T) {
	g := newThroughputGuard(profileFor("ghcr.io"))
	if err := g.check(2 * gib); err != nil {
		t.Errorf("the first layer was rejected on an invented speed: %v", err)
	}
}

// The alternative to refusing is not success — an upload timeout is not a rate
// limit, so nothing retries it — but finding out ten minutes later.
func TestThroughputGuardRefusesALayerThatCannotFinish(t *testing.T) {
	g := newThroughputGuard(profileFor("ghcr.io"))
	// 2 GiB in twenty minutes: about 1.8 MiB/s, well under the 4.3 MiB/s an
	// 8-minute budget needs.
	g.observe(2*gib, 20*time.Minute)

	err := g.check(2 * gib)
	if err == nil {
		t.Fatal("a layer projected at twenty minutes was allowed into a ten-minute budget")
	}
	for _, want := range []string{"measured", "need", "projected"} {
		if !strings.Contains(err.Error(), want) {
			t.Errorf("%q does not report what was %s", err, want)
		}
	}
}

func TestThroughputGuardAllowsALayerThatFits(t *testing.T) {
	g := newThroughputGuard(profileFor("ghcr.io"))
	g.observe(2*gib, time.Minute) // ~34 MiB/s
	if err := g.check(2 * gib); err != nil {
		t.Errorf("a layer projected at one minute was refused: %v", err)
	}
}

// Conservative on purpose: a layer is refused on what the link has actually done
// at its worst, not on an average a single good sample could hold up.
func TestThroughputGuardUsesTheSlowestRecentSample(t *testing.T) {
	g := newThroughputGuard(profileFor("ghcr.io"))
	g.observe(2*gib, time.Second)    // absurdly fast
	g.observe(2*gib, 20*time.Minute) // and then the link degrades
	if err := g.check(2 * gib); err == nil {
		t.Error("a fast early sample outvoted a slow recent one")
	}
}

// Recent, because a rate measured half an hour ago says little about the next
// ten minutes: enough good samples must be able to clear a bad one.
func TestThroughputGuardForgetsOldSamples(t *testing.T) {
	g := newThroughputGuard(profileFor("ghcr.io"))
	g.observe(2*gib, 20*time.Minute) // one bad layer
	for range throughputSamples {
		g.observe(2*gib, time.Minute)
	}
	if err := g.check(2 * gib); err != nil {
		t.Errorf("a recovered link stayed blocked by an old sample: %v", err)
	}
}

// A layer the registry already had returns at once, and would read as infinite
// bandwidth if it were counted.
func TestPresentLayersAreNotMeasured(t *testing.T) {
	path := writeArtifact(t, t.TempDir(), "sample.sqf", "abcdefghijklmnop")
	blobs := newMemoryBlobs()
	guard := newThroughputGuard(profileFor("ghcr.io"))
	ctx := withThroughputGuard(context.Background(), guard)

	if _, err := pushArtifactLayers(ctx, blobs, path, MediaTypeOverlayBlob, 4); err != nil {
		t.Fatal(err)
	}
	measured := len(guard.rates)

	// Every layer is present now, so the second push transfers nothing.
	if _, err := pushArtifactLayers(ctx, blobs, path, MediaTypeOverlayBlob, 4); err != nil {
		t.Fatal(err)
	}
	if len(guard.rates) != measured {
		t.Errorf("a push that transferred nothing recorded %d samples", len(guard.rates)-measured)
	}
}

// A nil guard is the common case — most registries document nothing — so its
// methods have to be free rather than guarded at every call site.
func TestNilThroughputGuardIsANoOp(t *testing.T) {
	var g *throughputGuard
	g.observe(gib, time.Second)
	if err := g.check(gib); err != nil {
		t.Errorf("a nil guard refused a layer: %v", err)
	}
	if err := throughputGuardFrom(context.Background()).check(gib); err != nil {
		t.Errorf("a push with no guard was blocked: %v", err)
	}
}
