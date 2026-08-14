package registry

import (
	"context"
	"testing"
	"time"
)

// A profile is a claim about a registry, so an unknown one claims nothing —
// zero means "make no claim", never "no limit applies".
func TestProfileFor(t *testing.T) {
	ghcr := profileFor("ghcr.io")
	if ghcr.MaxLayerSize != 10_000_000_000 || ghcr.UploadTimeout != 10*time.Minute || ghcr.MinMutationGap != time.Second {
		t.Errorf("ghcr.io profile = %+v", ghcr)
	}
	// Host matching is normalized, because the caller passes whatever the base
	// happened to be spelled as.
	if profileFor("GHCR.IO") != ghcr || profileFor("  ghcr.io ") != ghcr {
		t.Error("the host match is case- or space-sensitive")
	}
	for _, host := range []string{"docker.io", "registry.example.test", "localhost:5000", ""} {
		if got := profileFor(host); got != (transferProfile{}) {
			t.Errorf("profileFor(%q) = %+v, want no claims", host, got)
		}
	}
}

// The clamp is the profile's one effect on the plan, and it must land below the
// limit rather than on it: 10 GB is not a whole number of GiB.
func TestGHCRProfileClampsTheLayerPlan(t *testing.T) {
	profile := profileFor("ghcr.io")
	size, reason := planLayerSize(400*gib, profile.MaxLayerSize)
	if reason != layerSizeClamped {
		t.Fatalf("a 400 GiB artifact was not clamped: %s", reason)
	}
	if size != 9*gib {
		t.Errorf("clamped to %d, want 9 GiB", size)
	}
	if size >= profile.MaxLayerSize {
		t.Errorf("clamped to %d, which is not below the %d limit", size, profile.MaxLayerSize)
	}
	// An unknown registry makes no claim, so nothing clamps.
	if _, reason := planLayerSize(400*gib, profileFor("registry.example.test").MaxLayerSize); reason == layerSizeClamped {
		t.Error("an unknown registry clamped the plan on a limit nobody documented")
	}
}

// testPacer records what it would have slept instead of sleeping.
func testPacer(gap time.Duration) (*pacer, *[]time.Duration) {
	waits := &[]time.Duration{}
	p := newPacer(gap)
	p.sleep = func(_ context.Context, d time.Duration) error {
		*waits = append(*waits, d)
		return nil
	}
	return p, waits
}

func TestPacerSpacesConsecutiveWrites(t *testing.T) {
	p, waits := testPacer(time.Second)
	ctx := context.Background()

	// The first write waits for nothing — there is no previous one to space from.
	if err := p.wait(ctx); err != nil {
		t.Fatal(err)
	}
	if len(*waits) != 0 {
		t.Errorf("the first write waited %v", *waits)
	}

	if err := p.wait(ctx); err != nil {
		t.Fatal(err)
	}
	if len(*waits) != 1 || (*waits)[0] <= 0 || (*waits)[0] > time.Second {
		t.Errorf("waits = %v, want one pause inside the gap", *waits)
	}
}

// The clock runs from when the last write finished, so a layer that took four
// minutes pays nothing. The gap exists to stop a burst of small requests, not to
// add latency to work that was already slow.
func TestPacerDoesNotDelayAfterSlowWork(t *testing.T) {
	p, waits := testPacer(10 * time.Millisecond)
	ctx := context.Background()

	if err := p.wait(ctx); err != nil {
		t.Fatal(err)
	}
	time.Sleep(20 * time.Millisecond) // stand in for a slow upload
	if err := p.wait(ctx); err != nil {
		t.Fatal(err)
	}
	if len(*waits) != 0 {
		t.Errorf("waited %v after work that already outlasted the gap", *waits)
	}
}

// An unknown registry pays nothing at all.
func TestPacerWithNoGapIsFree(t *testing.T) {
	p, waits := testPacer(0)
	for range 5 {
		if err := p.wait(context.Background()); err != nil {
			t.Fatal(err)
		}
	}
	if len(*waits) != 0 {
		t.Errorf("a zero gap still waited: %v", *waits)
	}
	// And a push to a registry with no profile carries no pacer at all.
	if err := pacerFrom(context.Background()).wait(context.Background()); err != nil {
		t.Errorf("a nil pacer is not a no-op: %v", err)
	}
}

// Ctrl-C during a pause must exit in the moment it is pressed, like every other
// wait in a push.
func TestPacerStopsOnCancellation(t *testing.T) {
	p := newPacer(time.Hour)
	ctx, cancel := context.WithCancel(context.Background())

	if err := p.wait(ctx); err != nil {
		t.Fatal(err)
	}
	cancel()

	done := make(chan error, 1)
	go func() { done <- p.wait(ctx) }()
	select {
	case err := <-done:
		if err == nil {
			t.Error("a cancelled pause reported success")
		}
	case <-time.After(5 * time.Second):
		t.Fatal("the pacer ignored cancellation")
	}
}

// A presence hit creates nothing, so there is nothing to pace against: pushRange
// returns before the retry loop, which is the only place the pacer is consulted.
func TestExistingBlobIsNotPaced(t *testing.T) {
	path := writeArtifact(t, t.TempDir(), "sample.sqf", "abcdefghij")
	blobs := newMemoryBlobs()
	p, waits := testPacer(time.Hour)
	ctx := withPacer(context.Background(), p)

	if _, err := pushArtifactLayers(ctx, blobs, path, MediaTypeOverlayBlob, 4); err != nil {
		t.Fatal(err)
	}
	sent := len(*waits)

	// Every layer is present now, so the second push writes nothing.
	if _, err := pushArtifactLayers(ctx, blobs, path, MediaTypeOverlayBlob, 4); err != nil {
		t.Fatal(err)
	}
	if len(*waits) != sent {
		t.Errorf("a push that uploaded nothing paced %d times", len(*waits)-sent)
	}
}
