package registry

import (
	"testing"
)

const (
	gib = int64(1) << 30
	// ghcrMaxLayer is GHCR's documented 10 GB per-layer limit.
	ghcrMaxLayer = int64(10_000_000_000)
)

// The table from the plan, which is the whole argument for deriving the size:
// the request cost stops growing with the artifact instead of climbing past the
// ceiling that refused the first real push.
func TestPlanLayerSize(t *testing.T) {
	tests := []struct {
		why        string
		size       int64
		max        int64
		wantSize   int64
		wantReason layerSizeReason
		wantLayers int
	}{
		// Below 48 GiB the derivation would ask for less than the floor, so the
		// floor answers — which covers every artifact published today.
		{"a small artifact is one layer", 4 * gib, 0, 2 * gib, layerSizeFloor, 2},
		{"the 20 GiB fixture", 21474840576, 0, 2 * gib, layerSizeFloor, 11},
		{"a STAR human genome index", 26_000_000_000, 0, 2 * gib, layerSizeFloor, 13},
		{"the last size the floor covers", 48 * gib, 0, 2 * gib, layerSizeFloor, 24},

		// Above it the size follows the artifact, in whole GiB steps.
		{"one byte past the floor's reach", 48*gib + 1, 0, 3 * gib, layerSizeDerived, 17},
		{"60 GiB", 60 * gib, 0, 3 * gib, layerSizeDerived, 20},
		{"100 GiB", 100 * gib, 0, 5 * gib, layerSizeDerived, 20},
		{"200 GiB", 200 * gib, 0, 9 * gib, layerSizeDerived, 23},

		// A hard per-layer limit outranks the derivation, and 10 GB is not a
		// whole number of GiB, so it rounds down to 9 rather than up to 10.
		{"clamped at GHCR's per-layer limit", 400 * gib, ghcrMaxLayer, 9 * gib, layerSizeClamped, 45},
		{"a limit below the floor still wins", 4 * gib, gib, gib, layerSizeClamped, 4},
	}

	for _, tt := range tests {
		size, reason := planLayerSize(tt.size, tt.max)
		if size != tt.wantSize || reason != tt.wantReason {
			t.Errorf("%s: planLayerSize(%d, %d) = (%d, %s), want (%d, %s)",
				tt.why, tt.size, tt.max, size, reason, tt.wantSize, tt.wantReason)
			continue
		}
		if got := layerCount(tt.size, size); got != tt.wantLayers {
			t.Errorf("%s: layerCount = %d, want %d", tt.why, got, tt.wantLayers)
		}
	}
}

// The properties that hold for every artifact, checked across the whole range
// rather than at the few points the table names.
func TestPlanLayerSizeInvariants(t *testing.T) {
	for size := gib; size <= 512*gib; size += 7 * gib {
		got, _ := planLayerSize(size, ghcrMaxLayer)
		switch {
		case got%gib != 0:
			t.Fatalf("size %d: layer size %d is not a whole number of GiB", size, got)
		case got < minLayerSize:
			t.Fatalf("size %d: layer size %d is below the floor", size, got)
		case got > ghcrMaxLayer:
			t.Fatalf("size %d: layer size %d exceeds the hard limit", size, got)
		}
	}
}

// The reason the request cost is worth deriving: at a fixed size it climbs past
// the ceiling that refused the first push, and here it does not.
func TestPlanLayerSizeHoldsTheRequestCountDown(t *testing.T) {
	// Three requests per fresh layer, plus roughly ten fixed for the token,
	// tag probe, empty config blob, manifest, and tags.
	const fixedRequests = 10
	requests := func(size int64) int {
		layerSize, _ := planLayerSize(size, ghcrMaxLayer)
		return 3*layerCount(size, layerSize) + fixedRequests
	}

	// The push that failed: 41 layers of 512 MiB was about 123 requests, and
	// GHCR refused chunk 34.
	if got := 3*41 + fixedRequests; got <= 100 {
		t.Fatalf("the fixture's old plan cost %d requests, which would not have failed", got)
	}
	for _, size := range []int64{26_000_000_000, 60 * gib, 100 * gib, 200 * gib, 270 * gib} {
		if got := requests(size); got > 100 {
			t.Errorf("a %d-byte artifact plans %d requests, past the observed ceiling", size, got)
		}
	}
}

// A resumed push has to find the layers the first attempt committed, so the plan
// may not depend on anything but the artifact and its destination.
func TestPlanLayerSizeIsDeterministic(t *testing.T) {
	for range 100 {
		if got, _ := planLayerSize(137*gib, ghcrMaxLayer); got != 6*gib {
			t.Fatalf("planLayerSize = %d, want a stable 6 GiB", got)
		}
	}
}

// An empty artifact is one empty layer, not none: every artifact has a payload
// descriptor, so pull needs no special case for zero.
func TestLayerCountFloorsAtOne(t *testing.T) {
	for _, size := range []int64{0, -1} {
		if got := layerCount(size, 2*gib); got != 1 {
			t.Errorf("layerCount(%d) = %d, want 1", size, got)
		}
	}
	if got := layerCount(gib, 0); got != 1 {
		t.Errorf("layerCount with no layer size = %d, want 1", got)
	}
}

// A wrapped layer count is a plan that looks reasonable and is not, so the
// rounding saturates instead.
func TestLayerPlanArithmeticDoesNotOverflow(t *testing.T) {
	const maxInt64 = int64(^uint64(0) >> 1)
	if got := ceilDiv(maxInt64, 1); got != maxInt64 {
		t.Errorf("ceilDiv(max, 1) = %d", got)
	}
	if got := ceilTo(maxInt64, gib); got != maxInt64 {
		t.Errorf("ceilTo(max, 1GiB) = %d, want it saturated rather than wrapped", got)
	}
	if got, _ := planLayerSize(maxInt64, 0); got <= 0 {
		t.Errorf("planLayerSize(max) = %d, want a positive size", got)
	}
}
