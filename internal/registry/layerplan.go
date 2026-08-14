package registry

// How an artifact is cut into layers.
//
// The binding limit is a count of requests, not of bytes: a fresh layer costs
// three (HEAD to check presence, POST to open the session, PUT to commit), and
// GHCR refused a 20 GiB push at 512 MiB layers on chunk 34 of 41 — its ~101st
// request — with a secondary rate limit. So the layer *count* is held roughly
// constant and the size follows from the artifact; a fixed size makes the count
// grow with the artifact, turning every choice into a supported-size cliff.
//
// Nothing has to agree with the size chosen here: a pull takes the boundary from
// the manifest, so it is a per-push decision recorded nowhere.
const (
	// targetLayers is about 80% of what the observed ceiling allows, the margin
	// absorbing that the ceiling is one measurement and cannot be queried.
	targetLayers = 24
	// minLayerSize is the floor. Below it there is nothing left to optimize: a
	// 20 GiB push already costs 33 requests against the 123 that failed.
	minLayerSize = 2 << 30
	// layerGranularity keeps sizes to whole GiB so artifacts of similar size
	// share boundaries and deduplicate. A continuous size/targetLayers would give
	// two artifacts one byte apart entirely different boundaries.
	layerGranularity = 1 << 30
)

// layerSizeReason names why a layer size was chosen, so the upload plan can
// report it.
type layerSizeReason string

const (
	layerSizeFloor   layerSizeReason = "floor"
	layerSizeDerived layerSizeReason = "derived"
	layerSizeClamped layerSizeReason = "clamped"
)

// planLayerSize returns the layer size for an artifact of size bytes, and why.
// maxLayerSize is the destination's hard per-layer limit, or zero when none is
// known. There is no override, and deliberately no setting for one.
//
// Pure in its arguments — no clock, no measured throughput — so the same
// artifact plans identical boundaries twice, which is what lets a resumed push
// find the layers the first one committed.
func planLayerSize(size, maxLayerSize int64) (int64, layerSizeReason) {
	chosen, reason := ceilTo(ceilDiv(size, targetLayers), layerGranularity), layerSizeDerived
	if chosen <= minLayerSize {
		chosen, reason = minLayerSize, layerSizeFloor
	}
	// A hard limit outranks our floor: the registry will refuse what it refuses.
	if maxLayerSize > 0 && chosen > maxLayerSize {
		chosen, reason = floorTo(maxLayerSize, layerGranularity), layerSizeClamped
		if chosen == 0 {
			chosen = maxLayerSize
		}
	}
	return chosen, reason
}

// layerCount reports how many layers an artifact of size bytes cuts into. An
// empty artifact is one empty layer, so every artifact has a payload descriptor
// and pull needs no zero case.
func layerCount(size, layerSize int64) int {
	if size <= 0 || layerSize <= 0 {
		return 1
	}
	return int(ceilDiv(size, layerSize))
}

// ceilDiv divides, rounding up. Saturating rather than overflowing: a+b-1 wraps
// for a size near the int64 maximum, and a wrapped layer count is a plan that
// looks reasonable and is not.
func ceilDiv(a, b int64) int64 {
	if a <= 0 || b <= 0 {
		return 0
	}
	q := a / b
	if a%b != 0 {
		q++
	}
	return q
}

// ceilTo rounds up to a multiple of step, saturating rather than wrapping.
func ceilTo(v, step int64) int64 {
	if v <= 0 || step <= 0 {
		return 0
	}
	n := ceilDiv(v, step)
	if n > (1<<62)/step {
		return v
	}
	return n * step
}

// floorTo rounds down to a multiple of step, which may be zero when v is below
// one step. Callers decide what that means.
func floorTo(v, step int64) int64 {
	if v <= 0 || step <= 0 {
		return 0
	}
	return v / step * step
}
