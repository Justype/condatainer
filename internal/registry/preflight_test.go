package registry

import (
	"context"
	"encoding/json"
	"errors"
	"log/slog"
	"net/http"
	"net/http/httptest"
	"os"
	"path/filepath"
	"strings"
	"testing"
	"time"

	"github.com/Justype/condatainer/internal/logging"
)

// The probe's retryable cases are only observable once the budget cuts them off,
// so every test that exercises one pays it in wall time. The bound under test
// holds at any budget, and every server here is a local httptest, so the suite
// runs against a short one.
func TestMain(m *testing.M) {
	probeTimeout = 200 * time.Millisecond
	os.Exit(m.Run())
}

func preflightCtx(lines *[]string) context.Context {
	return logging.WithLogger(context.Background(), slog.New(recordingHandler{lines: lines}))
}

// The summary is the one place a publisher learns what the push decided, so it
// has to name the size, why that size, the count, and the expected cost.
func TestPreflightReportsThePlan(t *testing.T) {
	path := filepath.Join(t.TempDir(), "grch38-genome--star.sqf")
	if err := os.WriteFile(path, []byte("payload"), 0o600); err != nil {
		t.Fatal(err)
	}

	var lines []string
	plan, err := preflightUpload(preflightCtx(&lines), path, "ghcr.io/lab/cnt", []string{"star-2.7.11b"}, profileFor("ghcr.io"))
	if err != nil {
		t.Fatalf("preflightUpload: %v", err)
	}
	if plan.LayerSize != minLayerSize || plan.Reason != layerSizeFloor || plan.Layers != 1 {
		t.Errorf("plan = %+v, want one layer at the floor", plan)
	}
	if want := estimateRequests(1); plan.Requests != want {
		t.Errorf("Requests = %d, want %d", plan.Requests, want)
	}

	if len(lines) != 1 {
		t.Fatalf("logged %v, want one summary", lines)
	}
	for _, want := range []string{
		"upload plan", "layer-size=2.00 GB (floor)", "layers=1", "requests=~15", "registry=ghcr.io",
	} {
		if !strings.Contains(lines[0], want) {
			t.Errorf("summary %q does not mention %q", lines[0], want)
		}
	}
}

// A push past the observed ceiling is still correct — retry carries it — so it
// is a warning about what to expect, not a refusal. An operator told a long
// transfer will pause does not report it as a hang.
func TestPreflightWarnsPastTheObservedCeiling(t *testing.T) {
	// 400 GiB clamps to 9 GiB layers, 45 of them, 145 requests.
	if _, reason := planLayerSize(400*gib, ghcrMaxLayer); reason != layerSizeClamped {
		t.Fatalf("the fixture no longer clamps; the request estimate below is wrong")
	}
	layers := layerCount(400*gib, 9*gib)
	if got := estimateRequests(layers); got <= expectedRequestCeiling {
		t.Errorf("a 400 GiB push plans %d requests, which would not warn", got)
	}
}

// Every check here exists to fail before gigabytes have been digested and sent.
func TestPreflightRejectsBeforeAnyUpload(t *testing.T) {
	dir := t.TempDir()
	empty := filepath.Join(dir, "empty.sqf")
	if err := os.WriteFile(empty, nil, 0o600); err != nil {
		t.Fatal(err)
	}

	var lines []string
	ctx := preflightCtx(&lines)

	if _, err := preflightUpload(ctx, empty, "ghcr.io/lab/cnt", []string{"1.0"}, transferProfile{}); err == nil {
		t.Error("an empty artifact was accepted")
	}
	if _, err := preflightUpload(ctx, dir, "ghcr.io/lab/cnt", []string{"1.0"}, transferProfile{}); err == nil {
		t.Error("a directory was accepted as an artifact")
	}
	if _, err := preflightUpload(ctx, filepath.Join(dir, "absent.sqf"), "ghcr.io/lab/cnt", []string{"1.0"}, transferProfile{}); err == nil {
		t.Error("a missing artifact was accepted")
	}
	if len(lines) != 0 {
		t.Errorf("a refused push still logged a plan: %v", lines)
	}
}

// The manifest ceiling is the real bound on how small a layer may be, derived
// from a spec limit rather than a guess. It is why no lower bound is imposed
// anywhere else.
func TestManifestSizePreflight(t *testing.T) {
	if err := checkManifestFits(maxManifestBytes/descriptorBytes - 1); err != nil {
		t.Errorf("a manifest inside the ceiling was refused: %v", err)
	}
	err := checkManifestFits(maxManifestBytes/descriptorBytes + 1)
	if err == nil {
		t.Fatal("a manifest past the 4 MB ceiling was accepted")
	}
	if !strings.Contains(err.Error(), "layers") {
		t.Errorf("the refusal does not say what to change: %v", err)
	}
	// It cannot fire for any size the derivation produces: at the 2 GiB floor
	// this would need a 32 TiB artifact, and above the floor the count is held
	// near targetLayers. Checked rather than assumed.
	for _, size := range []int64{gib, 26_000_000_000, 200 * gib, 400 * gib} {
		layerSize, _ := planLayerSize(size, ghcrMaxLayer)
		if err := checkManifestFits(layerCount(size, layerSize)); err != nil {
			t.Errorf("a %d-byte artifact tripped the manifest ceiling: %v", size, err)
		}
	}
}

// A name checked at the first layer would miss the overrun: the chunk suffix
// lands on every title equally, so the last one is the longest.
func TestNamePreflight(t *testing.T) {
	long := strings.Repeat("a", 250) + ".sqf"
	if err := checkNames("ghcr.io/lab/cnt", "/images/"+long, []string{"1.0"}, 11); err == nil {
		t.Error("a layer title past the reference limit was accepted")
	}
	if err := checkNames("ghcr.io/lab/cnt", "/images/x.sqf", []string{strings.Repeat("v", 129)}, 1); err == nil {
		t.Error("a tag past the OCI limit of 128 was accepted")
	}
	if err := checkNames("ghcr.io/"+strings.Repeat("d", 250), "/images/x.sqf", []string{"1.0"}, 1); err == nil {
		t.Error("a reference past the length clients assume was accepted")
	}
	if err := checkNames("ghcr.io/lab/cnt", "/images/grch38-genome--star.sqf", []string{"star-2.7.11b", "latest"}, 13); err != nil {
		t.Errorf("an ordinary push was refused: %v", err)
	}
}

// probeAgainst runs the push preflight against a handler standing in for a
// registry, and reports what the registry saw.
func probeAgainst(t *testing.T, handler http.HandlerFunc) (*[]string, error) {
	t.Helper()
	seen := &[]string{}
	srv := httptest.NewServer(http.HandlerFunc(func(w http.ResponseWriter, r *http.Request) {
		*seen = append(*seen, r.Method+" "+r.URL.Path)
		handler(w, r)
	}))
	t.Cleanup(srv.Close)

	repository, err := newRepository(strings.TrimPrefix(srv.URL, "http://"), "lab/cnt")
	if err != nil {
		t.Fatal(err)
	}
	return seen, probePushAccess(context.Background(), repository)
}

// Two requests that answer what would otherwise stay unanswered until the first
// layer had been digested and sent.
func TestProbeOpensAndAbandonsASession(t *testing.T) {
	seen, err := probeAgainst(t, func(w http.ResponseWriter, r *http.Request) {
		if r.Method == http.MethodPost {
			w.Header().Set("Location", "/v2/lab/cnt/blobs/upload/probe-1")
			w.WriteHeader(http.StatusAccepted)
			return
		}
		w.WriteHeader(http.StatusAccepted)
	})
	if err != nil {
		t.Fatalf("probePushAccess: %v", err)
	}
	want := []string{"POST /v2/lab/cnt/blobs/uploads/", "DELETE /v2/lab/cnt/blobs/upload/probe-1"}
	if len(*seen) != 2 || (*seen)[0] != want[0] || (*seen)[1] != want[1] {
		t.Errorf("registry saw %v, want %v", *seen, want)
	}
}

// The two failures worth spending a round trip to find early. Both would
// otherwise surface after the first layer was hashed and sent.
func TestProbeReportsWhatTheRegistryWasClearAbout(t *testing.T) {
	tests := []struct {
		why      string
		status   int
		code     string
		want     error
		mentions string
	}{
		{"no push scope", http.StatusForbidden, "DENIED", ErrUnauthorized, "credential"},
		{"the repository does not exist", http.StatusNotFound, "NAME_UNKNOWN", ErrNotFound, "created before a push"},
		{"an unsupported manifest shape", http.StatusBadRequest, "UNSUPPORTED", ErrIncompatibleRegistry, ""},
	}
	for _, tt := range tests {
		_, err := probeAgainst(t, func(w http.ResponseWriter, _ *http.Request) {
			w.WriteHeader(tt.status)
			_ = json.NewEncoder(w).Encode(map[string]any{
				"errors": []map[string]string{{"code": tt.code, "message": "no"}},
			})
		})
		if !errors.Is(err, tt.want) {
			t.Errorf("%s: err = %v, want %v", tt.why, err, tt.want)
			continue
		}
		if tt.mentions != "" && !strings.Contains(err.Error(), tt.mentions) {
			t.Errorf("%s: %v does not say what to do about it", tt.why, err)
		}
	}
}

// A preflight that invents failure modes is worse than no preflight. Anything
// the registry was not clear about belongs to the push, which has retry and
// better messages.
func TestProbeDoesNotBlockOnAnInconclusiveAnswer(t *testing.T) {
	t.Parallel() // both wait out probeTimeout; overlap them
	// In parallel because the retryable statuses are retried by ORAS's transport
	// until probeTimeout cuts them off, which is the behaviour under test and
	// also the whole budget each.
	for _, status := range []int{
		http.StatusInternalServerError,
		http.StatusTooManyRequests,
		http.StatusMethodNotAllowed,
		http.StatusOK, // not the 202 the spec calls for, but not a refusal either
	} {
		t.Run(http.StatusText(status), func(t *testing.T) {
			t.Parallel()
			_, err := probeAgainst(t, func(w http.ResponseWriter, _ *http.Request) {
				w.WriteHeader(status)
			})
			if err != nil {
				t.Errorf("a %d from the probe blocked the push: %v", status, err)
			}
		})
	}
}

// The bound matters as much as the outcome: a preflight is only cheap if it
// cannot sit there being retried.
func TestProbeIsBounded(t *testing.T) {
	t.Parallel() // both wait out probeTimeout; overlap them
	start := time.Now()
	_, err := probeAgainst(t, func(w http.ResponseWriter, _ *http.Request) {
		w.WriteHeader(http.StatusInternalServerError)
	})
	if err != nil {
		t.Errorf("err = %v, want an inconclusive probe to be ignored", err)
	}
	if elapsed := time.Since(start); elapsed > 2*probeTimeout {
		t.Errorf("the probe took %v against a %v budget", elapsed, probeTimeout)
	}
}

// A registry that cannot be reached is the push's problem to report, not the
// preflight's to guess at.
func TestProbeDoesNotBlockOnAnUnreachableRegistry(t *testing.T) {
	repository, err := newRepository("127.0.0.1:1/lab", "cnt")
	if err != nil {
		t.Fatal(err)
	}
	if err := probePushAccess(context.Background(), repository); err != nil {
		t.Errorf("an unreachable registry blocked the push at preflight: %v", err)
	}
}

// Failing to tidy the probe's own session is not a reason to refuse a push that
// is otherwise fine.
func TestProbeSurvivesAFailedCancel(t *testing.T) {
	seen, err := probeAgainst(t, func(w http.ResponseWriter, r *http.Request) {
		if r.Method == http.MethodPost {
			w.Header().Set("Location", "/v2/lab/cnt/blobs/upload/probe-1")
			w.WriteHeader(http.StatusAccepted)
			return
		}
		w.WriteHeader(http.StatusNotFound)
	})
	if err != nil {
		t.Errorf("a failed cancel blocked the push: %v", err)
	}
	if len(*seen) != 2 {
		t.Errorf("registry saw %v, want the cancel to have been attempted", *seen)
	}
}

func TestRegistryHost(t *testing.T) {
	for in, want := range map[string]string{
		"ghcr.io/lab/cnt":       "ghcr.io",
		"oci://ghcr.io/lab/cnt": "ghcr.io",
		"localhost:5000/lab":    "localhost:5000",
		"ghcr.io":               "ghcr.io",
	} {
		if got := registryHost(in); got != want {
			t.Errorf("registryHost(%q) = %q, want %q", in, got, want)
		}
	}
}
