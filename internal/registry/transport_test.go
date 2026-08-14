package registry

import (
	"context"
	"encoding/json"
	"io"
	"net/http"
	"net/http/httptest"
	"strings"
	"testing"
	"time"

	"oras.land/oras-go/v2/registry/remote/errcode"
)

// roundTrip drives one request through the inspecting transport against srv,
// returning what the transport captured.
func roundTrip(t *testing.T, srv *httptest.Server, req *http.Request) (*http.Response, *transferFailure) {
	t.Helper()
	slot := &failureSlot{}
	transport := &inspectTransport{next: http.DefaultTransport}

	resp, err := transport.RoundTrip(req.WithContext(withFailureSlot(req.Context(), slot)))
	if err != nil {
		t.Fatalf("RoundTrip: %v", err)
	}
	t.Cleanup(func() { resp.Body.Close() }) //nolint:errcheck
	return resp, slot.get()
}

func newRequest(t *testing.T, method, url string, body io.Reader) *http.Request {
	t.Helper()
	req, err := http.NewRequestWithContext(context.Background(), method, url, body)
	if err != nil {
		t.Fatal(err)
	}
	return req
}

// The headers ORAS throws away. errcode.ErrorResponse keeps method, URL, status,
// and the parsed error document — no headers at all — so without capturing them
// here every wait falls back to the local schedule.
func TestTransportPreservesRateLimitHeaders(t *testing.T) {
	srv := httptest.NewServer(http.HandlerFunc(func(w http.ResponseWriter, _ *http.Request) {
		w.Header().Set("Retry-After", "120")
		w.Header().Set("X-RateLimit-Remaining", "0")
		w.Header().Set("X-RateLimit-Reset", "1700000000")
		w.WriteHeader(http.StatusTooManyRequests)
	}))
	defer srv.Close()

	_, failure := roundTrip(t, srv, newRequest(t, http.MethodGet, srv.URL+"/v2/lab/cnt/blobs/x", nil))
	if failure == nil {
		t.Fatal("nothing was captured from a 429")
	}
	if failure.RetryAfter != "120" || failure.Remaining != "0" || failure.Reset != "1700000000" {
		t.Errorf("captured %+v, want the three headers intact", failure)
	}
	if wait, ok := failure.delay(time.Now()); !ok || wait != 2*time.Minute {
		t.Errorf("delay = (%v, %v), want two minutes from Retry-After", wait, ok)
	}
}

// The body has to survive being read, or a specific registry error becomes
// "Forbidden": ORAS decodes the document itself and a consumed body decodes to
// nothing.
func TestTransportRestoresTheBodyForORAS(t *testing.T) {
	const message = "You have exceeded a secondary rate limit."
	srv := httptest.NewServer(http.HandlerFunc(func(w http.ResponseWriter, _ *http.Request) {
		w.WriteHeader(http.StatusForbidden)
		_ = json.NewEncoder(w).Encode(map[string]any{
			"errors": []map[string]string{{"code": "DENIED", "message": message}},
		})
	}))
	defer srv.Close()

	resp, _ := roundTrip(t, srv, newRequest(t, http.MethodGet, srv.URL+"/v2/lab/cnt/manifests/1.0", nil))

	// Exactly what ORAS does with it, so this proves the real consumer works.
	var doc struct{ Errors errcode.Errors }
	if err := json.NewDecoder(resp.Body).Decode(&doc); err != nil {
		t.Fatalf("the body ORAS would parse is gone: %v", err)
	}
	if len(doc.Errors) != 1 || doc.Errors[0].Message != message {
		t.Errorf("parsed %+v, want the secondary-limit message intact", doc.Errors)
	}
}

// A body larger than the captured prefix must still arrive whole: the prefix is
// put back in front of the rest, not in place of it.
func TestTransportRestoresABodyLargerThanTheCapture(t *testing.T) {
	body := strings.Repeat("x", maxCapturedBody*3)
	srv := httptest.NewServer(http.HandlerFunc(func(w http.ResponseWriter, _ *http.Request) {
		w.WriteHeader(http.StatusInternalServerError)
		_, _ = io.WriteString(w, body)
	}))
	defer srv.Close()

	resp, _ := roundTrip(t, srv, newRequest(t, http.MethodGet, srv.URL+"/v2/", nil))
	got, err := io.ReadAll(resp.Body)
	if err != nil {
		t.Fatal(err)
	}
	if string(got) != body {
		t.Errorf("body is %d bytes, want all %d", len(got), len(body))
	}
}

// The single largest reduction in wasted transfer: a rejected layer is rejected
// before its bytes move, not after.
func TestTransportAsksPermissionForALargeBody(t *testing.T) {
	var sawExpect, bodyRead bool
	srv := httptest.NewServer(http.HandlerFunc(func(w http.ResponseWriter, r *http.Request) {
		sawExpect = r.Header.Get("Expect") == "100-continue"
		n, _ := io.Copy(io.Discard, r.Body)
		bodyRead = n > 0
		w.WriteHeader(http.StatusCreated)
	}))
	defer srv.Close()

	req := newRequest(t, http.MethodPut, srv.URL+"/v2/lab/cnt/blobs/upload/7", strings.NewReader(strings.Repeat("p", expectContinueMinSize)))
	req.ContentLength = expectContinueMinSize
	roundTrip(t, srv, req)

	if !sawExpect {
		t.Error("a large blob upload did not ask permission before sending")
	}
	if !bodyRead {
		t.Error("the body never arrived")
	}
}

// Below the threshold the extra round trip costs more than it saves, and a
// manifest PUT is small.
func TestTransportDoesNotAskPermissionForASmallBody(t *testing.T) {
	var sawExpect bool
	srv := httptest.NewServer(http.HandlerFunc(func(w http.ResponseWriter, r *http.Request) {
		sawExpect = r.Header.Get("Expect") != ""
		w.WriteHeader(http.StatusCreated)
	}))
	defer srv.Close()

	req := newRequest(t, http.MethodPut, srv.URL+"/v2/lab/cnt/manifests/1.0", strings.NewReader("{}"))
	req.ContentLength = 2
	roundTrip(t, srv, req)

	if sawExpect {
		t.Error("a small manifest PUT paid for a 100-continue round trip")
	}
}

// ORAS never resumes an upload session, so a failed finalize leaves one that
// nothing will ever complete. Four retries would leave four.
func TestTransportCancelsAnAbandonedUploadSession(t *testing.T) {
	var deleted, deletedAuth string
	srv := httptest.NewServer(http.HandlerFunc(func(w http.ResponseWriter, r *http.Request) {
		if r.Method == http.MethodDelete {
			deleted, deletedAuth = r.URL.Path, r.Header.Get("Authorization")
			w.WriteHeader(http.StatusAccepted)
			return
		}
		w.WriteHeader(http.StatusForbidden)
	}))
	defer srv.Close()

	req := newRequest(t, http.MethodPut, srv.URL+"/v2/lab/cnt/blobs/upload/7.abc?digest=sha256%3Adead", nil)
	req.Header.Set("Authorization", "Bearer token-from-the-post")
	_, failure := roundTrip(t, srv, req)

	if failure == nil || failure.cancel == nil {
		t.Fatal("a failed blob finalize offered no way to abandon its session")
	}
	failure.cancel(context.Background())

	if deleted != "/v2/lab/cnt/blobs/upload/7.abc" {
		t.Errorf("cancelled %q, want the session without its digest query", deleted)
	}
	// The credential lives on the failed request and nowhere else by this point.
	if deletedAuth != "Bearer token-from-the-post" {
		t.Errorf("cancel sent %q, want the failed request's credential", deletedAuth)
	}
}

// Only a blob finalize has a session to abandon. A DELETE against a manifest or
// a GET would be a different request entirely.
func TestOnlyABlobFinalizeHasASessionToCancel(t *testing.T) {
	for _, tt := range []struct {
		method, url string
		want        bool
	}{
		{http.MethodPut, "https://ghcr.io/v2/lab/cnt/blobs/upload/7?digest=sha256%3Ax", true},
		{http.MethodPut, "https://ghcr.io/v2/lab/cnt/manifests/1.0", false},
		{http.MethodGet, "https://ghcr.io/v2/lab/cnt/blobs/upload/7", false},
		{http.MethodPost, "https://ghcr.io/v2/lab/cnt/blobs/uploads/", false},
	} {
		got := uploadSessionURL(newRequest(t, tt.method, tt.url, nil))
		if (got != "") != tt.want {
			t.Errorf("uploadSessionURL(%s %s) = %q, want a session: %v", tt.method, tt.url, got, tt.want)
		}
	}
}

// A success is left completely alone: no capture, no body rewrapping.
func TestTransportIgnoresSuccess(t *testing.T) {
	srv := httptest.NewServer(http.HandlerFunc(func(w http.ResponseWriter, _ *http.Request) {
		_, _ = io.WriteString(w, "ok")
	}))
	defer srv.Close()

	resp, failure := roundTrip(t, srv, newRequest(t, http.MethodGet, srv.URL+"/v2/", nil))
	if failure != nil {
		t.Errorf("a 200 was captured as a failure: %+v", failure)
	}
	body, _ := io.ReadAll(resp.Body)
	if string(body) != "ok" {
		t.Errorf("body = %q", body)
	}
}

// The whole chain, end to end: a registry sends Retry-After, ORAS builds an
// error that cannot carry it, and the wait still comes from the header rather
// than from the local schedule.
//
// Every other test here covers one link. This is the one that would catch the
// slot being scoped wrong, the transport not being installed in the auth client,
// or the annotation happening after the delay is chosen.
func TestRetryHonoursRetryAfterThroughORAS(t *testing.T) {
	var puts int
	srv := httptest.NewServer(http.HandlerFunc(func(w http.ResponseWriter, r *http.Request) {
		switch {
		case r.Method == http.MethodPost:
			w.Header().Set("Location", "/v2/lab/cnt/blobs/upload/session-1")
			w.WriteHeader(http.StatusAccepted)
		case r.Method == http.MethodDelete:
			w.WriteHeader(http.StatusAccepted)
		case r.Method == http.MethodPut:
			puts++
			if puts == 1 {
				_, _ = io.Copy(io.Discard, r.Body)
				w.Header().Set("Retry-After", "45")
				w.WriteHeader(http.StatusTooManyRequests)
				_, _ = io.WriteString(w, `{"errors":[{"code":"TOOMANYREQUESTS","message":"slow down"}]}`)
				return
			}
			_, _ = io.Copy(io.Discard, r.Body)
			w.WriteHeader(http.StatusCreated)
		default: // HEAD for the presence probe
			w.WriteHeader(http.StatusNotFound)
		}
	}))
	defer srv.Close()

	repository, err := newRepository(strings.TrimPrefix(srv.URL, "http://"), "lab/cnt")
	if err != nil {
		t.Fatal(err)
	}

	path := writeArtifact(t, t.TempDir(), "sample.sqf", "payload")
	policy := newTestPolicy()
	ctx := withRetryPolicy(context.Background(), policy.retryPolicy)

	if _, err := pushArtifactLayers(ctx, repository.Blobs(), path, MediaTypeOverlayBlob, 1<<20); err != nil {
		t.Fatalf("pushArtifactLayers: %v", err)
	}
	if puts != 2 {
		t.Errorf("PUT attempts = %d, want the refused one retried", puts)
	}
	if want := []time.Duration{45 * time.Second}; !equalDurations(*policy.waits, want) {
		t.Errorf("waits = %v, want %v — the header did not survive to the policy", *policy.waits, want)
	}
}

func TestTransferFailureDelay(t *testing.T) {
	now := time.Unix(1_700_000_000, 0)
	tests := []struct {
		why      string
		failure  transferFailure
		wantWait time.Duration
		wantOK   bool
	}{
		{"seconds", transferFailure{RetryAfter: "90"}, 90 * time.Second, true},
		{"an HTTP date", transferFailure{RetryAfter: now.Add(time.Minute).UTC().Format(http.TimeFormat)}, time.Minute, true},
		{"a date already past", transferFailure{RetryAfter: now.Add(-time.Hour).UTC().Format(http.TimeFormat)}, 0, true},
		{
			"a reset, once nothing remains",
			transferFailure{Remaining: "0", Reset: "1700000600"}, 10 * time.Minute, true,
		},
		// A reset on its own says when a window rolls over, not that anything is
		// being refused now.
		{"a reset with requests still left", transferFailure{Remaining: "42", Reset: "1700000600"}, 0, false},
		{"a reset with no remaining count", transferFailure{Reset: "1700000600"}, 0, false},
		{"nothing usable", transferFailure{}, 0, false},
		{"an unparseable Retry-After", transferFailure{RetryAfter: "soon"}, 0, false},
	}
	for _, tt := range tests {
		wait, ok := tt.failure.delay(now)
		if ok != tt.wantOK || wait != tt.wantWait {
			t.Errorf("%s: delay = (%v, %v), want (%v, %v)", tt.why, wait, ok, tt.wantWait, tt.wantOK)
		}
	}
}

// Providers key abuse heuristics and support triage on the User-Agent, and
// "oras-go" is not a useful answer when reporting a secondary rate limit.
func TestUserAgentIdentifiesCondaTainer(t *testing.T) {
	if got := userAgent(); !strings.HasPrefix(got, "condatainer/") {
		t.Errorf("userAgent() = %q", got)
	}
	var sent string
	srv := httptest.NewServer(http.HandlerFunc(func(w http.ResponseWriter, r *http.Request) {
		sent = r.UserAgent()
		w.WriteHeader(http.StatusOK)
	}))
	defer srv.Close()

	client := newAuthClient()
	resp, err := client.Do(newRequest(t, http.MethodGet, srv.URL+"/v2/", nil))
	if err != nil {
		t.Fatal(err)
	}
	defer resp.Body.Close() //nolint:errcheck
	if sent != userAgent() {
		t.Errorf("sent User-Agent %q, want %q", sent, userAgent())
	}
}
