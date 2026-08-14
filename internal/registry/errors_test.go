package registry

import (
	"context"
	"errors"
	"fmt"
	"net"
	"net/url"
	"testing"

	"oras.land/oras-go/v2/errdef"
	"oras.land/oras-go/v2/registry/remote/errcode"
)

func responseErr(status int) error {
	return &errcode.ErrorResponse{
		Method:     "GET",
		URL:        &url.URL{Scheme: "https", Host: "ghcr.io", Path: "/v2/lab/cnt/manifests/1.0"},
		StatusCode: status,
	}
}

func TestClassify(t *testing.T) {
	tests := []struct {
		why  string
		in   error
		want error // nil means "passed through unclassified"
	}{
		{"nothing to classify", nil, nil},
		{"oras says not found", errdef.ErrNotFound, ErrNotFound},
		{"a wrapped not found", fmt.Errorf("resolving: %w", errdef.ErrNotFound), ErrNotFound},
		{"404", responseErr(404), ErrNotFound},
		{"401", responseErr(401), ErrUnauthorized},
		{"403", responseErr(403), ErrUnauthorized},
		{"429", responseErr(429), ErrRateLimited},
		{"500", responseErr(500), ErrUnavailable},
		{"503", responseErr(503), ErrUnavailable},
		{"a dial failure", &net.DNSError{Err: "no such host", Name: "ghcr.io"}, ErrUnavailable},
		{"a timeout", context.DeadlineExceeded, ErrUnavailable},
		// The dangerous direction: anything unrecognized must not become the
		// category that authorizes a silent rebuild.
		{"a 400 nobody anticipated", responseErr(400), nil},
		{"an error from somewhere else", errors.New("mktemp: no space left"), nil},
	}

	for _, tt := range tests {
		got := classify(tt.in)
		if tt.want == nil {
			if got != nil && !errors.Is(got, tt.in) {
				t.Errorf("%s: classify(%v) = %v, want it passed through", tt.why, tt.in, got)
			}
			for _, category := range []error{ErrNotFound, ErrUnauthorized, ErrUnavailable, ErrRateLimited} {
				if errors.Is(got, category) {
					t.Errorf("%s: classify(%v) claimed %v", tt.why, tt.in, category)
				}
			}
			continue
		}
		if !errors.Is(got, tt.want) {
			t.Errorf("%s: classify(%v) = %v, want %v", tt.why, tt.in, got, tt.want)
		}
		// The category is added, never substituted: the message still has to say
		// what actually happened.
		if !errors.Is(got, tt.in) {
			t.Errorf("%s: classify dropped the original error: %v", tt.why, got)
		}
	}
}

// A cancelled context is the user pressing Ctrl-C. Reporting it as an unreachable
// registry would turn an interrupted pull into a forty-minute local rebuild.
func TestClassifyLeavesCancellationAlone(t *testing.T) {
	err := classify(fmt.Errorf("pulling: %w", context.Canceled))
	if errors.Is(err, ErrUnavailable) {
		t.Errorf("a cancelled pull was reported as an unavailable registry: %v", err)
	}
	if !errors.Is(err, context.Canceled) {
		t.Errorf("err = %v, want the cancellation intact", err)
	}
}

// codedErr builds the error shape ORAS produces from a registry's OCI error
// document: an ErrorResponse whose Unwrap exposes the coded errors.
func codedErr(status int, code, message string) error {
	return &errcode.ErrorResponse{
		Method:     "PUT",
		URL:        &url.URL{Scheme: "https", Host: "ghcr.io", Path: "/v2/lab/cnt/blobs/upload/7.8f99"},
		StatusCode: status,
		Errors:     errcode.Errors{{Code: code, Message: message}},
	}
}

// ghcrSecondaryLimitMessage is the message GHCR actually returned when a 20 GiB
// push in 41 layers of 512 MiB was refused on chunk 34. Recorded verbatim: the
// whole difficulty is that it arrives as a plain DENIED, so a paraphrase would
// test the paraphrase rather than the thing that broke.
const ghcrSecondaryLimitMessage = `permission_denied: Error from intermediary with HTTP status code 403 "Forbidden" - with-body: {
  "documentation_url": "https://docs.github.com/free-pro-team@latest/rest/overview/rate-limits-for-the-rest-api#about-secondary-rate-limits",
  "message": "You have exceeded a secondary rate limit. Please wait a few minutes before you try again. For more on scraping GitHub and how it may affect your rights, please review our Terms of Service (https://docs.github.com/en/site-policy/github-terms/github-terms-of-service) If you reach out to GitHub Support for help, please include the request ID D8B8:8C645:1EB456:232FB5:6A7DD4EF."
}`

// The regression test for the failure this whole change exists to fix. A
// secondary rate limit and a real permission denial are both 403 DENIED and
// differ only in the message, so they are asserted together.
func TestClassifyTellsARateLimitFromADenial(t *testing.T) {
	tests := []struct {
		why  string
		in   error
		want error
	}{
		{
			"the recorded GHCR secondary limit",
			codedErr(403, "DENIED", ghcrSecondaryLimitMessage),
			ErrRateLimited,
		},
		{
			"GitHub's older wording",
			codedErr(403, "DENIED", "You have triggered an abuse detection mechanism."),
			ErrRateLimited,
		},
		{
			"the OCI code, whatever the status",
			codedErr(403, "TOOMANYREQUESTS", "slow down"),
			ErrRateLimited,
		},
		// The one that must not move: an expired or under-scoped token is the
		// other reason a push gets a 403, and waiting four escalating backoffs
		// to arrive at "your token is wrong" would be the worse bug.
		{
			"an ordinary permission denial",
			codedErr(403, "DENIED", "permission_denied: write_package"),
			ErrUnauthorized,
		},
		{
			"a plain unauthorized",
			codedErr(401, "UNAUTHORIZED", "authentication required"),
			ErrUnauthorized,
		},
	}

	for _, tt := range tests {
		got := classify(tt.in)
		if !errors.Is(got, tt.want) {
			t.Errorf("%s: classify = %v, want %v", tt.why, got, tt.want)
		}
		// Reporting a rate limit as an ordinary denial is what sent the operator
		// to check a token that was fine; the reverse would retry a real denial.
		for _, other := range []error{ErrRateLimited, ErrUnauthorized} {
			if other != tt.want && errors.Is(got, other) {
				t.Errorf("%s: classify also claimed %v", tt.why, other)
			}
		}
	}
}

// ErrUnavailable is the category that authorizes a silent local rebuild, so a
// registry asking for a sixty-second pause must never land in it.
func TestRateLimitIsNotUnavailable(t *testing.T) {
	for _, in := range []error{
		responseErr(429),
		codedErr(403, "DENIED", ghcrSecondaryLimitMessage),
	} {
		if err := classify(in); errors.Is(err, ErrUnavailable) {
			t.Errorf("classify(%v) = %v, which would authorize a local rebuild", in, err)
		}
	}
}

// A registry that cannot store the OCI 1.1 shape is describing itself, not the
// artifact: no credential, retry, or rebuild changes it, only a destination.
func TestClassifyReportsAnUnsupportedManifestShape(t *testing.T) {
	for _, code := range []string{"MANIFEST_INVALID", "UNSUPPORTED"} {
		err := classify(codedErr(400, code, "unsupported config media type"))
		if !errors.Is(err, ErrIncompatibleRegistry) {
			t.Errorf("classify(%s) = %v, want ErrIncompatibleRegistry", code, err)
		}
		if errors.Is(err, ErrRateLimited) || errors.Is(err, ErrUnauthorized) {
			t.Errorf("classify(%s) = %v, which names the wrong next step", code, err)
		}
	}
}

// The categories are distinct values, so a caller branching on one never matches
// another by accident.
func TestSentinelsAreDistinct(t *testing.T) {
	all := []error{
		ErrNotFound, ErrUnsupportedPlatform, ErrUnavailable, ErrUnauthorized,
		ErrIncompatible, ErrMismatch, ErrNoAnnotations,
		ErrRateLimited, ErrIncompatibleRegistry,
	}
	for i, a := range all {
		for j, b := range all {
			if i != j && errors.Is(a, b) {
				t.Errorf("%v and %v are the same value", a, b)
			}
		}
	}
}
