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
		{"429", responseErr(429), ErrUnavailable},
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
			for _, category := range []error{ErrNotFound, ErrUnauthorized, ErrUnavailable} {
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

// The categories are distinct values, so a caller branching on one never matches
// another by accident.
func TestSentinelsAreDistinct(t *testing.T) {
	all := []error{
		ErrNotFound, ErrUnsupportedPlatform, ErrUnavailable, ErrUnauthorized,
		ErrIncompatible, ErrMismatch, ErrNoAnnotations,
	}
	for i, a := range all {
		for j, b := range all {
			if i != j && errors.Is(a, b) {
				t.Errorf("%v and %v are the same value", a, b)
			}
		}
	}
}
