package registry

import (
	"context"
	"errors"
	"fmt"
	"net"
	"net/http"

	"oras.land/oras-go/v2/errdef"
	"oras.land/oras-go/v2/registry/remote/errcode"
)

// Why a registry operation failed, as the categories a caller acts on.
//
// They are sentinels rather than one error type because nothing reads a field
// off them — a caller recognizes a case and branches. The split is by *what the
// caller does*, which is the only line worth drawing: an artifact that is absent
// or unreachable may be built locally instead, while a credential problem or a
// wrong artifact must be reported. Silently rebuilding a forty-minute index
// because a token expired is a worse outcome than an error.
var (
	// ErrNotFound reports that no such artifact is published. Fall back.
	ErrNotFound = errors.New("artifact is not published")
	// ErrUnsupportedPlatform reports that the artifact exists but carries no
	// payload for this architecture. Fall back, with a different message: this
	// one is somebody's missing push, not a missing artifact.
	ErrUnsupportedPlatform = errors.New("artifact is not published for this platform")
	// ErrUnavailable reports that the registry could not be reached or answered
	// that it was in trouble. Fall back, with a note.
	ErrUnavailable = errors.New("registry is unavailable")
	// ErrUnauthorized reports that the credential in hand does not open this
	// artifact. Report; never fall back, or an expired token silently becomes a
	// long rebuild.
	ErrUnauthorized = errors.New("not authorized for this registry")
	// ErrIncompatible reports metadata this build cannot read. Report.
	ErrIncompatible = errors.New("artifact needs a different CondaTainer version")
	// ErrMismatch reports a well-formed artifact that is not the one asked for.
	// Report, or keep looking; never install.
	ErrMismatch = errors.New("artifact is not the one requested")
	// ErrNoAnnotations reports a manifest carrying no CondaTainer metadata at
	// all, which means it was not published by CondaTainer.
	ErrNoAnnotations = errors.New("not a CondaTainer artifact")
	// ErrInvalidArtifact reports a published artifact that contradicts itself —
	// a payload whose regenerated keys are not the ones advertised, a manifest
	// whose layers are not the type it claims. Report; never install. One
	// sentinel covers them because a caller does the same thing about each: the
	// difference between them is a message, not a decision.
	ErrInvalidArtifact = errors.New("published artifact is not coherent")
)

// classify tags a transport error with the category a caller branches on,
// keeping the original wrapped so the message still says what actually happened.
//
// The default is to pass an unrecognized error through unclassified, so a caller
// treating anything it cannot place as fatal is right by default. Guessing
// [ErrUnavailable] here would be the dangerous direction — that is the category
// that authorizes a silent local rebuild.
func classify(err error) error {
	if err == nil {
		return nil
	}
	// A cancelled context is the user, not the registry.
	if errors.Is(err, context.Canceled) {
		return err
	}
	if errors.Is(err, errdef.ErrNotFound) {
		return fmt.Errorf("%w: %w", ErrNotFound, err)
	}

	var resp *errcode.ErrorResponse
	if errors.As(err, &resp) {
		switch {
		case resp.StatusCode == http.StatusUnauthorized, resp.StatusCode == http.StatusForbidden:
			return fmt.Errorf("%w: %w", ErrUnauthorized, err)
		case resp.StatusCode == http.StatusNotFound:
			return fmt.Errorf("%w: %w", ErrNotFound, err)
		case resp.StatusCode == http.StatusTooManyRequests, resp.StatusCode >= 500:
			return fmt.Errorf("%w: %w", ErrUnavailable, err)
		}
		return err
	}

	// Nothing answered at all: a dial failure, a DNS failure, or a timeout. Last,
	// because it is the least specific thing that can be said about a failure —
	// a registry that did answer has already been classified by what it said.
	var netErr net.Error
	if errors.As(err, &netErr) || errors.Is(err, context.DeadlineExceeded) {
		return fmt.Errorf("%w: %w", ErrUnavailable, err)
	}
	return err
}
