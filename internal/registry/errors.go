package registry

import (
	"context"
	"errors"
	"fmt"
	"net"
	"net/http"
	"strings"

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
	// ErrRateLimited reports that the registry is refusing requests for now and
	// will accept them again later. Wait and retry; never fall back. Kept apart
	// from ErrUnavailable on purpose, because that is the category that
	// authorizes a silent local rebuild, and rebuilding a forty-minute index
	// because a registry asked for a sixty-second pause is the worse outcome.
	ErrRateLimited = errors.New("registry is rate limiting this client")
	// ErrIncompatibleRegistry reports a registry that will not store the OCI 1.1
	// artifact shape CondaTainer publishes. Report; it is a property of the
	// destination, not of the artifact, and no retry or rebuild changes it.
	ErrIncompatibleRegistry = errors.New("registry does not accept this artifact format")
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

// categories is every verdict classify can reach, used to recognize an error it
// has already seen. Keep it complete: a category missing here is one that gets
// classified twice and reported twice.
var categories = []error{
	ErrNotFound, ErrUnsupportedPlatform, ErrUnavailable, ErrRateLimited,
	ErrUnauthorized, ErrIncompatibleRegistry, ErrIncompatible, ErrMismatch,
	ErrNoAnnotations, ErrInvalidArtifact,
}

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
	// Already classified deeper down. Retry runs inside the transport now, so a
	// blob or manifest error arrives here having passed through classify once
	// already; wrapping it again would double the message without adding a
	// category, and would let a second look reach a different verdict.
	for _, category := range categories {
		if errors.Is(err, category) {
			return err
		}
	}
	if errors.Is(err, errdef.ErrNotFound) {
		return fmt.Errorf("%w: %w", ErrNotFound, err)
	}

	var resp *errcode.ErrorResponse
	if errors.As(err, &resp) {
		switch {
		// Both ahead of the 401/403 branch, which would otherwise swallow them:
		// a secondary rate limit arrives as a 403, and a registry refusing the
		// manifest shape may too.
		case isRateLimited(err, resp):
			return fmt.Errorf("%w: %w", ErrRateLimited, err)
		case isUnsupportedFormat(err):
			return fmt.Errorf("%w: %w", ErrIncompatibleRegistry, err)
		case resp.StatusCode == http.StatusUnauthorized, resp.StatusCode == http.StatusForbidden:
			return fmt.Errorf("%w: %w", ErrUnauthorized, err)
		case resp.StatusCode == http.StatusNotFound:
			return fmt.Errorf("%w: %w", ErrNotFound, err)
		case resp.StatusCode >= 500:
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

// secondaryLimitPhrases are what GitHub's anti-abuse refusals say, lowercased.
// The current wording is "secondary rate limit"; "abuse detection mechanism" is
// the older one and still appears.
//
// A phrase rather than a status code because GHCR delivers this as a plain
// `403 DENIED` — indistinguishable by code from a real permission denial, and
// the message is the only thing that tells them apart.
var secondaryLimitPhrases = []string{"secondary rate limit", "abuse detection"}

// isRateLimited reports whether a registry said it is refusing requests for now.
//
// Positive evidence only. A bare 403, or a 403 whose code is DENIED and whose
// message says nothing about a limit, is an ordinary permission failure, and
// retrying one for four escalating waits before reporting it wastes ten minutes
// to arrive at "your token is wrong".
func isRateLimited(err error, resp *errcode.ErrorResponse) bool {
	if resp.StatusCode == http.StatusTooManyRequests {
		return true
	}
	var errs errcode.Errors
	if !errors.As(err, &errs) {
		return false
	}
	for _, e := range errs {
		if strings.EqualFold(e.Code, "TOOMANYREQUESTS") {
			return true
		}
		msg := strings.ToLower(e.Message)
		for _, phrase := range secondaryLimitPhrases {
			if strings.Contains(msg, phrase) {
				return true
			}
		}
	}
	return false
}

// isUnsupportedFormat reports whether a registry refused the OCI 1.1 artifact
// shape itself — an artifactType, or the empty config blob that goes with it —
// rather than refusing this particular artifact.
//
// The codes are enough: CondaTainer only ever pushes one manifest shape, so a
// registry calling it invalid or unsupported is describing what it can store.
// Reported apart from a generic failure because it names a different next step:
// no credential, retry, or rebuild helps, only a different destination.
func isUnsupportedFormat(err error) bool {
	var errs errcode.Errors
	if !errors.As(err, &errs) {
		return false
	}
	for _, e := range errs {
		switch strings.ToUpper(e.Code) {
		case errcode.ErrorCodeManifestInvalid, errcode.ErrorCodeUnsupported:
			return true
		}
	}
	return false
}
