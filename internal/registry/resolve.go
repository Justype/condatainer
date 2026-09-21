package registry

import (
	"context"
	"fmt"
	"runtime"
	"strings"
	"time"

	ocispec "github.com/opencontainers/image-spec/specs-go/v1"
	"oras.land/oras-go/v2"
	"oras.land/oras-go/v2/registry/remote"
)

// DefaultTimeout bounds a resolution: an unresolved reference is recoverable, a
// stalled build is not.
const DefaultTimeout = 20 * time.Second

// Resolve returns the digest ref names on this host's platform, as sha256:<hex>.
//
// Platform-specific, not the index digest a multi-arch tag gives by default: an
// index would give two hosts one digest for two different payloads.
func Resolve(ctx context.Context, ref string) (string, error) {
	return resolveOn(ctx, ref, &ocispec.Platform{OS: runtime.GOOS, Architecture: runtime.GOARCH})
}

// resolveOn is Resolve with the platform supplied, so a test can ask what
// another host would resolve.
func resolveOn(ctx context.Context, ref string, platform *ocispec.Platform) (string, error) {
	normalized, err := Normalize(ref)
	if err != nil {
		return "", err
	}

	repo, err := remote.NewRepository(normalized)
	if err != nil {
		return "", fmt.Errorf("%s is not a usable registry reference: %w", ref, err)
	}
	repo.Client = newAuthClient(repo.Reference.Registry + "/" + repo.Reference.Repository)

	ctx, cancel := context.WithTimeout(ctx, DefaultTimeout)
	defer cancel()

	desc, err := oras.Resolve(ctx, repo, repo.Reference.Reference, oras.ResolveOptions{TargetPlatform: platform})
	if err != nil {
		return "", fmt.Errorf("cannot resolve %s: %w", ref, err)
	}
	return desc.Digest.String(), nil
}

// Normalize expands a short reference the way a container runtime does, so
// `ubuntu:24.04` and `docker.io/library/ubuntu:24.04` resolve identically.
func Normalize(ref string) (string, error) {
	ref = strings.TrimSpace(ref)
	if ref == "" {
		return "", fmt.Errorf("empty registry reference")
	}
	if strings.Contains(ref, "://") {
		return "", fmt.Errorf("%q is a URI, not a registry reference", ref)
	}

	// A leading component is a registry host only if it looks like one, the same
	// heuristic Docker uses: `myorg/tool` is a Hub namespace, `reg.io/tool` is not.
	const dockerHub = "docker.io"
	head, _, hasSlash := strings.Cut(ref, "/")
	switch {
	case !hasSlash:
		return dockerHub + "/library/" + ref, nil
	case head == "localhost" || strings.ContainsAny(head, ".:"):
		return ref, nil
	default:
		return dockerHub + "/" + ref, nil
	}
}
