package registry

import (
	"context"
	"errors"
	"fmt"

	ocispec "github.com/opencontainers/image-spec/specs-go/v1"
)

// ResolveArtifact locates the artifact at "<base>/<repo>:<tag>" and returns the
// descriptor of the image manifest this machine can use, together with that
// manifest's annotations — without fetching the payload.
//
// That is the point of it: the annotations say what the artifact is, and they
// cost one manifest request against a payload measured in gigabytes. A caller
// decides whether to want it by passing them to [Check].
//
// tag may be a tag or a "sha256:…" digest, which is how an exact address is
// pulled. A multi-arch index is descended to this platform's child; an
// #ARCH:noarch artifact is published as a bare manifest and comes back as it is.
func ResolveArtifact(ctx context.Context, base, repo, tag string) (ocispec.Descriptor, map[string]string, error) {
	repository, err := newRepository(base, repo)
	if err != nil {
		return ocispec.Descriptor{}, nil, err
	}
	desc, annotations, err := resolvePlatformManifest(ctx, repository, tag)
	if err != nil {
		return desc, nil, fmt.Errorf("cannot resolve %s: %w", FullRef(base, repo, tag), classify(err))
	}
	return desc, annotations, nil
}

// Exists reports whether "<base>/<repo>:<tag>" is published.
//
// The error is returned rather than folded into false, because "nothing is
// published here" and "the registry would not say" are different answers and
// push acts on the first by overwriting nothing. Only [ErrNotFound] is reported
// as a plain false.
func Exists(ctx context.Context, base, repo, tag string) (bool, error) {
	repository, err := newRepository(base, repo)
	if err != nil {
		return false, err
	}
	if _, err := repository.Resolve(ctx, tag); err != nil {
		if err = classify(err); errors.Is(err, ErrNotFound) {
			return false, nil
		}
		return false, fmt.Errorf("cannot check %s: %w", FullRef(base, repo, tag), err)
	}
	return true, nil
}

// ListTags returns every tag published under "<base>/<repo>", in whatever order
// the registry pages them.
//
// A repository with nothing in it answers 404, which is reported as no tags
// rather than as an error: a name nobody has pushed yet is the ordinary state of
// a name, not a failure.
func ListTags(ctx context.Context, base, repo string) ([]string, error) {
	repository, err := newRepository(base, repo)
	if err != nil {
		return nil, err
	}
	var tags []string
	err = repository.Tags(ctx, "", func(page []string) error {
		tags = append(tags, page...)
		return nil
	})
	if err != nil {
		if err = classify(err); errors.Is(err, ErrNotFound) {
			return nil, nil
		}
		return nil, fmt.Errorf("cannot list tags for %s/%s: %w", TrimBaseScheme(base), repo, err)
	}
	return tags, nil
}
