package registry

import (
	"bytes"
	"context"
	"encoding/json"
	"errors"
	"fmt"
	"runtime"
	"slices"
	"sort"

	"github.com/opencontainers/go-digest"
	specs "github.com/opencontainers/image-spec/specs-go"
	ocispec "github.com/opencontainers/image-spec/specs-go/v1"
	"oras.land/oras-go/v2/content"
	"oras.land/oras-go/v2/errdef"
	"oras.land/oras-go/v2/registry/remote"
)

// platform returns the OCI platform of the running machine, used as an index
// child's platform descriptor.
func platform() (ocispec.Platform, bool) {
	switch runtime.GOARCH {
	case "amd64", "arm64":
		return ocispec.Platform{OS: "linux", Architecture: runtime.GOARCH}, true
	default:
		return ocispec.Platform{}, false
	}
}

// platformKey is a stable map key for an index child, os/arch[/variant].
func platformKey(p *ocispec.Platform) string {
	if p == nil {
		return ""
	}
	key := p.OS + "/" + p.Architecture
	if p.Variant != "" {
		key += "/" + p.Variant
	}
	return key
}

// publishIndex writes the image index assembled from entries and points every
// tag at it. An index already present is not an error — the digest is the same
// content, so a re-push is a no-op rather than a conflict.
func publishIndex(ctx context.Context, repository *remote.Repository, tags []string, entries map[string]ocispec.Descriptor, annotations map[string]string, artifactType string) error {
	indexBytes, err := json.Marshal(buildIndex(entries, annotations, artifactType))
	if err != nil {
		return fmt.Errorf("failed to encode index: %w", err)
	}
	indexDesc := ocispec.Descriptor{
		MediaType:    ocispec.MediaTypeImageIndex,
		ArtifactType: artifactType,
		Digest:       digest.FromBytes(indexBytes),
		Size:         int64(len(indexBytes)),
	}
	// Retried, so that a rate limit here is waited out rather than surfacing to
	// reconcileIndex, whose own retry means something else entirely: it would
	// report a throttled push as a concurrent publisher and send the reader
	// looking for someone who is not there.
	if err := retryPolicyFrom(ctx).run(ctx, mutation{verb: verbPublish, attrs: []any{"target", "index"}, do: func(ctx context.Context) error {
		err := repository.Push(ctx, indexDesc, bytes.NewReader(indexBytes))
		if err != nil && !errors.Is(err, errdef.ErrAlreadyExists) {
			return err
		}
		return nil
	}}); err != nil {
		return fmt.Errorf("failed to push index: %w", err)
	}
	for _, tag := range tags {
		if err := retryPolicyFrom(ctx).run(ctx, mutation{
			verb:  verbPublish,
			attrs: []any{"tag", tag},
			do:    func(ctx context.Context) error { return repository.Tag(ctx, indexDesc, tag) },
		}); err != nil {
			return fmt.Errorf("failed to add tag %s: %w", tag, err)
		}
	}
	return nil
}

// buildIndex assembles an OCI image index from per-platform children, sorted by
// platform key so the same set of children always encodes to the same digest.
func buildIndex(entries map[string]ocispec.Descriptor, annotations map[string]string, artifactType string) ocispec.Index {
	keys := make([]string, 0, len(entries))
	for k := range entries {
		keys = append(keys, k)
	}
	sort.Strings(keys)
	manifests := make([]ocispec.Descriptor, 0, len(entries))
	for _, k := range keys {
		manifests = append(manifests, entries[k])
	}
	return ocispec.Index{
		Versioned:    specs.Versioned{SchemaVersion: 2},
		MediaType:    ocispec.MediaTypeImageIndex,
		ArtifactType: artifactType,
		Manifests:    manifests,
		Annotations:  indexAnnotations(annotations),
	}
}

// indexAnnotations adapts a child's annotations for the index above it, so a
// registry UI shows a name and an identity without descending into a child.
//
// Only the compressor is dropped. Everything else an artifact publishes is a
// property of the artifact rather than of one build of it — the identity and
// equivalence keys are architecture-independent by construction, and the name,
// description, and source are the recipe's. The compressor is a property of the
// bytes, so the two children could disagree, and the consumer who needs it is
// resolving its own platform's child anyway.
//
// org.opencontainers.image.created is the pushing architecture's build time, not
// a claim about the others. That is the same looseness the date tag already has
// and is why the tag is an address rather than a fact.
func indexAnnotations(ann map[string]string) map[string]string {
	out := make(map[string]string, len(ann))
	for k, v := range ann {
		if k == AnnCompression {
			continue
		}
		out[k] = v
	}
	return out
}

// collectIndexEntries reads the children already published at any of tags, so a
// push preserves architectures it is not currently building. A tag that is absent
// or is not an index contributes nothing; this is a best-effort read, and an
// unreachable registry surfaces at the push that follows.
func collectIndexEntries(ctx context.Context, repo *remote.Repository, tags []string) map[string]ocispec.Descriptor {
	entries := map[string]ocispec.Descriptor{}
	for _, t := range tags {
		mergeIndexEntries(entries, readIndexEntries(ctx, repo, t))
	}
	return entries
}

func readIndexEntries(ctx context.Context, repo *remote.Repository, tag string) map[string]ocispec.Descriptor {
	entries := map[string]ocispec.Descriptor{}
	desc, err := repo.Resolve(ctx, tag)
	if err != nil || desc.MediaType != ocispec.MediaTypeImageIndex {
		return entries
	}
	data, err := content.FetchAll(ctx, repo, desc)
	if err != nil {
		return entries
	}
	var idx ocispec.Index
	if json.Unmarshal(data, &idx) != nil {
		return entries
	}
	for _, manifest := range idx.Manifests {
		key := platformKey(manifest.Platform)
		if key == "" {
			// A child with no platform cannot be keyed by one, but dropping it
			// would silently unpublish it on the next push.
			key = manifest.Digest.String()
		}
		entries[key] = manifest
	}
	return entries
}

// mergeIndexEntries adds the entries of src that dst has no platform for. Keeping
// dst's is what makes the local architecture win over the published one: the
// pusher put its own child in before merging what the registry had.
func mergeIndexEntries(dst, src map[string]ocispec.Descriptor) {
	for key, desc := range src {
		if _, exists := dst[key]; !exists {
			dst[key] = desc
		}
	}
}

// containsIndexEntries reports whether have carries every platform in want at the
// same digest. A platform present at a different digest is a failure, not a
// match — it is somebody else's push landing on top of this one.
func containsIndexEntries(have, want map[string]ocispec.Descriptor) bool {
	for key, desc := range want {
		got, ok := have[key]
		if !ok || got.Digest != desc.Digest {
			return false
		}
	}
	return true
}

// verifyTaggedIndexes re-reads every tag and reports whether all of them now
// carry want, along with the union of what was actually observed.
//
// Every alias is checked, not just the canonical tag: a registry offers no
// compare-and-swap on tags, so a concurrent push can land between tagging the
// date and tagging "latest" and leave the two pointing at different indexes. The
// union feeds the next reconciliation attempt.
func verifyTaggedIndexes(ctx context.Context, repo *remote.Repository, tags []string, want map[string]ocispec.Descriptor) (observed map[string]ocispec.Descriptor, complete bool) {
	observed = map[string]ocispec.Descriptor{}
	complete = true
	for _, tag := range tags {
		atTag := readIndexEntries(ctx, repo, tag)
		if !containsIndexEntries(atTag, want) {
			complete = false
		}
		mergeIndexEntries(observed, atTag)
	}
	return observed, complete
}

// resolvePlatformManifest resolves tag to the image manifest this machine can
// use, returning its descriptor and manifest-level annotations without fetching
// the payload — which is the point, since the payload is measured in gigabytes
// and the annotations are what decide whether to want it.
//
// An image index is descended to the child matching linux/<arch>; a plain image
// manifest is returned as it is, which is how an #ARCH:noarch artifact arrives.
func resolvePlatformManifest(ctx context.Context, repository *remote.Repository, tag string) (ocispec.Descriptor, map[string]string, error) {
	desc, err := repository.Resolve(ctx, tag)
	if err != nil {
		return ocispec.Descriptor{}, nil, err
	}
	if desc.MediaType != ocispec.MediaTypeImageIndex {
		manifest, err := fetchManifest(ctx, repository, desc)
		if err != nil {
			return desc, nil, err
		}
		desc.ArtifactType = manifest.ArtifactType
		return desc, manifest.Annotations, nil
	}

	data, err := content.FetchAll(ctx, repository, desc)
	if err != nil {
		return desc, nil, fmt.Errorf("failed to fetch index: %w", err)
	}
	var idx ocispec.Index
	if err := json.Unmarshal(data, &idx); err != nil {
		return desc, nil, fmt.Errorf("failed to parse index: %w", err)
	}
	plat, ok := platform()
	if !ok {
		return desc, nil, fmt.Errorf("%w: %s is not a distribution architecture",
			ErrUnsupportedPlatform, runtime.GOARCH)
	}
	for _, child := range idx.Manifests {
		if child.Platform == nil || child.Platform.OS != plat.OS || child.Platform.Architecture != plat.Architecture {
			continue
		}
		manifest, err := fetchManifest(ctx, repository, child)
		if err != nil {
			return child, nil, err
		}
		child.ArtifactType = manifest.ArtifactType
		return child, manifest.Annotations, nil
	}
	return desc, nil, fmt.Errorf("%w: index carries %v, not %s",
		ErrUnsupportedPlatform, indexPlatforms(idx), platformKey(&plat))
}

// indexPlatforms lists what an index does carry, so a refusal names the
// architectures somebody did publish instead of only the one that is missing.
func indexPlatforms(idx ocispec.Index) []string {
	var keys []string
	for _, child := range idx.Manifests {
		if key := platformKey(child.Platform); key != "" {
			keys = append(keys, key)
		}
	}
	slices.Sort(keys)
	return slices.Compact(keys)
}

// fetchManifest reads and decodes an image manifest.
func fetchManifest(ctx context.Context, repository *remote.Repository, desc ocispec.Descriptor) (ocispec.Manifest, error) {
	data, err := content.FetchAll(ctx, repository, desc)
	if err != nil {
		return ocispec.Manifest{}, fmt.Errorf("failed to fetch manifest: %w", err)
	}
	var manifest ocispec.Manifest
	if err := json.Unmarshal(data, &manifest); err != nil {
		return ocispec.Manifest{}, fmt.Errorf("invalid image manifest: %w", err)
	}
	return manifest, nil
}
