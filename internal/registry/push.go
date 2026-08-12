package registry

import (
	"context"
	"errors"
	"fmt"

	ocispec "github.com/opencontainers/image-spec/specs-go/v1"
	"oras.land/oras-go/v2"
	"oras.land/oras-go/v2/registry/remote"

	"github.com/Justype/condatainer/catalog"
	"github.com/Justype/condatainer/internal/artifact/compare"
	"github.com/Justype/condatainer/internal/artifact/meta"
	"github.com/Justype/condatainer/internal/image/sif"
	"github.com/Justype/condatainer/internal/image/squashfs"
	"github.com/Justype/condatainer/internal/logging"
	"github.com/Justype/condatainer/internal/utils"
)

// Visibility declares who can pull from an endpoint, and so decides what may be
// pushed to it.
//
// It states a fact about a registry, and CondaTainer derives the permitted set
// from it. Config must never enumerate types directly: `types: [app]` beside a
// public endpoint would erase the rule with no error, whereas a wrong visibility
// is a claim someone has to write down and defend. Nothing here verifies the
// registry really is private — it is a declaration, not enforcement.
type Visibility string

const (
	// Public is the default. Anyone can pull, so only what is ours to
	// redistribute may go here.
	Public Visibility = "public"
	// Internal is a registry only the organization can pull from.
	Internal Visibility = "internal"
)

// maxIndexReconcileAttempts bounds the retry when another architecture is
// publishing at the same time. Registries offer no compare-and-swap on tags, so
// this converges rather than locks; CI should still publish sequentially.
const maxIndexReconcileAttempts = 3

// Accepts reports why an artifact may not be published to an endpoint of this
// visibility, or nil.
//
// A public endpoint takes `base`, `os`, and `data` — a container root, packages
// from a public distribution, and public reference data with the indexes built
// from it. It does not take an `app`: an app installs someone else's software,
// and publishing a copy redistributes their binaries under our name, which is
// theirs to permit and not ours to assume. A Conda build is refused for a second
// and independent reason — its solve inputs are embedded, so a consumer rebuilds
// from a few kilobytes instead of downloading gigabytes.
func (v Visibility) Accepts(m meta.Manifest) error {
	if v == Internal {
		return nil
	}
	if m.BuildType == "conda" {
		return fmt.Errorf("a Conda build is not published to a %s registry: its solve inputs travel with it, so %s rebuilds from kilobytes", v, m.Name)
	}
	switch m.Type {
	case catalog.TypeBase, catalog.TypeOS, catalog.TypeData:
		return nil
	}
	return fmt.Errorf("a %s artifact is not published to a %s registry: %s installs software that is not ours to redistribute", m.Type, v, m.Name)
}

// PublishRequest is one artifact and where it is going.
type PublishRequest struct {
	// Path is the local artifact. Its embedded manifest decides everything
	// about the destination, so renaming the file cannot smuggle it anywhere.
	Path string
	// Base is the registry base, owner and prefix included.
	Base string
	// Visibility is what the endpoint declared.
	Visibility Visibility
	// Force allows replacing an existing versioned tag, which is otherwise
	// immutable.
	Force bool
}

// Publish validates an artifact against its destination, then pushes it.
//
// Everything published is derived from the artifact's own bytes: the type that
// decides whether the endpoint will take it, the keys, the name, and every
// annotation. Nothing is accepted from the caller except where it is going.
func Publish(ctx context.Context, req PublishRequest) error {
	log := logging.FromContext(ctx)
	if req.Visibility == "" {
		req.Visibility = Public
	}
	if _, _, err := imageTypes(req.Path); err != nil {
		return err
	}

	m, err := meta.ReadManifest(req.Path)
	if err != nil {
		return err
	}
	if err := req.Visibility.Accepts(m); err != nil {
		return err
	}
	if m.Keys.Identity.Empty() || m.Keys.Equiv.Empty() {
		return fmt.Errorf("%s records no scheme-backed identity, so nothing could pin what was pulled", m.Name)
	}
	// Held to its own claim before anyone else has to trust it: a published
	// artifact whose files do not reproduce its recorded keys is one that pull
	// would refuse on arrival, and finding that out here is cheaper for everyone.
	got, err := compare.Read(req.Path)
	if err != nil {
		return fmt.Errorf("%w: %w", ErrInvalidArtifact, err)
	}
	if err := checkRegeneratedKeys(got, m.Keys.Identity, m.Keys.Equiv, "recorded"); err != nil {
		return err
	}

	repo, tags, err := PushReference(m)
	if err != nil {
		return err
	}
	if err := checkTagIsFree(ctx, req, m, repo, tags[0]); err != nil {
		return err
	}

	annotations := Annotations(m, compressionOf(req.Path))
	log.Info("publishing", "artifact", m.Name, "reference", FullRef(req.Base, repo, tags[0]))
	return Push(ctx, req.Path, req.Base, repo, tags, annotations, m.Platform.Arch == meta.ArchNone)
}

// checkTagIsFree refuses to displace already-published content without Force.
//
// Two things are deliberately not refusals. A version-less artifact is addressed
// by build date and a rolling tag, so re-publishing it is the point of the
// scheme. And a tag that already resolves to an index this architecture is not
// in yet is an *addition*: one tag serving several architectures is what the
// index is for, so refusing there would break the multi-arch flow the reconcile
// loop exists to support.
//
// What is refused is replacing a payload someone may already have pulled: a
// child for this platform, or any content under a tag that carries no index at
// all.
func checkTagIsFree(ctx context.Context, req PublishRequest, m meta.Manifest, repo, tag string) error {
	if req.Force || isVersionLess(m.Type, m.Name) {
		return nil
	}
	repository, err := newRepository(req.Base, repo)
	if err != nil {
		return err
	}
	desc, err := repository.Resolve(ctx, tag)
	if err != nil {
		if err = classify(err); errors.Is(err, ErrNotFound) {
			return nil
		}
		return err
	}

	occupant := "it"
	if desc.MediaType == ocispec.MediaTypeImageIndex && m.Platform.Arch != meta.ArchNone {
		plat, ok := platform()
		if !ok {
			return fmt.Errorf("%w: this build cannot publish from its own architecture", ErrUnsupportedPlatform)
		}
		if _, taken := readIndexEntries(ctx, repository, tag)[platformKey(&plat)]; !taken {
			return nil
		}
		occupant = platformKey(&plat)
	}
	return fmt.Errorf("%s already publishes %s; a versioned tag is immutable, so replacing it needs --force",
		FullRef(req.Base, repo, tag), occupant)
}

// Push uploads artifactPath to "<base>/<repo>" under every tag, the first of
// which is canonical.
//
// annotations stays a parameter rather than being derived here: the transport
// carries what it is given and reads none of it, which is what keeps a change to
// the metadata from being a change to the upload. [Publish] is where policy
// lives.
//
// When archIndependent is false the manifest becomes a child of an OCI image
// index carrying this architecture's platform descriptor, merged with whatever
// is already published, so running push on each architecture builds one
// multi-arch tag. An #ARCH:noarch artifact is tagged directly: the payload is
// identical everywhere, and an index over one child would only imply otherwise.
func Push(ctx context.Context, artifactPath, base, repo string, tags []string, annotations map[string]string, archIndependent bool) error {
	if len(tags) == 0 {
		return fmt.Errorf("no tags to push %s under", artifactPath)
	}
	artifactType, layerType, err := imageTypes(artifactPath)
	if err != nil {
		return err
	}
	var plat ocispec.Platform
	if !archIndependent {
		var ok bool
		if plat, ok = platform(); !ok {
			return fmt.Errorf("%w: this build cannot publish from its own architecture", ErrUnsupportedPlatform)
		}
	}

	repository, err := newRepository(base, repo)
	if err != nil {
		return err
	}

	// Blobs first: a manifest may not reference what the registry does not yet
	// hold. Nothing is staged on the way — see pushArtifactLayers.
	layers, err := pushArtifactLayers(ctx, repository.Blobs(), artifactPath, layerType)
	if err != nil {
		return classify(err)
	}

	// Packed straight against the repository. There is no local store to copy
	// from, because the layers were never written to one.
	manifestDesc, err := oras.PackManifest(ctx, repository, oras.PackManifestVersion1_1, artifactType, oras.PackManifestOptions{
		Layers:              layers,
		ManifestAnnotations: annotations,
	})
	if err != nil {
		return fmt.Errorf("failed to publish manifest: %w", classify(err))
	}
	manifestDesc.ArtifactType = artifactType

	if archIndependent {
		for _, tag := range tags {
			if err := repository.Tag(ctx, manifestDesc, tag); err != nil {
				return fmt.Errorf("failed to tag %s: %w", FullRef(base, repo, tag), classify(err))
			}
		}
		return nil
	}
	manifestDesc.Platform = &plat
	return reconcileIndex(ctx, repository, tags, manifestDesc, annotations, artifactType, base, repo)
}

// reconcileIndex points every tag at an index containing this architecture's
// child, without dropping architectures somebody else published.
//
// Re-read and retry rather than write once: tagging is not atomic across tags
// and a registry offers no compare-and-swap, so a concurrent publisher can land
// between tagging the date and tagging "latest". Each attempt folds in what was
// actually observed, so the two pushes converge instead of overwriting.
func reconcileIndex(ctx context.Context, repository *remote.Repository, tags []string, child ocispec.Descriptor, annotations map[string]string, artifactType, base, repo string) error {
	entries := collectIndexEntries(ctx, repository, tags)
	// Inserted before the merge, so the child being pushed wins over the one
	// already published for this platform.
	entries[platformKey(child.Platform)] = child

	for range maxIndexReconcileAttempts {
		if err := publishIndex(ctx, repository, tags, entries, annotations, artifactType); err != nil {
			return classify(err)
		}
		observed, complete := verifyTaggedIndexes(ctx, repository, tags, entries)
		if complete {
			return nil
		}
		mergeIndexEntries(entries, observed)
	}
	return fmt.Errorf("%s kept changing under a concurrent push; publish architectures sequentially and retry",
		FullRef(base, repo, tags[0]))
}

// compressionOf reports the SquashFS compressor of an image, or "" when it
// cannot be read. Best-effort by design: the annotation tells a puller whether
// its kernel can mount the payload, and not knowing is a reason to omit the
// claim rather than to refuse the push.
func compressionOf(path string) string {
	var offset int64
	if utils.IsSif(path) {
		part, err := sif.PrimarySystemPartition(path)
		if err != nil {
			return ""
		}
		offset = part.Offset
	}
	stats, err := squashfs.GetSquashFSStatsAt(path, offset)
	if err != nil || stats == nil {
		return ""
	}
	return stats.Compression
}
