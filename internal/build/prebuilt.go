package build

import (
	"context"
	"errors"
	"fmt"

	"github.com/Justype/condatainer/internal/artifact/key"
	"github.com/Justype/condatainer/internal/artifact/meta"
	"github.com/Justype/condatainer/internal/logging"
	"github.com/Justype/condatainer/internal/registry"
)

// prebuiltResult says whether build should stop because an artifact was
// installed. A false result means the caller may continue with the local build.
type prebuiltResult bool

var (
	resolvePrebuilt = registry.ResolveArtifact
	pullPrebuilt    = registry.PullLocked
)

// tryPrebuilt tries the selected recipe source's ordered pull endpoints. The
// caller already holds the target's producer lock, so pull uses the locked
// transport entry point and installs atomically into the final pathname.
func (b *BuildObject) tryPrebuilt(ctx context.Context) (prebuiltResult, error) {
	if b.catalogSource == nil || b.catalogSource.DescriptorErr != nil || len(b.catalogSource.Desc.OCI.Pull) == 0 {
		return false, nil
	}
	if b.catalogSource.Desc.OCI.Audience == "public" && b.spec.Image.Type == "app" {
		return false, nil
	}

	want, err := b.prebuiltEquivalence(ctx)
	if err != nil {
		return false, fmt.Errorf("cannot derive expected equivalence for %s: %w", b.spec.Image.Name, err)
	}
	repo, tag, err := registry.PullReference(b.spec.Image.Type, b.spec.Image.Name)
	if err != nil {
		return false, err
	}

	log := logging.FromContext(ctx)
	var lastUnavailable error
	for _, endpoint := range b.catalogSource.Desc.OCI.Pull {
		desc, annotations, err := resolvePrebuilt(ctx, endpoint, repo, tag)
		if err != nil {
			switch {
			case errors.Is(err, registry.ErrNotFound), errors.Is(err, registry.ErrUnsupportedPlatform):
				continue
			case errors.Is(err, registry.ErrUnavailable):
				lastUnavailable = err
				continue
			default:
				return false, fmt.Errorf("cannot use prebuilt %s from %s: %w", b.spec.Image.Name, endpoint, err)
			}
		}
		if err := registry.Check(annotations, registry.Want{Name: b.spec.Image.Name}); err != nil {
			return false, fmt.Errorf("cannot use prebuilt %s from %s: %w", b.spec.Image.Name, endpoint, err)
		}
		if got := registry.Equiv(annotations); got != want {
			return false, fmt.Errorf("%w: prebuilt equivalence %s, selected recipe derives %s",
				registry.ErrInvalidArtifact, describePrebuiltKey(got), describePrebuiltKey(want))
		}
		// Said before the download rather than after it: everything above is
		// metadata, and the gigabytes start here. Without this the operator
		// watches a long transfer with nothing saying what is being fetched or
		// that it has already been checked against the local recipe.
		log.Info("prebuilt found and verified", "artifact", b.spec.Image.Name,
			"endpoint", endpoint, "equivalence", describePrebuiltKey(want))
		if err := pullPrebuilt(ctx, endpoint, repo, desc, annotations, b.tgt.Path); err != nil {
			switch {
			case errors.Is(err, registry.ErrNotFound), errors.Is(err, registry.ErrUnsupportedPlatform):
				continue
			case errors.Is(err, registry.ErrUnavailable):
				lastUnavailable = err
				continue
			default:
				return false, fmt.Errorf("cannot pull prebuilt %s from %s: %w", b.spec.Image.Name, endpoint, err)
			}
		}
		invalidateInstalledOverlays()
		log.Info("prebuilt image ready", "kind", "success", "path", b.tgt.Path, "endpoint", endpoint)
		return true, nil
	}
	if lastUnavailable != nil {
		log.Warn("registry endpoints unavailable; building selected recipe locally",
			"name", b.spec.Image.Name, "err", lastUnavailable)
	}
	return false, nil
}

// prebuiltEquivalence derives the key from the selected local recipe and the
// dependencies currently installed. Identity is intentionally ignored: a data
// artifact built against equivalent dependency instances is substitutable even
// when their exact identities differ.
func (b *BuildObject) prebuiltEquivalence(ctx context.Context) (meta.KeyRef, error) {
	recipe, ok := b.spec.Source.RecipeFile()
	if !ok {
		return meta.KeyRef{}, errors.New("selected source has no recipe")
	}
	manifest := b.Manifest()
	manifest.Dependencies, _ = key.Manifest(key.Artifact{
		Name:         b.spec.Image.Name,
		Type:         b.spec.Image.Type,
		Env:          b.spec.Image.Env,
		Recipe:       recipe.Data,
		Placeholders: b.spec.Source.Placeholders,
		Deps:         b.dependencyKeys(ctx),
		From:         b.spec.Source.UpstreamDigest(),
	})
	derived, err := key.Generate(manifest, b.keySources())
	if err != nil {
		return meta.KeyRef{}, err
	}
	return derived.Equiv.Ref, nil
}

func describePrebuiltKey(k meta.KeyRef) string {
	if k.Empty() {
		return "absent"
	}
	return k.Scheme + " " + k.Digest()
}
