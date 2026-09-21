package build

import (
	"context"
	"errors"
	"fmt"

	ocispec "github.com/opencontainers/image-spec/specs-go/v1"

	"github.com/Justype/condatainer/catalog"
	"github.com/Justype/condatainer/internal/artifact/key"
	"github.com/Justype/condatainer/internal/artifact/meta"
	"github.com/Justype/condatainer/internal/config"
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

// prebuiltChoice is what planning decided about a node's prebuilt artifact.
type prebuiltChoice int

const (
	prebuiltUndecided prebuiltChoice = iota // settled when the build runs
	prebuiltNone                            // build from the recipe
	prebuiltPull                            // pull the planned candidate
)

// prebuiltCandidate is a published artifact that matches the local recipe.
// next is the index of the pull endpoint after the one it came from.
type prebuiltCandidate struct {
	endpoint    string
	repo        string
	next        int
	equiv       meta.KeyRef
	desc        ocispec.Descriptor
	annotations map[string]string
}

// prebuiltPlan is a node's prebuilt decision, made before anything runs.
type prebuiltPlan struct {
	choice    prebuiltChoice
	candidate prebuiltCandidate
}

// plannedDep is what a dependency that is not installed yet will contribute to
// a dependent's equivalence: its name and type, and for data its equivalence key.
type plannedDep struct {
	Name  string
	Type  catalog.Type
	Equiv meta.KeyRef
}

// plannedAs reports what this node will contribute as a dependency. A data
// dependency contributes its equivalence key, so it is unknown until derived.
func (b *BuildObject) plannedAs() (plannedDep, bool) {
	if b.spec.Image.Type == catalog.TypeData && b.plannedEquiv.Empty() {
		return plannedDep{}, false
	}
	return plannedDep{Name: b.spec.Image.Name, Type: b.spec.Image.Type, Equiv: b.plannedEquiv}, true
}

// prebuiltMiss says why no candidate was found, most serious first.
type prebuiltMiss int

const (
	missNone        prebuiltMiss = iota // nothing published, or none for this platform
	missMismatch                        // published from a different recipe
	missUnavailable                     // an endpoint could not be reached
)

// prebuiltEligible reports whether a prebuilt can apply at all. A Conda build
// has none to pull.
func (b *BuildObject) prebuiltEligible() bool {
	return !config.Global.Build.SkipPrebuilt && b.buildType != BuildTypeConda &&
		b.catalogSource != nil && b.catalogSource.DescriptorErr == nil &&
		len(b.catalogSource.Desc.OCI.Pull) > 0
}

// findPrebuilt reads the pull endpoints from start, in order and by metadata
// alone, and returns the first artifact made from the local recipe. Every other
// outcome is a miss; a mismatch is warned about here.
func (b *BuildObject) findPrebuilt(ctx context.Context, start int) (prebuiltCandidate, bool, prebuiltMiss, error) {
	want, err := b.prebuiltEquivalence(ctx)
	if err != nil {
		return prebuiltCandidate{}, false, missNone, fmt.Errorf("cannot derive expected equivalence for %s: %w", b.spec.Image.Name, err)
	}
	repo, tag, err := registry.PullReference(b.spec.Image.Type, b.spec.Image.Name)
	if err != nil {
		return prebuiltCandidate{}, false, missNone, err
	}

	log := logging.FromContext(ctx)
	miss := missNone
	endpoints := b.catalogSource.Desc.OCI.Pull
	for i := start; i < len(endpoints); i++ {
		endpoint := endpoints[i]
		desc, annotations, err := resolvePrebuilt(ctx, endpoint, repo, tag)
		if err != nil {
			switch {
			case errors.Is(err, registry.ErrNotFound), errors.Is(err, registry.ErrUnsupportedPlatform),
				closedDoor(ctx, err, endpoint):
				continue
			case errors.Is(err, registry.ErrUnavailable):
				log.Debug("prebuilt unavailable", "name", b.spec.Image.Name, "err", err)
				miss = missUnavailable
				continue
			default:
				return prebuiltCandidate{}, false, miss, fmt.Errorf("cannot use prebuilt %s from %s: %w", b.spec.Image.Name, endpoint, err)
			}
		}
		if err := registry.Check(annotations, registry.Want{Name: b.spec.Image.Name}); err != nil {
			return prebuiltCandidate{}, false, miss, fmt.Errorf("cannot use prebuilt %s from %s: %w", b.spec.Image.Name, endpoint, err)
		}
		if got := registry.Equiv(annotations); got != want {
			log.Warn(fmt.Sprintf("Prebuilt %s at %s was made from a different recipe (prebuilt %s, this recipe %s)",
				b.spec.Image.Name, endpoint, describePrebuiltKey(got), describePrebuiltKey(want)))
			miss = max(miss, missMismatch)
			continue
		}
		return prebuiltCandidate{endpoint: endpoint, repo: repo, next: i + 1, equiv: want,
			desc: desc, annotations: annotations}, true, miss, nil
	}
	return prebuiltCandidate{}, false, miss, nil
}

// planPrebuilt settles, before anything runs, whether this node pulls a
// prebuilt artifact or builds from the recipe, so the plan can say which. known
// says every dependency's contribution is available, installed or planned; a
// node that cannot derive its equivalence yet stays undecided until it builds.
func (b *BuildObject) planPrebuilt(ctx context.Context, known bool) error {
	if known {
		if equiv, err := b.prebuiltEquivalence(ctx); err == nil {
			b.plannedEquiv = equiv
		}
	}
	switch {
	case !b.prebuiltEligible():
		b.prebuilt = prebuiltPlan{choice: prebuiltNone}
	case known:
		cand, found, miss, err := b.findPrebuilt(ctx, 0)
		if err != nil {
			return err
		}
		if !found {
			if miss == missUnavailable {
				logging.FromContext(ctx).Warn(fmt.Sprintf("Cannot fetch prebuilt %s, building from the recipe", b.spec.Image.Name))
			}
			b.prebuilt = prebuiltPlan{choice: prebuiltNone}
			return nil
		}
		b.prebuilt = prebuiltPlan{choice: prebuiltPull, candidate: cand}
	}
	return nil
}

// tryPrebuilt installs the planned prebuilt, or looks for one when planning left
// the node undecided, trying the selected recipe source's ordered pull
// endpoints. The caller already holds the target's producer lock, so pull uses
// the locked transport entry point and installs atomically into the final
// pathname. A prebuilt made from a different recipe is skipped with a warning.
func (b *BuildObject) tryPrebuilt(ctx context.Context) (prebuiltResult, error) {
	if !b.prebuiltEligible() || b.prebuilt.choice == prebuiltNone {
		return false, nil
	}

	log := logging.FromContext(ctx)
	cand, found := b.prebuilt.candidate, b.prebuilt.choice == prebuiltPull
	miss, unavailable := missNone, false
	for {
		if !found {
			var err error
			if cand, found, miss, err = b.findPrebuilt(ctx, cand.next); err != nil {
				return false, err
			}
			if !found {
				break
			}
		}
		// Said before the download rather than after it: everything above is
		// metadata, and the gigabytes start here. Without this the operator
		// watches a long transfer with nothing saying what is being fetched or
		// that it has already been checked against the local recipe.
		log.Info("prebuilt found and verified", "artifact", b.spec.Image.Name,
			"endpoint", cand.endpoint, "equivalence", describePrebuiltKey(cand.equiv))
		if err := pullPrebuilt(ctx, cand.endpoint, cand.repo, cand.desc, cand.annotations, b.tgt.Path,
			registry.KindFor(b.spec.Image.Type)); err != nil {
			switch {
			case errors.Is(err, registry.ErrNotFound), errors.Is(err, registry.ErrUnsupportedPlatform),
				closedDoor(ctx, err, cand.endpoint):
			case errors.Is(err, registry.ErrUnavailable):
				log.Debug("prebuilt unavailable", "name", b.spec.Image.Name, "err", err)
				unavailable = true
			default:
				return false, fmt.Errorf("cannot pull prebuilt %s from %s: %w", b.spec.Image.Name, cand.endpoint, err)
			}
			found = false
			continue
		}
		invalidateInstalledOverlays()
		log.Info("prebuilt image ready", "kind", "success", "path", b.tgt.Path, "endpoint", cand.endpoint)
		return true, nil
	}
	switch {
	case unavailable || miss == missUnavailable:
		log.Warn(fmt.Sprintf("Cannot fetch prebuilt %s, building from the recipe", b.spec.Image.Name))
	case miss == missMismatch:
		log.Info(fmt.Sprintf("Building %s from the recipe", b.spec.Image.Name))
	default:
		log.Info(fmt.Sprintf("No prebuilt %s available, building from the recipe", b.spec.Image.Name))
	}
	return false, nil
}

// closedDoor reports a refusal that leaves nothing to report and nothing to fix:
// the endpoint turned away a request carrying no credential at all.
//
// ErrUnauthorized is otherwise fatal on purpose — an expired token must not
// silently become a long rebuild. That reasoning needs a credential to have
// expired. With none in hand the same status means only that this endpoint is
// not one this user can pull from, which is the outcome ErrNotFound already
// falls back on, and a collection may legitimately list an endpoint that only
// some of its readers can reach.
func closedDoor(ctx context.Context, err error, endpoint string) bool {
	return errors.Is(err, registry.ErrUnauthorized) && !registry.HasCredential(ctx, endpoint)
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
