package build

import (
	"context"
	"errors"
	"fmt"
	"strings"

	"github.com/Justype/condatainer/catalog"
	"github.com/Justype/condatainer/internal/artifact/capsule"
	"github.com/Justype/condatainer/internal/artifact/key"
	"github.com/Justype/condatainer/internal/artifact/meta"
	"github.com/Justype/condatainer/internal/logging"
	"github.com/Justype/condatainer/internal/runtime/container"
	"github.com/Justype/condatainer/internal/store"
)

// generateKeys derives both keys for what this build produced.
//
// A locked rebuild derives with the schemes its lock recorded rather than the
// current pair. A scheme version that shipped since the lock was written would
// otherwise stamp a key nothing can relate to the recorded one — comparison
// rejects two different schemes before looking at content — so every existing
// lock would fail to restore for a reason that says nothing about the artifact.
func (b *BuildObject) generateKeys() (key.Derived, error) {
	manifest := b.Manifest()
	if b.locked && !b.lockedKeys.Identity.Empty() {
		manifest.Keys = b.lockedKeys
		return key.Regenerate(manifest, b.keySources())
	}
	return key.Generate(manifest, b.keySources())
}

// deriveKeys computes the scheme-backed identity and equivalence immediately
// before staging, when every source and dependency input is known.
func (b *BuildObject) deriveKeys(ctx context.Context) error {
	if b.spec.Source.Conda != nil {
		b.keysFromSources(ctx)
		return nil
	}
	return b.deriveRecipeKeys(ctx)
}

// keysFromSources derives Conda keys from the two captured canonical exports. A
// failed capture leaves both keys absent, preserving the existing unrecorded
// behavior rather than making a half-claim.
func (b *BuildObject) keysFromSources(ctx context.Context) {
	b.keys = meta.Keys{}
	derived, err := b.generateKeys()
	if err != nil {
		logging.FromContext(ctx).Warn("conda environment recorded without keys",
			"name", b.spec.Image.Name, "reason", err)
		return
	}
	b.keys = derived.Keys()
}

// recipeKeyArtifact assembles the key inputs of a recipe-backed build.
//
// One spelling, used both to derive the real keys after the build and to predict
// the identity before it. A second construction of the same inputs is a second
// preimage, and two preimages that drift produce two identities for one artifact.
func (b *BuildObject) recipeKeyArtifact(ctx context.Context) (key.Artifact, bool) {
	recipe, ok := b.spec.Source.RecipeFile()
	if !ok {
		return key.Artifact{}, false
	}
	return key.Artifact{
		Name:         b.spec.Image.Name,
		Type:         b.spec.Image.Type,
		Env:          b.spec.Image.Env,
		Recipe:       recipe.Data,
		Placeholders: b.spec.Source.Placeholders,
		Deps:         b.dependencyKeys(ctx),
		From:         b.spec.Source.UpstreamDigest(),
		Fetched:      b.spec.Source.Fetched,
	}, true
}

// deriveRecipeKeys freezes dependency adjacency, derives both keys from the
// manifest and recipe, and composes the rebuild-source capsule.
func (b *BuildObject) deriveRecipeKeys(ctx context.Context) error {
	artifact, ok := b.recipeKeyArtifact(ctx)
	if !ok {
		return nil
	}

	b.dependencies, b.provenanceComplete = key.Manifest(artifact)
	b.keys = meta.Keys{}
	derived, err := b.generateKeys()
	if err != nil {
		return fmt.Errorf("refusing to pack %s: %w", b.spec.Image.Name, err)
	}
	b.keys = derived.Keys()
	return b.composeCapsule(ctx)
}

func (b *BuildObject) keySources() key.Sources {
	sources := make(key.Sources)
	for _, file := range b.embedded {
		if file.isSource {
			sources[file.Name] = file.Data
		}
	}
	return sources
}

// composeCapsule embeds the closure this artifact was built from, by union from
// its dependencies' own images. Only data has dependencies, so only data has one.
//
// Completeness is the manifest's, corrected by what the capsule found: an
// unrecorded dependency anywhere below makes everything above it incomplete.
func (b *BuildObject) composeCapsule(ctx context.Context) error {
	if len(b.dependencies) == 0 {
		return nil
	}

	deps := make([]capsule.Dep, 0, len(b.dependencies))
	for _, dep := range b.dependencies {
		entry := capsule.Dep{Name: dep.Name, Identity: dep.Identity}
		if !dep.Identity.Empty() {
			// The path dependencyKeys already read this dependency's manifest
			// from. Resolving the name again could land on a different image
			// than the one whose keys were just recorded, and a locked rebuild
			// has no name resolution to fall back on at all.
			entry.ImagePath = b.depImagePaths[dep.Name]
			if entry.ImagePath == "" {
				logging.FromContext(ctx).Warn("dependency vanished before its provenance was read", "dep", dep.Name)
				continue
			}
		}
		deps = append(deps, entry)
	}

	complete, err := capsule.Compose(b.ws.MetaDir, deps)
	if err != nil {
		return err
	}
	if b.provenanceComplete != nil && *b.provenanceComplete {
		b.provenanceComplete = &complete
	}
	return nil
}

// dependencyKeys reads what each direct dependency's image says about itself.
//
// A dependency that cannot be read contributes its name and type alone, marked
// unrecorded. A manifest from the removed file-backed format has no schemes and
// follows this same path.
//
// An edge records the name the dependency's own manifest declares, never the
// string that located it, and reads the image the build mounts: for a constrained
// #DEP: the installed version that satisfies it, not the preferred one. A #DEP: may be an overlay path, and a path is
// machine-local: recorded as a name it would reach the manifest, the capsule
// entry directory, and any published copy, and it could never match the child
// it points at.
func (b *BuildObject) dependencyKeys(ctx context.Context) []key.Dep {
	if len(b.spec.Dependencies) == 0 {
		return nil
	}
	log := logging.FromContext(ctx)

	b.depImagePaths = make(map[string]string, len(b.spec.Dependencies))
	out := make([]key.Dep, 0, len(b.spec.Dependencies))
	for _, raw := range b.spec.Dependencies {
		// A path is recognized before it is parsed, never after. ParseDep
		// succeeds on an overlay path and normalizes it into a name — turning
		// .../hello--1.0.sqf into ".../hello/1.0.sqf", which names nothing, reads
		// no manifest, and records an edge with no keys at all. A locked rebuild
		// supplies every dependency as a path, so this is its whole edge set.
		requested, constrained := raw, raw
		if !catalog.IsPathDep(raw) {
			parsed, err := catalog.ParseDep(raw)
			if err != nil {
				log.Warn("skipping an unparsable dependency", "dep", raw, "err", err)
				continue
			}
			requested, constrained = parsed.NameVersion(), parsed.String()
		}
		dep := key.Dep{Name: requested, Type: catalog.TypeApp}

		// Read from what the build mounts: a constraint picks the highest installed
		// version it admits, so the edge names that one and carries its keys.
		path, manifest, err := readDependencyManifest(constrained)
		if err != nil {
			if planned, ok := b.plannedDeps[requested]; ok {
				out = append(out, key.Dep{Name: planned.Name, Type: planned.Type, Equiv: planned.Equiv})
				continue
			}
			log.Debug("dependency carries no scheme-backed keys", "dep", requested, "err", err)
			out = append(out, dep)
			continue
		}
		if manifest.Name != "" {
			dep.Name = manifest.Name
		}
		dep.Type = manifest.Type
		dep.Identity = manifest.Keys.Identity
		dep.Equiv = manifest.Keys.Equiv
		b.depImagePaths[dep.Name] = path
		out = append(out, dep)
	}
	return out
}

// readDependencyManifest resolves a dependency to the image providing it and
// reads what that image records about itself, reporting both. The path is what
// the capsule is then composed from, so provenance is read out of the same
// image the keys came from.
func readDependencyManifest(nameVersion string) (string, meta.Manifest, error) {
	paths, err := container.ResolveOverlayPaths([]string{nameVersion})
	if err != nil {
		return "", meta.Manifest{}, err
	}
	if len(paths) == 0 {
		return "", meta.Manifest{}, fmt.Errorf("no image provides %s", nameVersion)
	}
	path := strings.TrimSuffix(strings.TrimSuffix(paths[0], ":ro"), ":rw")
	manifest, err := meta.ReadManifest(path)
	return path, manifest, err
}

// ErrNoPrediction reports that this build's identity cannot be known before it
// runs. It is never a failure: the caller builds, and the identity is derived
// from what was produced, as it always was.
var ErrNoPrediction = errors.New("identity cannot be predicted before the build")

// PredictIdentity computes the identity this build will produce, without
// producing it.
//
// It matters because --store files by identity, so "is this already installed"
// is an identity question, not a name question — and answering it after the
// build has run costs the whole build. Every input to a recipe-backed key is
// known beforehand: the recipe text is embedded at resolution, placeholders and
// environment are fixed, dependency edges settle once dependencies are resolved,
// and #SOURCE: digests once the sources are fetched. Call it after both, or the
// answer describes a different artifact.
//
// A prediction is only ever acted on when it *matches* something installed, so
// the direction of any error matters: an identity that cannot be predicted, or
// is predicted wrongly, costs a build that the publication step then adopts.
func (b *BuildObject) PredictIdentity(ctx context.Context) (meta.KeyRef, error) {
	if b.spec.Source.Conda != nil {
		return b.predictCondaIdentity(ctx)
	}
	artifact, ok := b.recipeKeyArtifact(ctx)
	if !ok {
		return meta.KeyRef{}, ErrNoPrediction
	}

	// A copy, never b's own fields: prediction must leave the build exactly as it
	// found it, and deriveRecipeKeys sets these for real once the build is done.
	manifest := b.Manifest()
	manifest.Dependencies, manifest.ProvenanceComplete = key.Manifest(artifact)
	manifest.Keys = meta.Keys{}

	derived, err := key.Generate(manifest, b.keySources())
	if err != nil {
		return meta.KeyRef{}, fmt.Errorf("%w: %v", ErrNoPrediction, err)
	}
	return derived.Identity.Ref, nil
}

// InstalledAt reports where an identity of this build's name is already
// installed, flat or in a store, across every readable images directory.
func (b *BuildObject) InstalledAt(identity meta.KeyRef) (string, bool) {
	if identity.Empty() {
		return "", false
	}
	candidate, _, err := store.ResolveIdentity(b.spec.Image.Name,
		store.IdentityQuery{Scheme: identity.Scheme, SHA256: identity.SHA256}, nil)
	if err != nil {
		return "", false
	}
	return candidate.Path, true
}

// skipIfInstalled reports whether this build's product is already installed,
// deciding by identity rather than by name.
//
// Only under storeOverflow. Without it the bare name is the address, and the
// name check has already answered; with it the name is expected to be taken and
// says nothing about whether this exact build exists.
func (b *BuildObject) skipIfInstalled(ctx context.Context) bool {
	if !b.storeOverflow {
		return false
	}
	identity, err := b.PredictIdentity(ctx)
	if err != nil {
		logging.FromContext(ctx).Debug("building without a predicted identity",
			"name", b.spec.Image.Name, "reason", err)
		return false
	}
	path, found := b.InstalledAt(identity)
	if !found {
		return false
	}
	logging.FromContext(ctx).Info("this exact build is already installed, skipping",
		"kind", "note", "name", b.spec.Image.Name, "identity", identity.Digest(), "path", path)
	return true
}
