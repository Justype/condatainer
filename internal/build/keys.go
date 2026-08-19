package build

import (
	"context"
	"fmt"
	"strings"

	"github.com/Justype/condatainer/catalog"
	"github.com/Justype/condatainer/internal/artifact/capsule"
	"github.com/Justype/condatainer/internal/artifact/key"
	"github.com/Justype/condatainer/internal/artifact/meta"
	"github.com/Justype/condatainer/internal/logging"
	"github.com/Justype/condatainer/internal/runtime/container"
	"github.com/Justype/condatainer/internal/utils"
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

// deriveRecipeKeys freezes dependency adjacency, derives both keys from the
// manifest and recipe, and composes the rebuild-source capsule.
func (b *BuildObject) deriveRecipeKeys(ctx context.Context) error {
	recipe, ok := b.spec.Source.RecipeFile()
	if !ok {
		return nil
	}

	artifact := key.Artifact{
		Name:         b.spec.Image.Name,
		Type:         b.spec.Image.Type,
		Env:          b.spec.Image.Env,
		Recipe:       recipe.Data,
		Placeholders: b.spec.Source.Placeholders,
		Deps:         b.dependencyKeys(ctx),
		From:         b.spec.Source.UpstreamDigest(),
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
// string that located it. A #DEP: may be an overlay path, and a path is
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
		requested := raw
		if parsed, err := catalog.ParseDep(raw); err == nil {
			requested = parsed.NameVersion()
		} else if !utils.IsOverlay(raw) {
			log.Warn("skipping an unparsable dependency", "dep", raw, "err", err)
			continue
		}
		dep := key.Dep{Name: requested, Type: catalog.TypeApp}

		path, manifest, err := readDependencyManifest(requested)
		if err != nil {
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
