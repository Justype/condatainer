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
)

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
	derived, err := key.Generate(b.Manifest(), b.keySources())
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
	derived, err := key.Generate(b.Manifest(), b.keySources())
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
			paths, err := container.ResolveOverlayPaths([]string{dep.Name})
			if err != nil || len(paths) == 0 {
				logging.FromContext(ctx).Warn("dependency vanished before its provenance was read", "dep", dep.Name)
				continue
			}
			entry.ImagePath = strings.TrimSuffix(strings.TrimSuffix(paths[0], ":ro"), ":rw")
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
func (b *BuildObject) dependencyKeys(ctx context.Context) []key.Dep {
	if len(b.spec.Dependencies) == 0 {
		return nil
	}
	log := logging.FromContext(ctx)

	out := make([]key.Dep, 0, len(b.spec.Dependencies))
	for _, raw := range b.spec.Dependencies {
		parsed, err := catalog.ParseDep(raw)
		if err != nil {
			log.Warn("skipping an unparsable dependency", "dep", raw, "err", err)
			continue
		}
		dep := key.Dep{Name: parsed.NameVersion(), Type: catalog.TypeApp}

		manifest, err := readDependencyManifest(dep.Name)
		if err != nil {
			log.Debug("dependency carries no scheme-backed keys", "dep", dep.Name, "err", err)
			out = append(out, dep)
			continue
		}
		dep.Type = manifest.Type
		dep.Identity = manifest.Keys.Identity
		dep.Equiv = manifest.Keys.Equiv
		out = append(out, dep)
	}
	return out
}

// readDependencyManifest resolves a dependency to the image providing it and
// reads what that image records about itself.
func readDependencyManifest(nameVersion string) (meta.Manifest, error) {
	paths, err := container.ResolveOverlayPaths([]string{nameVersion})
	if err != nil {
		return meta.Manifest{}, err
	}
	if len(paths) == 0 {
		return meta.Manifest{}, fmt.Errorf("no image provides %s", nameVersion)
	}
	path := strings.TrimSuffix(strings.TrimSuffix(paths[0], ":ro"), ":rw")
	return meta.ReadManifest(path)
}
