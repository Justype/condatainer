package build

import (
	"context"
	"fmt"
	"strings"

	"github.com/Justype/condatainer/catalog"
	"github.com/Justype/condatainer/internal/artifact/capsule"
	"github.com/Justype/condatainer/internal/artifact/key"
	"github.com/Justype/condatainer/internal/artifact/meta"
	"github.com/Justype/condatainer/internal/artifact/record"
	"github.com/Justype/condatainer/internal/conda"
	"github.com/Justype/condatainer/internal/logging"
	"github.com/Justype/condatainer/internal/runtime/container"
)

// recordKeys computes what the image is addressed and compared by, and embeds
// whatever backs it. Every backend calls this immediately before staging, which
// is the last point at which everything a key depends on is known.
//
// A Conda app gets keys pointing at the exports it already embeds, since a
// record for it would hold only a format tag, a type line and one digest.
// Everything else gets identity.record and equiv.record — a base included, which
// nothing may depend on but anything may identify.
func (b *BuildObject) recordKeys(ctx context.Context) error {
	if b.spec.Source.Conda != nil {
		b.keysFromExports(ctx)
		return nil
	}
	return b.recordRecipeKeys(ctx)
}

// keysFromExports names the Conda exports as the keys themselves. A capture that
// failed leaves no keys rather than a claim about a file that is not there.
func (b *BuildObject) keysFromExports(ctx context.Context) {
	var keys meta.Keys
	for _, file := range b.embedded {
		switch file.Name {
		case conda.ExplicitFileName:
			keys.Identity = meta.KeyRef{SHA256: record.Sum(file.Data), File: file.Name}
		case conda.EnvironmentFileName:
			keys.Equiv = meta.KeyRef{SHA256: record.Sum(file.Data), File: file.Name}
		}
	}
	if keys.Identity.File == "" || keys.Equiv.File == "" {
		logging.FromContext(ctx).Warn("conda environment recorded without keys",
			"name", b.spec.Image.Name, "reason", "an export was not captured")
		return
	}
	b.keys = keys
}

// recordRecipeKeys builds both records for a recipe build and embeds them.
func (b *BuildObject) recordRecipeKeys(ctx context.Context) error {
	recipe, ok := b.spec.Source.RecipeFile()
	if !ok {
		return nil // nothing was resolved from a recipe; nothing to key
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

	identity, equiv, err := key.Records(artifact)
	if err != nil {
		return fmt.Errorf("refusing to pack %s: %w", b.spec.Image.Name, err)
	}
	identityBytes, err := record.Marshal(identity)
	if err != nil {
		return err
	}
	equivBytes, err := record.Marshal(equiv)
	if err != nil {
		return err
	}

	b.embedRecord(meta.IdentityFileName, identityBytes)
	b.embedRecord(meta.EquivFileName, equivBytes)
	b.keys = meta.Keys{
		Identity: meta.KeyRef{SHA256: record.Sum(identityBytes), File: meta.IdentityFileName},
		Equiv:    meta.KeyRef{SHA256: record.Sum(equivBytes), File: meta.EquivFileName},
	}
	b.dependencies, b.provenanceComplete = key.Manifest(artifact)
	return b.composeCapsule(ctx)
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
		if dep.Identity != "" {
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
// unrecorded: every image built before this format is that case, and failing a
// build over it would strand every site that has not rebuilt yet.
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
			log.Debug("dependency carries no records", "dep", dep.Name, "err", err)
			out = append(out, dep)
			continue
		}
		dep.Type = manifest.Type
		dep.Identity = manifest.Keys.Identity.Digest()
		dep.Equiv = manifest.Keys.Equiv.Digest()
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
