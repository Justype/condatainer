package key

import (
	"bytes"
	"fmt"
	"os"
	"path/filepath"

	"github.com/Justype/condatainer/catalog"
	"github.com/Justype/condatainer/internal/artifact/meta"
	"github.com/Justype/condatainer/internal/conda"
)

// Scheme identifies one immutable key derivation algorithm.
type Scheme string

const (
	ScriptIdentityV1     Scheme = "script-identity-v1"
	ScriptEquivV1        Scheme = "script-equiv-v1"
	DefinitionIdentityV1 Scheme = "definition-identity-v1"
	DefinitionEquivV1    Scheme = "definition-equiv-v1"
	CondaExplicitV1      Scheme = "conda-explicit-v1"
	CondaEnvironmentV1   Scheme = "conda-environment-v1"
)

// Sources are the rebuild-source files embedded beside manifest.json.
type Sources map[string][]byte

// Value is one derived key and the canonical bytes it hashes.
type Value struct {
	Ref      meta.KeyRef
	Preimage []byte
}

// Derived is both keys regenerated from a manifest and its rebuild sources.
// IdentityModel is present for recipe-backed artifacts so comparison can report
// semantic differences and validate runtime environment contributions.
type Derived struct {
	Identity      Value
	Equiv         Value
	IdentityModel *Model
}

// Keys returns the manifest representation of these derived values.
func (d Derived) Keys() meta.Keys {
	return meta.Keys{Identity: d.Identity.Ref, Equiv: d.Equiv.Ref}
}

// Generate derives keys using the current scheme pair for the manifest's build
// type. It is used by new builds.
func Generate(m meta.Manifest, sources Sources) (Derived, error) {
	identity, equiv, err := latest(m.BuildType)
	if err != nil {
		return Derived{}, err
	}
	return derive(m, sources, identity, equiv)
}

// Regenerate derives both keys using the schemes m names, without checking the
// digests recorded beside them.
//
// This is what a rebuild against an existing record must use. Generate picks the
// current scheme pair, so after a scheme version ships it would stamp a key that
// cannot be related to the old one at all: comparison refuses to compare across
// schemes before it looks at any content, so the result would read as a
// difference in the scheme rather than as anything about the artifact.
func Regenerate(m meta.Manifest, sources Sources) (Derived, error) {
	if err := meta.ValidateManifest(m); err != nil {
		return Derived{}, err
	}
	if IsSnapshot(m) {
		return Derived{}, fmt.Errorf("%w: %s is a frozen environment; its identity is hashed from the packed payload, which only the artifact holds", ErrUnkeyed, m.Name)
	}
	if m.Keys.Identity.Empty() {
		return Derived{}, fmt.Errorf("artifact records no keys")
	}
	return derive(m, sources, Scheme(m.Keys.Identity.Scheme), Scheme(m.Keys.Equiv.Scheme))
}

// Verify regenerates the scheme-backed keys named by m and checks their stored
// digests. Manifests from the old file-backed format have no schemes and fail.
//
// A snapshot's keys are hashed from the packed payload, which sources do not
// hold, so its recorded keys are returned as they stand.
func Verify(m meta.Manifest, sources Sources) (Derived, error) {
	if IsSnapshot(m) {
		return snapshotDerived(m)
	}
	identity := Scheme(m.Keys.Identity.Scheme)
	equiv := Scheme(m.Keys.Equiv.Scheme)
	d, err := Regenerate(m, sources)
	if err != nil {
		return Derived{}, err
	}
	if d.Identity.Ref.SHA256 != m.Keys.Identity.SHA256 {
		return Derived{}, fmt.Errorf("%s derives identity %s, manifest says %s",
			identity, short(d.Identity.Ref.SHA256), short(m.Keys.Identity.SHA256))
	}
	if d.Equiv.Ref.SHA256 != m.Keys.Equiv.SHA256 {
		return Derived{}, fmt.Errorf("%s derives equivalence %s, manifest says %s",
			equiv, short(d.Equiv.Ref.SHA256), short(m.Keys.Equiv.SHA256))
	}
	return d, nil
}

// VerifyDir reads the source files required by m from dir and verifies its keys.
func VerifyDir(dir string, m meta.Manifest) (Derived, error) {
	if IsSnapshot(m) {
		return Verify(m, nil)
	}
	sources, err := ReadSources(dir, m)
	if err != nil {
		return Derived{}, err
	}
	return Verify(m, sources)
}

// snapshotDerived returns a snapshot's recorded keys as derived values, with no
// preimage, nothing having been computed.
func snapshotDerived(m meta.Manifest) (Derived, error) {
	if m.Keys.Identity.Empty() {
		return Derived{}, fmt.Errorf("%w: %s records no identity", ErrUnkeyed, m.Name)
	}
	return Derived{
		Identity: Value{Ref: m.Keys.Identity},
		Equiv:    Value{Ref: m.Keys.Equiv},
	}, nil
}

// ReadSources reads the fixed source set required by a manifest's build type.
func ReadSources(dir string, m meta.Manifest) (Sources, error) {
	var names []string
	switch m.BuildType {
	case meta.BuildTypeScript, meta.BuildTypeDef:
		names = []string{meta.RecipeFileName}
	case meta.BuildTypeConda:
		names = []string{conda.ExplicitFileName, conda.EnvironmentFileName}
	case meta.BuildTypeSnapshot:
		return nil, fmt.Errorf("%w: %s is a frozen environment", ErrUnkeyed, m.Name)
	default:
		return nil, fmt.Errorf("unknown build type %q", m.BuildType)
	}

	sources := make(Sources, len(names))
	for _, name := range names {
		data, err := os.ReadFile(filepath.Join(dir, name))
		if err != nil {
			return nil, fmt.Errorf("manifest requires %s: %w", name, err)
		}
		sources[name] = data
	}
	return sources, nil
}

func latest(buildType meta.BuildType) (Scheme, Scheme, error) {
	switch buildType {
	case meta.BuildTypeScript:
		return ScriptIdentityV1, ScriptEquivV1, nil
	case meta.BuildTypeDef:
		return DefinitionIdentityV1, DefinitionEquivV1, nil
	case meta.BuildTypeConda:
		return CondaExplicitV1, CondaEnvironmentV1, nil
	case meta.BuildTypeSnapshot:
		return "", "", ErrUnkeyed
	default:
		return "", "", fmt.Errorf("unknown build type %q", buildType)
	}
}

func derive(m meta.Manifest, sources Sources, identity, equiv Scheme) (Derived, error) {
	switch identity {
	case ScriptIdentityV1:
		if err := requireSchemePair(m, meta.BuildTypeScript, ScriptIdentityV1, ScriptEquivV1, equiv); err != nil {
			return Derived{}, err
		}
		if err := validateSourceNames(m); err != nil {
			return Derived{}, err
		}
		return deriveScriptV1(m, sources)

	case DefinitionIdentityV1:
		if err := requireSchemePair(m, meta.BuildTypeDef, DefinitionIdentityV1, DefinitionEquivV1, equiv); err != nil {
			return Derived{}, err
		}
		if err := validateSourceNames(m); err != nil {
			return Derived{}, err
		}
		return deriveDefinitionV1(m, sources)

	case CondaExplicitV1:
		if err := requireSchemePair(m, meta.BuildTypeConda, CondaExplicitV1, CondaEnvironmentV1, equiv); err != nil {
			return Derived{}, err
		}
		if err := validateSourceNames(m); err != nil {
			return Derived{}, err
		}
		return deriveCondaV1(m, sources)

	default:
		return Derived{}, fmt.Errorf("unknown identity scheme %q", identity)
	}
}

func validateSourceNames(m meta.Manifest) error {
	var required []string
	switch m.BuildType {
	case meta.BuildTypeScript, meta.BuildTypeDef:
		required = []string{meta.RecipeFileName}
	case meta.BuildTypeConda:
		required = []string{conda.ExplicitFileName, conda.EnvironmentFileName}
	default:
		return fmt.Errorf("unknown build type %q", m.BuildType)
	}
	if len(m.Source.Files) != len(required) {
		return fmt.Errorf("build type %s requires source files %v, manifest names %v",
			m.BuildType, required, m.Source.Files)
	}
	seen := make(map[string]bool, len(m.Source.Files))
	for _, name := range m.Source.Files {
		seen[name] = true
	}
	for _, name := range required {
		if !seen[name] {
			return fmt.Errorf("build type %s requires source file %s", m.BuildType, name)
		}
	}
	return nil
}

func requireSchemePair(m meta.Manifest, buildType meta.BuildType, identity, wantEquiv, gotEquiv Scheme) error {
	if m.BuildType != buildType {
		return fmt.Errorf("identity scheme %s requires build type %s, got %s", identity, buildType, m.BuildType)
	}
	if gotEquiv != wantEquiv {
		return fmt.Errorf("identity scheme %s requires equivalence scheme %s, got %s",
			identity, wantEquiv, gotEquiv)
	}
	return nil
}

func deriveScriptV1(m meta.Manifest, sources Sources) (Derived, error) {
	artifact, err := recipeArtifact(m, sources)
	if err != nil {
		return Derived{}, err
	}
	identity, identityModel, err := deriveScriptIdentityV1(artifact)
	if err != nil {
		return Derived{}, err
	}
	equiv, _, err := deriveScriptEquivV1(artifact)
	if err != nil {
		return Derived{}, err
	}
	return Derived{Identity: identity, Equiv: equiv, IdentityModel: &identityModel}, nil
}

func deriveDefinitionV1(m meta.Manifest, sources Sources) (Derived, error) {
	artifact, err := recipeArtifact(m, sources)
	if err != nil {
		return Derived{}, err
	}
	identity, identityModel, err := deriveDefinitionIdentityV1(artifact)
	if err != nil {
		return Derived{}, err
	}
	equiv, _, err := deriveDefinitionEquivV1(artifact)
	if err != nil {
		return Derived{}, err
	}
	return Derived{Identity: identity, Equiv: equiv, IdentityModel: &identityModel}, nil
}

func deriveCondaV1(m meta.Manifest, sources Sources) (Derived, error) {
	explicit, ok := sources[conda.ExplicitFileName]
	if !ok {
		return Derived{}, fmt.Errorf("missing %s", conda.ExplicitFileName)
	}
	environment, ok := sources[conda.EnvironmentFileName]
	if !ok {
		return Derived{}, fmt.Errorf("missing %s", conda.EnvironmentFileName)
	}
	identity, err := deriveCondaExplicitV1(m, explicit)
	if err != nil {
		return Derived{}, err
	}
	equiv, err := deriveCondaEnvironmentV1(m, environment)
	if err != nil {
		return Derived{}, err
	}
	return Derived{Identity: identity, Equiv: equiv}, nil
}

// recipeArtifact reconstructs neutral recipe inputs from the stored source and
// manifest. It does not decide which inputs enter a key; each scheme does that.
func recipeArtifact(m meta.Manifest, sources Sources) (Artifact, error) {
	recipeData, ok := sources[meta.RecipeFileName]
	if !ok {
		return Artifact{}, fmt.Errorf("missing %s", meta.RecipeFileName)
	}
	recipe, err := catalog.ParseRecipe(meta.RecipeFileName, bytes.NewReader(recipeData))
	if err != nil {
		return Artifact{}, fmt.Errorf("cannot parse %s: %w", meta.RecipeFileName, err)
	}
	if len(m.Source.Placeholders) > 0 {
		recipe, err = catalog.Expand(recipe, m.Source.Placeholders)
		if err != nil {
			return Artifact{}, fmt.Errorf("cannot apply recorded placeholders: %w", err)
		}
	}

	artifact := Artifact{
		Name:         m.Name,
		Type:         m.Type,
		Recipe:       recipeData,
		Placeholders: m.Source.Placeholders,
		From:         upstreamDigest(m),
	}
	for _, env := range recipe.Env {
		artifact.Env = append(artifact.Env, meta.EnvVar{
			Key: env.Key, Value: env.Value(nil), Note: env.Note,
		})
	}
	for _, dep := range m.Dependencies {
		if err := validateRole(dep); err != nil {
			return Artifact{}, err
		}
		artifact.Deps = append(artifact.Deps, Dep{
			Name: dep.Name, Type: dep.Type, Identity: dep.Identity, Equiv: dep.Equiv, Role: dep.Role,
		})
	}
	return artifact, nil
}

func upstreamDigest(m meta.Manifest) string {
	if m.BuildType == meta.BuildTypeDef && m.Build.From != nil {
		return m.Build.From.Digest
	}
	return ""
}

func validateRole(dep meta.Dependency) error {
	switch dep.Type {
	case catalog.TypeData:
		if dep.Role != meta.RoleData {
			return fmt.Errorf("data dependency %s has role %q", dep.Name, dep.Role)
		}
	case catalog.TypeApp, catalog.TypeOS:
		if dep.Role != meta.RoleApp && dep.Role != meta.RoleHistory {
			return fmt.Errorf("%s dependency %s has role %q", dep.Type, dep.Name, dep.Role)
		}
	case catalog.TypeBase:
		return fmt.Errorf("base %s cannot be a dependency", dep.Name)
	default:
		return fmt.Errorf("dependency %s has unknown type %q", dep.Name, dep.Type)
	}
	return nil
}

func value(scheme Scheme, preimage []byte) Value {
	return Value{
		Ref:      meta.KeyRef{Scheme: string(scheme), SHA256: Sum(preimage)},
		Preimage: preimage,
	}
}

func short(digest string) string {
	if len(digest) > 12 {
		return digest[:12]
	}
	return digest
}
