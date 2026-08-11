package key

import (
	"fmt"

	"github.com/Justype/condatainer/catalog"
	"github.com/Justype/condatainer/internal/artifact/meta"
)

// ScriptIdentityV1 answers: "Which exact script build is this?"
//
// Preimage, in canonical order:
//   - artifact type
//   - every #ENV value
//   - comment-stripped recipe digest
//   - every selected placeholder
//   - every direct dependency's name and exact identity
//
// An unkeyed or half-keyed dependency contributes its name plus "unrecorded".
// The artifact name and dependency equivalence keys never enter this scheme.
func deriveScriptIdentityV1(a Artifact) (Value, Model, error) {
	if a.Type != catalog.TypeApp && a.Type != catalog.TypeData {
		return Value{}, Model{}, fmt.Errorf("%s cannot derive type %s", ScriptIdentityV1, a.Type)
	}
	if a.Type != catalog.TypeData && len(a.Deps) > 0 {
		return Value{}, Model{}, fmt.Errorf("%s %s cannot have dependencies", ScriptIdentityV1, a.Type)
	}

	model := Model{
		Kind:         KindIdentity,
		Type:         a.Type,
		Env:          Env(a.Env),
		Recipe:       RecipeDigest(a.Recipe),
		Placeholders: Placeholders(a.Placeholders),
		Deps:         scriptIdentityDependenciesV1(a.Deps),
	}
	derived, err := valueFromModel(ScriptIdentityV1, model)
	return derived, model, err
}

func scriptIdentityDependenciesV1(deps []Dep) []DependencyValue {
	if len(deps) == 0 {
		return nil
	}
	out := make([]DependencyValue, 0, len(deps))
	for _, dep := range deps {
		identity := dep.Identity
		if !dep.Recorded() {
			identity = meta.Unrecorded
		}
		out = append(out, DependencyValue{
			Type:   dep.Type,
			Fields: []string{dep.Name, identity},
		})
	}
	return out
}
