package key

import (
	"fmt"

	"github.com/Justype/condatainer/catalog"
	"github.com/Justype/condatainer/internal/artifact/meta"
)

// ScriptEquivV1 answers: "Can this script-built artifact substitute for the
// requested artifact?"
//
// Preimage, in canonical order:
//   - artifact type
//   - every #ENV value
//   - comment-stripped recipe digest
//   - every selected placeholder
//   - named app/OS dependencies as name/version only
//   - data dependencies by equivalence key only
//
// History-only app/OS dependencies are omitted. An unkeyed data dependency uses
// "unrecorded" plus its name. Exact dependency identities never enter.
func deriveScriptEquivV1(a Artifact) (Value, Model, error) {
	if a.Type != catalog.TypeApp && a.Type != catalog.TypeData {
		return Value{}, Model{}, fmt.Errorf("%s cannot derive type %s", ScriptEquivV1, a.Type)
	}
	if a.Type != catalog.TypeData && len(a.Deps) > 0 {
		return Value{}, Model{}, fmt.Errorf("%s %s cannot have dependencies", ScriptEquivV1, a.Type)
	}

	model := Model{
		Kind:         KindEquiv,
		Type:         a.Type,
		Env:          Env(a.Env),
		Recipe:       RecipeDigest(a.Recipe),
		Placeholders: Placeholders(a.Placeholders),
		Deps:         scriptEquivDependenciesV1(a.Name, a.Deps),
	}
	derived, err := valueFromModel(ScriptEquivV1, model)
	return derived, model, err
}

func scriptEquivDependenciesV1(name string, deps []Dep) []DependencyValue {
	var out []DependencyValue
	for _, dep := range deps {
		role := dep.Role
		if role == "" {
			role = Role(name, dep)
		}
		switch role {
		case meta.RoleData:
			if dep.Recorded() {
				out = append(out, DependencyValue{
					Type:   catalog.TypeData,
					Fields: []string{dep.Equiv.Scheme, dep.Equiv.Digest()},
				})
			} else {
				out = append(out, DependencyValue{
					Type: catalog.TypeData, Fields: []string{meta.Unrecorded, dep.Name},
				})
			}
		case meta.RoleApp:
			out = append(out, DependencyValue{
				Type: dep.Type, Fields: []string{dep.Name},
			})
		case meta.RoleHistory:
			// Mounted during the build but irrelevant to substitution.
		}
	}
	return out
}

// Role classifies a dependency according to ScriptEquivV1.
//
// Data always contributes by content. An app or OS contributes by name/version
// only when its complete component sequence occurs in the artifact name;
// otherwise it is build history and is omitted from equivalence.
func Role(name string, dep Dep) string {
	switch {
	case dep.Type == catalog.TypeData:
		return meta.RoleData
	case catalog.HasComponents(name, dep.Name):
		return meta.RoleApp
	default:
		return meta.RoleHistory
	}
}
