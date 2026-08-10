package key

import (
	"fmt"

	"github.com/Justype/condatainer/catalog"
	"github.com/Justype/condatainer/internal/artifact/meta"
	"github.com/Justype/condatainer/internal/artifact/record"
)

// Artifact is everything a recipe build knows about itself when its records are
// written. It is the input to Records and to nothing else, so the rules below are
// stated once and read the same at build time and at comparison time.
type Artifact struct {
	// Name is the resolved full name. It never enters a record — it decides
	// which dependencies count toward equivalence, and nothing more.
	Name string
	Type catalog.Type
	Env  []meta.EnvVar
	// Recipe is the recipe as stored at /.cnt/recipe, tokens and all.
	Recipe       []byte
	Placeholders map[string]string
	Deps         []Dep
	// From is the digest of the upstream image a definition bootstrapped from, or
	// meta.Unrecorded when the registry could not be reached. Empty when there is
	// no upstream at all — a script recipe, or a bootstrap from scratch or disk.
	From string
}

// Dep is one direct build dependency and what its image says about itself. An
// empty Identity means the image carried no records — see meta.Unrecorded.
type Dep struct {
	Name     string // name/version, normalized
	Type     catalog.Type
	Identity string
	Equiv    string
}

// Recorded reports whether the dependency's image carried keys of its own.
//
// Both are required. A dependency pinned in one record and not the other would
// be a half-claim, and the pair is all-or-nothing in practice: a recipe build
// writes both records, a Conda app names both exports or neither, and a base has
// neither. An image offering one is treated as unrecorded, which is a weaker
// claim honestly stated rather than a record that cannot be written.
func (d Dep) Recorded() bool { return d.Identity != "" && d.Equiv != "" }

// Records builds an artifact's identity and equivalence records.
//
// Both carry the same type, env, recipe and placeholder lines. They diverge only
// at dependencies, which data alone has, and at the upstream digest, which a
// definition alone has — the same idea twice: identity pins what was actually
// used, equivalence keeps the looser contract the recipe asked for.
func Records(a Artifact) (identity, equiv record.Record, err error) {
	if len(a.Deps) > 0 && a.Type != catalog.TypeData {
		return record.Record{}, record.Record{}, fmt.Errorf("key: %s is type %s and cannot have dependencies", a.Name, a.Type)
	}

	shared := record.Record{
		Type:         a.Type,
		Env:          Env(a.Env),
		Recipe:       RecipeDigest(a.Recipe),
		Placeholders: Placeholders(a.Placeholders),
	}

	identity = shared
	identity.Kind = record.KindIdentity
	identity.From = a.From
	identity.Deps = identityDeps(a.Deps)

	equiv = shared
	equiv.Kind = record.KindEquiv
	equiv.Deps = equivDeps(a.Name, a.Deps)
	return identity, equiv, nil
}

// identityDeps records every direct dependency, pinned to the exact build that
// was mounted. Anything mounted while the recipe ran could have shaped the
// payload, so nothing is projected away here.
func identityDeps(deps []Dep) []record.Dep {
	if len(deps) == 0 {
		return nil
	}
	out := make([]record.Dep, 0, len(deps))
	for _, d := range deps {
		digest := d.Identity
		if !d.Recorded() {
			digest = meta.Unrecorded
		}
		out = append(out, record.Dep{Type: d.Type, Fields: []string{d.Name, digest}})
	}
	return out
}

// equivDeps records only what decides whether one artifact substitutes for
// another. Each line carries the thing that actually determines the result and
// drops the rest:
//
//   - a data dependency contributes its equivalence and no name: the content is
//     the contract, and equivalence composes transitively without walking,
//     because that dependency's own equiv already covers its important inputs.
//   - an app or OS dependency the artifact's name mentions contributes its
//     name/version and no digest: the version is the contract, so a Conda-built
//     and a script-built STAR of one version are interchangeable producers.
//   - an app or OS dependency the name does not mention contributes nothing. It
//     was mounted and is recorded in identity, but it does not decide
//     substitution.
func equivDeps(name string, deps []Dep) []record.Dep {
	var out []record.Dep
	for _, d := range deps {
		switch d.Type {
		case catalog.TypeData:
			if d.Recorded() {
				out = append(out, record.Dep{Type: d.Type, Fields: []string{d.Equiv}})
			} else {
				// No digest to stand for it, so say which one it was instead.
				out = append(out, record.Dep{Type: d.Type, Fields: []string{meta.Unrecorded, d.Name}})
			}
		default:
			// The naming convention is the contract: a tool whose version
			// changes the result belongs in the artifact's name.
			if catalog.HasComponents(name, d.Name) {
				out = append(out, record.Dep{Type: d.Type, Fields: []string{d.Name}})
			}
		}
	}
	return out
}

// Role reports how a dependency reached the artifact's equivalence, for the
// manifest to explain what the frozen record already contains.
func Role(name string, d Dep) string {
	switch {
	case d.Type == catalog.TypeData:
		return meta.RoleData
	case catalog.HasComponents(name, d.Name):
		return meta.RoleApp
	default:
		return meta.RoleHistory
	}
}

// Manifest renders the dependency list and completeness flag an artifact records
// beside its records. Every direct dependency appears, whatever its role: the
// list is what a diff reads to say *which* one moved.
func Manifest(a Artifact) (deps []meta.Dependency, complete *bool) {
	if len(a.Deps) == 0 {
		return nil, nil
	}
	whole := true
	for _, d := range a.Deps {
		entry := meta.Dependency{
			Name:     d.Name,
			Type:     d.Type,
			Identity: d.Identity,
			Equiv:    d.Equiv,
			Role:     Role(a.Name, d),
		}
		if !d.Recorded() {
			entry.Records = meta.Unrecorded
			whole = false
		}
		deps = append(deps, entry)
	}
	return deps, &whole
}
