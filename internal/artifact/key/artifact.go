package key

import (
	"github.com/Justype/condatainer/catalog"
	"github.com/Justype/condatainer/internal/artifact/meta"
)

// Artifact contains the resolved inputs from which recipe-backed schemes derive
// keys. It is not a canonical model: each scheme explicitly chooses which of
// these fields contribute to its own preimage.
type Artifact struct {
	// Name never enters a preimage. ScriptEquivV1 uses it only to classify
	// dependency roles.
	Name         string
	Type         catalog.Type
	Env          []meta.EnvVar
	Recipe       []byte
	Placeholders map[string]string
	Deps         []Dep
	// From is the resolved upstream image digest for a definition, or
	// meta.Unrecorded when resolution failed. It is empty without an upstream.
	From string
}

// Dep is one direct build dependency and the keys advertised by its image.
type Dep struct {
	Name     string
	Type     catalog.Type
	Identity string
	Equiv    string
	// Role is empty while a build is deciding policy and frozen in a manifest.
	Role string
}

// Recorded reports whether both keys needed to describe a dependency exist.
func (d Dep) Recorded() bool { return d.Identity != "" && d.Equiv != "" }

// Manifest freezes every direct dependency and the role selected by the current
// script equivalence scheme. This is explanatory provenance; scheme files
// independently decide which values enter their canonical preimages.
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
