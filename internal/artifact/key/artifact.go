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
	// Fetched are the #SOURCE: inputs the build downloaded, by the name the
	// recipe gave each and the digest of the bytes that arrived. A recipe that
	// fetches without declaring a source contributes none, and its keys are what
	// they were before sources existed.
	Fetched []meta.SourceFile
}

// sourceValues renders fetched inputs as preimage lines. Only identity schemes
// call it: a re-cut upstream file makes a different build, but not one that
// stops substituting for what a recipe asks for.
func sourceValues(sources []meta.SourceFile) []SourceValue {
	if len(sources) == 0 {
		return nil
	}
	out := make([]SourceValue, 0, len(sources))
	for _, s := range sources {
		out = append(out, SourceValue{Name: s.Name, Digest: s.Digest()})
	}
	return out
}

// Dep is one direct build dependency and the keys advertised by its image.
type Dep struct {
	Name     string
	Type     catalog.Type
	Identity meta.KeyRef
	Equiv    meta.KeyRef
	// Role is empty while a build is deciding policy and frozen in a manifest.
	Role string
}

// Recorded reports whether both keys needed to describe a dependency exist,
// each complete: a scheme without a digest, or the reverse, describes nothing.
func (d Dep) Recorded() bool { return d.Identity.Digest() != "" && d.Equiv.Digest() != "" }

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
