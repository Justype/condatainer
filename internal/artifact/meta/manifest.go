package meta

import (
	"encoding/json"
	"errors"
	"fmt"
	"strings"

	"github.com/Justype/condatainer/catalog"
	"github.com/Justype/condatainer/internal/image/tool"
)

// Path is where the manifest lives inside every image.
const Path = "/" + DirName + "/" + FileName

// FileName is the manifest's basename, for callers staging it into a directory.
const FileName = "manifest.json"

// Manifest is what an image records about what it is and where it came from. It
// is read on demand — by info, comparison, and restore — never on the mount
// path, so it is free to grow. There is no runtime block: that lives in
// runtime.json and nowhere else.
type Manifest struct {
	SchemaVersion int          `json:"schema_version"`
	Name          string       `json:"name"`
	Type          catalog.Type `json:"type"`       // base, os, app, data
	BuildType     string       `json:"build_type"` // conda, script, def
	Description   string       `json:"description,omitempty"`
	URL           string       `json:"url,omitempty"`
	Platform      Platform     `json:"platform"`
	Source        Source       `json:"source,omitzero"`
	Keys          Keys         `json:"keys,omitzero"`
	Dependencies  []Dependency `json:"dependencies,omitempty"`
	// ProvenanceComplete reports whether every dependency carried records of its
	// own. It is nil when the question does not arise — an artifact with no
	// dependencies is neither complete nor incomplete.
	ProvenanceComplete *bool `json:"provenance_complete,omitempty"`
	Build              Build `json:"build,omitzero"`
}

// Keys names the files an artifact's identity and equivalence are the digests
// of. Every build type with a source has them, a base included; a Conda app
// names its exports rather than records.
type Keys struct {
	Identity KeyRef `json:"identity,omitzero"`
	Equiv    KeyRef `json:"equiv,omitzero"`
}

// KeyRef is one key: the file it hashes, relative to /.cnt, and that hash.
//
// Naming the file is what lets a reader verify a key without knowing which build
// type produced it — a recipe build points at identity.record, a Conda app at
// explicit.txt, and `sha256sum` on the named file reproduces the value either way.
type KeyRef struct {
	SHA256 string `json:"sha256"`
	File   string `json:"file"`
}

// Digest renders the key the way a record writes one inline, sha256:<hex>, or
// empty when there is no key.
func (k KeyRef) Digest() string {
	if k.SHA256 == "" {
		return ""
	}
	return "sha256:" + k.SHA256
}

// Dependency is one direct build dependency, as the artifact recorded it. The
// list is adjacency only: following each manifest's own through the capsule
// reconstructs the graph without a second representation of it.
type Dependency struct {
	Name     string       `json:"name"`
	Type     catalog.Type `json:"type"`
	Identity string       `json:"identity,omitempty"`
	Equiv    string       `json:"equiv,omitempty"`
	// Records is "unrecorded" when the image that satisfied this dependency
	// carried no keys of its own, which is every image built before this format.
	Records string `json:"records,omitempty"`
	// Role explains what the frozen equiv.record already contains — data, app,
	// or history. It is never an input to recomputing it: a reader whose rules
	// disagree warns and trusts the record, because recomputing equivalence
	// under newer rules would silently rewrite the past.
	Role string `json:"role"`
}

// Roles a dependency can play in its dependent's equivalence.
const (
	RoleData    = "data"    // contributes its equivalence
	RoleApp     = "app"     // named in the artifact's name; contributes name/version
	RoleHistory = "history" // mounted, but decides nothing about substitution
)

// Unrecorded marks a dependency satisfied by an image carrying no records.
const Unrecorded = "unrecorded"

// The record files a recipe build embeds.
const (
	IdentityFileName = "identity.record"
	EquivFileName    = "equiv.record"
)

// Build is what the build knew that the source does not say.
type Build struct {
	// Channels are the Conda channels in the priority order the solve used. The
	// embedded environment.yml carries the channels that actually provided
	// packages; this is the order they were offered in, which the export cannot
	// show because Micromamba alphabetizes on the way out.
	Channels []string `json:"channels,omitempty"`
	// From is the upstream image a definition bootstrapped from. Nil when there
	// is no upstream.
	From *From `json:"from,omitzero"`
}

// From is what a definition bootstrapped from: the Bootstrap and From directives
// as written, and the digest they named at build time.
//
// Only the digest reaches the identity record. The other two are the recipe's own
// text rather than the mirror that served the pull, which never enters an image.
type From struct {
	// Bootstrap is the Bootstrap: directive. Without it Ref is ambiguous: docker
	// and library serve different images under one name.
	Bootstrap string `json:"bootstrap"`
	// Ref is the From: reference verbatim, e.g. "ubuntu:24.04".
	Ref string `json:"ref"`
	// Digest is the platform-specific manifest digest Ref resolved to, or
	// Unrecorded when the registry could not be reached.
	Digest string `json:"digest"`
}

// URI renders the bootstrap as a source URI, e.g. "docker://ubuntu:24.04".
func (f From) URI() string {
	if f.Bootstrap == "" {
		return f.Ref
	}
	return f.Bootstrap + "://" + f.Ref
}

// Source describes what the image was built from, for a reader that has the
// embedded files in front of it and needs to know how to use them.
type Source struct {
	// Files are the embedded sources, relative to /.cnt — "recipe" for a recipe
	// build, the Conda exports for a Conda one.
	Files []string `json:"files,omitempty"`
	// Placeholders are the selected #PH: values. This is their only stored copy:
	// the embedded recipe keeps its {placeholder} tokens, so without these
	// nothing says which variant of a template this is.
	Placeholders map[string]string `json:"placeholders,omitempty"`
	// TargetTemplate is the #TARGET: the name was rendered from, empty for a
	// recipe that is not a template.
	TargetTemplate string `json:"target_template,omitempty"`
	// RequiresInput reports that the recipe declared #INPUT: prompts, so a
	// rebuild needs a human. The answers themselves are never recorded.
	RequiresInput bool `json:"requires_input,omitempty"`
}

// Normalize fills in what a manifest is allowed to leave out: an absent or
// unrecognized type means app, and an absent OS means linux.
func (m *Manifest) Normalize() {
	m.Type = normalizeType(m.Type)
	if m.Platform.OS == "" {
		m.Platform.OS = "linux"
	}
}

// ValidateManifest reports whether m describes a usable image: a known schema, a
// name, a known type, and an architecture. Call Normalize first; the payload is
// not checked.
func ValidateManifest(m Manifest) error {
	if m.SchemaVersion != SchemaVersion {
		return fmt.Errorf("%w: manifest is %d (this build reads %d)", ErrUnsupportedSchema, m.SchemaVersion, SchemaVersion)
	}
	if strings.TrimSpace(m.Name) == "" {
		return fmt.Errorf("%w: manifest name is empty", ErrInvalid)
	}
	if m.Platform.Arch == "" {
		return fmt.Errorf("%w: %q records no architecture", ErrInvalid, m.Name)
	}

	switch m.Type {
	case catalog.TypeBase, catalog.TypeOS, catalog.TypeApp, catalog.TypeData:
	default:
		return fmt.Errorf("%w: unknown type %q", ErrInvalid, m.Type)
	}
	return nil
}

// MarshalManifest renders a manifest the way StageManifest writes it.
func MarshalManifest(m Manifest) ([]byte, error) { return marshalJSON("manifest", m) }

// StageManifest writes the manifest into dir, ready to be packed into an image
// at Path.
func StageManifest(dir string, m Manifest) error {
	data, err := MarshalManifest(m)
	if err != nil {
		return err
	}
	return stageFile(dir, FileName, data)
}

// ReadManifest returns the manifest embedded in an image. Only a genuinely
// absent manifest is ErrNoManifest; a host failure keeps its own cause. Reads
// are uncached: unlike the runtime document, nothing asks for this per exec.
func ReadManifest(imagePath string) (Manifest, error) {
	data, err := readRaw(imagePath, Path)
	if err != nil {
		if errors.Is(err, tool.ErrFileNotFound) {
			return Manifest{}, fmt.Errorf("%w: %s", ErrNoManifest, imagePath)
		}
		return Manifest{}, err
	}

	var m Manifest
	if err := json.Unmarshal(data, &m); err != nil {
		return Manifest{}, fmt.Errorf("%w: %s: %w", ErrInvalid, imagePath, err)
	}
	m.Normalize()
	if err := ValidateManifest(m); err != nil {
		return Manifest{}, fmt.Errorf("%s: %w", imagePath, err)
	}
	return m, nil
}
