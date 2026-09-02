package meta

import (
	"crypto/sha256"
	"encoding/json"
	"errors"
	"fmt"
	"strings"
	"time"

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
	Type          catalog.Type `json:"type"`
	BuildType     BuildType    `json:"build_type"`
	Description   string       `json:"description,omitempty"`
	URL           string       `json:"url,omitempty"`
	// License is the recipe's #LICENSE:, an SPDX expression kept verbatim and
	// never parsed. It is published as org.opencontainers.image.licenses.
	License string `json:"license,omitempty"`
	// Redistribute is the recipe's #REDISTRIBUTE: answer, and nil when it did
	// not answer — three states, because "the author said this may not be
	// republished" and "nobody has been asked" default differently.
	//
	// It decides publication and nothing else: it is read at push time, enters
	// no key preimage, and changing it moves neither identity nor equivalence.
	Redistribute *bool        `json:"redistribute,omitempty"`
	Platform     Platform     `json:"platform"`
	Source       Source       `json:"source,omitzero"`
	Keys         Keys         `json:"keys,omitzero"`
	Dependencies []Dependency `json:"dependencies,omitempty"`
	// ProvenanceComplete reports whether every dependency carried records of its
	// own. It is nil when the question does not arise — an artifact with no
	// dependencies is neither complete nor incomplete.
	ProvenanceComplete *bool `json:"provenance_complete,omitempty"`
	Build              Build `json:"build,omitzero"`
	// Snapshot is what a frozen environment captured, and is nil for anything
	// built from a recipe.
	Snapshot *Snapshot `json:"snapshot,omitzero"`
}

// Snapshot describes a frozen writable overlay. Nothing here gates anything:
// two environments are kept apart by their prefix, which every one of them
// records as EnvPrefix. The sizes are recorded because unfreeze needs them
// before it can create the image, and recovering them from the archive means
// decompressing its metadata entry by entry.
type Snapshot struct {
	// Convention is the whiteout convention the source overlay used, before
	// translation: chardev, whfile, mixed, or none. Diagnostics — a file that
	// came back from the dead reads differently on each driver.
	Convention string `json:"convention,omitempty"`
	// Whiteouts is how many deletions the artifact carries, after translation.
	// An opaque directory expands to one per hidden base entry, so this counts
	// nodes rather than user-visible deletions.
	Whiteouts int `json:"whiteouts,omitempty"`
	// PayloadMB is the space the payload needs as a filesystem, rounded up to
	// blocks: a 14-byte file still occupies one, and summing apparent sizes
	// understates a conda prefix by orders of magnitude.
	PayloadMB int `json:"payload_mb,omitempty"`
	// Entries is how many things the payload holds, which sizes the inode table.
	Entries int `json:"entries,omitempty"`
}

// Keys holds the versioned derivation scheme and expected digest for both keys.
type Keys struct {
	Identity KeyRef `json:"identity,omitzero"`
	Equiv    KeyRef `json:"equiv,omitzero"`
}

// KeyRef is one derived key: its immutable scheme and expected SHA-256.
type KeyRef struct {
	Scheme string `json:"scheme"`
	SHA256 string `json:"sha256"`
}

// Digest renders sha256:<hex>, or empty when the reference is absent.
func (k KeyRef) Digest() string {
	if k.Scheme == "" || k.SHA256 == "" {
		return ""
	}
	return "sha256:" + k.SHA256
}

// Empty reports whether a key reference is completely absent.
func (k KeyRef) Empty() bool { return k.Scheme == "" && k.SHA256 == "" }

// Dependency is one direct build dependency, as the artifact recorded it. The
// list is adjacency only: following each manifest's own through the capsule
// reconstructs the graph without a second representation of it.
type Dependency struct {
	Name string       `json:"name"`
	Type catalog.Type `json:"type"`
	// Identity and Equiv are complete keys, scheme and SHA-256 both, so an edge
	// is held to the same contract as the artifact it points at.
	Identity KeyRef `json:"identity,omitzero"`
	Equiv    KeyRef `json:"equiv,omitzero"`
	// Records is "unrecorded" when the image that satisfied this dependency
	// carried no keys of its own, which is every image built before this format.
	Records string `json:"records,omitempty"`
	// Role freezes how this dependency contributes to equivalence. A reader trusts
	// it rather than applying current policy and silently rewriting the past.
	Role string `json:"role"`
}

// Roles a dependency can play in its dependent's equivalence.
const (
	RoleData    = "data"    // contributes its equivalence
	RoleApp     = "app"     // named in the artifact's name; contributes name/version
	RoleHistory = "history" // mounted, but decides nothing about substitution
)

// Unrecorded marks a dependency satisfied by an image carrying no scheme-backed keys.
const Unrecorded = "unrecorded"

// Build is what the build knew that the recipe does not say.
type Build struct {
	// Tools identify the implementations that performed the build. They are
	// diagnostic provenance only: key schemes select their own inputs and do not
	// implicitly hash this block.
	Tools BuildTools `json:"tools,omitzero"`
	// Channels are the Conda channels in the priority order the solve used. The
	// embedded environment.yml carries the channels that actually provided
	// packages; this is the order they were offered in, which the export cannot
	// show because Micromamba alphabetizes on the way out.
	Channels []string `json:"channels,omitempty"`
	// From is the upstream image a definition bootstrapped from. Nil when there
	// is no upstream.
	From *From `json:"from,omitzero"`
	// Source is the repository of the collection that supplied the recipe, empty
	// when it declares none. Never the local handle, which names nothing outside
	// one installation's config.
	Source string `json:"source,omitempty"`
	// Created is when the build finished. The SquashFS superblock time is not
	// usable instead: a reproducible build pins it, and a SIF has none.
	Created time.Time `json:"created,omitzero"`
}

// BuildTools are the tools Condatainer directly used for a build. Apptainer
// names the compatible tool family; Tool.Name distinguishes an actual
// Apptainer binary from the supported Singularity fallback.
type BuildTools struct {
	Condatainer Tool `json:"condatainer,omitzero"`
	Apptainer   Tool `json:"apptainer,omitzero"`
	Micromamba  Tool `json:"micromamba,omitzero"`
}

// Empty reports whether no build tool was recorded.
func (t BuildTools) Empty() bool {
	return t.Condatainer.Empty() && t.Apptainer.Empty() && t.Micromamba.Empty()
}

// Tool is one directly used build implementation. Name is omitted when the
// enclosing field already identifies it; it is set for Apptainer because the
// configured binary may instead be Singularity.
type Tool struct {
	Name    string `json:"name,omitempty"`
	Version string `json:"version"`
}

// Empty reports whether a tool is completely absent.
func (t Tool) Empty() bool { return t.Name == "" && t.Version == "" }

// From is what a definition bootstrapped from: the Bootstrap and From directives
// as written, and the digest they named at build time.
//
// Only the digest reaches the identity scheme. The other two are the recipe's own
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

// Normalize fills in what a manifest is allowed to leave out: an absent type
// means app, and an absent OS means linux.
func (m *Manifest) Normalize() {
	m.Type = normalizeType(m.Type)
	if m.Platform.OS == "" {
		m.Platform.OS = "linux"
	}
}

// BuildType is how an artifact was produced, beside catalog.Type, which is what
// the artifact is. internal/build aliases this type and its constants.
type BuildType string

const (
	// BuildTypeConda is a Micromamba solve. It embeds no recipe, only its exports.
	BuildTypeConda BuildType = "conda"
	// BuildTypeScript is a recipe run as a shell script.
	BuildTypeScript BuildType = "script"
	// BuildTypeDef is an Apptainer definition file.
	BuildTypeDef BuildType = "def"
	// BuildTypeSnapshot is a writable overlay packed by `overlay freeze`.
	// ValidateManifest requires it to accompany catalog.TypeEnv and the name
	// EnvName.
	BuildTypeSnapshot BuildType = "snapshot"
)

// String renders the build type, or "unknown" when it is unset.
func (b BuildType) String() string {
	if b == "" {
		return "unknown"
	}
	return string(b)
}

// EnvPrefix is where an environment's conda prefix lives inside the container,
// for a writable .img and the frozen artifact made from one alike. Every
// environment records this or nothing, which is what stops two mounting together.
const EnvPrefix = "/cnt_env"

// EnvName is what every frozen environment is called. It is fixed: an
// environment is addressed by its path, so a name would distinguish nothing.
const EnvName = "env"

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
	case catalog.TypeBase, catalog.TypeOS, catalog.TypeApp, catalog.TypeData, catalog.TypeEnv:
	default:
		return fmt.Errorf("%w: unknown type %q", ErrInvalid, m.Type)
	}
	switch m.BuildType {
	case BuildTypeConda, BuildTypeScript, BuildTypeDef, BuildTypeSnapshot:
	default:
		return fmt.Errorf("%w: %s records unknown build type %q", ErrInvalid, m.Name, m.BuildType)
	}

	// A snapshot and the env type imply each other. Neither is derivable from the
	// other by a reader, and an artifact carrying one without the other would be
	// read correctly by whichever rule looked at the field it had — the publishing
	// rule branches on the build type, the mount rules on the type — so the two
	// must not be able to disagree.
	if (m.BuildType == BuildTypeSnapshot) != (m.Type == catalog.TypeEnv) {
		return fmt.Errorf("%w: %s is type %q with build type %q; a frozen environment is %q and %q, and nothing else is either",
			ErrInvalid, m.Name, m.Type, m.BuildType, catalog.TypeEnv, BuildTypeSnapshot)
	}
	if m.BuildType == BuildTypeSnapshot && m.Name != EnvName {
		return fmt.Errorf("%w: a frozen environment is always called %q, not %q",
			ErrInvalid, EnvName, m.Name)
	}
	// A snapshot's two keys are one value. It has no inputs to abstract away, so
	// "is this the same environment" and "can this substitute for it" cannot come
	// apart, and a manifest claiming otherwise describes a distinction that does
	// not exist.
	if m.BuildType == BuildTypeSnapshot && m.Keys.Identity != m.Keys.Equiv {
		return fmt.Errorf("%w: %s is a frozen environment; its identity and equivalence are the same value",
			ErrInvalid, m.Name)
	}

	identityEmpty := m.Keys.Identity.Empty()
	equivEmpty := m.Keys.Equiv.Empty()
	if identityEmpty != equivEmpty {
		return fmt.Errorf("%w: identity and equivalence keys must both be present or absent", ErrInvalid)
	}
	if !identityEmpty {
		refs := []struct {
			kind string
			ref  KeyRef
		}{
			{kind: "identity", ref: m.Keys.Identity},
			{kind: "equiv", ref: m.Keys.Equiv},
		}
		for _, item := range refs {
			if strings.TrimSpace(item.ref.Scheme) == "" {
				return fmt.Errorf("%w: %s key has no scheme", ErrInvalid, item.kind)
			}
			if !validSHA256(item.ref.SHA256) {
				return fmt.Errorf("%w: %s key has invalid sha256 %q", ErrInvalid, item.kind, item.ref.SHA256)
			}
		}
	}
	if err := validateBuildTools(m.BuildType, m.Build.Tools); err != nil {
		return err
	}
	return nil
}

func validateBuildTools(buildType BuildType, tools BuildTools) error {
	if tools.Empty() {
		return nil
	}
	if tools.Condatainer.Version == "" {
		return fmt.Errorf("%w: build tools have no Condatainer version", ErrInvalid)
	}
	if tools.Apptainer.Name == "" || tools.Apptainer.Version == "" {
		return fmt.Errorf("%w: build tools have incomplete Apptainer information", ErrInvalid)
	}

	switch buildType {
	case BuildTypeConda:
		if tools.Micromamba.Version == "" {
			return fmt.Errorf("%w: Conda build tools have no Micromamba version", ErrInvalid)
		}
	case BuildTypeDef, BuildTypeScript, BuildTypeSnapshot:
		// A freeze packs what is already there and never solves, so Micromamba is
		// not one of its tools even though the payload is a conda environment.
		if !tools.Micromamba.Empty() {
			return fmt.Errorf("%w: %s build records Micromamba as a direct build tool", ErrInvalid, buildType)
		}
	default:
		return fmt.Errorf("%w: build tools accompany unknown build type %q", ErrInvalid, buildType)
	}
	return nil
}

func validSHA256(s string) bool {
	if len(s) != 2*sha256.Size {
		return false
	}
	for _, c := range s {
		if (c < '0' || c > '9') && (c < 'a' || c > 'f') {
			return false
		}
	}
	return true
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
	return DecodeManifest(data, imagePath)
}

// DecodeManifest turns a manifest document into a validated, normalized
// Manifest. source names the document's origin in errors. Split from the
// archive read so what the bytes mean is decided in one place, whether they
// came from an image, a staging directory, or a registry blob.
func DecodeManifest(data []byte, source string) (Manifest, error) {
	var m Manifest
	if err := json.Unmarshal(data, &m); err != nil {
		return Manifest{}, fmt.Errorf("%w: %s: %w", ErrInvalid, source, err)
	}
	m.Normalize()
	if err := ValidateManifest(m); err != nil {
		return Manifest{}, fmt.Errorf("%s: %w", source, err)
	}
	return m, nil
}
