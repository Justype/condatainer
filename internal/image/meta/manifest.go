// Package meta is the metadata CondaTainer embeds in an image and reads back out
// of it: one manifest.json describing what the image is, where its payload sits,
// and what it contributes to the environment. It is trusted, not verified.
package meta

import (
	"errors"
	"fmt"
	"strings"

	"github.com/Justype/condatainer/catalog"
)

// Path is where the manifest lives inside every image.
const Path = "/.cnt/manifest.json"

// FileName is the manifest's basename, for callers staging it into a directory.
const FileName = "manifest.json"

// DirName is the directory holding the manifest, relative to the image root.
const DirName = ".cnt"

// SchemaVersion is the schema this build reads and writes. An unknown version
// is reported as ErrUnsupportedSchema and treated like a missing manifest.
// See the README's Manifests.
const SchemaVersion = 1

// Manifest is what an image records about itself.
type Manifest struct {
	SchemaVersion int          `json:"schema_version"`
	Name          string       `json:"name"`
	Type          catalog.Type `json:"type"`       // base, os, app, data
	BuildType     string       `json:"build_type"` // conda, script, def
	Description   string       `json:"description,omitempty"`
	URL           string       `json:"url,omitempty"`
	Runtime       Runtime      `json:"runtime"`
}

// Runtime is what loading the image does to the environment.
type Runtime struct {
	// Prefix is the image's install prefix, /cnt/<name>, for app and data: where
	// its payload sits inside the container, and what the recipe wrote to as
	// $CNT_PREFIX. Empty for base and os, whose files apply at the root.
	Prefix string `json:"prefix,omitempty"`
	// Env keeps {prefix} intact; it is substituted when the image is loaded,
	// because the install prefix is not known when the image is built.
	Env []EnvVar `json:"env,omitempty"`
}

// EnvVar is one #ENV: contribution.
type EnvVar struct {
	Key   string `json:"key"`
	Value string `json:"value"`
	Note  string `json:"note,omitempty"`
}

// Errors a caller distinguishes. Everything else is a plain validation failure.
var (
	// ErrNoManifest reports that the image contains no manifest at all. It is
	// expected for an image built before this format and for a plain Apptainer
	// SIF, so it is the one read outcome callers routinely tolerate.
	ErrNoManifest = errors.New("image has no CondaTainer manifest")
	// ErrUnsupportedSchema reports a schema_version this build does not know.
	ErrUnsupportedSchema = errors.New("unsupported manifest schema version")
	// ErrInvalid reports a manifest that is present and decodes but does not
	// describe a usable image.
	ErrInvalid = errors.New("invalid manifest")
)

// Resolved renders the value with {prefix} replaced by the image's install prefix.
func (v EnvVar) Resolved(prefix string) string {
	return strings.ReplaceAll(v.Value, "{prefix}", prefix)
}

// Prefix returns the install prefix an image of this name and type gets. It is
// empty for base and os, which apply at the container root.
func Prefix(name string, typ catalog.Type) string {
	if typ == catalog.TypeBase || typ == catalog.TypeOS {
		return ""
	}
	return "/cnt/" + strings.Trim(name, "/")
}

// Normalize fills in what a manifest is allowed to leave out: an absent or
// unrecognized type means app.
func (m *Manifest) Normalize() {
	switch m.Type {
	case catalog.TypeBase, catalog.TypeOS, catalog.TypeApp, catalog.TypeData:
	default:
		m.Type = catalog.TypeApp
	}
}

// Validate reports whether m describes a usable image: a known schema, a name
// and prefix where they are load-bearing, and an environment that applies
// without collisions. Call Normalize first; the payload is not checked.
func Validate(m Manifest) error {
	if m.SchemaVersion != SchemaVersion {
		return fmt.Errorf("%w: %d (this build reads %d)", ErrUnsupportedSchema, m.SchemaVersion, SchemaVersion)
	}
	if strings.TrimSpace(m.Name) == "" {
		return fmt.Errorf("%w: name is empty", ErrInvalid)
	}

	switch m.Type {
	case catalog.TypeApp, catalog.TypeData:
		if m.Runtime.Prefix == "" {
			return fmt.Errorf("%w: %s image %q has no runtime prefix", ErrInvalid, m.Type, m.Name)
		}
		if !strings.HasPrefix(m.Runtime.Prefix, "/") {
			return fmt.Errorf("%w: prefix %q is not absolute", ErrInvalid, m.Runtime.Prefix)
		}
	case catalog.TypeBase, catalog.TypeOS:
		// No prefix: these apply at the container root. A stray one is ignored
		// rather than rejected, since nothing reads it for these types.
	default:
		return fmt.Errorf("%w: unknown type %q", ErrInvalid, m.Type)
	}

	seen := make(map[string]bool, len(m.Runtime.Env))
	for _, env := range m.Runtime.Env {
		if !validEnvKey(env.Key) {
			return fmt.Errorf("%w: %q is not a valid environment variable name", ErrInvalid, env.Key)
		}
		if seen[env.Key] {
			return fmt.Errorf("%w: %s is set more than once", ErrInvalid, env.Key)
		}
		seen[env.Key] = true
	}
	return nil
}

// validEnvKey reports whether s can be exported as a shell variable.
func validEnvKey(s string) bool {
	if s == "" {
		return false
	}
	for i, r := range s {
		switch {
		case r == '_':
		case r >= 'A' && r <= 'Z', r >= 'a' && r <= 'z':
		case r >= '0' && r <= '9' && i > 0:
		default:
			return false
		}
	}
	return true
}
