// Package meta is the metadata CondaTainer embeds in an image and reads back out
// of it, split by how often it is read: runtime.json is the mount-time contract,
// read on every exec; manifest.json describes what the image is and where it came
// from, read on demand. Both are trusted, not verified.
package meta

import (
	"encoding/json"
	"errors"
	"fmt"
	"os"
	"path/filepath"
	"runtime"
	"strings"

	"github.com/Justype/condatainer/catalog"
	"github.com/Justype/condatainer/internal/image/sif"
	"github.com/Justype/condatainer/internal/image/squashfs"
	"github.com/Justype/condatainer/internal/image/tool"
	"github.com/Justype/condatainer/internal/utils"
)

// DirName is the directory holding CondaTainer metadata, relative to the image root.
const DirName = ".cnt"

// RecipeFileName is the recipe a build embeds beside its metadata, byte for byte
// as it was fetched. RecipePath is where it lands inside the image.
const (
	RecipeFileName = "recipe"
	RecipePath     = "/" + DirName + "/" + RecipeFileName
)

// SchemaVersion is the schema this build reads and writes, for both documents. An
// unknown version is reported as ErrUnsupportedSchema and treated like a missing
// document. See the README's Manifests.
const SchemaVersion = 1

// Errors a caller distinguishes. Everything else is a plain validation failure.
var (
	// ErrNoRuntime reports that the image carries no runtime.json. It is expected
	// for any image built before the runtime split, which mounts and contributes
	// nothing, so it is the one read outcome callers routinely tolerate.
	ErrNoRuntime = errors.New("image has no CondaTainer runtime metadata")
	// ErrNoManifest reports that the image contains no manifest at all.
	ErrNoManifest = errors.New("image has no CondaTainer manifest")
	// ErrUnsupportedSchema reports a schema_version this build does not know.
	ErrUnsupportedSchema = errors.New("unsupported metadata schema version")
	// ErrInvalid reports metadata that is present and decodes but does not
	// describe a usable image.
	ErrInvalid = errors.New("invalid image metadata")
)

// Platform is where an image was built and where it may run.
type Platform struct {
	// OS is recorded as linux and never compared: everything runs in a Linux
	// container, so it is a label rather than a gate.
	OS string `json:"os"`
	// Arch is the build architecture in uname form, or the literal "noarch" for
	// an artifact whose recipe declared it is architecture-independent.
	Arch string `json:"arch"`
}

// ArchNone is the Arch value of an artifact that runs anywhere.
const ArchNone = "noarch"

// NativePlatform is the platform of the machine this process runs on.
func NativePlatform() Platform {
	return Platform{OS: "linux", Arch: NativeArch()}
}

// NativeArch reports the running architecture in uname form, which is what
// recipes, Conda subdirs and users all say.
func NativeArch() string {
	switch runtime.GOARCH {
	case "amd64":
		return "x86_64"
	case "arm64":
		return "aarch64"
	default:
		return runtime.GOARCH
	}
}

// EnvVar is one #ENV: contribution.
type EnvVar struct {
	Key   string `json:"key"`
	Value string `json:"value"`
	Note  string `json:"note,omitempty"`
}

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

// normalizeType maps an absent or unrecognized type to app.
func normalizeType(typ catalog.Type) catalog.Type {
	switch typ {
	case catalog.TypeBase, catalog.TypeOS, catalog.TypeApp, catalog.TypeData:
		return typ
	default:
		return catalog.TypeApp
	}
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

// validateEnv reports whether the environment applies without collisions.
func validateEnv(env []EnvVar) error {
	seen := make(map[string]bool, len(env))
	for _, e := range env {
		if !validEnvKey(e.Key) {
			return fmt.Errorf("%w: %q is not a valid environment variable name", ErrInvalid, e.Key)
		}
		if seen[e.Key] {
			return fmt.Errorf("%w: %s is set more than once", ErrInvalid, e.Key)
		}
		seen[e.Key] = true
	}
	return nil
}

// marshalJSON renders a metadata document the way staging writes it: indented,
// newline-terminated, and byte-identical for the same input.
func marshalJSON(what string, v any) ([]byte, error) {
	data, err := json.MarshalIndent(v, "", "  ")
	if err != nil {
		return nil, fmt.Errorf("failed to encode %s: %w", what, err)
	}
	return append(data, '\n'), nil
}

// StageBytes writes a metadata file into dir under its own name, ready to be
// packed into an image at /.cnt/<name>. It is how a build stages the sources it
// embeds verbatim — the recipe, a Conda export — beside the two documents.
func StageBytes(dir, name string, data []byte) error { return stageFile(dir, name, data) }

// stageFile writes data into dir/name, ready to be packed into an image. The
// write goes through a temp file and an atomic rename.
func stageFile(dir, name string, data []byte) error {
	if err := utils.MkdirAllShared(dir); err != nil {
		return fmt.Errorf("failed to create metadata dir %s: %w", dir, err)
	}

	final := filepath.Join(dir, name)
	tmp := final + ".tmp"
	f, err := utils.CreateFileWritable(tmp)
	if err != nil {
		return fmt.Errorf("failed to create %s: %w", tmp, err)
	}
	defer os.Remove(tmp) // no-op after a successful rename

	if _, err := f.Write(data); err != nil {
		f.Close() //nolint:errcheck
		return fmt.Errorf("failed to write %s: %w", tmp, err)
	}
	if err := f.Close(); err != nil {
		return fmt.Errorf("failed to close %s: %w", tmp, err)
	}
	if err := os.Rename(tmp, final); err != nil {
		return fmt.Errorf("failed to install %s: %w", final, err)
	}
	utils.ShareWithParentGroup(final)
	return nil
}

// readRaw pulls one metadata file's bytes out of whichever image format this is.
// A writable .img is not handled: its environment comes from its .env sidecar.
func readRaw(imagePath, inImagePath string) ([]byte, error) {
	switch {
	case utils.IsSqf(imagePath):
		return squashfs.CatFile(imagePath, inImagePath, 0)
	case utils.IsSif(imagePath):
		return sif.ReadFile(imagePath, inImagePath)
	default:
		return nil, fmt.Errorf("%w: %s is not a CondaTainer image", tool.ErrCorrupt, imagePath)
	}
}
