package meta

import (
	"encoding/json"
	"errors"
	"fmt"
	"log/slog"
	"os"
	"path/filepath"
	"strings"

	"github.com/Justype/condatainer/catalog"
	"github.com/Justype/condatainer/internal/artifactcache"
	"github.com/Justype/condatainer/internal/image/tool"
)

var runtimeCache = artifactcache.Default()

// RuntimePath is where the runtime document lives inside every image.
const RuntimePath = "/" + DirName + "/" + RuntimeFileName

// RuntimeFileName is the runtime document's basename, for callers staging it
// into a directory.
const RuntimeFileName = "runtime.json"

// Runtime is everything container setup needs, and nothing else. It is the only
// metadata read on the mount path, so it stays small and stops growing: anything
// added for provenance belongs in the manifest, which no exec pays for.
type Runtime struct {
	SchemaVersion int          `json:"schema_version"`
	Name          string       `json:"name"`
	Type          catalog.Type `json:"type"` // base, os, app, data
	Description   string       `json:"description,omitempty"`
	Platform      Platform     `json:"platform"`
	// Prefix is where the payload sits inside the container, /cnt/<name> for app
	// and data, and what the recipe wrote to as $CNT_PREFIX. Empty for base and
	// os, whose files apply at the root. Readers take it as authoritative rather
	// than reconstructing it from the name.
	Prefix string `json:"prefix,omitempty"`
	// Env keeps {prefix} intact; it is substituted when the image is loaded,
	// because the install prefix is not known when the image is built.
	Env []EnvVar `json:"env,omitempty"`
}

// Normalize fills in what a runtime document is allowed to leave out: an absent
// or unrecognized type means app, and an absent OS means linux.
func (r *Runtime) Normalize() {
	r.Type = normalizeType(r.Type)
	if r.Platform.OS == "" {
		r.Platform.OS = "linux"
	}
}

// ValidateRuntime reports whether r describes a usable mount: a known schema, a
// name and prefix where they are load-bearing, an architecture, and an
// environment that applies without collisions. Call Normalize first.
func ValidateRuntime(r Runtime) error {
	if r.SchemaVersion != SchemaVersion {
		return fmt.Errorf("%w: runtime is %d (this build reads %d)", ErrUnsupportedSchema, r.SchemaVersion, SchemaVersion)
	}
	if strings.TrimSpace(r.Name) == "" {
		return fmt.Errorf("%w: runtime name is empty", ErrInvalid)
	}
	if r.Platform.Arch == "" {
		return fmt.Errorf("%w: %q records no architecture", ErrInvalid, r.Name)
	}

	switch r.Type {
	case catalog.TypeApp, catalog.TypeData:
		if r.Prefix == "" {
			return fmt.Errorf("%w: %s image %q has no runtime prefix", ErrInvalid, r.Type, r.Name)
		}
		if !strings.HasPrefix(r.Prefix, "/") {
			return fmt.Errorf("%w: prefix %q is not absolute", ErrInvalid, r.Prefix)
		}
	case catalog.TypeBase, catalog.TypeOS:
		// No prefix: these apply at the container root. A stray one is ignored
		// rather than rejected, since nothing reads it for these types.
	default:
		return fmt.Errorf("%w: unknown type %q", ErrInvalid, r.Type)
	}

	return validateEnv(r.Env)
}

// MarshalRuntime renders a runtime document the way StageRuntime writes it.
func MarshalRuntime(r Runtime) ([]byte, error) { return marshalJSON("runtime", r) }

// StageRuntime writes the runtime document into dir, ready to be packed into an
// image at RuntimePath.
func StageRuntime(dir string, r Runtime) error {
	data, err := MarshalRuntime(r)
	if err != nil {
		return err
	}
	return stageFile(dir, RuntimeFileName, data)
}

// ReadRuntime returns the runtime document embedded in an image. Only a
// genuinely absent one is ErrNoRuntime; a host failure keeps its own cause.
// Cached by path, size and mtime, negative verdicts included. See the README's
// Manifests.
func ReadRuntime(imagePath string) (Runtime, error) {
	return ReadRuntimeWithCache(imagePath, runtimeCache)
}

// ReadRuntimeWithCache is ReadRuntime using records. A batch lets callers that
// inspect many images persist all misses once at the end of their scan.
func ReadRuntimeWithCache(imagePath string, records artifactcache.Access) (Runtime, error) {
	abs, err := filepath.Abs(imagePath)
	if err != nil {
		abs = imagePath
	}
	fi, err := os.Lstat(abs)
	if err != nil {
		records.Forget(abs)
		return Runtime{}, fmt.Errorf("%w: %s: %w", tool.ErrUnreadable, imagePath, err)
	}

	if record, ok := records.Lookup(abs, fi); ok && record.RuntimeKnown {
		if len(record.Runtime) == 0 {
			return Runtime{}, fmt.Errorf("%w: %s", ErrNoRuntime, imagePath)
		}
		var rt Runtime
		if err := json.Unmarshal(record.Runtime, &rt); err != nil {
			records.Forget(abs)
		} else {
			rt.Normalize()
			if err := ValidateRuntime(rt); err == nil {
				return rt, nil
			}
			records.Forget(abs)
		}
	}

	rt, err := readRuntimeUncached(abs)
	switch {
	case err == nil:
		data, marshalErr := json.Marshal(rt)
		if marshalErr == nil {
			records.Merge(abs, fi, func(record *artifactcache.Record) {
				record.RuntimeKnown = true
				record.Runtime = data
			})
		}
	case errors.Is(err, ErrNoRuntime):
		records.Merge(abs, fi, func(record *artifactcache.Record) {
			record.RuntimeKnown = true
			record.Runtime = nil
		})
	}
	return rt, err
}

// readRuntimeUncached does the archive read and decode, bypassing the cache.
func readRuntimeUncached(imagePath string) (Runtime, error) {
	data, err := readRaw(imagePath, RuntimePath)
	if err != nil {
		if errors.Is(err, tool.ErrFileNotFound) {
			return Runtime{}, fmt.Errorf("%w: %s", ErrNoRuntime, imagePath)
		}
		return Runtime{}, err
	}

	var rt Runtime
	if err := json.Unmarshal(data, &rt); err != nil {
		return Runtime{}, fmt.Errorf("%w: %s: %w", ErrInvalid, imagePath, err)
	}
	rt.Normalize()
	if err := ValidateRuntime(rt); err != nil {
		return Runtime{}, fmt.Errorf("%s: %w", imagePath, err)
	}
	return rt, nil
}

// CheckBase reports whether an image may serve as a container root. Runtime
// metadata that reads has to declare type base; a missing document passes, an
// unreadable one warns and passes. See the README's Manifests.
func CheckBase(imagePath string) error {
	rt, err := ReadRuntime(imagePath)
	switch {
	case errors.Is(err, ErrNoRuntime):
		return nil
	case err != nil:
		slog.Default().Warn("could not read base image metadata", "path", imagePath, "err", err)
		return nil
	case rt.Type != catalog.TypeBase:
		return fmt.Errorf("%s is not a base image: its metadata says type %s", imagePath, rt.Type)
	}
	return nil
}
