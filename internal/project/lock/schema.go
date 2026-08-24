// Package lock reads and writes cnt-lock/, the tracked record of which exact
// artifact identity satisfies each dependency a project declares.
//
// What lives here is what Git should carry: the selection map and, beside it,
// the vendored manifests and rebuild sources. What does not is anything
// machine-local — no payload, no absolute path, no hostname, no restore result.
// A checkout is therefore a complete rebuild specification on a machine that has
// never run CondaTainer.
package lock

import (
	"bytes"
	"encoding/json"
	"errors"
	"fmt"
	"path"
	"sort"
	"strings"
)

// SchemaVersion is the lock format this build reads and writes.
const SchemaVersion = 1

// DirName is the tracked directory at a project root.
const DirName = "cnt-lock"

// FileName is the lock document inside DirName.
const FileName = "lock.json"

// ProvenanceDir holds one directory per vendored artifact, inside DirName.
const ProvenanceDir = "provenance"

// PathPrefix marks a selection key that addresses an overlay by project path
// rather than by name. A path request has no catalog.Dep to render, and the
// prefix keeps the two kinds of key from ever colliding.
const PathPrefix = "path:"

var (
	// ErrSchema reports a lock this build cannot read.
	ErrSchema = errors.New("unsupported lock schema")
	// ErrInvalid reports a structurally invalid lock.
	ErrInvalid = errors.New("invalid lock")
)

// Lock is cnt-lock/lock.json.
//
// Everything derivable is absent by design: requests come from rescanning the
// scripts, and identity, type, platform, sources and dependency edges come from
// the vendored manifests. Only the mapping and the non-derivable acquisition
// locations are stored.
type Lock struct {
	SchemaVersion int `json:"schema_version"`
	// Selections maps a canonical request to the artifact directory that
	// satisfies it. The value is an object rather than a bare path so
	// selection-specific policy can be added without copying artifact facts
	// into the mapping.
	Selections map[string]Selection `json:"selections"`
	// Remotes records where an artifact can be fetched from, keyed by the same
	// relative artifact path selections use. It is absent until something
	// records a location: `project lock select --from`, and later `project
	// push`, append to it. A closure-only artifact may have remotes too, which
	// is why this is keyed by artifact rather than nested under a selection.
	Remotes map[string][]Remote `json:"remotes,omitempty"`
}

// Selection is which artifact satisfies one request.
type Selection struct {
	// Artifact is the slash-separated path of the vendored artifact directory,
	// relative to the lock directory.
	Artifact string `json:"artifact"`
}

// Remote is one exact place an artifact can be fetched from. It is a location,
// never artifact metadata: nothing here is compared against the payload, which
// carries its own keys and is verified on arrival.
type Remote struct {
	// Repository is the complete OCI repository coordinate, with no scheme,
	// tag, or digest.
	Repository string `json:"repository"`
	// ManifestDigest is the selected platform manifest digest — not a mutable
	// tag, and not merely the multi-platform index digest.
	ManifestDigest string `json:"manifest_digest"`
}

// New returns an empty lock at the current schema version.
func New() *Lock {
	return &Lock{SchemaVersion: SchemaVersion, Selections: map[string]Selection{}}
}

// EntryPath renders the relative artifact directory for one entry name.
func EntryPath(entry string) string { return ProvenanceDir + "/" + entry }

// Unmarshal parses a lock and rejects anything this build does not understand.
//
// Unknown fields are an error rather than a silent drop: a lock written by a
// newer CondaTainer may mean something this one would quietly discard on the
// next write.
func Unmarshal(data []byte) (*Lock, error) {
	decoder := json.NewDecoder(bytes.NewReader(data))
	decoder.DisallowUnknownFields()

	var lock Lock
	if err := decoder.Decode(&lock); err != nil {
		return nil, fmt.Errorf("%w: %v", ErrInvalid, err)
	}
	if decoder.More() {
		return nil, fmt.Errorf("%w: trailing content after the lock document", ErrInvalid)
	}
	if lock.SchemaVersion != SchemaVersion {
		return nil, fmt.Errorf("%w: lock is %d (this build reads %d)", ErrSchema, lock.SchemaVersion, SchemaVersion)
	}
	if lock.Selections == nil {
		lock.Selections = map[string]Selection{}
	}
	if err := lock.Validate(); err != nil {
		return nil, err
	}
	return &lock, nil
}

// Marshal renders the lock as the bytes to publish: fixed field order, keys
// sorted, two-space indentation, one final newline. Two runs over equal content
// produce equal bytes, so a lock only appears in a diff when it changed.
func (l *Lock) Marshal() ([]byte, error) {
	if err := l.Validate(); err != nil {
		return nil, err
	}
	out := *l
	out.SchemaVersion = SchemaVersion
	if len(out.Remotes) == 0 {
		out.Remotes = nil
	}
	// encoding/json sorts map keys, and struct fields keep declaration order.
	data, err := json.MarshalIndent(out, "", "  ")
	if err != nil {
		return nil, err
	}
	return append(data, '\n'), nil
}

// Validate checks the lock's internal shape: key spelling, artifact paths, and
// remote references. It reads nothing from disk — a checkout with no images and
// no configuration validates exactly the same.
func (l *Lock) Validate() error {
	for request, selection := range l.Selections {
		if strings.TrimSpace(request) == "" {
			return fmt.Errorf("%w: a selection key is empty", ErrInvalid)
		}
		if err := validEntryPath(selection.Artifact); err != nil {
			return fmt.Errorf("%w: selection %q: %v", ErrInvalid, request, err)
		}
		if destination, ok := strings.CutPrefix(request, PathPrefix); ok {
			if err := validProjectPath(destination); err != nil {
				return fmt.Errorf("%w: selection %q: %v", ErrInvalid, request, err)
			}
		}
	}
	for artifact, remotes := range l.Remotes {
		if err := validEntryPath(artifact); err != nil {
			return fmt.Errorf("%w: remote key: %v", ErrInvalid, err)
		}
		seen := make(map[Remote]bool, len(remotes))
		for _, remote := range remotes {
			if err := remote.validate(); err != nil {
				return fmt.Errorf("%w: remote for %q: %v", ErrInvalid, artifact, err)
			}
			if seen[remote] {
				return fmt.Errorf("%w: remote for %q is listed twice: %s@%s",
					ErrInvalid, artifact, remote.Repository, remote.ManifestDigest)
			}
			seen[remote] = true
		}
	}
	return nil
}

// validEntryPath requires a relative, slash-separated path inside the
// provenance directory. Every lock path is hostile input: it names a directory
// a later step will read, so traversal and absolute paths are refused before
// anything touches the filesystem.
func validEntryPath(p string) error {
	switch {
	case p == "":
		return errors.New("artifact path is empty")
	case strings.ContainsRune(p, '\\'):
		return fmt.Errorf("artifact path %q is not slash-separated", p)
	case path.IsAbs(p), strings.HasPrefix(p, "/"):
		return fmt.Errorf("artifact path %q is absolute", p)
	case p != path.Clean(p):
		return fmt.Errorf("artifact path %q is not clean", p)
	case !strings.HasPrefix(p, ProvenanceDir+"/"):
		return fmt.Errorf("artifact path %q is outside %s/", p, ProvenanceDir)
	case strings.Contains(p, "/../"), strings.HasSuffix(p, "/.."):
		return fmt.Errorf("artifact path %q escapes the lock directory", p)
	}
	if len(strings.Split(p, "/")) != 2 {
		return fmt.Errorf("artifact path %q is not %s/<entry>", p, ProvenanceDir)
	}
	return nil
}

// validProjectPath requires a relative, slash-separated path inside the project
// for a `path:` selection.
//
// That path is not somewhere to read from — it is where restore *writes* the
// artifact. An absolute or escaping one would have restore create a file
// outside the checkout entirely, so it is refused here, where the lock is
// parsed, rather than trusted by everything downstream.
func validProjectPath(p string) error {
	switch {
	case p == "":
		return errors.New("project path is empty")
	case strings.ContainsRune(p, '\\'):
		return fmt.Errorf("project path %q is not slash-separated", p)
	case path.IsAbs(p), strings.HasPrefix(p, "/"):
		return fmt.Errorf("project path %q is absolute", p)
	case p != path.Clean(p):
		return fmt.Errorf("project path %q is not clean", p)
	case p == "..", strings.HasPrefix(p, "../"), strings.Contains(p, "/../"), strings.HasSuffix(p, "/.."):
		return fmt.Errorf("project path %q escapes the project", p)
	case strings.HasPrefix(p, DirName+"/"):
		return fmt.Errorf("project path %q is inside %s, which holds no payload", p, DirName)
	case !strings.HasSuffix(p, ".sqf"):
		return fmt.Errorf("project path %q is not a .sqf", p)
	}
	return nil
}

func (o Remote) validate() error {
	if strings.TrimSpace(o.Repository) == "" {
		return errors.New("repository is empty")
	}
	if strings.Contains(o.Repository, "://") {
		return fmt.Errorf("repository %q carries a scheme", o.Repository)
	}
	if i := strings.IndexAny(o.Repository, "@:"); i >= 0 && !isPortColon(o.Repository, i) {
		return fmt.Errorf("repository %q carries a tag or digest", o.Repository)
	}
	digest, ok := strings.CutPrefix(o.ManifestDigest, "sha256:")
	if !ok || len(digest) != 64 || strings.TrimLeft(digest, "0123456789abcdef") != "" {
		return fmt.Errorf("manifest digest %q is not sha256:<64 hex>", o.ManifestDigest)
	}
	return nil
}

// isPortColon reports whether the colon at i separates a registry host from its
// port, as in localhost:5000/name, rather than introducing a tag.
func isPortColon(repository string, i int) bool {
	if repository[i] != ':' {
		return false
	}
	return strings.Contains(repository[i:], "/")
}

// AddRemote records one fetch location for an artifact, keeping the list
// deduplicated and in insertion order. Order is retry priority, so an existing
// entry stays where it is rather than moving to the front.
//
// This is the seam `project push` writes through: publishing appends the
// repository and platform digest it produced, for selected and closure-only
// artifacts alike.
func (l *Lock) AddRemote(artifact string, remote Remote) error {
	if err := validEntryPath(artifact); err != nil {
		return fmt.Errorf("%w: %v", ErrInvalid, err)
	}
	if err := remote.validate(); err != nil {
		return fmt.Errorf("%w: %v", ErrInvalid, err)
	}
	for _, have := range l.Remotes[artifact] {
		if have == remote {
			return nil
		}
	}
	if l.Remotes == nil {
		l.Remotes = map[string][]Remote{}
	}
	l.Remotes[artifact] = append(l.Remotes[artifact], remote)
	return nil
}

// Requests returns every selection key, sorted, so callers report in a stable
// order without each sorting for themselves.
func (l *Lock) Requests() []string {
	out := make([]string, 0, len(l.Selections))
	for request := range l.Selections {
		out = append(out, request)
	}
	sort.Strings(out)
	return out
}

// SelectedArtifacts is the set of artifact paths the selections point at
// directly. These are the *roots* of reachability, not reachability itself:
// each one's vendored manifest names dependency edges that pull further
// artifact directories into the closure, and walking those needs to read the
// manifests. Pruning against this set alone would delete the closure.
//
// Remotes are deliberately not roots. An remote says where an artifact can be
// fetched, which is meaningless once nothing selects it.
func (l *Lock) SelectedArtifacts() map[string]bool {
	out := make(map[string]bool, len(l.Selections))
	for _, selection := range l.Selections {
		out[selection.Artifact] = true
	}
	return out
}
