// Package lock reads and writes cnt-lock/, the tracked record of which exact
// artifact identity satisfies each dependency a project declares.
//
// What lives here is what Git should carry: the pins and, beside them,
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
	"regexp"
	"sort"
	"strings"

	"github.com/Justype/condatainer/catalog"
)

// SchemaVersion is the lock format this build reads and writes.
const SchemaVersion = 1

// DirName is the tracked directory at a project root.
const DirName = "cnt-lock"

// FileName is the lock document inside DirName.
const FileName = "lock.json"

// ProvenanceDir holds one directory per vendored artifact, inside DirName.
const ProvenanceDir = "provenance"

// PathPrefix marks a pin key that addresses an overlay by project path
// rather than by name. A path request has no catalog.Dep to render, and the
// prefix keeps the two kinds of key from ever colliding.
const PathPrefix = "path:"

// ociRepoSegment is one path segment of an OCI repository, matching the
// distribution spec's grammar. Duplicated from internal/registry rather than
// imported: this package must stay free of the transport so that project
// validation is checkout-local by construction.
var ociRepoSegment = regexp.MustCompile(`^[a-z0-9]+(?:[._-][a-z0-9]+)*$`)

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
	// Source is the project's own code repository, the same key and the same
	// meaning a collection's source.json gives it. It becomes
	// org.opencontainers.image.source on everything the project publishes, which
	// is what links a GHCR package to the repository.
	//
	// Recorded rather than read from `git remote origin` at push time: two
	// collaborators with different remote spellings would otherwise publish two
	// different annotations for one project.
	Source string `json:"source,omitempty"`
	// OCI is where this project publishes. Absent until `project registry set`
	// records one.
	OCI OCI `json:"oci,omitzero"`
	// Pins maps a canonical request to the artifact directory that
	// satisfies it. The value is an object rather than a bare path so
	// pin-specific policy can be added without copying artifact facts
	// into the mapping.
	Pins map[string]PinEntry `json:"pins"`
	// Remotes records where an artifact can be fetched from, keyed by the same
	// relative artifact path pins use. It is absent until something
	// records a location: `project pin`, and later `project push`, append to it. A closure-only artifact may have remotes too, which
	// is why this is keyed by artifact rather than nested under a pin.
	Remotes map[string][]Remote `json:"remotes,omitempty"`
}

// OCI is a project's own publishing destination. The field names mirror
// catalog.Descriptor's OCI block, because a collection declares these about
// itself and a project declares the same things about itself.
//
// There is no Pull: a project has no name-addressed pull list. Where its
// artifacts can be fetched from is per-artifact and lives in Remotes.
type OCI struct {
	// Push is the complete repository coordinate, registry host included and
	// with no scheme, tag, or digest.
	Push string `json:"push,omitempty"`
	// Audience is who can pull from the endpoint, and CondaTainer derives from
	// it what may be published there. Empty means public, which is the
	// restrictive answer — see registry.Audience.
	Audience string `json:"audience,omitempty"`
}

// Empty reports that no destination has been recorded.
func (o OCI) Empty() bool { return strings.TrimSpace(o.Push) == "" }

// PinEntry is which artifact satisfies one request.
type PinEntry struct {
	// Artifact is the slash-separated path of the vendored artifact directory,
	// relative to the lock directory.
	Artifact string `json:"artifact"`
	// Manual marks a pin the project recorded rather than one derived from a
	// #DEP:. Pins are otherwise a projection of the scan, and Reconcile deletes
	// every key the scan no longer produces; a manual one is removed only when
	// someone says so.
	//
	// It exists because a project depends on artifacts no script declares: a
	// helper's #REQUIRED_OVERLAYS: (build-essential, an RStudio image), and a
	// frozen environment. Recorded rather than derived, since "not currently
	// declared" is what the sweep tests and deriving it would make the sweep a
	// no-op.
	Manual bool `json:"manual,omitempty"`
}

// Remote is one exact place an artifact can be fetched from. It is a location,
// never artifact metadata: nothing here is compared against the payload, which
// carries its own keys and is verified on arrival.
type Remote struct {
	// Repository is the complete OCI repository coordinate, with no scheme,
	// tag, or digest.
	Repository string `json:"repository"`
	// ManifestDigest is the pinned platform manifest digest — not a mutable
	// tag, and not merely the multi-platform index digest.
	ManifestDigest string `json:"manifest_digest"`
}

// New returns an empty lock at the current schema version.
func New() *Lock {
	return &Lock{SchemaVersion: SchemaVersion, Pins: map[string]PinEntry{}}
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
	if lock.Pins == nil {
		lock.Pins = map[string]PinEntry{}
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
	out.Source = strings.TrimSpace(out.Source)
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
	if err := validSourceURL(l.Source); err != nil {
		return fmt.Errorf("%w: source: %v", ErrInvalid, err)
	}
	if err := l.OCI.validate(); err != nil {
		return fmt.Errorf("%w: oci: %v", ErrInvalid, err)
	}
	for request, pin := range l.Pins {
		if strings.TrimSpace(request) == "" {
			return fmt.Errorf("%w: a pin key is empty", ErrInvalid)
		}
		if err := validatePin(request, pin); err != nil {
			return fmt.Errorf("%w: pin %q: %v", ErrInvalid, request, err)
		}
	}
	for artifact, remotes := range l.Remotes {
		if err := validEntryPath(artifact); err != nil {
			return fmt.Errorf("%w: remote key: %v", ErrInvalid, err)
		}
		if err := validateRemotes(artifact, remotes); err != nil {
			return fmt.Errorf("%w: %v", ErrInvalid, err)
		}
	}
	return nil
}

// validatePin checks one entry: a vendored artifact directory, and a project
// path when the key names one.
func validatePin(request string, pin PinEntry) error {
	if destination, ok := strings.CutPrefix(request, PathPrefix); ok {
		if err := validProjectPath(destination); err != nil {
			return err
		}
	}
	return validEntryPath(pin.Artifact)
}

// validateRemotes checks one artifact's fetch locations.
func validateRemotes(artifact string, remotes []Remote) error {
	seen := make(map[Remote]bool, len(remotes))
	for _, remote := range remotes {
		if err := remote.validate(); err != nil {
			return fmt.Errorf("remote for %q: %v", artifact, err)
		}
		if seen[remote] {
			return fmt.Errorf("remote for %q is listed twice: %s@%s",
				artifact, remote.Repository, remote.ManifestDigest)
		}
		seen[remote] = true
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
// for a `path:` pin.
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

// validate checks a recorded publishing destination.
//
// The coordinate is held to what a registry will actually accept rather than
// merely to "not obviously wrong", because it is written into a tracked file and
// every later push composes tags onto it. Nothing is down-cased: two names
// differing only in case would collide into one repository, and silently
// publishing one over the other is worse than refusing both.
func (o OCI) validate() error {
	push := strings.TrimSpace(o.Push)
	switch strings.ToLower(strings.TrimSpace(o.Audience)) {
	case "", "public", "restricted":
	default:
		return fmt.Errorf("audience %q is not public or restricted", o.Audience)
	}
	if push == "" {
		if strings.TrimSpace(o.Audience) != "" {
			return errors.New("an audience was recorded without a push destination")
		}
		return nil
	}
	if strings.Contains(push, "://") {
		return fmt.Errorf("push %q carries a scheme", push)
	}
	if strings.HasSuffix(push, "/") || strings.HasPrefix(push, "/") {
		return fmt.Errorf("push %q is not a bare registry/repository coordinate", push)
	}
	host, path, found := strings.Cut(push, "/")
	if !found || host == "" || path == "" {
		return fmt.Errorf("push %q is not a registry/repository coordinate", push)
	}
	if i := strings.IndexAny(path, "@:"); i >= 0 {
		return fmt.Errorf("push %q carries a tag or digest", push)
	}
	for _, segment := range strings.Split(path, "/") {
		if !ociRepoSegment.MatchString(segment) {
			return fmt.Errorf("push repository segment %q must already be lowercase and OCI-safe", segment)
		}
	}
	return nil
}

// validSourceURL keeps junk out of an annotation every published artifact
// carries. The rule is catalog's, because a project's source.json and a lock's
// fill the same annotation and two definitions would drift.
func validSourceURL(raw string) error {
	return catalog.ValidSourceURL(raw)
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
// repository and platform digest it produced, for pinned and closure-only
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

// Requests returns every pin key, sorted, so callers report in a stable
// order without each sorting for themselves.
func (l *Lock) Requests() []string {
	out := make([]string, 0, len(l.Pins))
	for request := range l.Pins {
		out = append(out, request)
	}
	sort.Strings(out)
	return out
}

// PinnedArtifacts is the set of artifact paths the pins point at
// directly. These are the *roots* of reachability, not reachability itself:
// each one's vendored manifest names dependency edges that pull further
// artifact directories into the closure, and walking those needs to read the
// manifests. Pruning against this set alone would delete the closure.
//
// Remotes are deliberately not roots. An remote says where an artifact can be
// fetched, which is meaningless once nothing pins it.
func (l *Lock) PinnedArtifacts() map[string]bool {
	out := make(map[string]bool, len(l.Pins))
	for _, pin := range l.Pins {
		out[pin.Artifact] = true
	}
	return out
}
