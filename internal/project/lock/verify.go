package lock

import (
	"encoding/json"
	"fmt"
	"os"
	"path/filepath"
	"sort"
	"strings"

	"github.com/Justype/condatainer/internal/artifact/capsule"
	"github.com/Justype/condatainer/internal/artifact/key"
	"github.com/Justype/condatainer/internal/artifact/meta"
)

// MaxSourceBytes bounds any file read out of a lock. A vendored source is a
// recipe or a Conda export; anything at this size is not one, and a lock is
// untrusted input that must not be able to exhaust memory.
const MaxSourceBytes = 8 << 20

// Entry is one verified artifact directory in cnt-lock/artifacts/.
type Entry struct {
	// Path is the artifact directory relative to the lock directory.
	Path     string
	Manifest meta.Manifest
	// Identity and Equiv were regenerated from the vendored sources, never
	// taken on trust from the manifest.
	Identity meta.KeyRef
	Equiv    meta.KeyRef
}

// Verified is a checkout's artifact closure, keyed by relative artifact path.
type Verified struct {
	Entries map[string]*Entry
	// Roots are the artifact paths selections point at, sorted.
	Roots []string
	// Reachable is every entry the roots reach through dependency edges. An
	// entry outside it is unreferenced and may be pruned.
	Reachable map[string]bool
}

// Problem is one validation failure, attributed to what it is about.
type Problem struct {
	// Artifact is the relative artifact path, empty for a project-level issue.
	Artifact string
	// Request is the selection key, empty when the problem is not about one.
	Request string
	Reason  string
}

func (p Problem) String() string {
	switch {
	case p.Artifact != "" && p.Request != "":
		return fmt.Sprintf("%s (%s): %s", p.Request, p.Artifact, p.Reason)
	case p.Artifact != "":
		return fmt.Sprintf("%s: %s", p.Artifact, p.Reason)
	case p.Request != "":
		return fmt.Sprintf("%s: %s", p.Request, p.Reason)
	}
	return p.Reason
}

// Verify validates a lock against the vendored artifacts beside it, using only
// the checkout.
//
// Nothing here reads an installed overlay, a store, the catalog, build
// configuration, or the network, and no payload is needed. That is the CI
// contract: a clean checkout on a machine without CondaTainer images either is
// a complete, internally consistent rebuild specification or is not, and this
// says which.
//
// Problems are collected rather than returned at the first failure, because a
// caller fixing a lock wants the whole list.
func Verify(root string, l *Lock) (*Verified, []Problem) {
	dir := Dir(root)
	out := &Verified{Entries: map[string]*Entry{}}
	var problems []Problem

	present, err := listArtifactDirs(root)
	if err != nil {
		return out, []Problem{{Reason: err.Error()}}
	}
	for _, relative := range present {
		entry, entryProblems := readEntry(dir, relative)
		problems = append(problems, entryProblems...)
		if entry != nil {
			out.Entries[relative] = entry
		}
	}

	// Selections must point at something that verified.
	for _, request := range l.Requests() {
		artifact := l.Selections[request].Artifact
		entry, ok := out.Entries[artifact]
		if !ok {
			problems = append(problems, Problem{Request: request, Artifact: artifact,
				Reason: "selection points at a missing or invalid artifact directory"})
			continue
		}
		if reason := satisfies(request, entry.Manifest.Name); reason != "" {
			problems = append(problems, Problem{Request: request, Artifact: artifact, Reason: reason})
		}
	}
	out.Roots = sortedKeys(l.SelectedArtifacts())

	// Walk dependency edges from the roots. Reachability is the closure, not the
	// selection set: pruning against the roots alone would delete the closure.
	out.Reachable = map[string]bool{}
	var walk func(artifact string, trail []string)
	walk = func(artifact string, trail []string) {
		if out.Reachable[artifact] {
			return
		}
		out.Reachable[artifact] = true
		entry, ok := out.Entries[artifact]
		if !ok {
			return
		}
		for _, dep := range entry.Manifest.Dependencies {
			child, reason := resolveEdge(out.Entries, dep)
			if reason != "" {
				problems = append(problems, Problem{Artifact: artifact, Reason: reason})
				continue
			}
			walk(child, append(trail, artifact))
		}
	}
	for _, artifact := range out.Roots {
		walk(artifact, nil)
	}

	for _, relative := range present {
		if !out.Reachable[relative] {
			problems = append(problems, Problem{Artifact: relative,
				Reason: "artifact directory is not reachable from any selection"})
		}
	}

	// An origin must name an artifact that is actually here.
	for _, artifact := range sortedKeys(originKeys(l)) {
		if _, ok := out.Entries[artifact]; !ok {
			problems = append(problems, Problem{Artifact: artifact,
				Reason: "origin refers to an artifact that is not vendored"})
		}
	}

	sort.Slice(problems, func(i, j int) bool { return problems[i].String() < problems[j].String() })
	return out, problems
}

// resolveEdge finds the entry a dependency edge points at, by name and complete
// identity. An edge whose provenance was never recorded cannot be followed, and
// a lock that contains one is not a rebuild specification.
func resolveEdge(entries map[string]*Entry, dep meta.Dependency) (string, string) {
	if dep.Identity.Empty() || dep.Equiv.Empty() {
		return "", fmt.Sprintf("dependency %s carries no complete keys, so its provenance cannot be followed", dep.Name)
	}
	want := ArtifactPath(capsule.EntryName(dep.Name, dep.Identity.Digest()))
	entry, ok := entries[want]
	if !ok {
		return "", fmt.Sprintf("dependency %s@%s is not vendored at %s", dep.Name, short(dep.Identity.SHA256), want)
	}
	if entry.Manifest.Name != dep.Name {
		return "", fmt.Sprintf("dependency edge names %s but %s records %s", dep.Name, want, entry.Manifest.Name)
	}
	if entry.Identity != dep.Identity {
		return "", fmt.Sprintf("dependency %s expects identity %s but %s regenerates %s",
			dep.Name, short(dep.Identity.SHA256), want, short(entry.Identity.SHA256))
	}
	if entry.Equiv != dep.Equiv {
		return "", fmt.Sprintf("dependency %s expects equivalence %s but %s regenerates %s",
			dep.Name, short(dep.Equiv.SHA256), want, short(entry.Equiv.SHA256))
	}
	return want, ""
}

// readEntry loads and fully verifies one vendored artifact directory.
func readEntry(lockDir, relative string) (*Entry, []Problem) {
	fail := func(format string, args ...any) []Problem {
		return []Problem{{Artifact: relative, Reason: fmt.Sprintf(format, args...)}}
	}
	dir := filepath.Join(lockDir, relative)

	manifestBytes, err := readBounded(filepath.Join(dir, meta.FileName))
	if err != nil {
		return nil, fail("cannot read %s: %v", meta.FileName, err)
	}
	var manifest meta.Manifest
	if err := json.Unmarshal(manifestBytes, &manifest); err != nil {
		return nil, fail("cannot parse %s: %v", meta.FileName, err)
	}
	manifest.Normalize()
	if err := meta.ValidateManifest(manifest); err != nil {
		return nil, fail("invalid manifest: %v", err)
	}

	// The directory name is addressing, never identity: the complete key is
	// regenerated and the name and prefix are checked against it.
	if want := ArtifactPath(capsule.EntryName(manifest.Name, manifest.Keys.Identity.Digest())); want != relative {
		return nil, fail("directory name does not match its manifest; expected %s", want)
	}

	sources, err := readSources(dir, manifest)
	if err != nil {
		return nil, fail("%v", err)
	}
	derived, err := key.Verify(manifest, sources)
	if err != nil {
		return nil, fail("sources do not regenerate the recorded keys: %v", err)
	}

	return &Entry{Path: relative, Manifest: manifest,
		Identity: derived.Identity.Ref, Equiv: derived.Equiv.Ref}, nil
}

// Sources reads one verified entry's vendored build inputs, keyed as
// manifest.source.files names them. A rebuild starts from exactly these bytes.
func Sources(root string, entry *Entry) (map[string][]byte, error) {
	sources, err := readSources(filepath.Join(Dir(root), entry.Path), entry.Manifest)
	if err != nil {
		return nil, fmt.Errorf("%s: %w", entry.Path, err)
	}
	return sources, nil
}

// readSources reads exactly the source set the manifest's build type requires.
// A recipe build must carry its recipe byte for byte, and a Conda build both of
// its exports; a URL or a catalog reference is never a substitute.
func readSources(dir string, manifest meta.Manifest) (key.Sources, error) {
	sources := key.Sources{}
	for _, name := range manifest.Source.Files {
		if name != filepath.Base(name) || name == "." || name == ".." {
			return nil, fmt.Errorf("source file %q is not a plain name", name)
		}
		data, err := readBounded(filepath.Join(dir, name))
		if err != nil {
			return nil, fmt.Errorf("cannot read source %s: %w", name, err)
		}
		sources[name] = data
	}
	return sources, nil
}

// readBounded reads a lock file, refusing a symlink and anything oversized.
// A vendored source is untrusted input: it must not escape the checkout through
// a link and must not be able to exhaust memory.
func readBounded(path string) ([]byte, error) {
	info, err := os.Lstat(path)
	if err != nil {
		return nil, err
	}
	if !info.Mode().IsRegular() {
		return nil, fmt.Errorf("%s is not a regular file", filepath.Base(path))
	}
	if info.Size() > MaxSourceBytes {
		return nil, fmt.Errorf("%s is %d bytes, over the %d limit", filepath.Base(path), info.Size(), MaxSourceBytes)
	}
	return os.ReadFile(path)
}

// listArtifactDirs returns every artifact directory present, relative to the
// lock directory and sorted. A non-directory or symlink is refused rather than
// skipped: it is in the tracked tree and something put it there.
func listArtifactDirs(root string) ([]string, error) {
	entries, err := os.ReadDir(ArtifactsPath(root))
	if os.IsNotExist(err) {
		return nil, nil
	}
	if err != nil {
		return nil, err
	}
	var out []string
	for _, entry := range entries {
		if entry.Type()&os.ModeSymlink != 0 {
			return nil, fmt.Errorf("%s/%s is a symlink", ArtifactsDir, entry.Name())
		}
		// A dotted directory is a write in progress, not an artifact: skipping it
		// keeps an interrupted write invisible instead of a spurious problem.
		if !entry.IsDir() || strings.HasPrefix(entry.Name(), ".") {
			continue
		}
		out = append(out, ArtifactPath(entry.Name()))
	}
	sort.Strings(out)
	return out, nil
}

// satisfies reports why a concrete artifact name does not answer a request, or
// "". A request names one exact artifact, so the two names must be equal.
func satisfies(request, name string) string {
	if strings.HasPrefix(request, PathPrefix) {
		return ""
	}
	dep, err := parseRequest(request)
	if err != nil {
		return fmt.Sprintf("selection key is not a usable dependency: %v", err)
	}
	if reason := ConstraintReason(dep, request); reason != "" {
		return reason
	}
	if dep.NameVersion() != name {
		return fmt.Sprintf("request names %s but the artifact is %s", dep.NameVersion(), name)
	}
	return ""
}

func originKeys(l *Lock) map[string]bool {
	out := make(map[string]bool, len(l.Origins))
	for artifact := range l.Origins {
		out[artifact] = true
	}
	return out
}

func sortedKeys(set map[string]bool) []string {
	out := make([]string, 0, len(set))
	for key := range set {
		out = append(out, key)
	}
	sort.Strings(out)
	return out
}

func short(sha string) string {
	if len(sha) > capsule.IdentityChars {
		return sha[:capsule.IdentityChars]
	}
	return sha
}
