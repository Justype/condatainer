package project

import (
	"fmt"
	"os"
	"path/filepath"

	"github.com/Justype/condatainer/catalog"
	"github.com/Justype/condatainer/internal/artifact/compare"
	"github.com/Justype/condatainer/internal/artifact/meta"
	"github.com/Justype/condatainer/internal/project/lock"
	"github.com/Justype/condatainer/internal/store"
	"github.com/Justype/condatainer/internal/utils"
)

// Match is which of an artifact's two keys a locally available copy has to
// agree with.
//
// The two keys exist to answer different questions, and this is where the
// question gets asked: identity is "which exact build is this?", equivalence is
// "can this substitute for the requested artifact?".
type Match string

const (
	// MatchEquivalent accepts anything that can substitute for what was locked.
	// The default, because it is what the equivalence key exists to decide.
	//
	// It is also the only workable rule for ordinary analysis. A Conda solve
	// that lands on the same versions but different build strings or package
	// URLs writes a different explicit.txt and an identical environment.yml, so
	// it is equivalent and not identical. Requiring identity would fail that,
	// and would additionally require every data artifact to be rebuilt against
	// the exact tool identities it was first built with, since a script
	// identity hashes its dependencies' identities.
	MatchEquivalent Match = "equivalent"
	// MatchIdentity accepts only the exact build the lock names.
	//
	// Deliberately stringent: a Conda artifact has to replay its explicit.txt to
	// the byte, and a data artifact has to be rebuilt in the same environment.
	// That is the point when bit-level provenance is the requirement, and too
	// much to ask when it is not.
	MatchIdentity Match = "identity"
)

// Normalize fills in the default, so a zero value is the ordinary mode rather
// than the strict one.
func (m Match) Normalize() Match {
	if m == MatchIdentity {
		return MatchIdentity
	}
	return MatchEquivalent
}

// LookupFunc finds a local artifact answering to keys under a match mode.
type LookupFunc func(name string, keys meta.Keys, match Match, dirs []string) (store.Candidate, bool)

// LookupAtFunc reports whether the artifact at one exact path satisfies keys.
type LookupAtFunc func(path, name string, keys meta.Keys, match Match) (store.Candidate, bool)

// LookupLocal finds the locked artifact, or under MatchEquivalent something
// that can stand in for it.
//
// The exact identity is tried first in both modes, so an available exact copy
// is never passed over for a substitute and a caller can report which it got.
// Any equivalent candidate is as good as any other by definition, so the sorted
// first is taken for a stable answer rather than for being better.
func LookupLocal(name string, keys meta.Keys, match Match, dirs []string) (store.Candidate, bool) {
	exact := store.IdentityQuery{Scheme: keys.Identity.Scheme, SHA256: keys.Identity.SHA256}
	if candidate, _, err := store.ResolveIdentity(name, exact, dirs); err == nil {
		return candidate, true
	}
	if match.Normalize() == MatchIdentity {
		return store.Candidate{}, false
	}
	report, err := store.Equivalent(name,
		store.IdentityQuery{Scheme: keys.Equiv.Scheme, SHA256: keys.Equiv.SHA256}, dirs)
	if err != nil || len(report.Candidates) == 0 {
		return store.Candidate{}, false
	}
	store.SortReport(&report)
	return report.Candidates[0], true
}

// LookupInput finds a local artifact that may satisfy one dependency edge of a
// rebuild.
//
// This is a different question from LookupLocal, which asks whether a *pinned*
// artifact is here. A build dependency is not pinned: what may stand in for it is
// whatever leaves the dependent's equivalence key unchanged, and the role the
// edge records is exactly the scheme's statement about that.
//
//	RoleData    the dependency's equivalence key is in the preimage → same equiv
//	RoleApp     its name/version is in the preimage         → same version, any build
//	RoleHistory nothing is in the preimage                  → any version of that name
//
// The exact recorded identity is preferred in every case, because it yields the
// exact identity for the dependent too, which is the result that satisfies every
// other project that locked it. Under MatchIdentity nothing else is accepted:
// ScriptIdentityV1 hashes dependency identities, so any substitution would
// produce a dependent the mode rejects anyway.
func LookupInput(dep meta.Dependency, match Match, dirs []string) (store.Candidate, bool) {
	exact := store.IdentityQuery{Scheme: dep.Identity.Scheme, SHA256: dep.Identity.SHA256}
	if candidate, _, err := store.ResolveIdentity(dep.Name, exact, dirs); err == nil {
		return candidate, true
	}
	if match.Normalize() == MatchIdentity {
		return store.Candidate{}, false
	}

	if dep.Role == meta.RoleData {
		report, err := store.Equivalent(dep.Name,
			store.IdentityQuery{Scheme: dep.Equiv.Scheme, SHA256: dep.Equiv.SHA256}, dirs)
		if err != nil || len(report.Candidates) == 0 {
			return store.Candidate{}, false
		}
		return report.Candidates[0], true
	}

	// Any build of the recorded name/version, which is all an app contributes and
	// is closer to the record than another version would be.
	name := catalog.Normalize(dep.Name)
	report := store.Scan(store.ScanOptions{Dirs: dirs, Name: name})
	if len(report.Candidates) > 0 {
		return report.Candidates[0], true
	}
	if dep.Role != meta.RoleHistory {
		return store.Candidate{}, false
	}
	return newestOf(toolName(name), dirs)
}

// toolName drops the version component, so "samtools/1.21" matches any samtools.
func toolName(nameVersion string) string {
	parsed, err := catalog.ParseDep(nameVersion)
	if err != nil {
		return nameVersion
	}
	return parsed.Name
}

// newestOf picks the highest installed version of one tool.
//
// Scan visits the image roots in order, so an earlier root's copy wins a version
// tie and the nearest one is preferred. Newest rather than first because any
// version is correct by the scheme, and the newest is what a person reaching for
// the tool would expect.
func newestOf(tool string, dirs []string) (store.Candidate, bool) {
	var best store.Candidate
	var bestVersion string
	for _, candidate := range store.Scan(store.ScanOptions{Dirs: dirs}).Candidates {
		parsed, err := catalog.ParseDep(candidate.Name)
		if err != nil || parsed.Name != tool || parsed.Version == "" {
			continue
		}
		if bestVersion == "" || higher(parsed.Version, bestVersion) {
			best, bestVersion = candidate, parsed.Version
		}
	}
	return best, bestVersion != ""
}

// higher reports whether version a sorts above b. Equal versions are not higher,
// which is what leaves the nearest root's copy in place on a tie.
func higher(a, b string) bool {
	return a != b && utils.SortVersionsDescending([]string{a, b})[0] == a
}

// LookupAt reports whether the artifact a project path declares is there and
// satisfies the match mode.
//
// The path is the whole lookup. A project-path artifact is project-owned, so a
// copy of it in an images root is not a substitute — mounting that would leave
// the declared path empty and the script still broken.
//
// The filename is read as an address and never as a name. That rule holds only
// here: in an images root or the store the filename *is* how an artifact is
// found, so a scan there requires it to encode the artifact's own name. A path
// artifact is found by the path the project declared, and its name comes from
// the script's #TARGET: — deliberately unrelated to what the file is called — so
// decoding the filename would reject every artifact that used one.
func LookupAt(path, name string, keys meta.Keys, match Match) (store.Candidate, bool) {
	info, err := os.Lstat(path)
	if err != nil || !info.Mode().IsRegular() {
		return store.Candidate{}, false
	}
	artifact, err := compare.Read(path)
	if err != nil || artifact.Name != name {
		return store.Candidate{}, false
	}
	candidate := store.Candidate{
		Name: artifact.Name, Path: path, Root: filepath.Dir(path),
		Layout: store.LayoutFlat, Size: info.Size(),
		Identity: artifact.IdentityRef(), Equiv: artifact.EquivRef(),
	}
	if candidate.Identity == keys.Identity {
		return candidate, true
	}
	if match.Normalize() == MatchEquivalent && candidate.Equiv == keys.Equiv {
		return candidate, true
	}
	return store.Candidate{}, false
}

// Mount is one declaration resolved to something the runtime can mount.
type Mount struct {
	// Request is the pin key this answers.
	Request string
	// Name is the artifact name, empty for an unpinnable declaration.
	Name string
	// Path is absolute. Inside a project a relative declaration is relative to
	// the project root, and the runtime resolves a relative overlay against the
	// process working directory — which a scheduler chooses. Resolving here is
	// what keeps the two from disagreeing.
	Path string
	// Identity is what the lock records; Found is what is actually there, set
	// only when an equivalent artifact stands in.
	Identity string
	Found    string
	// Unpinned marks a declaration that carries no pin by design and is
	// mounted as the literal path it names.
	Unpinned bool
}

// Unresolved is one declaration that cannot be mounted, and why.
type Unresolved struct {
	Request string
	Reason  string
}

// Resolution is everything one script's declarations resolve to.
type Resolution struct {
	Root   string
	Mounts []Mount
	// Unresolved is non-empty when the project cannot run as locked.
	Unresolved []Unresolved
}

// Complete reports whether every declaration resolved.
func (r *Resolution) Complete() bool { return len(r.Unresolved) == 0 }

// ResolveOptions tunes resolution.
type ResolveOptions struct {
	// Match is which key a local copy must agree with. Empty means equivalent.
	Match Match
	// SearchDirs overrides the configured image roots.
	SearchDirs []string
	// lookup and lookupAt are injected for tests.
	lookup   LookupFunc
	lookupAt LookupAtFunc
}

func (o ResolveOptions) resolver() LookupFunc {
	if o.lookup != nil {
		return o.lookup
	}
	return LookupLocal
}

func (o ResolveOptions) pathResolver() LookupAtFunc {
	if o.lookupAt != nil {
		return o.lookupAt
	}
	return LookupAt
}

// Resolve turns one script's declarations into absolute paths to mount, using
// the project's lock and nothing else.
//
// This is the single entry point execution and restore share, so they cannot
// disagree about what a lock means. It acquires nothing: a declaration whose
// artifact is absent is reported unresolved, never fetched or built. Compute
// nodes routinely lack the network, credentials, build tools and writable
// images directories that acquiring would need.
//
// Inside a project an unsatisfiable declaration is always an error. There is no
// falling back to whatever currently answers to the name — that is the failure
// the lock exists to prevent.
func Resolve(root string, l *lock.Lock, requests []lock.Request, opts ResolveOptions) (*Resolution, error) {
	root, err := filepath.Abs(root)
	if err != nil {
		return nil, err
	}
	verified, problems := lock.Verify(root, l)
	if len(problems) > 0 {
		resolution := &Resolution{Root: root}
		for _, problem := range problems {
			resolution.Unresolved = append(resolution.Unresolved,
				Unresolved{Request: problem.Request, Reason: problem.Reason})
		}
		return resolution, nil
	}

	match := opts.Match.Normalize()
	resolve, resolveAt := opts.resolver(), opts.pathResolver()
	resolution := &Resolution{Root: root}

	for _, request := range requests {
		// An unpinned declaration is mounted as written.
		//
		// Both halves matter. A kind that cannot be pinned has nowhere for a
		// pin to point, and G13 requires it to carry the marker — the
		// scanner has already reported one that does not. The marker on a
		// *pinnable* kind is different: it is an author deliberately opting out
		// of pinning something that could have been pinned, which is how a
		// per-project scratch overlay stays out of the lock. Honouring it here
		// is what keeps `project validate` and `run` agreeing, since validate
		// already stops asking for a pin once the marker is present.
		if request.Unpinned || !request.Kind.Pinnable() {
			resolution.Mounts = append(resolution.Mounts, Mount{
				Request: request.Key, Path: literalPath(root, request), Unpinned: true,
			})
			continue
		}
		pin, ok := l.Pins[request.Key]
		if !ok {
			resolution.Unresolved = append(resolution.Unresolved, Unresolved{Request: request.Key,
				Reason: "declared but not pinned; run `condatainer project pin` to pin it"})
			continue
		}
		entry, ok := verified.Entries[pin.Artifact]
		if !ok {
			resolution.Unresolved = append(resolution.Unresolved, Unresolved{Request: request.Key,
				Reason: "the pinned artifact is not vendored"})
			continue
		}
		keys := meta.Keys{Identity: entry.Identity, Equiv: entry.Equiv}
		mount := Mount{Request: request.Key, Name: entry.Manifest.Name, Identity: entry.Identity.Digest()}

		var candidate store.Candidate
		if request.Kind == lock.KindPath {
			at := filepath.Join(root, filepath.FromSlash(request.Path))
			candidate, ok = resolveAt(at, entry.Manifest.Name, keys, match)
		} else {
			candidate, ok = resolve(entry.Manifest.Name, keys, match, opts.SearchDirs)
		}
		if !ok {
			resolution.Unresolved = append(resolution.Unresolved, Unresolved{Request: request.Key,
				Reason: fmt.Sprintf("%s is locked at %s but no %s copy is available here",
					entry.Manifest.Name, short(entry.Identity.Digest()), match)})
			continue
		}
		mount.Path = candidate.Path
		if candidate.Identity != entry.Identity {
			mount.Found = candidate.Identity.Digest()
		}
		resolution.Mounts = append(resolution.Mounts, mount)
	}
	return resolution, nil
}

// literalPath is where an unpinnable declaration points. A project-relative one
// is anchored on the root like every other path in a project; an external one
// already says where it is.
func literalPath(root string, request lock.Request) string {
	if filepath.IsAbs(request.Path) {
		return request.Path
	}
	return filepath.Join(root, filepath.FromSlash(request.Path))
}

func short(digest string) string {
	if len(digest) > 19 {
		return digest[:19]
	}
	return digest
}
