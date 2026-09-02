package lock

import (
	"errors"
	"fmt"
	"os"
	"path/filepath"
	"sort"
	"strings"

	"github.com/Justype/condatainer/internal/artifact/capsule"
	"github.com/Justype/condatainer/internal/artifact/compare"
	"github.com/Justype/condatainer/internal/artifact/meta"
	"github.com/Justype/condatainer/internal/image"
	"github.com/Justype/condatainer/internal/store"
	"github.com/Justype/condatainer/internal/utils"
)

// ErrNoCandidate reports that nothing local answers a pin.
var ErrNoCandidate = errors.New("no artifact matches")

// PinOptions tunes candidate discovery.
type PinOptions struct {
	// SearchDirs overrides the configured image roots.
	SearchDirs []string
}

// Pinned is what one pin resolved to.
type Pinned struct {
	// Request is the canonical pin key.
	Request string
	// Artifact is the vendored directory, relative to the lock directory.
	Artifact string
	// Path is the image the records came out of. Recorded nowhere: it is
	// machine-local, and reported only so a user can see what was read.
	Path     string
	Name     string
	Identity meta.KeyRef
	// Vendored are every artifact directory this pin wrote or confirmed,
	// the pinned artifact plus its transitive closure, sorted.
	Vendored []string
	// Others are the installed identities answering the same name that this did
	// not take, in search order. Empty unless PinAll filled them in.
	Others []store.Candidate
	// Manual records the pin as one the project stated rather than one a #DEP:
	// produced, so a rescan cannot sweep it. Set by the caller, which is the only
	// layer that knows whether anything declares the request.
	Manual bool
}

// Pin resolves one request to an exact local artifact, vendors its complete
// source closure, and returns what to record. It does not publish: the caller
// updates the lock and publishes last, so a failure leaves the old lock intact.
//
// target is an identity — `scheme@sha256:<hex>`, a full SHA, or an unambiguous
// prefix. It is optional for a name request, which otherwise takes what the
// search order would mount, and refused for a path request, which already says
// which file it means.
//
// A file is never named here. Nothing downstream could use one: restore and push
// both find a pinned artifact through the image roots and the store, so
// pointing at a file elsewhere would record an identity only a rebuild could
// satisfy.
func Pin(root, request, target string, opts PinOptions) (*Pinned, error) {
	request, err := CanonicalRequest(request)
	if err != nil {
		return nil, err
	}
	if reason := pinnableRequest(request); reason != "" {
		return nil, fmt.Errorf("%w: %s", ErrInvalid, reason)
	}
	candidate, err := resolveCandidate(root, request, target, opts)
	if err != nil {
		return nil, err
	}
	if reason := satisfies(request, candidate.Name); reason != "" {
		return nil, fmt.Errorf("%w: %s", ErrInvalid, reason)
	}

	staging, err := os.MkdirTemp("", "cnt-select-")
	if err != nil {
		return nil, err
	}
	defer os.RemoveAll(staging) //nolint:errcheck

	// One extraction, not one read per file: every archive read spawns a process.
	if err := image.ExtractDir(candidate.Path, "/"+meta.DirName, staging); err != nil {
		return nil, fmt.Errorf("cannot read metadata from %s: %w", candidate.Path, err)
	}
	metaDir := filepath.Join(staging, meta.DirName)

	vendored, err := vendorClosure(root, metaDir, candidate)
	if err != nil {
		return nil, err
	}
	sort.Strings(vendored)

	return &Pinned{
		Request:  request,
		Artifact: EntryPath(capsule.EntryName(candidate.Name, candidate.Identity.Digest())),
		Path:     candidate.Path,
		Name:     candidate.Name,
		Identity: candidate.Identity,
		Vendored: vendored,
	}, nil
}

// CanonicalRequest normalizes what a caller typed into the key a scan produces.
//
// It runs the same grammar a #DEP: does, so `overlays/tool.sqf` on the command
// line and in a script cannot mean different things, and an extension is what
// separates a path from a name. The lock's own `path:` key is taken as written.
// Nothing on disk is read, so a request naming a file that is gone still
// normalizes.
func CanonicalRequest(request string) (string, error) {
	request = strings.TrimSpace(request)
	if request == "" {
		return "", fmt.Errorf("%w: empty request", ErrInvalid)
	}
	if strings.HasPrefix(request, PathPrefix) {
		return request, nil
	}
	parsed, reason := ParseDeclaration(request)
	if reason != "" {
		return "", fmt.Errorf("%w: %s", ErrInvalid, reason)
	}
	return parsed.Key, nil
}

// resolveCandidate finds the one artifact a request pins.
//
// A path request with no identity reads the file it names, which is the only way
// to select one: the project may call that file anything, while identity lookup
// applies the flat-name rule and would refuse overlays/combined.sqf for not
// being testdata--combined--1.0.sqf.
func resolveCandidate(root, request, target string, opts PinOptions) (store.Candidate, error) {
	target = strings.TrimSpace(target)
	if destination, isPath := strings.CutPrefix(request, PathPrefix); isPath {
		if target != "" {
			return store.Candidate{}, fmt.Errorf(
				"%w: %s already names the file it means, so it takes no identity", ErrInvalid, request)
		}
		return candidateFromPath(filepath.Join(root, filepath.FromSlash(destination)))
	}
	if target == "" {
		candidate, _, err := resolveByName(request, opts)
		return candidate, err
	}
	if looksLikePath(target) {
		return store.Candidate{}, fmt.Errorf(
			"%w: %s is a file, and a pin names an identity. Install it into an images root to select it by identity, or keep it in the project and declare it as a path",
			ErrInvalid, target)
	}
	query, err := store.ParseIdentityQuery(target)
	if err != nil {
		return store.Candidate{}, err
	}
	name, err := requestName(request)
	if err != nil {
		return store.Candidate{}, err
	}
	// Every copy in every readable root is a candidate, not just the nearest
	// one: a pin names an exact artifact, and nearest-first name
	// resolution would hide the copy the user meant.
	candidate, _, err := store.ResolveIdentity(name, query, opts.SearchDirs)
	if err != nil {
		if errors.Is(err, store.ErrNotFound) {
			return store.Candidate{}, fmt.Errorf("%w: %s at %s", ErrNoCandidate, name, target)
		}
		return store.Candidate{}, err
	}
	return candidate, nil
}

// resolveByName finds the artifact a name request means with no identity given,
// and reports the other installed identities for that name.
//
// Search order is the configured one, so a project pins what the same name would
// mount outside a project. Within one root the flat name--version file wins: the
// store's destination rule puts exactly one identity there and the rest under
// store/, which makes "the obvious one" a fact rather than a guess.
func resolveByName(request string, opts PinOptions) (store.Candidate, []store.Candidate, error) {
	name, err := requestName(request)
	if err != nil {
		return store.Candidate{}, nil, err
	}
	report := store.Scan(store.ScanOptions{Dirs: opts.SearchDirs, Name: name})
	if len(report.Candidates) == 0 {
		return store.Candidate{}, nil, fmt.Errorf("%w: %s is not installed", ErrNoCandidate, name)
	}

	chosen := report.Candidates[0]
	for _, candidate := range report.Candidates {
		if candidate.Layout == store.LayoutFlat {
			chosen = candidate
			break
		}
	}
	var others []store.Candidate
	seen := map[meta.KeyRef]bool{chosen.Identity: true}
	for _, candidate := range report.Candidates {
		if seen[candidate.Identity] {
			continue
		}
		seen[candidate.Identity] = true
		others = append(others, candidate)
	}
	return chosen, others, nil
}

// looksLikePath reports whether a target addresses a file rather than an
// identity. An identity never contains a separator — `scheme@sha256:<hex>` and a
// bare digest are both flat — so anything that does is a path, as is anything
// carrying an image extension. Routing those here means a .sif is refused for
// being a .sif rather than for failing to parse as a digest.
func looksLikePath(target string) bool {
	return strings.ContainsRune(target, filepath.Separator) ||
		utils.IsOverlay(target) || utils.IsSif(target)
}

// candidateFromPath verifies the .sqf a path request names. A writable .img has
// no identity to pin and a .sif is a container root, so neither can be pinned.
func candidateFromPath(target string) (store.Candidate, error) {
	if !strings.HasSuffix(target, ".sqf") {
		return store.Candidate{}, fmt.Errorf("%w: only .sqf can be pinned, not %s", ErrInvalid, filepath.Base(target))
	}
	absolute, err := filepath.Abs(target)
	if err != nil {
		return store.Candidate{}, err
	}
	// Read the file, never scan its directory. A scan applies the flat-name
	// rule — the filename must encode the artifact name — which is right where
	// the filename *is* the address and wrong here: a path pin exists so a
	// project can call the file whatever it likes, and overlays/combined.sqf
	// would be refused for not being testdata--combined--1.0.sqf. This is the
	// same rule project.LookupAt applies when it later resolves the pin,
	// and the two have to agree or a lock could be written and never resolved.
	info, err := os.Lstat(absolute)
	if err != nil {
		return store.Candidate{}, fmt.Errorf("%w: %s: %v", ErrInvalid, absolute, err)
	}
	if !info.Mode().IsRegular() {
		return store.Candidate{}, fmt.Errorf("%w: %s is not a regular file", ErrInvalid, absolute)
	}
	artifact, err := compare.Read(absolute)
	if err != nil {
		return store.Candidate{}, fmt.Errorf("%w: %s: %v", ErrInvalid, absolute, err)
	}
	identity, equiv := artifact.IdentityRef(), artifact.EquivRef()
	if identity.Empty() || equiv.Empty() {
		return store.Candidate{}, fmt.Errorf("%w: %s carries no verifiable keys", ErrNoCandidate, absolute)
	}
	return store.Candidate{
		Name: artifact.Name, Path: absolute, Root: filepath.Dir(absolute),
		Layout: store.LayoutFlat, Size: info.Size(),
		Identity: identity, Equiv: equiv,
	}, nil
}

// pinnableRequest reports why a pin key cannot be pinned, or "". It says what
// closes the gap, since neither alternative is obvious from the refusal.
func pinnableRequest(request string) string {
	destination, isPath := strings.CutPrefix(request, PathPrefix)
	if !isPath {
		return ""
	}
	if utils.IsImg(destination) {
		return fmt.Sprintf("%s is writable and has no identity to pin; freeze it with `condatainer overlay freeze %s overlays/<name>.sqf` and pin that",
			destination, destination)
	}
	if err := validProjectPath(destination); err != nil {
		return fmt.Sprintf("%v; restore cannot own this path, so copy it under the project and pin that", err)
	}
	return ""
}

// requestName is the artifact name a request addresses, for candidate lookup.
func requestName(request string) (string, error) {
	dep, err := parseRequest(request)
	if err != nil {
		return "", fmt.Errorf("%w: %v", ErrInvalid, err)
	}
	if why := ConstraintReason(dep, request); why != "" {
		return "", fmt.Errorf("%w: %s", ErrInvalid, why)
	}
	return dep.NameVersion(), nil
}

// vendorClosure writes the pinned artifact and every entry of its embedded
// capsule into cnt-lock/provenance/.
//
// The capsule is copied across, never re-derived: it is already the transitive
// union its builds composed, so there is no recursive resolution, no catalog
// access, and no network. Entry directory names are `capsule.EntryName`, which
// is what an image's own provenance uses too, so the closure maps across
// one-to-one.
func vendorClosure(root, metaDir string, candidate store.Candidate) ([]string, error) {
	var vendored []string

	manifest, err := readStagedManifest(metaDir)
	if err != nil {
		return nil, err
	}
	if manifest.Name != candidate.Name {
		return nil, fmt.Errorf("%w: %s records the name %s", ErrInvalid, candidate.Path, manifest.Name)
	}
	if manifest.Keys.Identity != candidate.Identity {
		return nil, fmt.Errorf("%w: %s records an identity its sources do not derive", ErrInvalid, candidate.Path)
	}
	// A pin is only as good as its provenance: an artifact built over an
	// unrecorded dependency cannot be rebuilt from the checkout, whatever its
	// own sources say.
	if manifest.ProvenanceComplete != nil && !*manifest.ProvenanceComplete {
		return nil, fmt.Errorf("%w: %s was built over a dependency that carried no records",
			ErrInvalid, candidate.Name)
	}

	files, err := stagedFiles(metaDir, capsule.FileNames(manifest))
	if err != nil {
		return nil, err
	}
	relative, err := StageEntry(root, capsule.EntryName(manifest.Name, manifest.Keys.Identity.Digest()), files)
	if err != nil {
		return nil, err
	}
	vendored = append(vendored, relative)

	capsuleDir := filepath.Join(metaDir, capsule.DirName)
	entries, err := capsule.Entries(capsuleDir)
	if err != nil {
		return nil, err
	}
	for _, entry := range entries {
		body, err := stagedFiles(filepath.Join(capsuleDir, entry.Dir), entry.Files)
		if err != nil {
			return nil, err
		}
		relative, err := StageEntry(root, entry.Dir, body)
		if err != nil {
			return nil, err
		}
		vendored = append(vendored, relative)
	}
	return vendored, nil
}

func readStagedManifest(metaDir string) (meta.Manifest, error) {
	manifest, err := capsule.ReadManifest(metaDir, readBounded)
	if err != nil {
		return meta.Manifest{}, fmt.Errorf("%w: %v", ErrInvalid, err)
	}
	if manifest.Keys.Identity.Empty() || manifest.Keys.Equiv.Empty() {
		return meta.Manifest{}, fmt.Errorf("%w: the artifact records no complete keys", ErrInvalid)
	}
	return manifest, nil
}

// stagedFiles reads a fixed file set out of an extracted directory. runtime.json
// is never among them: it is not part of any hashed preimage, so a lock that
// carried it would be carrying something it cannot verify. The name set comes
// from capsule.FileNames, so the lock stages exactly what an image carries.
func stagedFiles(dir string, names []string) (map[string][]byte, error) {
	out := make(map[string][]byte, len(names))
	for _, name := range names {
		if name != filepath.Base(name) || name == "." || name == ".." {
			return nil, fmt.Errorf("%w: source file %q is not a plain name", ErrInvalid, name)
		}
		if name == meta.RuntimeFileName {
			continue
		}
		data, err := readBounded(filepath.Join(dir, name))
		if err != nil {
			return nil, fmt.Errorf("%w: cannot read %s: %v", ErrInvalid, name, err)
		}
		out[name] = data
	}
	return out, nil
}

// Apply records a resolved pin and publishes the lock, verifying the whole
// closure first. Publication is last, so a lock is never left pointing at
// something that does not validate.
func Apply(root string, l *Lock, pinned *Pinned) error {
	previous, had := l.Pins[pinned.Request]
	entry := PinEntry{Artifact: pinned.Artifact, Manual: pinned.Manual}
	// Re-pinning never turns a manual pin back into a sweepable one: the artifact
	// changed, not the reason it is in the lock.
	if had && previous.Manual {
		entry.Manual = true
	}
	l.Pins[pinned.Request] = entry

	if _, problems := Verify(root, l); len(problems) > 0 {
		if had {
			l.Pins[pinned.Request] = previous
		} else {
			delete(l.Pins, pinned.Request)
		}
		return fmt.Errorf("%w: pinning %s leaves the project invalid:\n  %s",
			ErrInvalid, pinned.Request, joinProblems(problems))
	}
	return Publish(root, l)
}

func joinProblems(problems []Problem) string {
	out := make([]string, 0, len(problems))
	for _, problem := range problems {
		out = append(out, problem.String())
	}
	return strings.Join(out, "\n  ")
}

// Reconcile rescans a project and drops pins nothing requests any more,
// reporting what remains needPin. It never invents a pin: choosing an
// artifact is an explicit act.
func Reconcile(root string, l *Lock, result *ScanResult) (needPin []Request) {
	requested := make(map[string]bool, len(result.Requests))
	for _, request := range result.Requests {
		requested[request.Key] = true
	}
	for _, key := range l.Requests() {
		// A manual pin is recorded by the project, not derived from its scripts,
		// so a scan cannot speak to whether it belongs.
		if l.Pins[key].Manual {
			continue
		}
		if !requested[key] {
			delete(l.Pins, key)
		}
	}

	verified, _ := Verify(root, l)
	for _, request := range result.Requests {
		// An unpinnable request is not needPin: there is nowhere for restore
		// to put an answer, and the scanner has already reported it.
		if !request.Kind.Pinnable() {
			continue
		}
		pin, ok := l.Pins[request.Key]
		if !ok {
			needPin = append(needPin, request)
			continue
		}
		if _, valid := verified.Entries[pin.Artifact]; !valid {
			delete(l.Pins, request.Key)
			needPin = append(needPin, request)
		}
	}
	return needPin
}

// PinAll resolves every request in unpinned and records what it found, so a lock
// says which exact artifact satisfies each declaration rather than only which
// declarations exist.
//
// Each request is applied and published on its own, so an artifact that is not
// installed costs its own line in the report and nothing else: what was pinned
// stays pinned, and re-running after `create` finishes the job.
func PinAll(root string, l *Lock, needPin []Request, opts PinOptions) (pinned []*Pinned, failed []error) {
	for _, request := range needPin {
		var others []store.Candidate
		if !strings.HasPrefix(request.Key, PathPrefix) {
			// Reported rather than refused: the flat name is what this name
			// resolves to everywhere else, so pinning it is the answer that
			// agrees with the rest of the system.
			if _, rest, err := resolveByName(request.Key, opts); err == nil {
				others = rest
			}
		}
		entry, err := Pin(root, request.Key, "", opts)
		if err != nil {
			failed = append(failed, err)
			continue
		}
		if err := Apply(root, l, entry); err != nil {
			failed = append(failed, err)
			continue
		}
		entry.Others = others
		pinned = append(pinned, entry)
	}
	return pinned, failed
}
