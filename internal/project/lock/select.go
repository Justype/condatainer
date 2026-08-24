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

// ErrNoCandidate reports that nothing local answers a selection.
var ErrNoCandidate = errors.New("no artifact matches")

// SelectOptions tunes candidate discovery.
type SelectOptions struct {
	// SearchDirs overrides the configured image roots.
	SearchDirs []string
}

// Selected is what one selection resolved to.
type Selected struct {
	// Request is the canonical selection key.
	Request string
	// Artifact is the vendored directory, relative to the lock directory.
	Artifact string
	// Path is the image the records came out of. Recorded nowhere: it is
	// machine-local, and reported only so a user can see what was read.
	Path     string
	Name     string
	Identity meta.KeyRef
	// Vendored are every artifact directory this selection wrote or confirmed,
	// the selected artifact plus its transitive closure, sorted.
	Vendored []string
}

// Select resolves one request to an exact local artifact, vendors its complete
// source closure, and returns what to record. It does not publish: the caller
// updates the lock and publishes last, so a failure leaves the old lock intact.
//
// target is an identity — `scheme@sha256:<hex>`, a full SHA, or an unambiguous
// prefix — or a path to an immutable .sqf. A path may point outside the
// project; what gets recorded is the identity, never the path.
func Select(root, request, target string, opts SelectOptions) (*Selected, error) {
	if strings.TrimSpace(request) == "" {
		return nil, fmt.Errorf("%w: empty request", ErrInvalid)
	}
	if reason := pinnableRequest(request); reason != "" {
		return nil, fmt.Errorf("%w: %s", ErrInvalid, reason)
	}
	candidate, err := resolveCandidate(request, target, opts)
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

	return &Selected{
		Request:  request,
		Artifact: EntryPath(capsule.EntryName(candidate.Name, candidate.Identity.Digest())),
		Path:     candidate.Path,
		Name:     candidate.Name,
		Identity: candidate.Identity,
		Vendored: vendored,
	}, nil
}

// resolveCandidate finds the one artifact a target names, by path or identity.
func resolveCandidate(request, target string, opts SelectOptions) (store.Candidate, error) {
	if looksLikePath(target) {
		return candidateFromPath(target)
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
	// one: a selection names an exact artifact, and nearest-first name
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

// looksLikePath reports whether a target addresses a file rather than an
// identity. An identity never contains a separator — `scheme@sha256:<hex>` and a
// bare digest are both flat — so anything that does is a path, as is anything
// carrying an image extension. Routing those here means a .sif is refused for
// being a .sif rather than for failing to parse as a digest.
func looksLikePath(target string) bool {
	return strings.ContainsRune(target, filepath.Separator) ||
		utils.IsOverlay(target) || utils.IsSif(target)
}

// candidateFromPath verifies an explicit .sqf, which may live outside any images
// root. A writable .img has no identity to pin and a .sif is a container root,
// so neither can be selected.
func candidateFromPath(target string) (store.Candidate, error) {
	if !strings.HasSuffix(target, ".sqf") {
		return store.Candidate{}, fmt.Errorf("%w: only .sqf can be selected, not %s", ErrInvalid, filepath.Base(target))
	}
	absolute, err := filepath.Abs(target)
	if err != nil {
		return store.Candidate{}, err
	}
	// Read the file, never scan its directory. A scan applies the flat-name
	// rule — the filename must encode the artifact name — which is right where
	// the filename *is* the address and wrong here: a path selection exists so a
	// project can call the file whatever it likes, and overlays/combined.sqf
	// would be refused for not being testdata--combined--1.0.sqf. This is the
	// same rule project.LookupAt applies when it later resolves the selection,
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

// pinnableRequest reports why a selection key cannot be pinned, or "".
//
// It says what to do instead, because the alternative is not obvious: an
// unpinnable dependency is declared unpinnable rather than left unselected.
func pinnableRequest(request string) string {
	destination, isPath := strings.CutPrefix(request, PathPrefix)
	if !isPath {
		return ""
	}
	if utils.IsImg(destination) {
		return fmt.Sprintf("%s is writable and has no identity to pin; declare it `## %s` instead",
			destination, UnpinnedMarker)
	}
	if err := validProjectPath(destination); err != nil {
		return fmt.Sprintf("%v; restore cannot own this path, so declare it `## %s` instead", err, UnpinnedMarker)
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

// vendorClosure writes the selected artifact and every entry of its embedded
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
	// A selection is only as good as its provenance: an artifact built over an
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

// Apply records a resolved selection and publishes the lock, verifying the whole
// closure first. Publication is last, so a lock is never left pointing at
// something that does not validate.
func Apply(root string, l *Lock, selected *Selected) error {
	previous, had := l.Selections[selected.Request]
	l.Selections[selected.Request] = Selection{Artifact: selected.Artifact}

	if _, problems := Verify(root, l); len(problems) > 0 {
		if had {
			l.Selections[selected.Request] = previous
		} else {
			delete(l.Selections, selected.Request)
		}
		return fmt.Errorf("%w: selecting %s leaves the project invalid:\n  %s",
			ErrInvalid, selected.Request, joinProblems(problems))
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

// Reconcile rescans a project and drops selections nothing requests any more,
// reporting what remains unselected. It never invents a selection: choosing an
// artifact is an explicit act.
func Reconcile(root string, l *Lock, result *ScanResult) (unselected []Request) {
	requested := make(map[string]bool, len(result.Requests))
	for _, request := range result.Requests {
		requested[request.Key] = true
	}
	for _, key := range l.Requests() {
		if !requested[key] {
			delete(l.Selections, key)
		}
	}

	verified, _ := Verify(root, l)
	for _, request := range result.Requests {
		// An unpinnable request is not unselected: there is nowhere for restore
		// to put an answer, and G13's marker is how that is declared.
		if !request.Kind.Pinnable() {
			continue
		}
		selection, ok := l.Selections[request.Key]
		if !ok {
			unselected = append(unselected, request)
			continue
		}
		if _, valid := verified.Entries[selection.Artifact]; !valid {
			delete(l.Selections, request.Key)
			unselected = append(unselected, request)
		}
	}
	return unselected
}
