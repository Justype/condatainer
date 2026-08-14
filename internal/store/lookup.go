package store

import (
	"errors"
	"fmt"
	"os"
	"path/filepath"
	"sort"
	"strings"

	"github.com/Justype/condatainer/catalog"
	"github.com/Justype/condatainer/internal/artifact/compare"
	"github.com/Justype/condatainer/internal/artifact/meta"
	"github.com/Justype/condatainer/internal/artifactcache"
	"github.com/Justype/condatainer/internal/config"
	"github.com/Justype/condatainer/internal/image"
)

var (
	ErrNotFound  = errors.New("artifact identity not found")
	ErrAmbiguous = errors.New("artifact identity is ambiguous")
)

// Layout identifies the namespace containing a candidate.
type Layout string

const (
	LayoutFlat   Layout = "flat"
	LayoutStored Layout = "store"
)

// Candidate is a verified local artifact. Identity and Equiv were regenerated
// from its embedded sources and are not inferred from its filename.
type Candidate struct {
	Name     string      `json:"name"`
	Path     string      `json:"path"`
	Root     string      `json:"root"`
	Layout   Layout      `json:"layout"`
	Size     int64       `json:"size"`
	Identity meta.KeyRef `json:"identity"`
	Equiv    meta.KeyRef `json:"equiv"`
}

// Issue is an inventory entry that could not be trusted.
type Issue struct {
	Path  string `json:"path"`
	Error string `json:"error"`
}

// ScanOptions selects namespaces and an optional normalized artifact name.
// With neither namespace selected, both are scanned.
type ScanOptions struct {
	Dirs   []string
	Name   string
	Flat   bool
	Stored bool
	// Uncached forces key regeneration for explicit validation commands.
	Uncached bool
}

// Report separates verified candidates from untrusted entries.
type Report struct {
	Candidates []Candidate `json:"candidates"`
	Issues     []Issue     `json:"issues,omitempty"`
}

type artifactReader func(string) (compare.Artifact, error)

// Scan regenerates keys for every selected candidate. Directory errors and bad
// artifacts are returned as issues so one corrupt lower-priority entry cannot
// hide a valid one.
func Scan(opts ScanOptions) Report {
	if opts.Uncached {
		return scan(opts, compare.Read, nil)
	}
	batch := artifactcache.Default().NewBatch()
	report := scan(opts, compare.Read, batch)
	batch.Flush()
	return report
}

func scan(opts ScanOptions, read artifactReader, cache artifactcache.Access) Report {
	dirs := opts.Dirs
	if dirs == nil {
		dirs = config.GetImageSearchPaths()
	}
	if !opts.Flat && !opts.Stored {
		opts.Flat, opts.Stored = true, true
	}
	wanted := catalog.Normalize(opts.Name)
	var report Report
	for _, root := range dirs {
		if opts.Flat {
			report.scanFlat(root, wanted, read, cache)
		}
		if opts.Stored {
			report.scanStored(root, wanted, read, cache)
		}
	}
	return report
}

func (r *Report) scanFlat(root, wanted string, read artifactReader, cache artifactcache.Access) {
	entries, err := os.ReadDir(root)
	if err != nil {
		if !os.IsNotExist(err) {
			r.addIssue(root, err)
		}
		return
	}
	for _, entry := range entries {
		if entry.IsDir() || entry.Type()&os.ModeSymlink != 0 || filepath.Ext(entry.Name()) != ".sqf" {
			continue
		}
		name := image.DecodeArtifactName(strings.TrimSuffix(entry.Name(), ".sqf"))
		if wanted != "" && name != wanted {
			continue
		}
		r.verify(filepath.Join(root, entry.Name()), root, LayoutFlat, name, "", read, cache)
	}
}

func (r *Report) scanStored(root, wanted string, read artifactReader, cache artifactcache.Access) {
	dir := filepath.Join(root, DirName)
	entries, err := os.ReadDir(dir)
	if err != nil {
		if !os.IsNotExist(err) {
			r.addIssue(dir, err)
		}
		return
	}
	for _, entry := range entries {
		if entry.IsDir() {
			continue
		}
		if entry.Type()&os.ModeSymlink != 0 {
			if filepath.Ext(entry.Name()) == ".sqf" {
				r.addIssue(filepath.Join(dir, entry.Name()), errors.New("store entries may not be symbolic links"))
			}
			continue
		}
		parsed, err := ParseFilename(entry.Name())
		if err != nil {
			if filepath.Ext(entry.Name()) == ".sqf" {
				r.addIssue(filepath.Join(dir, entry.Name()), err)
			}
			continue
		}
		if wanted != "" && parsed.Name != wanted {
			continue
		}
		r.verify(filepath.Join(dir, entry.Name()), root, LayoutStored, parsed.Name, parsed.Prefix, read, cache)
	}
}

func (r *Report) verify(path, root string, layout Layout, addressedName, prefix string, read artifactReader, cache artifactcache.Access) {
	info, err := os.Lstat(path)
	if err != nil {
		r.addIssue(path, err)
		return
	}
	if !info.Mode().IsRegular() {
		r.addIssue(path, errors.New("not a regular file"))
		return
	}
	if cache != nil {
		if record, ok := cache.Lookup(path, info); ok && record.KeysVerified && record.Name != "" {
			identity := meta.KeyRef{Scheme: record.Identity.Scheme, SHA256: record.Identity.SHA256}
			equiv := meta.KeyRef{Scheme: record.Equiv.Scheme, SHA256: record.Equiv.SHA256}
			if cachedAddressValid(record.Name, identity, equiv, addressedName, prefix) {
				r.Candidates = append(r.Candidates, Candidate{
					Name: record.Name, Path: path, Root: root, Layout: layout,
					Size: info.Size(), Identity: identity, Equiv: equiv,
				})
				return
			}
		}
	}
	artifact, err := read(path)
	if err != nil {
		r.addIssue(path, err)
		return
	}
	identity := keyRef(artifact.IdentityScheme, artifact.Identity)
	equiv := keyRef(artifact.EquivScheme, artifact.Equiv)
	if _, err := ParseKeyRef(FormatKeyRef(identity)); err != nil {
		r.addIssue(path, fmt.Errorf("artifact carries no valid complete identity: %w", err))
		return
	}
	if _, err := ParseKeyRef(FormatKeyRef(equiv)); err != nil {
		r.addIssue(path, fmt.Errorf("artifact carries no valid complete equivalence key: %w", err))
		return
	}
	if artifact.Name != addressedName {
		r.addIssue(path, fmt.Errorf("filename names %q but artifact names %q", addressedName, artifact.Name))
		return
	}
	if prefix != "" && !strings.HasPrefix(identity.SHA256, prefix) {
		r.addIssue(path, fmt.Errorf("filename prefix %s does not match identity %s", prefix, identity.Digest()))
		return
	}
	r.Candidates = append(r.Candidates, Candidate{
		Name: artifact.Name, Path: path, Root: root, Layout: layout,
		Size: info.Size(), Identity: identity, Equiv: equiv,
	})
	if cache != nil {
		cache.Merge(path, info, func(record *artifactcache.Record) {
			record.Name = artifact.Name
			record.Identity = artifactcache.Key{Scheme: identity.Scheme, SHA256: identity.SHA256}
			record.Equiv = artifactcache.Key{Scheme: equiv.Scheme, SHA256: equiv.SHA256}
			record.KeysVerified = true
			record.Layout = string(layout)
			record.Root = root
			record.Size = info.Size()
		})
	}
}

func cachedAddressValid(name string, identity, equiv meta.KeyRef, addressedName, prefix string) bool {
	if name != addressedName {
		return false
	}
	if _, err := ParseKeyRef(FormatKeyRef(identity)); err != nil {
		return false
	}
	if _, err := ParseKeyRef(FormatKeyRef(equiv)); err != nil {
		return false
	}
	return prefix == "" || strings.HasPrefix(identity.SHA256, prefix)
}

func (r *Report) addIssue(path string, err error) {
	r.Issues = append(r.Issues, Issue{Path: path, Error: err.Error()})
}

func keyRef(scheme, digest string) meta.KeyRef {
	return meta.KeyRef{Scheme: scheme, SHA256: strings.TrimPrefix(digest, "sha256:")}
}

// ResolveIdentity returns the nearest verified candidate matching q. Multiple
// paths carrying the same complete identity are copies, not ambiguity.
func ResolveIdentity(name string, q IdentityQuery, dirs []string) (Candidate, Report, error) {
	report := Scan(ScanOptions{Dirs: dirs, Name: name})
	candidate, err := resolveReport(name, q, report)
	return candidate, report, err
}

func resolveReport(name string, q IdentityQuery, report Report) (Candidate, error) {
	var matches []Candidate
	identities := make(map[meta.KeyRef]struct{})
	for _, candidate := range report.Candidates {
		if q.Matches(candidate.Identity) {
			matches = append(matches, candidate)
			identities[candidate.Identity] = struct{}{}
		}
	}
	if len(matches) == 0 {
		return Candidate{}, fmt.Errorf("%w: %s %s", ErrNotFound, catalog.Normalize(name), queryString(q))
	}
	if len(identities) > 1 {
		return Candidate{}, fmt.Errorf("%w: %s %s matches %d complete identities", ErrAmbiguous, catalog.Normalize(name), queryString(q), len(identities))
	}
	return matches[0], nil
}

// Equivalent returns verified candidates with the requested equivalence key.
func Equivalent(name string, q IdentityQuery, dirs []string) (Report, error) {
	report := Scan(ScanOptions{Dirs: dirs, Name: name})
	filtered := report.Candidates[:0]
	for _, candidate := range report.Candidates {
		if q.Matches(candidate.Equiv) {
			filtered = append(filtered, candidate)
		}
	}
	report.Candidates = filtered
	if len(filtered) == 0 {
		return report, fmt.Errorf("%w: no equivalent candidate for %s", ErrNotFound, catalog.Normalize(name))
	}
	return report, nil
}

// SortReport gives CLI and tests stable output independent of ReadDir details.
func SortReport(report *Report) {
	sort.SliceStable(report.Candidates, func(i, j int) bool { return report.Candidates[i].Path < report.Candidates[j].Path })
	sort.SliceStable(report.Issues, func(i, j int) bool { return report.Issues[i].Path < report.Issues[j].Path })
}

func queryString(q IdentityQuery) string {
	if q.Scheme != "" {
		return q.Scheme + "@sha256:" + q.SHA256
	}
	return q.SHA256
}
