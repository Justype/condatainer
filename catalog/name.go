package catalog

import (
	"errors"
	"strings"
)

// ErrEmptyDep is returned by ParseDep for a blank dependency.
var ErrEmptyDep = errors.New("catalog: empty dependency")

// Dep is one #DEP: edge. The preferred version is the implicit upper bound, so
// Min "2.7.0" with Version "2.7.11b" admits [2.7.0, 2.7.11b].
type Dep struct {
	Name    string // module name, e.g. "star" or "grch38/genome"
	Version string // preferred version; empty when the dep names none
	Min     string // from >= or >; empty when unconstrained
	Op      string // ">=" or ">"; empty when unconstrained
}

// String renders the dep back to its #DEP: form.
func (d Dep) String() string {
	return d.NameVersion() + d.Op + d.Min
}

// NameVersion is the module path without the constraint: what names an artifact,
// as against what a dep will accept.
func (d Dep) NameVersion() string {
	if d.Version == "" {
		return d.Name
	}
	return d.Name + "/" + d.Version
}

// Normalize folds the spellings of one module into name/version form,
// converting =, @ and -- to / and trimming space. Constraint suffixes survive.
//
// Exported because it decides what a name is: index keys use this form, so a
// second implementation elsewhere would make one string resolve two ways.
func Normalize(nameVersion string) string {
	nv, op, min := splitConstraint(strings.TrimSpace(nameVersion))
	nv = strings.ReplaceAll(nv, "--", "/")
	nv = strings.ReplaceAll(nv, "=", "/")
	nv = strings.ReplaceAll(nv, "@", "/")
	return nv + op + min
}

// ParseDep parses a #DEP: value such as "star/2.7.11b>=2.7.0". The last path
// component is the version, so a template dep keeps its token there.
func ParseDep(raw string) (Dep, error) {
	nv, op, min := splitConstraint(Normalize(raw))
	if nv == "" {
		return Dep{}, ErrEmptyDep
	}
	d := Dep{Name: nv, Op: op, Min: min}
	if i := strings.LastIndex(nv, "/"); i >= 0 {
		d.Name, d.Version = nv[:i], nv[i+1:]
	}
	return d, nil
}

// splitConstraint separates a spec from its >= or > suffix.
func splitConstraint(raw string) (nameVersion, op, min string) {
	for _, sep := range []string{">=", ">"} {
		if before, after, ok := strings.Cut(raw, sep); ok {
			return strings.TrimSpace(before), sep, strings.TrimSpace(after)
		}
	}
	return raw, "", ""
}

// Satisfies reports whether version falls in the range the dep admits: at or
// above Min per Op, never above Version. Unconstrained admits any version
// sharing Version as a dot-component prefix — "17" admits "17.0.18" the same
// family match Conda's own MatchSpec syntax gives a bare "=17".
func (d Dep) Satisfies(version string) bool {
	if d.Op == "" {
		return d.Version == "" || versionHasPrefix(version, d.Version)
	}
	switch cmp := CompareVersions(version, d.Min); d.Op {
	case ">=":
		if cmp < 0 {
			return false
		}
	case ">":
		if cmp <= 0 {
			return false
		}
	default:
		return false
	}
	return d.Version == "" || CompareVersions(version, d.Version) <= 0
}

// versionHasPrefix reports whether every dot-separated component of prefix
// equals the corresponding component of version, in order. A full version is
// trivially its own prefix, so this subsumes equality rather than replacing
// it — only a genuinely shorter spec ("17", "17.0") gains new matches.
func versionHasPrefix(version, prefix string) bool {
	vParts := strings.Split(version, ".")
	pParts := strings.Split(prefix, ".")
	if len(pParts) > len(vParts) {
		return false
	}
	for i, p := range pParts {
		if vParts[i] != p {
			return false
		}
	}
	return true
}
