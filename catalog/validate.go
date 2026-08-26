package catalog

import (
	"errors"
	"fmt"
	"path/filepath"
	"slices"
	"strings"
)

// ErrInvalidRecipe reports a recipe that declares something its type may not.
var ErrInvalidRecipe = errors.New("catalog: invalid recipe")

// IsPathDep reports whether a #DEP: value names an overlay file rather than a
// catalog name/version.
//
// Exported because it decides which grammar a dependency is read with, and this
// package owns that grammar: Normalize and ParseDep apply to a name/version and
// to nothing else. A second answer elsewhere would let one string be a name in
// one place and a path in another. The extensions are utils.IsOverlay's set,
// which is what a running script's declaration is already read with.
func IsPathDep(value string) bool {
	switch strings.ToLower(filepath.Ext(strings.TrimSpace(value))) {
	case ".sqf", ".sqsh", ".squashfs", ".img", ".ext3":
		return true
	}
	return false
}

// ValidateDeps reports why a **build** of this type may not declare these
// dependencies, or nil.
//
// Only builds answer to these rules. A running script may declare an overlay
// path or an external `.sqf`, because it mounts what it names and records
// nothing; a build's declaration becomes an edge in an artifact that has to mean
// the same thing on another machine. Accordingly the only callers are
// Recipe.Validate and build.FromExternalSource — never `run` or `check`.
//
// Separate from Validate because an external build (`create -p -f <script>`)
// reaches the same rules without ever being a catalog Recipe — its type comes
// from #TYPE: and its name from #TARGET:, so it has no recipe path to be parsed
// from. One definition, both callers.
//
//   - Only data may declare #DEP:. An app is prebuilt and self-contained, an OS
//     is self-contained by definition, and a base *is* the build environment. A
//     recipe that genuinely needs a compiler is an os artifact providing that
//     toolchain, not an app depending on one.
//   - A build dependency is a name/version, never an overlay path. A built
//     artifact records each edge as a name plus a complete identity, and a path
//     supplies neither: nothing could re-resolve it on another machine, and no
//     key could be regenerated from it. Keeping a project's vendored source
//     closure total depends on it — see internal/project/README.md.
func ValidateDeps(name string, typ Type, deps []string) error {
	if len(deps) == 0 {
		return nil
	}
	var errs []error

	if typ != TypeData {
		errs = append(errs, fmt.Errorf("%w: %s is type %s and may not declare #DEP: (%s); only data has build dependencies",
			ErrInvalidRecipe, name, typ, strings.Join(deps, ", ")))
	}

	var paths []string
	for _, dep := range deps {
		if IsPathDep(dep) {
			paths = append(paths, dep)
		}
	}
	if len(paths) > 0 {
		errs = append(errs, fmt.Errorf("%w: %s declares #DEP: %s; a build dependency must be a name/version, not an overlay path",
			ErrInvalidRecipe, name, strings.Join(paths, ", ")))
	}

	return errors.Join(errs...)
}

// Validate reports headers a recipe of this type may not declare. It is called
// once the recipe is fetched and about to be built, not while indexing: one bad
// recipe must not take a whole collection out of a listing.
//
// The #DEP: rules live in ValidateDeps, which an external build shares. The one
// rule that is only a recipe's is #ARCH:, which only app and data may declare,
// and only as native or noarch: an OS and a base are root filesystems and are
// always architecture-specific.
//
// #REDISTRIBUTE: is rejected rather than ignored when it is neither yes nor no.
// It decides whether a payload may be published, so a typo that silently read as
// "unanswered" would fall back to the type default and publish the artifact the
// author was trying to hold back.
func (r *Recipe) Validate() error {
	var errs []error

	if err := ValidateDeps(r.Name, r.Type, r.Deps); err != nil {
		errs = append(errs, err)
	}

	if r.Arch != "" {
		switch r.Type {
		case TypeApp, TypeData:
			if r.Arch != ArchNative && r.Arch != ArchNoarch {
				errs = append(errs, fmt.Errorf("%w: %s declares #ARCH:%s; the values are %s and %s",
					ErrInvalidRecipe, r.Name, r.Arch, ArchNative, ArchNoarch))
			}
		default:
			errs = append(errs, fmt.Errorf("%w: %s is type %s and may not declare #ARCH:; a root filesystem is always architecture-specific",
				ErrInvalidRecipe, r.Name, r.Type))
		}
	}

	if r.Redistribute != "" && r.Redistribute != "yes" && r.Redistribute != "no" {
		errs = append(errs, fmt.Errorf("%w: %s declares #REDISTRIBUTE:%s; the values are yes and no",
			ErrInvalidRecipe, r.Name, r.Redistribute))
	}

	return errors.Join(errs...)
}

// Lint reports declarations that are legal but probably not what the author
// meant. A lint blocks nothing and changes no key; it is returned rather than
// printed, because this package never writes to a terminal.
//
// The one rule so far is the near-miss on the naming convention that decides
// which tools define a dataset. An app or OS dependency shapes what a dataset may
// substitute for only when its name/version appears as a clean run of components
// in the dataset's own name — so a #TARGET: rendering .../star{star_version}/...
// produces the single component star2.7.11b, and a version the author clearly
// meant to be load-bearing quietly stops counting. The fix is to the #TARGET:,
// which is why this names the dependency rather than guessing at a repair.
func (r *Recipe) Lint() []string {
	var out []string
	for _, raw := range r.Deps {
		dep, err := ParseDep(raw)
		if err != nil || dep.Version == "" {
			continue
		}
		nameVersion := dep.NameVersion()
		if HasComponents(r.Name, nameVersion) {
			continue
		}
		// The tool's own component, not the whole dep path: an OS dependency
		// carries its distro, and a dataset that mentions ".../pytorch/2.9/..."
		// without it is exactly the near-miss worth reporting.
		tool := dep.Name
		if i := strings.LastIndex(tool, "/"); i >= 0 {
			tool = tool[i+1:]
		}
		if strings.Contains(r.Name, tool) && strings.Contains(r.Name, dep.Version) {
			out = append(out, fmt.Sprintf("%s: dependency %s is in the name but not as whole components, so it does not count toward equivalence; fix the #TARGET: or the name",
				r.Name, nameVersion))
		}
	}
	return out
}

// HasComponents reports whether dep's slash-separated components occur as a
// contiguous run of name's.
//
// Components are compared as exact strings and versions are never compared
// semantically, since 2.7.11b and 2.7.11 name different builds. Substring
// matching is never used: star/2.7.11b must not be found inside a component such
// as star2.7.11b-old.
func HasComponents(name, dep string) bool {
	n := strings.Split(Normalize(name), "/")
	d := strings.Split(Normalize(dep), "/")
	if len(d) == 0 || len(d) > len(n) {
		return false
	}
	for i := 0; i+len(d) <= len(n); i++ {
		if slices.Equal(n[i:i+len(d)], d) {
			return true
		}
	}
	return false
}
