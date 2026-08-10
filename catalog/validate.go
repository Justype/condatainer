package catalog

import (
	"errors"
	"fmt"
	"slices"
	"strings"
)

// ErrInvalidRecipe reports a recipe that declares something its type may not.
var ErrInvalidRecipe = errors.New("catalog: invalid recipe")

// Validate reports headers a recipe of this type may not declare. It is called
// once the recipe is fetched and about to be built, not while indexing: one bad
// recipe must not take a whole collection out of a listing.
//
// Two rules, each narrow:
//
//   - Only data may declare #DEP:. An app is prebuilt and self-contained, an OS
//     is self-contained by definition, and a base *is* the build environment. A
//     recipe that genuinely needs a compiler is an os artifact providing that
//     toolchain, not an app depending on one.
//   - Only app and data may declare #ARCH:, and only as native or noarch. An OS
//     and a base are root filesystems and are always architecture-specific.
func (r *Recipe) Validate() error {
	var errs []error

	if len(r.Deps) > 0 && r.Type != TypeData {
		errs = append(errs, fmt.Errorf("%w: %s is type %s and may not declare #DEP: (%s); "+
			"only data has build dependencies", ErrInvalidRecipe, r.Name, r.Type, strings.Join(r.Deps, ", ")))
	}

	if r.Arch != "" {
		switch r.Type {
		case TypeApp, TypeData:
			if r.Arch != ArchNative && r.Arch != ArchNoarch {
				errs = append(errs, fmt.Errorf("%w: %s declares #ARCH:%s; the values are %s and %s",
					ErrInvalidRecipe, r.Name, r.Arch, ArchNative, ArchNoarch))
			}
		default:
			errs = append(errs, fmt.Errorf("%w: %s is type %s and may not declare #ARCH:; "+
				"a root filesystem is always architecture-specific", ErrInvalidRecipe, r.Name, r.Type))
		}
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
			out = append(out, fmt.Sprintf(
				"%s: dependency %s is mentioned in the name but not as whole components, "+
					"so it does not count toward equivalence — fix the #TARGET: or the name",
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
