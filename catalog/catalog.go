// Package catalog reads recipe collections: what a source offers, what a recipe
// declares, and in what order a dependency graph has to be built.
//
// Every decision that follows is the caller's, so nothing here reads process
// state, prints, or knows what an overlay is.
package catalog

import "strings"

// Type is what a payload is, and what a runtime does with it. It is what a
// recipe declares as #TYPE: and what an image records as its manifest `type`.
//
// Go reserves `type`, so a local or parameter holding one is spelled `typ`.
type Type string

const (
	// TypeBase is retired: no recipe derives it any more (DeriveType always
	// returns TypeOS for a .def, regardless of name), and nothing selects a
	// container root by type-equals-base. It stays defined only so that an
	// artifact manifest written before this retirement, which may still
	// literally record "type": "base", remains a syntactically known type
	// wherever validation reads one back — it is never produced, and never
	// specially selected, again.
	TypeBase Type = "base"
	TypeOS   Type = "os"   // adds to a root that already exists; any os can be the container root, chosen per invocation
	TypeApp  Type = "app"  // contributes to PATH
	TypeData Type = "data" // does not
	// TypeEnv is what `overlay freeze` captures from a writable overlay. No
	// recipe may declare it: DeriveType never returns it and Recipe.Validate
	// rejects it. It stacks at the container root like an os, but is its own
	// type because type carries publishing defaults and an undeclared os may be
	// published publicly. Written "environment" in anything a user reads.
	TypeEnv Type = "env"
)

// DeriveType reports the type of a recipe from its path and headers. A .def is
// always os; everything else takes #TYPE: when declared, else data at two or
// more name components. A template counts the components of target, its
// module path. Only "app" and "data" are accepted.
func DeriveType(name, target string, isDef bool, declared string) Type {
	if isDef {
		return TypeOS
	}
	switch Type(declared) {
	case TypeApp, TypeData:
		return Type(declared)
	}
	path := target
	if path == "" {
		path = name
	}
	if strings.Count(path, "/") >= 2 {
		return TypeData
	}
	return TypeApp
}
