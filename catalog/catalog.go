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
	TypeBase Type = "base" // produces the container root
	TypeOS   Type = "os"   // adds to a root that already exists
	TypeApp  Type = "app"  // contributes to PATH
	TypeData Type = "data" // does not
)

// DeriveType reports the type of a recipe from its path and headers. A .def is
// base when its name ends in /base, os otherwise; everything else takes #TYPE:
// when declared, else data at two or more name components. A template counts the
// components of target, its module path. Only "app" and "data" are accepted.
func DeriveType(name, target string, isDef bool, declared string) Type {
	if isDef {
		if strings.HasSuffix(name, "/base") {
			return TypeBase
		}
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
