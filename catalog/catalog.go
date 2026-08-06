// Package catalog reads recipe collections: what a source offers, what a recipe
// declares, and in what order a dependency graph has to be built.
//
// Every decision that follows is the caller's, so nothing here reads process
// state, prints, or knows what an overlay is.
package catalog

import "strings"

// Kind is what a payload is, and what a runtime does with it.
type Kind string

const (
	KindBase Kind = "base" // produces the container root
	KindOS   Kind = "os"   // adds to a root that already exists
	KindApp  Kind = "app"  // contributes to PATH
	KindData Kind = "data" // does not
)

// DeriveKind reports the kind of a recipe from its path and headers.
//
// A .def is base when its name ends in /base, os otherwise. Everything else is
// app or data: from #TYPE: when declared, from the slash count otherwise, where
// two or more components means data. A template counts slashes in target, its
// module path — grch38/star-gencode has one, its target three.
func DeriveKind(name, target string, isDef bool, declared string) Kind {
	if isDef {
		if strings.HasSuffix(name, "/base") {
			return KindBase
		}
		return KindOS
	}
	switch Kind(declared) {
	case KindApp, KindData:
		return Kind(declared)
	}
	path := target
	if path == "" {
		path = name
	}
	if strings.Count(path, "/") >= 2 {
		return KindData
	}
	return KindApp
}
