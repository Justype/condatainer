package catalog

import (
	"context"
	"encoding/json"
	"errors"
	"fmt"
	"strings"
	"sync"
)

const (
	recipesDir     = "recipes"
	indexDir       = "index"
	descriptorFile = "source.json"
)

// ErrNotProvided reports that no source offers a name. It is an outcome rather
// than a failure: the caller has a fallback this package must not know about.
var ErrNotProvided = errors.New("catalog: no source provides this")

// Spec is one configured source: a local handle and a base. This is how config
// gets in without being imported.
type Spec struct{ Name, Base string }

// Descriptor is a source's source.json.
type Descriptor struct {
	Schema      int    `json:"schema"`
	Repository  string `json:"repository"`
	DefaultBase string `json:"default_base"`
}

// Source is one collection of recipes and helpers.
type Source struct {
	Name string // the local handle from config
	Base string // HTTP base URL or filesystem path
	Desc Descriptor

	// Stale is set when an index was served from cache after a fetch failed.
	// What that means — a note, or a stop — is the caller's.
	Stale bool
	// Err is set when the source could not be reached and had nothing cached.
	// It does not remove the source from the catalog, so list can say so.
	Err error

	b backend

	mu     sync.Mutex
	loaded bool
	cached map[string]*Entry
}

// Entries reads what this source offers, once. A walk or an index fetch is
// repeated for every lookup otherwise, and a collection does not change under a
// running command.
func (s *Source) Entries(ctx context.Context) (map[string]*Entry, error) {
	s.mu.Lock()
	defer s.mu.Unlock()
	if s.loaded {
		return s.cached, s.Err
	}
	s.loaded = true
	s.cached, s.Err = s.b.entries(ctx)
	return s.cached, s.Err
}

// Catalog is the configured sources in order. Earlier entries shadow later
// ones, matching image search, the config chain, and PATH.
type Catalog []*Source

// backend is how a source is read. Open picks one from the shape of Base and
// nothing switches on it afterwards.
type backend interface {
	entries(ctx context.Context) (map[string]*Entry, error)
	read(ctx context.Context, path string) ([]byte, error)
}

// Open prepares the configured sources. A source that cannot be reached keeps
// its place with Err set, since dropping it would silently promote the next
// source's recipes under first-wins ordering.
func Open(ctx context.Context, specs []Spec, cache Cache) (Catalog, error) {
	if len(specs) == 0 {
		return nil, errors.New("catalog: no sources configured")
	}
	cat := make(Catalog, 0, len(specs))
	for _, spec := range specs {
		if spec.Base == "" {
			return nil, fmt.Errorf("catalog: source %q has no base", spec.Name)
		}
		s := &Source{Name: spec.Name, Base: spec.Base}
		if isURL(spec.Base) {
			s.b = &httpBackend{src: s, cache: cache}
		} else {
			s.b = &dirBackend{src: s}
		}
		s.loadDescriptor(ctx)
		cat = append(cat, s)
	}
	return cat, nil
}

// loadDescriptor reads source.json. A source without one still works — a plain
// directory of recipes is the simplest collection there is — so only the
// provenance and default base are lost.
func (s *Source) loadDescriptor(ctx context.Context) {
	data, err := s.b.read(ctx, descriptorFile)
	if err != nil {
		return
	}
	_ = json.Unmarshal(data, &s.Desc)
}

// Entries returns every entry the catalog offers, merged by name with the first
// source winning.
//
// An unreachable source contributes nothing and is skipped, the same way Lookup
// skips it — including when that leaves the result empty, which is a catalog
// that offers nothing rather than a failure. Why each source came back empty is
// on the source, as Err and Stale, so a caller reports it once for all of them
// instead of receiving it joined onto every call.
func (c Catalog) Entries(ctx context.Context) map[string]*Entry {
	out := map[string]*Entry{}
	for _, s := range c {
		found, err := s.Entries(ctx)
		if err != nil {
			continue
		}
		for name, e := range found {
			if _, taken := out[name]; !taken {
				out[name] = e
			}
		}
	}
	return out
}

// DefaultBase is the base recipe named by the first source declaring one.
func (c Catalog) DefaultBase() string {
	for _, s := range c {
		if s.Desc.DefaultBase != "" {
			return s.Desc.DefaultBase
		}
	}
	return ""
}

func isURL(base string) bool {
	return strings.HasPrefix(base, "http://") || strings.HasPrefix(base, "https://")
}
