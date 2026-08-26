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
	Schema int `json:"schema"`
	// Source is the collection's own code repository, recorded into every
	// artifact it builds as manifest build.source and published as
	// org.opencontainers.image.source.
	//
	// Named for what it becomes rather than what it looks like: `repository`
	// means an OCI repository coordinate everywhere else in the tree
	// (lock.Remote, oras' Reference), and one word cannot carry both.
	Source      string `json:"source"`
	DefaultBase string `json:"default_base"`
	OCI         OCI    `json:"oci,omitzero"`
}

// OCI is how artifacts belonging to one source are published and fetched.
// Push is singular because replication is an explicit publishing operation;
// Pull is ordered so a site-local mirror can precede an external registry.
type OCI struct {
	Push     string   `json:"push,omitempty"`
	Pull     []string `json:"pull,omitempty"`
	Audience string   `json:"audience,omitempty"`
}

// Source is one collection of recipes and helpers.
type Source struct {
	Name string // the local handle from config
	Base string // HTTP base URL or filesystem path
	Desc Descriptor
	// DescriptorErr reports that source.json existed but was not usable. Recipes
	// remain available: a broken optional descriptor must not erase a collection,
	// but callers must not trust its registry or provenance defaults.
	DescriptorErr error

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
	desc, err := ParseDescriptor(data)
	if err != nil {
		s.DescriptorErr = err
		return
	}
	s.Desc = desc
}

// ParseDescriptor decodes and validates source.json. Schema zero is accepted
// for backward-compatible local descriptors; a non-zero unsupported schema is
// rejected. OCI endpoints are normalized to registry/repository roots.
func ParseDescriptor(data []byte) (Descriptor, error) {
	var desc Descriptor
	if err := json.Unmarshal(data, &desc); err != nil {
		return Descriptor{}, fmt.Errorf("catalog: decode source descriptor: %w", err)
	}
	if desc.Schema != 0 && desc.Schema != 1 {
		return Descriptor{}, fmt.Errorf("catalog: unsupported source descriptor schema %d", desc.Schema)
	}
	if err := ValidSourceURL(desc.Source); err != nil {
		return Descriptor{}, fmt.Errorf("catalog: source descriptor: %w", err)
	}
	if err := desc.OCI.normalize(); err != nil {
		return Descriptor{}, err
	}
	return desc, nil
}

// ValidSourceURL reports whether raw is usable as
// org.opencontainers.image.source. Empty is valid and means the annotation is
// omitted.
//
// Shape only: whether the repository exists is not a question anything here can
// answer offline. The scheme is the part that matters, because a registry links
// a package to a repository by exact URL match and anything else links nothing.
func ValidSourceURL(raw string) error {
	source := strings.TrimSpace(raw)
	if source == "" {
		return nil
	}
	if strings.ContainsAny(source, " \t\r\n") {
		return fmt.Errorf("source %q contains whitespace", source)
	}
	if !strings.HasPrefix(source, "https://") && !strings.HasPrefix(source, "http://") {
		return fmt.Errorf("source %q is not an http(s) URL", source)
	}
	return nil
}

func (o *OCI) normalize() error {
	o.Audience = strings.ToLower(strings.TrimSpace(o.Audience))
	if o.Audience == "" {
		o.Audience = "public"
	}
	if o.Audience != "public" && o.Audience != "restricted" {
		return fmt.Errorf("catalog: OCI audience must be public or restricted, got %q", o.Audience)
	}

	declared := strings.TrimSpace(o.Push) != "" || len(o.Pull) != 0
	if !declared {
		return nil
	}
	if strings.TrimSpace(o.Push) == "" {
		return errors.New("catalog: OCI endpoints require push")
	}
	if len(o.Pull) == 0 {
		return errors.New("catalog: OCI endpoints require at least one pull endpoint")
	}

	var err error
	if o.Push, err = normalizeOCIEndpoint(o.Push); err != nil {
		return fmt.Errorf("catalog: OCI push endpoint: %w", err)
	}
	for i := range o.Pull {
		o.Pull[i], err = normalizeOCIEndpoint(o.Pull[i])
		if err != nil {
			return fmt.Errorf("catalog: OCI pull endpoint %d: %w", i+1, err)
		}
	}
	return nil
}

func normalizeOCIEndpoint(endpoint string) (string, error) {
	endpoint = strings.TrimSpace(endpoint)
	endpoint = strings.TrimPrefix(endpoint, "oci://")
	endpoint = strings.TrimRight(endpoint, "/")
	if endpoint == "" || strings.Contains(endpoint, "://") ||
		strings.ContainsAny(endpoint, " \t\r\n") || !strings.Contains(endpoint, "/") {
		return "", fmt.Errorf("%q must be a registry/repository root", endpoint)
	}
	return endpoint, nil
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
