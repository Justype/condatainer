package config

import (
	"context"
	"os"
	"path/filepath"
	"strings"
	"sync"

	"github.com/Justype/condatainer/catalog"
	"github.com/Justype/condatainer/internal/logging"
)

// defaultSource is the collection every resolution ends at.
//
// It ships as a default *value*, not as a fallback the resolver reaches for:
// naming it here means a site replaces it by writing its own `cnt` entry, rather
// than patching a constant (build-recipes.md §2.5). Nothing else in the code
// knows this URL, and no code path prefers it over a configured entry.
var defaultSource = catalog.Spec{
	Name: "cnt",
	Base: "https://raw.githubusercontent.com/condatainer/recipes/main",
}

// layerSources reads the `sources` key from every config layer and concatenates
// them strongest first, so a user entry shadows a site entry of the same name.
//
// Each entry is a single-key mapping — `- lab: /shared/lab/recipes` — because
// both the order and the handle are load-bearing and a plain map gives neither.
// CNT_SOURCES overrides the lot: "cnt=https://…|lab=/shared/lab".
func layerSources() []catalog.Spec {
	if ev := os.Getenv("CNT_SOURCES"); ev != "" {
		return withDefaultSource(parseSourceSpecs(strings.Split(ev, "|")))
	}
	var out []catalog.Spec
	seen := map[string]bool{}
	for _, v := range configLayers {
		if !v.InConfig("sources") {
			continue
		}
		for _, spec := range decodeSourceList(v.Get("sources")) {
			if !seen[spec.Name] {
				seen[spec.Name] = true
				out = append(out, spec)
			}
		}
	}
	return withDefaultSource(out)
}

// withDefaultSource appends the default collection when nothing already answers
// to its handle, so a fresh install resolves recipes unconfigured and a site
// adding its own collection does not silently lose the public one.
//
// Appended, never prepended: every configured entry outranks it. Redefining
// `cnt` replaces it outright, which is how a site points the handle somewhere
// else without rewriting the `#DEP:` lines that name it.
func withDefaultSource(specs []catalog.Spec) []catalog.Spec {
	for _, s := range specs {
		if s.Name == defaultSource.Name {
			return specs
		}
	}
	return append(specs, defaultSource)
}

// decodeSourceList turns viper's view of the YAML sequence into specs. Entries
// may also be written "name=base", which is what the env form uses.
func decodeSourceList(raw any) []catalog.Spec {
	items, ok := raw.([]any)
	if !ok {
		return nil
	}
	var out []catalog.Spec
	for _, item := range items {
		switch v := item.(type) {
		case string:
			out = append(out, parseSourceSpecs([]string{v})...)
		case map[string]any:
			for name, base := range v {
				if s, ok := base.(string); ok {
					out = append(out, catalog.Spec{Name: name, Base: strings.TrimRight(s, "/")})
				}
			}
		case map[any]any:
			for name, base := range v {
				n, okN := name.(string)
				s, okS := base.(string)
				if okN && okS {
					out = append(out, catalog.Spec{Name: n, Base: strings.TrimRight(s, "/")})
				}
			}
		}
	}
	return out
}

// parseSourceSpecs parses "name=base" pairs.
func parseSourceSpecs(pairs []string) []catalog.Spec {
	var out []catalog.Spec
	for _, pair := range pairs {
		name, base, ok := strings.Cut(strings.TrimSpace(pair), "=")
		if !ok || name == "" || base == "" {
			continue
		}
		out = append(out, catalog.Spec{Name: name, Base: strings.TrimRight(base, "/")})
	}
	return out
}

var (
	catalogOnce sync.Once
	catalogVal  catalog.Catalog
	catalogErr  error
)

// OpenCatalog opens the configured sources, once per process.
//
// Config owns which sources exist and where the cache lives; the catalog owns
// everything under that directory.
func OpenCatalog(ctx context.Context) (catalog.Catalog, error) {
	catalogOnce.Do(func() {
		// Defensive: layerSources always leaves defaultSource in. An empty
		// catalog provides nothing, which callers already handle.
		if len(Global.Sources) == 0 {
			return
		}
		cache := catalog.Cache{Dir: CatalogCacheDir(), TTL: Global.MetadataCacheTTL}
		catalogVal, catalogErr = catalog.Open(ctx, Global.Sources, cache)
	})
	return catalogVal, catalogErr
}

var warnSourcesOnce sync.Once

// WarnUnreachableSources reports sources that could not be read, once per
// process however many names get resolved.
//
// Not fatal: the remaining sources still answer, and a compute node with no
// route out is ordinary. But silence is not an option either — sources are
// first-wins, so an unreachable one promotes the next source's recipe, or falls
// through to conda, and the build would otherwise look normal.
//
// Call after the catalog has been consulted: Err is set when a source is first
// read, not when it is opened.
func WarnUnreachableSources(ctx context.Context, cat catalog.Catalog) {
	warnSourcesOnce.Do(func() {
		log := logging.FromContext(ctx)
		for _, s := range cat {
			switch {
			case s.Err != nil:
				log.Warn("source unreachable, skipping it", "source", s.Name, "err", s.Err)
			case s.Stale:
				log.Warn("source not refreshed, using the cached index", "source", s.Name)
			}
		}
	})
}

// BaseRecipeName returns the module name of the base recipe, e.g. "ubuntu24/base"
// for recipes/ubuntu24/base.def. Empty when no base is configured.
func BaseRecipeName() string {
	if base := ResolvedBase(); base != "" {
		return base + "/base"
	}
	return ""
}

// BaseRecipeNameFrom is BaseRecipeName with the default_base fallback, for
// callers that already hold an open catalog. Config `base` still wins.
func BaseRecipeNameFrom(cat catalog.Catalog) string {
	if name := BaseRecipeName(); name != "" {
		return name
	}
	if def := cat.DefaultBase(); def != "" {
		return def + "/base"
	}
	return ""
}

// EnsureBase records the base in config the first time one is needed, taking it
// from the first source declaring a default_base.
//
// Once written it is never revised. Changing the base rebuilds the container
// root and every os overlay stacked on it, so following an upstream bump would
// invalidate a whole set of images on an ordinary update. A later default is
// something the user opts into with `config set base`.
//
// Returns the resolved base, "" when nothing supplies one. A failed write is not
// an error: the base still resolves for this run.
func EnsureBase(cat catalog.Catalog) string {
	if Global.Base != "" {
		return Global.Base
	}
	def := cat.DefaultBase()
	if def == "" {
		return ""
	}
	Global.Base = def
	if path, _, err := ResolveWritableConfigPath(""); err == nil {
		_ = SetConfigKey(path, "base", def)
	}
	return def
}

// SourceDefaultBase returns the base the first source recommends, which may
// differ from the recorded one after an upstream change.
func SourceDefaultBase(cat catalog.Catalog) string { return cat.DefaultBase() }

// ResolvedBase returns the configured base, e.g. "ubuntu24".
//
// Config only — it never opens the catalog. This is the bare-name prefix for
// installed overlays (`build-essential` -> `ubuntu24/build-essential`), so it is
// called on offline paths like list and info, where reaching for a source
// descriptor would mean a network fetch to expand a local name.
//
// The default_base fallback belongs where a catalog is already open: base image
// resolution (internal/build).
func ResolvedBase() string { return Global.Base }

// CatalogCacheDir is where fetched index and recipe bytes are kept.
func CatalogCacheDir() string {
	dir, err := GetWritableCacheDir()
	if err != nil {
		return ""
	}
	return filepath.Join(dir, "catalog")
}

// RefreshCatalogCache discards everything cached for every source, so the next
// read refetches. Removing the whole directory also drops entries for sources no
// longer configured, which is what pruning used to do separately.
func RefreshCatalogCache() error {
	dir := CatalogCacheDir()
	if dir == "" {
		return nil
	}
	if err := os.RemoveAll(dir); err != nil {
		return err
	}
	ResetCatalog()
	return nil
}

// ResetCatalog drops the memoized catalog. Tests reconfigure sources between
// cases; nothing in a command run needs it.
func ResetCatalog() {
	catalogOnce = sync.Once{}
	warnSourcesOnce = sync.Once{}
	catalogVal, catalogErr = nil, nil
}
