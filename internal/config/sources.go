package config

import (
	"context"
	"fmt"
	"os"
	"path/filepath"
	"strings"
	"sync"

	"github.com/Justype/condatainer/catalog"
	"github.com/Justype/condatainer/internal/logging"
)

// defaultSource is the collection every resolution ends at. It is a default
// value rather than a fallback the resolver reaches for, so a site replaces it
// by writing its own `cnt` entry. See the README's Recipe sources.
var defaultSource = catalog.Spec{
	Name: "cnt",
	Base: "https://raw.githubusercontent.com/condatainer/cnt/main",
}

// layerSources reads the `sources` key from every config layer and concatenates
// them strongest first, so a user entry shadows a site entry of the same name.
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
// to its handle, so a fresh install resolves recipes unconfigured. Appended,
// never prepended: every configured entry outranks it.
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

// SelectSources restricts recipe resolution to the named configured source
// handles, in the order given. An empty selection leaves the configured source
// list unchanged. Duplicate handles are ignored after their first occurrence.
func SelectSources(names []string) error {
	if len(names) == 0 {
		return nil
	}

	configured := make(map[string]catalog.Spec, len(Global.Sources))
	available := make([]string, 0, len(Global.Sources))
	for _, source := range Global.Sources {
		configured[source.Name] = source
		available = append(available, source.Name)
	}

	selected := make([]catalog.Spec, 0, len(names))
	seen := make(map[string]bool, len(names))
	for _, raw := range names {
		name := strings.TrimSpace(raw)
		source, ok := configured[name]
		if !ok || name == "" {
			return fmt.Errorf("unknown source %q (configured: %s)",
				name, strings.Join(available, ", "))
		}
		if seen[name] {
			continue
		}
		seen[name] = true
		selected = append(selected, source)
	}

	Global.Sources = selected
	ResetCatalog()
	return nil
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
// process however many names get resolved. Call it after the catalog has been
// consulted. See the README's Recipe sources.
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
// for recipes/ubuntu24/base.def. Empty when no default distro is configured.
func BaseRecipeName() string {
	if distro := ResolvedDefaultDistro(); distro != "" {
		return distro + "/base"
	}
	return ""
}

// BaseRecipeNameFrom is BaseRecipeName with the default_distro fallback, for
// callers that already hold an open catalog. Config `default_distro` still wins.
func BaseRecipeNameFrom(cat catalog.Catalog) string {
	if name := BaseRecipeName(); name != "" {
		return name
	}
	if def := cat.DefaultDistro(); def != "" {
		return def + "/base"
	}
	return ""
}

// EnsureDefaultDistro records the default distro in config the first time one
// is needed, taking it from the first source declaring a default_distro, and
// never revises it. Returns the resolved distro, "" when nothing supplies
// one; a failed write is not an error.
func EnsureDefaultDistro(cat catalog.Catalog) string {
	if Global.DefaultDistro != "" {
		return Global.DefaultDistro
	}
	def := cat.DefaultDistro()
	if def == "" {
		return ""
	}
	Global.DefaultDistro = def
	if path, _, err := ResolveWritableConfigPath(""); err == nil {
		_ = SetConfigKey(path, "default_distro", def)
	}
	return def
}

// SourceDefaultDistro returns the default distro the first source recommends,
// which may differ from the recorded one after an upstream change.
func SourceDefaultDistro(cat catalog.Catalog) string { return cat.DefaultDistro() }

// ResolvedDefaultDistro returns the configured default distro, e.g. "ubuntu24"
// — the bare-name prefix for installed overlays. Config only: it never opens
// the catalog, so offline paths like list and info stay offline.
func ResolvedDefaultDistro() string { return Global.DefaultDistro }

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
