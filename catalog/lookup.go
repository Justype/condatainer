package catalog

import (
	"context"
	"fmt"
	"maps"
	"slices"
	"strings"
)

// Match is a resolved name: the entry, the source it came from, and — when the
// name was a full template target — the values that produce it.
type Match struct {
	Entry  *Entry
	Source *Source
	Vars   map[string]string // empty unless the name matched a template target
}

// Lookup resolves a name against the catalog, earlier sources winning. A name is
// either an index key or a full target match, and not found is an outcome rather
// than an error. See the README's What may surprise you.
func (c Catalog) Lookup(ctx context.Context, name string) (*Match, bool, error) {
	name = Normalize(name)
	if name == "" {
		return nil, false, nil
	}

	for _, s := range c {
		entries, err := s.Entries(ctx)
		if err != nil {
			continue // recorded on the source; a dead source is not a dead lookup
		}
		if e, ok := entries[name]; ok {
			return &Match{Entry: e, Source: s}, true, nil
		}
		if m := matchTemplates(name, entries, s); m != nil {
			return m, true, nil
		}
	}
	return nil, false, nil
}

// matchTemplates tries every template in one source. Two matching within a
// single source is a collision the index generator can detect; across sources
// the earlier one wins by having been tried first.
func matchTemplates(name string, entries map[string]*Entry, s *Source) *Match {
	for _, key := range slices.Sorted(maps.Keys(entries)) {
		e := entries[key]
		if !e.IsTemplate {
			continue
		}
		if vars, ok := NewTemplate(e.TargetTemplate).Match(name, e.PH); ok {
			return &Match{Entry: e, Source: s, Vars: vars}
		}
	}
	return nil
}

// Versions lists the versions of a name, newest first, unioned across sources.
//
// It serves completion rather than creation: what could be typed next. A name
// with no version axis has no answer — grch38/star-gencode expands over three
// placeholders and none of them is "the version" — so it comes back empty and a
// caller wanting the family uses Enumerate.
func (c Catalog) Versions(ctx context.Context, name string) ([]string, error) {
	name = strings.TrimSuffix(Normalize(name), "/")
	if name == "" {
		return nil, nil
	}

	var out []string
	seen := map[string]bool{}
	add := func(v string) {
		if v != "" && v != "*" && !seen[v] {
			seen[v] = true
			out = append(out, v)
		}
	}

	for _, s := range c {
		entries, err := s.Entries(ctx)
		if err != nil {
			continue
		}
		for key, e := range entries {
			// A plain recipe has one index key per version.
			if rest, ok := strings.CutPrefix(key, name+"/"); ok && !strings.Contains(rest, "/") {
				add(rest)
			}
			if axis, ok := versionAxis(name, e); ok {
				for _, v := range e.PH[axis] {
					add(v)
				}
			}
		}
	}

	slices.SortStableFunc(out, func(a, b string) int { return CompareVersions(b, a) })
	return out, nil
}

// versionAxis reports the placeholder that is a name's version: the target must
// be exactly the name followed by one placeholder, as in ubuntu24/r/{version}.
func versionAxis(name string, e *Entry) (string, bool) {
	if !e.IsTemplate {
		return "", false
	}
	t := NewTemplate(e.TargetTemplate)
	if len(t.parts) != 2 || t.parts[0].IsToken || !t.parts[1].IsToken {
		return "", false
	}
	if t.parts[0].Text != name+"/" {
		return "", false
	}
	return t.parts[1].Text, true
}

// Open resolves a name, fetches its recipe and expands it.
//
// vars must complete whatever the name left open. A full target name supplies
// its own, so vars may be nil; a bare template name supplies none, and the
// caller passes what it chose.
func (c Catalog) Open(ctx context.Context, name string, vars map[string]string) (*Recipe, error) {
	m, found, err := c.Lookup(ctx, name)
	if err != nil {
		return nil, err
	}
	if !found {
		return nil, fmt.Errorf("%w: %s", ErrNotProvided, Normalize(name))
	}

	data, err := c.ReadPath(ctx, m.Source, m.Entry.Path)
	if err != nil {
		return nil, err
	}
	rec, err := ParseRecipe(m.Entry.Path, strings.NewReader(string(data)))
	if err != nil {
		return nil, err
	}
	if !rec.IsTemplate {
		return rec, nil
	}

	// What the name recovered comes first; what the caller chose overrides it.
	merged := make(map[string]string, len(m.Vars)+len(vars))
	maps.Copy(merged, m.Vars)
	maps.Copy(merged, vars)
	return Expand(rec, merged)
}

// ReadPath moves bytes from one source, reusing the fetch and cache path
// without this package growing an API for what they contain.
func (c Catalog) ReadPath(ctx context.Context, s *Source, path string) ([]byte, error) {
	return s.b.read(ctx, path)
}
