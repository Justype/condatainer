package catalog

import (
	"encoding/json"
	"maps"
	"net/http"
	"net/http/httptest"
	"os"
	"path/filepath"
	"slices"
	"testing"
	"time"
)

// fakeSource writes a minimal collection: a descriptor, two recipes, and the
// index an HTTP source would be served from.
func fakeSource(t *testing.T) string {
	t.Helper()
	root := t.TempDir()

	write := func(rel, body string) {
		p := filepath.Join(root, filepath.FromSlash(rel))
		if err := os.MkdirAll(filepath.Dir(p), 0o755); err != nil {
			t.Fatal(err)
		}
		if err := os.WriteFile(p, []byte(body), 0o644); err != nil {
			t.Fatal(err)
		}
	}

	write("source.json", `{"schema":1,"repository":"https://example.invalid/r","default_base":"ubuntu24"}`)
	write("recipes/cellranger/9.0.1", "#DESC:cellranger\n#URL:https://example.invalid\n")
	write("recipes/ubuntu24/base.def", "#DESC:base\n\nBootstrap: docker\n")
	write("recipes/grch38/star-gencode", starRecipe)
	write("recipes/README.md", "not a recipe")

	// The index an HTTP backend reads, built from the same files.
	cat, err := Open(t.Context(), []Spec{{Name: "local", Base: root}}, Cache{})
	if err != nil {
		t.Fatal(err)
	}
	entries := cat.Entries(t.Context())
	data, err := json.Marshal(entries)
	if err != nil {
		t.Fatal(err)
	}
	write("index/recipes.json", string(data))
	return root
}

func TestOpenDirSource(t *testing.T) {
	root := fakeSource(t)
	cat, err := Open(t.Context(), []Spec{{Name: "local", Base: root}}, Cache{})
	if err != nil {
		t.Fatal(err)
	}

	if got := cat.DefaultBase(); got != "ubuntu24" {
		t.Errorf("DefaultBase = %q, want ubuntu24", got)
	}
	if cat[0].Desc.Repository != "https://example.invalid/r" {
		t.Errorf("descriptor not loaded: %+v", cat[0].Desc)
	}

	entries := cat.Entries(t.Context())
	want := []string{"cellranger/9.0.1", "grch38/star-gencode", "ubuntu24/base"}
	if got := slices.Sorted(maps.Keys(entries)); !slices.Equal(got, want) {
		t.Errorf("entries = %v, want %v (README.md skipped)", got, want)
	}
	if e := entries["ubuntu24/base"]; e.Type != TypeBase || e.Path != "recipes/ubuntu24/base.def" {
		t.Errorf("base entry = %+v", e)
	}
	if e := entries["grch38/star-gencode"]; !e.IsTemplate || e.Type != TypeData {
		t.Errorf("template entry = %+v", e)
	}
}

// Both backends must produce the same Entry for the same recipe, which is the
// one place the index and the walk can drift apart.
func TestBackendsAgree(t *testing.T) {
	root := fakeSource(t)
	srv := httptest.NewServer(http.FileServer(http.Dir(root)))
	defer srv.Close()

	local, err := Open(t.Context(), []Spec{{Name: "local", Base: root}}, Cache{})
	if err != nil {
		t.Fatal(err)
	}
	remote, err := Open(t.Context(), []Spec{{Name: "remote", Base: srv.URL}},
		Cache{Dir: t.TempDir(), TTL: time.Minute})
	if err != nil {
		t.Fatal(err)
	}
	if remote[0].Desc.DefaultBase != "ubuntu24" {
		t.Errorf("descriptor not fetched over HTTP: %+v", remote[0].Desc)
	}

	lEntries := local.Entries(t.Context())
	rEntries := remote.Entries(t.Context())
	if len(lEntries) != len(rEntries) {
		t.Fatalf("walk found %d entries, index %d", len(lEntries), len(rEntries))
	}
	for name, l := range lEntries {
		r, ok := rEntries[name]
		if !ok {
			t.Errorf("%s: missing from the index", name)
			continue
		}
		if l.Type != r.Type || l.Path != r.Path || l.Description != r.Description || l.URL != r.URL ||
			l.IsTemplate != r.IsTemplate || l.TargetTemplate != r.TargetTemplate ||
			!slices.Equal(l.Deps, r.Deps) {
			t.Errorf("%s:\n walk  %+v\n index %+v", name, l, r)
		}
		for k, v := range l.PH {
			if !slices.Equal(r.PH[k], v) {
				t.Errorf("%s: ph[%s] walk=%v index=%v", name, k, v, r.PH[k])
			}
		}
	}
}

func TestOpenRejectsBadSpec(t *testing.T) {
	if _, err := Open(t.Context(), nil, Cache{}); err == nil {
		t.Error("Open with no sources should fail")
	}
	if _, err := Open(t.Context(), []Spec{{Name: "x"}}, Cache{}); err == nil {
		t.Error("Open with an empty base should fail")
	}
}

// An unreachable source keeps its place, so first-wins ordering cannot silently
// promote the next source's recipes.
func TestUnreachableSourceKeepsItsPlace(t *testing.T) {
	root := fakeSource(t)
	cat, err := Open(t.Context(), []Spec{
		{Name: "dead", Base: filepath.Join(root, "does-not-exist")},
		{Name: "local", Base: root},
	}, Cache{})
	if err != nil {
		t.Fatal(err)
	}
	if len(cat) != 2 || cat[0].Name != "dead" {
		t.Fatalf("catalog = %v", cat)
	}
	entries := cat.Entries(t.Context())
	if cat[0].Err == nil {
		t.Error("unreachable source should carry its error")
	}
	if _, ok := entries["cellranger/9.0.1"]; !ok {
		t.Error("the reachable source should still contribute")
	}
}
