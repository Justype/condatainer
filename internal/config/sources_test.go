package config

import (
	"strings"
	"testing"

	"github.com/Justype/condatainer/catalog"
	"github.com/spf13/viper"
)

func TestDecodeSourceList(t *testing.T) {
	// What viper hands back for a YAML sequence of single-key mappings.
	raw := []any{
		map[string]any{"cnt": "https://example.invalid/recipes/"},
		map[string]any{"lab": "/shared/lab/recipes"},
	}
	got := decodeSourceList(raw)
	if len(got) != 2 {
		t.Fatalf("got %d specs: %+v", len(got), got)
	}
	// Order comes from the sequence, the handle from the key, and a trailing
	// slash is trimmed so a base never doubles up when a path is joined.
	if got[0].Name != "cnt" || got[0].Base != "https://example.invalid/recipes" {
		t.Errorf("first = %+v", got[0])
	}
	if got[1].Name != "lab" || got[1].Base != "/shared/lab/recipes" {
		t.Errorf("second = %+v", got[1])
	}

	if got := decodeSourceList([]any{"lab=/shared/lab"}); len(got) != 1 ||
		got[0].Name != "lab" || got[0].Base != "/shared/lab" {
		t.Errorf("name=base form = %+v", got)
	}
	if got := decodeSourceList("not a list"); got != nil {
		t.Errorf("non-list = %+v, want nil", got)
	}
}

func TestParseSourceSpecs(t *testing.T) {
	got := parseSourceSpecs([]string{"cnt=https://a.invalid/", "lab=/l", "", "broken", "=/x", "n="})
	if len(got) != 2 {
		t.Fatalf("got %+v, want the two well-formed pairs", got)
	}
	if got[0].Base != "https://a.invalid" || got[1].Name != "lab" {
		t.Errorf("specs = %+v", got)
	}
}

func TestSelectSources(t *testing.T) {
	original := Global.Sources
	t.Cleanup(func() {
		Global.Sources = original
		ResetCatalog()
	})

	configured := []catalog.Spec{
		{Name: "site", Base: "/site"},
		{Name: "lab", Base: "/lab"},
		{Name: "cnt", Base: "https://example.invalid"},
	}

	t.Run("empty keeps configured order", func(t *testing.T) {
		Global.Sources = append([]catalog.Spec(nil), configured...)
		if err := SelectSources(nil); err != nil {
			t.Fatal(err)
		}
		if got := Global.Sources; len(got) != 3 || got[0].Name != "site" || got[2].Name != "cnt" {
			t.Fatalf("sources = %+v", got)
		}
	})

	t.Run("flag order wins and duplicates collapse", func(t *testing.T) {
		Global.Sources = append([]catalog.Spec(nil), configured...)
		if err := SelectSources([]string{"cnt", "site", "cnt"}); err != nil {
			t.Fatal(err)
		}
		if got := Global.Sources; len(got) != 2 || got[0].Name != "cnt" || got[1].Name != "site" {
			t.Fatalf("sources = %+v, want cnt then site", got)
		}
	})

	t.Run("unknown is an error and changes nothing", func(t *testing.T) {
		Global.Sources = append([]catalog.Spec(nil), configured...)
		err := SelectSources([]string{"lab", "missing"})
		if err == nil || !strings.Contains(err.Error(), "unknown source") {
			t.Fatalf("error = %v", err)
		}
		if got := Global.Sources; len(got) != 3 || got[0].Name != "site" {
			t.Fatalf("sources changed after error: %+v", got)
		}
	})
}

func TestLayerSourcesEnvWins(t *testing.T) {
	t.Setenv("CNT_SOURCES", "lab=/shared/lab|cnt=https://a.invalid")
	got := layerSources()
	if len(got) != 2 || got[0].Name != "lab" || got[1].Name != "cnt" {
		t.Errorf("layerSources = %+v", got)
	}
}

// Layers concatenate strongest first, and a name already taken is not repeated:
// a user entry shadows a site entry of the same handle.
func TestLayerSourcesMergesLayers(t *testing.T) {
	t.Setenv("CNT_SOURCES", "")

	user := viper.New()
	user.SetConfigType("yaml")
	if err := user.ReadConfig(stringReader(`
sources:
  - lab: /user/lab
`)); err != nil {
		t.Fatal(err)
	}
	site := viper.New()
	site.SetConfigType("yaml")
	if err := site.ReadConfig(stringReader(`
sources:
  - lab: /site/lab
  - cnt: https://site.invalid
`)); err != nil {
		t.Fatal(err)
	}

	saved := configLayers
	t.Cleanup(func() { configLayers = saved })
	configLayers = []*viper.Viper{user, site}

	got := layerSources()
	if len(got) != 2 {
		t.Fatalf("got %+v, want lab and cnt", got)
	}
	if got[0].Name != "lab" || got[0].Base != "/user/lab" {
		t.Errorf("first = %+v, want the user's lab to win", got[0])
	}
	if got[1].Name != "cnt" {
		t.Errorf("second = %+v, want the site's cnt to survive", got[1])
	}
}

// The default collection is always reachable, however the list is written, and
// always last so every configured entry outranks it.
func TestLayerSourcesAppendsDefault(t *testing.T) {
	useLayers := func(t *testing.T, yaml string) {
		t.Helper()
		v := viper.New()
		v.SetConfigType("yaml")
		if err := v.ReadConfig(stringReader(yaml)); err != nil {
			t.Fatal(err)
		}
		saved := configLayers
		t.Cleanup(func() { configLayers = saved })
		configLayers = []*viper.Viper{v}
	}

	// No `sources` key at all: a fresh install resolves recipes unconfigured.
	t.Run("absent", func(t *testing.T) {
		t.Setenv("CNT_SOURCES", "")
		useLayers(t, "base: ubuntu24\n")
		got := layerSources()
		if len(got) != 1 || got[0] != defaultSource {
			t.Fatalf("layerSources = %+v, want just the default", got)
		}
	})

	// An explicit list keeps its own order and gains the default at the end.
	t.Run("appended last", func(t *testing.T) {
		t.Setenv("CNT_SOURCES", "")
		useLayers(t, "sources:\n  - lab: /shared/lab\n")
		got := layerSources()
		if len(got) != 2 || got[0].Name != "lab" || got[1] != defaultSource {
			t.Fatalf("layerSources = %+v, want lab then the default", got)
		}
	})

	// Redefining the handle replaces the default rather than duplicating it:
	// that is how a site points `cnt` at its own collection.
	t.Run("redefined wins", func(t *testing.T) {
		t.Setenv("CNT_SOURCES", "")
		useLayers(t, "sources:\n  - cnt: /site/recipes\n")
		got := layerSources()
		if len(got) != 1 || got[0].Base != "/site/recipes" {
			t.Fatalf("layerSources = %+v, want only the site's cnt", got)
		}
	})

	// The env form takes the same treatment; it overrides the list, not the default.
	t.Run("env", func(t *testing.T) {
		t.Setenv("CNT_SOURCES", "lab=/shared/lab")
		got := layerSources()
		if len(got) != 2 || got[0].Name != "lab" || got[1] != defaultSource {
			t.Fatalf("layerSources = %+v, want lab then the default", got)
		}
	})
}

func stringReader(s string) *strings.Reader { return strings.NewReader(s) }
