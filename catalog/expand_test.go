package catalog

import (
	"slices"
	"strings"
	"testing"
)

func TestExpand(t *testing.T) {
	r := parse(t, "grch38/star-gencode", starRecipe)
	vars := map[string]string{"star_version": "2.7.11b", "gencode_version": "47", "read_length": "101"}

	got, err := Expand(r, vars)
	if err != nil {
		t.Fatal(err)
	}

	if want := "grch38/star/2.7.11b/gencode47-101"; got.Name != want {
		t.Errorf("Name = %q, want %q", got.Name, want)
	}
	if got.IsTemplate || got.TargetTemplate != "" {
		t.Error("expanded recipe is still marked a template")
	}
	if want := "STAR 2.7.11b index for GENCODE 47"; got.Description != want {
		t.Errorf("Description = %q, want %q", got.Description, want)
	}
	if want := []string{"star/2.7.11b", "samtools/1.23.1>=1.10"}; !slices.Equal(got.Deps, want) {
		t.Errorf("Deps = %v, want %v", got.Deps, want)
	}
	// PH records the values chosen, one per placeholder.
	for k, v := range vars {
		if !slices.Equal(got.PH[k], []string{v}) {
			t.Errorf("PH[%q] = %v, want [%s]", k, got.PH[k], v)
		}
	}

	// The body is substituted, since the recipe hash is taken over the expanded
	// text — but a shell ${VAR} is not a placeholder and must survive.
	body := string(got.Text)
	if strings.Contains(body, "{star_version}") {
		t.Error("Text still contains an unexpanded placeholder")
	}
	if !strings.Contains(body, `"$CNT_PREFIX/index"`) || !strings.Contains(body, `"$NCPUS"`) {
		t.Error("Text lost a shell variable")
	}

	// {prefix} has no var and survives to the manifest.
	if v := got.Env[0].Value(nil); v != "{prefix}/index" {
		t.Errorf("Env value = %q, want {prefix} intact", v)
	}

	// The original is untouched.
	if r.Name != "grch38/star-gencode" || !r.IsTemplate {
		t.Error("Expand mutated its input")
	}
}

func TestExpandRequiresEveryVar(t *testing.T) {
	r := parse(t, "grch38/star-gencode", starRecipe)
	if _, err := Expand(r, map[string]string{"star_version": "2.7.11b"}); err == nil {
		t.Error("Expand with an incomplete var set should fail")
	}
	plain := parse(t, "cellranger/9.0.1", "#DESC:cellranger\n")
	if _, err := Expand(plain, nil); err == nil {
		t.Error("Expand on a non-template should fail")
	}
}
