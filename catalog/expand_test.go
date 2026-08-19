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
	if got.IsTemplate {
		t.Error("expanded recipe is still marked a template")
	}
	// The target survives so a manifest can record which template produced this.
	if got.TargetTemplate != r.TargetTemplate {
		t.Errorf("TargetTemplate = %q, want it kept", got.TargetTemplate)
	}
	if want := "STAR 2.7.11b index for GENCODE 47"; got.Description != want {
		t.Errorf("Description = %q, want %q", got.Description, want)
	}
	// The fixture's heredoc declares one too; annotations are read wherever they
	// are written.
	if want := []string{"star/2.7.11b", "samtools/1.23.1>=1.10", "not/a/dep"}; !slices.Equal(got.Deps, want) {
		t.Errorf("Deps = %v, want %v", got.Deps, want)
	}
	// PH records the values chosen, one per placeholder.
	for k, v := range vars {
		if !slices.Equal(got.PH[k], []string{v}) {
			t.Errorf("PH[%q] = %v, want [%s]", k, got.PH[k], v)
		}
	}

	// Text is the template, tokens intact: it is what the artifact embeds and
	// what a rebuild starts from. Rendered is the copy the build runs.
	if !strings.Contains(string(got.Text), "{star_version}") {
		t.Error("Text was expanded; the template is what gets embedded")
	}
	rendered := string(got.Rendered)
	if strings.Contains(rendered, "{star_version}") {
		t.Error("Rendered still contains an unexpanded placeholder")
	}
	if string(got.Script()) != rendered {
		t.Error("Script() did not return the rendered copy")
	}
	// A shell ${VAR} is not a placeholder and must survive substitution.
	if !strings.Contains(rendered, `"$CNT_PREFIX/index"`) || !strings.Contains(rendered, `"$NCPUS"`) {
		t.Error("Rendered lost a shell variable")
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

// A recipe that is not a template runs the bytes it was fetched as, so Script()
// and Text are the same thing and nothing has to know which case it is in.
func TestScriptFallsBackToText(t *testing.T) {
	plain := parse(t, "cellranger/9.0.1", "#DESC:cellranger\necho build\n")
	if string(plain.Script()) != string(plain.Text) {
		t.Errorf("Script() = %q, want the recipe itself", plain.Script())
	}
}
