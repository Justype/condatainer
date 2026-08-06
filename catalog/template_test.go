package catalog

import (
	"slices"
	"strings"
	"testing"
)

const starTarget = "grch38/star/{star_version}/gencode{gencode_version}-{read_length}"

var starPH = map[string][]string{
	"star_version":    {"2.7.11b", "2.7.11a", "2.7.9a"},
	"gencode_version": {"49", "48", "47"},
	"read_length":     {"101", "151", "*"},
}

func TestTemplateNames(t *testing.T) {
	got := NewTemplate(starTarget).Names()
	want := []string{"star_version", "gencode_version", "read_length"}
	if !slices.Equal(got, want) {
		t.Errorf("Names() = %v, want %v (target order)", got, want)
	}
}

func TestTemplateFill(t *testing.T) {
	tmpl := NewTemplate(starTarget)
	vars := map[string]string{"star_version": "2.7.11b", "gencode_version": "47", "read_length": "101"}
	got, err := tmpl.Fill(vars)
	if err != nil {
		t.Fatal(err)
	}
	if want := "grch38/star/2.7.11b/gencode47-101"; got != want {
		t.Errorf("Fill = %q, want %q", got, want)
	}

	delete(vars, "read_length")
	if _, err := tmpl.Fill(vars); err == nil {
		t.Error("Fill with a missing var should fail, not leave a token standing")
	}
}

func TestTemplateMatch(t *testing.T) {
	tmpl := NewTemplate(starTarget)
	tests := []struct {
		name string
		want map[string]string
	}{
		{"grch38/star/2.7.11b/gencode47-101", map[string]string{
			"star_version": "2.7.11b", "gencode_version": "47", "read_length": "101"}},
		// 2.7.11a and 2.7.11b share a prefix; alternation order must not matter.
		{"grch38/star/2.7.11a/gencode49-151", map[string]string{
			"star_version": "2.7.11a", "gencode_version": "49", "read_length": "151"}},
		// read_length is open-ended, so a value outside the set still matches.
		{"grch38/star/2.7.9a/gencode48-75", map[string]string{
			"star_version": "2.7.9a", "gencode_version": "48", "read_length": "75"}},
	}
	for _, tt := range tests {
		got, ok := tmpl.Match(tt.name, starPH)
		if !ok {
			t.Errorf("Match(%q) = false", tt.name)
			continue
		}
		for k, v := range tt.want {
			if got[k] != v {
				t.Errorf("Match(%q)[%q] = %q, want %q", tt.name, k, got[k], v)
			}
		}
	}

	// A closed set is its own pattern: an undeclared star version is not found.
	if _, ok := tmpl.Match("grch38/star/9.9.9/gencode47-101", starPH); ok {
		t.Error("Match accepted a star_version outside the #PH: set")
	}
	// A partially filled name determines nothing, so there is no partial form.
	for _, name := range []string{
		"grch38/star/2.7.11b",
		"grch38/star-gencode",
		"grch38/star/2.7.11b/gencode47-101/extra",
		"other/star/2.7.11b/gencode47-101",
	} {
		if _, ok := tmpl.Match(name, starPH); ok {
			t.Errorf("Match(%q) = true, want not found", name)
		}
	}
	// * must never swallow a path separator.
	if _, ok := tmpl.Match("grch38/star/2.7.9a/gencode48-75/deeper", starPH); ok {
		t.Error("open-ended value matched across a /")
	}
}

func TestTemplateEnumerate(t *testing.T) {
	got := NewTemplate(starTarget).Enumerate(starPH)
	// 3 star versions x 3 gencode versions x 2 concrete read lengths; * cannot
	// be enumerated.
	if len(got) != 18 {
		t.Fatalf("Enumerate produced %d variants, want 18", len(got))
	}
	if got[0].Name != "grch38/star/2.7.11b/gencode49-101" {
		t.Errorf("first variant = %q, want the newest of each axis", got[0].Name)
	}
	// Vars come back with the name, so nothing has to re-Match to recover them.
	for _, v := range got {
		filled, err := NewTemplate(starTarget).Fill(v.Vars)
		if err != nil || filled != v.Name {
			t.Errorf("variant %q does not round-trip through its vars (%v)", v.Name, v.Vars)
		}
	}
}

func TestTemplateRoundTrip(t *testing.T) {
	tmpl := NewTemplate(starTarget)
	for _, v := range tmpl.Enumerate(starPH) {
		vars, ok := tmpl.Match(v.Name, starPH)
		if !ok {
			t.Errorf("Match(%q) = false for an enumerated variant", v.Name)
			continue
		}
		for k, want := range v.Vars {
			if vars[k] != want {
				t.Errorf("%s: Match gave %s=%q, Enumerate gave %q", v.Name, k, vars[k], want)
			}
		}
	}
}

func TestValidateTemplate(t *testing.T) {
	tests := []struct {
		name   string
		target string
		ph     map[string][]string
		want   string // substring of the expected problem; empty means sound
	}{
		{name: "sound", target: starTarget, ph: starPH},
		{name: "no placeholders at all", target: "", ph: nil},
		{
			name:   "declared but not in target",
			target: "grch38/star/{star_version}",
			ph:     map[string][]string{"star_version": {"2.7.11b"}, "unused": {"1"}},
			want:   "#PH:unused is declared",
		},
		{
			name:   "in target but not declared",
			target: "grch38/star/{star_version}/{missing}",
			ph:     map[string][]string{"star_version": {"2.7.11b"}},
			want:   "there is no #PH:missing",
		},
		{
			name:   "PH without TARGET",
			target: "",
			ph:     map[string][]string{"version": {"1"}},
			want:   "there is no #TARGET:",
		},
		{
			name:   "adjacent open-ended",
			target: "x/{a}{b}",
			ph:     map[string][]string{"a": {"*"}, "b": {"*"}},
			want:   "adjacent and both are open-ended",
		},
		{
			name:   "adjacent but closed is fine",
			target: "x/{a}{b}",
			ph:     map[string][]string{"a": {"1"}, "b": {"2"}},
		},
	}
	for _, tt := range tests {
		got := ValidateTemplate(tt.target, tt.ph)
		if tt.want == "" {
			if got != nil {
				t.Errorf("%s: got problems %v, want none", tt.name, got)
			}
			continue
		}
		if !slices.ContainsFunc(got, func(p string) bool { return strings.Contains(p, tt.want) }) {
			t.Errorf("%s: got %v, want one containing %q", tt.name, got, tt.want)
		}
	}
}
