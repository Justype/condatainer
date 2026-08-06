package cmd

import "testing"

// MatchesOrAliasWithText matches against a name and a description: AND terms may be
// split across the two, while anchored patterns still match the name on its own.
func TestMatchesOrAliasWithText(t *testing.T) {
	const (
		name        = "grch38/cellranger/2024-A"
		description = "Cell Ranger GRCh38 2024-A index"
	)

	cases := []struct {
		desc        string
		terms       []string
		name        string
		description string
		want        bool
	}{
		{"term in name only", []string{"cellranger"}, name, description, true},
		{"term in description only", []string{"index"}, name, description, true},
		{"AND split across both fields", []string{"cellranger", "index"}, name, description, true},
		{"AND with one term absent", []string{"cellranger", "bowtie"}, name, description, false},
		{"AND both terms in description", []string{"cell", "ranger"}, name, description, true},
		{"no term matches", []string{"salmon"}, name, description, false},
		{"empty description falls back to name", []string{"index"}, name, "", false},
		{"empty description still matches name", []string{"cellranger"}, name, "", true},
		// "*server" compiles to ^.*server$, which a name+description concatenation
		// would no longer satisfy.
		{"anchored wildcard matches name", []string{"*server"}, "ubuntu24/code-server", "Code Server IDE in the browser", true},
	}

	for _, tc := range cases {
		t.Run(tc.desc, func(t *testing.T) {
			q := NewSearchQuery(normalizeFilters(tc.terms), nil)
			if got := q.MatchesOrAliasWithText(tc.name, "", tc.description); got != tc.want {
				t.Errorf("MatchesOrAliasWithText(%q, \"\", %q) with terms %v = %v, want %v",
					tc.name, tc.description, tc.terms, got, tc.want)
			}
		})
	}
}

// An empty query matches everything, with or without a description.
func TestMatchesOrAliasWithText_EmptyQuery(t *testing.T) {
	q := NewSearchQuery(nil, nil)
	if !q.MatchesOrAliasWithText("anything", "", "") {
		t.Error("empty query should match every entry")
	}
}
