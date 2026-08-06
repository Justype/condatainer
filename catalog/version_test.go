package catalog

import (
	"slices"
	"testing"
)

func TestCompareVersions(t *testing.T) {
	tests := []struct {
		a, b string
		want int
	}{
		// The case a numeric-only comparator gets wrong: it stops at the first
		// non-numeric segment, making every 2.7.x STAR release equal.
		{"2.7.11b", "2.7.11a", 1},
		{"2.7.11a", "2.7.9a", 1},
		{"2.7.10b", "2.7.10a", 1},
		{"2.7.9a", "2.7.8a", 1},
		{"2.7.11b", "2.7.11b", 0},

		{"1.23.1", "1.10", 1},
		{"1.10", "1.9", 1},
		{"10", "9", 1},
		{"1.16", "1.16.0", -1}, // a prefix sorts below what extends it
		{"1.16rc1", "1.16", 1},
		{"", "1.0", -1},
		{"", "", 0},
		{"4.6.1", "4.5.3", 1},
		{"2026.06.0-242", "2025.12.1-100", 1},
		{"01", "1", 0}, // leading zeros are not significant

		// A pre-release sorts below its own release, so "newest" never picks an
		// alpha over the thing it precedes.
		{"1.1.0", "1.1.0-alpha", 1},
		{"1.1.0", "1.1.0-rc2", 1},
		{"1.1.0-beta", "1.1.0-alpha", 1},
		{"1.1.0-rc1", "1.1.0-beta", 1},
		{"1.1.0-rc2", "1.1.0-rc1", 1},
		{"1.1.0-alpha", "1.0.9", 1}, // still a 1.1.0
		// Labels order alphabetically, which is what alpha/beta/rc rely on and
		// what puts -dev after -beta rather than before -alpha.
		{"1.1.0-dev", "1.1.0-alpha", 1},
		{"1.1.0-rc1", "1.1.0-dev", 1},
		// A digit after the dash is an ordinary component, not a pre-release.
		{"2026.06.0-242", "2026.06.0", 1},
		{"2026.06.0-242", "2026.06.0-99", 1},
		// An uppercase label is still a pre-release.
		{"1.1.0", "1.1.0-RC1", 1},
		{"1.1.0", "1.1.0-Alpha", 1},
		{"1.1.0-RC2", "1.1.0-RC1", 1},
	}
	for _, tt := range tests {
		if got := CompareVersions(tt.a, tt.b); got != tt.want {
			t.Errorf("CompareVersions(%q, %q) = %d, want %d", tt.a, tt.b, got, tt.want)
		}
		if got := CompareVersions(tt.b, tt.a); got != -tt.want {
			t.Errorf("CompareVersions(%q, %q) = %d, want %d (antisymmetry)", tt.b, tt.a, got, -tt.want)
		}
	}
}

// TestCompareVersionsMatchesGenerator pins the comparator to the ordering the
// index generator produces, since ph[name][0] is the newest value and the
// resolver's "no preferred version means the newest" rule reads it.
func TestCompareVersionsMatchesGenerator(t *testing.T) {
	tests := []struct {
		name   string
		sorted []string // as recipe/scripts/generate_recipe_index.py emits them
	}{
		{"star", []string{"2.7.11b", "2.7.11a", "2.7.10b", "2.7.10a", "2.7.9a", "2.7.8a"}},
		{"r", []string{"4.6.1", "4.6.0", "4.5.3", "4.4.10", "4.4.3", "3.6.3", "3.1.3"}},
		{"gencode", []string{"49", "48", "47", "46", "10", "9"}},
	}
	for _, tt := range tests {
		got := slices.Clone(tt.sorted)
		slices.SortFunc(got, func(a, b string) int { return CompareVersions(b, a) })
		if !slices.Equal(got, tt.sorted) {
			t.Errorf("%s: sorted descending = %v, want %v", tt.name, got, tt.sorted)
		}
	}
}

// TestPreReleaseLadder pins the full ordering of one release's pre-releases.
func TestPreReleaseLadder(t *testing.T) {
	ladder := []string{"1.1.0-alpha", "1.1.0-beta", "1.1.0-rc1", "1.1.0-rc2", "1.1.0"}
	shuffled := []string{"1.1.0-rc2", "1.1.0", "1.1.0-alpha", "1.1.0-rc1", "1.1.0-beta"}
	slices.SortFunc(shuffled, CompareVersions)
	if !slices.Equal(shuffled, ladder) {
		t.Errorf("ascending = %v, want %v", shuffled, ladder)
	}
}
