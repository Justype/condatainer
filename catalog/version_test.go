package catalog

import (
	"os"
	"slices"
	"strconv"
	"strings"
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

// TestSharedVersionCases pins CompareVersions to the same table the recipe
// collection's generator is checked against.
//
// The index a collection publishes is consumed verbatim over HTTP while a local
// checkout is parsed and sorted here, so a disagreement between the two would
// make ph[0] — the default version offered to users — depend on how the
// collection was read. Skips when no collection is checked out beside the repo.
func TestSharedVersionCases(t *testing.T) {
	const cases = "../recipe/scripts/version_cases.txt"
	data, err := os.ReadFile(cases)
	if err != nil {
		t.Skip("no recipe collection checked out")
	}
	checked := 0
	for i, line := range strings.Split(string(data), "\n") {
		line = strings.TrimSpace(line)
		if line == "" || strings.HasPrefix(line, "#") {
			continue
		}
		parts := strings.Split(line, "|")
		if len(parts) != 3 {
			t.Errorf("%s:%d: malformed case %q", cases, i+1, line)
			continue
		}
		want, err := strconv.Atoi(parts[2])
		if err != nil {
			t.Errorf("%s:%d: %v", cases, i+1, err)
			continue
		}
		if got := CompareVersions(parts[0], parts[1]); got != want {
			t.Errorf("%s:%d: CompareVersions(%q, %q) = %d, want %d",
				cases, i+1, parts[0], parts[1], got, want)
		}
		checked++
	}
	t.Logf("%d shared cases agree", checked)
}
