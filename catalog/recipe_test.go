package catalog

import (
	"slices"
	"strings"
	"testing"
)

const starRecipe = `#!/usr/bin/env bash
# Build a STAR index. This comment must not end the header block.
#DESC:STAR {star_version} index for GENCODE {gencode_version}
#URL:https://github.com/alexdobin/STAR

#TARGET:grch38/star/{star_version}/gencode{gencode_version}-{read_length}
#PH:star_version:2.7.11b,2.7.11a,2.7.9a
#PH:gencode_version:47-49
#PH:read_length:101|151|*

#DEP:star/{star_version}
#DEP:samtools/1.23.1>=1.10
#ENV:STAR_INDEX={prefix}/index   ## pass to --genomeDir
#INPUT:paste the download link

#SBATCH --cpus-per-task=8
#SBATCH --mem=64G

set -euo pipefail
cat > "$CNT_TMP/job.sh" <<'EOF'
#SBATCH --this-is-not-a-directive
#DEP:not/a/dep
EOF
STAR --runThreadN "$NCPUS" --genomeDir "$CNT_PREFIX/index"
`

func parse(t *testing.T, path, body string) *Recipe {
	t.Helper()
	r, err := ParseRecipe(path, strings.NewReader(body))
	if err != nil {
		t.Fatalf("ParseRecipe(%q): %v", path, err)
	}
	return r
}

func TestParseRecipe(t *testing.T) {
	r := parse(t, "grch38/star-gencode", starRecipe)

	if r.Name != "grch38/star-gencode" || r.Path != "grch38/star-gencode" {
		t.Errorf("name/path = %q/%q", r.Name, r.Path)
	}
	if r.Type != TypeData {
		t.Errorf("Type = %q, want data (from the target's slash count)", r.Type)
	}
	if !r.IsTemplate {
		t.Error("IsTemplate = false")
	}
	if r.Description != "STAR {star_version} index for GENCODE {gencode_version}" {
		t.Errorf("Description = %q", r.Description)
	}
	if r.URL != "https://github.com/alexdobin/STAR" {
		t.Errorf("URL = %q", r.URL)
	}

	wantDeps := []string{"star/{star_version}", "samtools/1.23.1>=1.10"}
	if !slices.Equal(r.Deps, wantDeps) {
		t.Errorf("Deps = %v, want %v", r.Deps, wantDeps)
	}
	wantDirectives := []string{"#SBATCH --cpus-per-task=8", "#SBATCH --mem=64G"}
	if !slices.Equal(r.Directives, wantDirectives) {
		t.Errorf("Directives = %v, want %v", r.Directives, wantDirectives)
	}

	if len(r.Env) != 1 || r.Env[0].Key != "STAR_INDEX" || r.Env[0].Note != "pass to --genomeDir" {
		t.Fatalf("Env = %+v", r.Env)
	}
	if got := r.Env[0].Value(nil); got != "{prefix}/index" {
		t.Errorf("Env value = %q, want the token intact", got)
	}
	if len(r.Inputs) != 1 || r.Inputs[0] != "paste the download link" {
		t.Errorf("Inputs = %+v", r.Inputs)
	}
}

// A plain comment must not end the block, and nothing below the block is a
// header — including a heredoc that writes a job script.
func TestParseRecipeHeaderBlock(t *testing.T) {
	r := parse(t, "grch38/star-gencode", starRecipe)
	if len(r.Deps) != 2 {
		t.Errorf("Deps = %v, want the heredoc's #DEP: ignored", r.Deps)
	}
	if len(r.Directives) != 2 {
		t.Errorf("Directives = %v, want the heredoc's #SBATCH ignored", r.Directives)
	}

	// The block ends at the first line that is neither comment nor blank.
	below := parse(t, "x/y", "#DESC:kept\nBootstrap: docker\n#URL:below the block\n")
	if below.Description != "kept" || below.URL != "" {
		t.Errorf("description=%q url=%q, want the header below Bootstrap ignored", below.Description, below.URL)
	}

	// A .def has no shebang and needs no exception.
	def := parse(t, "ubuntu24/base.def", "#DESC:base\n\nBootstrap: docker\n")
	if def.Name != "ubuntu24/base" || def.Type != TypeBase || def.Description != "base" {
		t.Errorf("def = %q/%q/%q", def.Name, def.Type, def.Description)
	}
}

func TestParseRecipeDescriptionUsesDescOnly(t *testing.T) {
	rec := parse(t, "x/y", "#DESCRIPTION:legacy\n#DESC:short\n")
	if rec.Description != "short" {
		t.Fatalf("description = %q, want DESC value", rec.Description)
	}

	legacy := parse(t, "x/y", "#DESCRIPTION:legacy\n")
	if legacy.Description != "" {
		t.Fatalf("legacy DESCRIPTION parsed as %q", legacy.Description)
	}
}

func TestParseValues(t *testing.T) {
	tests := []struct {
		name string
		raw  string
		want []string
	}{
		{"comma sorts newest first", "3.1.3,4.0.0,3.6.2", []string{"4.0.0", "3.6.2", "3.1.3"}},
		{"range expands and sorts", "47-49", []string{"49", "48", "47"}},
		{"pipe keeps written order", "101|151", []string{"101", "151"}},
		{"star is always last", "101,151,*", []string{"151", "101", "*"}},
		{"pipe with star", "101|151|*", []string{"101", "151", "*"}},
		{"duplicates collapse", "1.0,1.0,2.0", []string{"2.0", "1.0"}},
		{"reversed range", "49-47", []string{"49", "48", "47"}},
		{"single value", "2026.06.0-242", []string{"2026.06.0-242"}},
		// A range is digits-only, so a literal carrying a dash survives whole.
		{"dashed literal is not a range", "2024-A", []string{"2024-A"}},
		{"range expands in pipe form too", "a | 1-3", []string{"a", "1", "2", "3"}},
		{"surrounding space trimmed", "kasm | turbo", []string{"kasm", "turbo"}},
		{"empty yields no values", "  ,  ", nil},
	}
	for _, tt := range tests {
		if got := ParseValues(tt.raw); !slices.Equal(got, tt.want) {
			t.Errorf("%s: ParseValues(%q) = %v, want %v", tt.name, tt.raw, got, tt.want)
		}
	}
}

func TestParseRecipePH(t *testing.T) {
	r := parse(t, "grch38/star-gencode", starRecipe)
	want := map[string][]string{
		"star_version":    {"2.7.11b", "2.7.11a", "2.7.9a"},
		"gencode_version": {"49", "48", "47"},
		"read_length":     {"101", "151", "*"},
	}
	for key, vals := range want {
		if !slices.Equal(r.PH[key], vals) {
			t.Errorf("PH[%q] = %v, want %v", key, r.PH[key], vals)
		}
	}
}
