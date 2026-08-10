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

// The boundary decides what ParseRecipe reads, so it is pinned byte-exactly:
// everything about it that could plausibly drift gets a case.
func TestHeaderBoundary(t *testing.T) {
	tests := []struct {
		name string
		text string
		want string // the body, verbatim
	}{
		{"shebang is header", "#!/bin/bash\necho hi\n", "echo hi\n"},
		{"a blank line before the first command is header", "#DESC:x\n\n\necho hi\n", "echo hi\n"},
		{"indented comments are header", "  # note\n\t#DESC:x\necho hi\n", "echo hi\n"},
		{"a whitespace-only line is header", "#DESC:x\n   \necho hi\n", "echo hi\n"},
		{"an indented command starts the body", "#DESC:x\n  echo hi\n", "  echo hi\n"},
		{"a later comment stays in the body", "#DESC:x\necho hi\n# trailing\n", "echo hi\n# trailing\n"},
		{"a # inside a string does not move it", "#DESC:x\necho \"a # b\"\n", "echo \"a # b\"\n"},
		{"a ${x#y} expansion does not move it", "#DESC:x\necho \"${v#pre}\"\n", "echo \"${v#pre}\"\n"},
		{"a heredoc's directives stay in the body", "#DESC:x\ncat <<'EOF'\n#DEP:not/a/dep\nEOF\n", "cat <<'EOF'\n#DEP:not/a/dep\nEOF\n"},
		{"header only", "#DESC:x\n#URL:y\n", ""},
		{"empty", "", ""},
		{"no header", "echo hi\n", "echo hi\n"},
		{"no trailing newline", "#DESC:x\necho hi", "echo hi"},
		{"CRLF", "#DESC:x\r\n\r\necho hi\r\n", "echo hi\r\n"},
	}
	for _, tt := range tests {
		t.Run(tt.name, func(t *testing.T) {
			text := []byte(tt.text)
			if got := string(text[headerBoundary(text):]); got != tt.want {
				t.Errorf("body = %q, want %q", got, tt.want)
			}
		})
	}
}

// Header edits must not move the boundary, or a heredoc's #DEP: would stop being
// inert the moment someone reworded a description.
func TestBoundaryIsStableUnderHeaderEdits(t *testing.T) {
	body := func(text string) string {
		b := []byte(text)
		return string(b[headerBoundary(b):])
	}
	grown := strings.Replace(starRecipe,
		"#PH:star_version:2.7.11b,2.7.11a,2.7.9a",
		"#PH:star_version:2.7.11b,2.7.11a,2.7.9a,2.7.8a\n#DESC:reworded\n#SBATCH --time=4:00:00", 1)
	if body(grown) != body(starRecipe) {
		t.Error("growing a #PH: menu, rewording #DESC:, or adding a directive moved the boundary")
	}
}

// Both keys hash the comment-stripped recipe, so this is the one derivation the
// whole design rests on.
func TestStripComments(t *testing.T) {
	tests := []struct {
		name string
		body string
		want string
	}{
		{"the header goes with the rest", "#!/bin/bash\n#DESC:x\necho a\n", "echo a\n"},
		{"a whole-line comment goes", "echo a\n# note\necho b\n", "echo a\necho b\n"},
		{"an indented one goes", "echo a\n   \t# note\necho b\n", "echo a\necho b\n"},
		{"a heredoc's comment goes", "cat <<'EOF'\n#!/bin/bash\nexec x\nEOF\n", "cat <<'EOF'\nexec x\nEOF\n"},
		{"blank lines stay", "echo a\n\necho b\n", "echo a\n\necho b\n"},
		{"a trailing comment stays", "make -j2  # parallel\n", "make -j2  # parallel\n"},
		{"a # inside a string stays", "echo \"a # b\"\n", "echo \"a # b\"\n"},
		{"a ${x#y} expansion stays", "echo \"${v#pre}\"\n", "echo \"${v#pre}\"\n"},
		{"no trailing newline", "echo a\n# note", "echo a\n"},
		{"CRLF", "echo a\r\n# note\r\necho b\r\n", "echo a\r\necho b\r\n"},
		{"empty", "", ""},
		{"comments only", "# a\n# b\n", ""},
	}
	for _, tt := range tests {
		t.Run(tt.name, func(t *testing.T) {
			if got := string(StripComments([]byte(tt.body))); got != tt.want {
				t.Errorf("StripComments = %q, want %q", got, tt.want)
			}
		})
	}
}

// Nothing that only describes a recipe may reach a key. Every header item that
// should is carried by its own record line instead, so the preimage has to be
// blind to all of them — including the ones tooling rewrites unprompted.
func TestRecipePreimageIgnoresDescription(t *testing.T) {
	preimage := func(text string) string { return string(StripComments([]byte(text))) }
	base := preimage(starRecipe)

	cosmetic := map[string]string{
		"a #PH: menu autoupdate grew": strings.Replace(starRecipe,
			"#PH:star_version:2.7.11b,2.7.11a,2.7.9a",
			"#PH:star_version:2.7.12,2.7.11b,2.7.11a,2.7.9a", 1),
		"a reworded #DESC:": strings.Replace(starRecipe,
			"#DESC:STAR {star_version} index", "#DESC:STAR {star_version} genome index", 1),
		"a changed directive": strings.Replace(starRecipe,
			"#SBATCH --mem=64G", "#SBATCH --mem=128G", 1),
		"a reworded body comment": strings.Replace(starRecipe,
			"# Build a STAR index. This comment must not end the header block.",
			"# Builds the index.", 1),
	}
	for what, edited := range cosmetic {
		if edited == starRecipe {
			t.Fatalf("%s: the fixture did not change", what)
		}
		if got := preimage(edited); got != base {
			t.Errorf("%s moved the preimage:\n%q\nvs\n%q", what, got, base)
		}
	}

	if preimage(strings.Replace(starRecipe, "--runThreadN", "--runThreadN 1 --sjdbOverhang", 1)) == base {
		t.Error("a code change did not move the preimage")
	}
}
