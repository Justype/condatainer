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

	// The heredoc's lines count too: position carries no meaning, so a recipe
	// that writes a job script also declares whatever that script declares.
	wantDeps := []string{"star/{star_version}", "samtools/1.23.1>=1.10", "not/a/dep"}
	if !slices.Equal(r.Deps, wantDeps) {
		t.Errorf("Deps = %v, want %v", r.Deps, wantDeps)
	}
	wantDirectives := []string{"#SBATCH --cpus-per-task=8", "#SBATCH --mem=64G", "#SBATCH --this-is-not-a-directive"}
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

// Position never disqualifies an annotation. A recipe that writes a job script
// therefore adopts that script's declarations as its own — the cost of the rule,
// recorded here so it is a decision rather than a surprise.
func TestParseRecipeReadsAnnotationsAnywhere(t *testing.T) {
	r := parse(t, "grch38/star-gencode", starRecipe)
	if !slices.Contains(r.Deps, "not/a/dep") {
		t.Errorf("Deps = %v, want the heredoc's #DEP: included", r.Deps)
	}
	if !slices.Contains(r.Directives, "#SBATCH --this-is-not-a-directive") {
		t.Errorf("Directives = %v, want the heredoc's #SBATCH included", r.Directives)
	}

	// Nothing ends the annotations, so a key after an executable line is read.
	below := parse(t, "x/y", "#DESC:kept\nBootstrap: docker\n#URL:below\n")
	if below.Description != "kept" || below.URL != "below" {
		t.Errorf("description=%q url=%q, want both read", below.Description, below.URL)
	}

	// A .def has no shebang and needs no exception.
	def := parse(t, "ubuntu24/base.def", "#DESC:base\n\nBootstrap: docker\n")
	if def.Name != "ubuntu24/base" || def.Type != TypeBase || def.Description != "base" {
		t.Errorf("def = %q/%q/%q", def.Name, def.Type, def.Description)
	}
}

// Annotations are found wherever they are written, so what qualifies as one is
// pinned exactly: position never disqualifies a line, and shape always can.
func TestScanAnnotations(t *testing.T) {
	tests := []struct {
		name string
		text string
		want []Annotation
	}{
		{"a leading annotation", "#DESC: hello\n", []Annotation{{Key: "#DESC", Value: "hello", Line: 1}}},
		{"one below the first command", "echo hi\n#DEP: star/2.7\n",
			[]Annotation{{Key: "#DEP", Value: "star/2.7", Line: 2}}},
		{"one inside a heredoc", "cat <<'EOF'\n#DEP: star/2.7\nEOF\n",
			[]Annotation{{Key: "#DEP", Value: "star/2.7", Line: 2}}},
		{"indented", "  \t#DEP: star/2.7\n", []Annotation{{Key: "#DEP", Value: "star/2.7", Line: 1}}},
		{"a note is kept apart from the value", "#ENV: PATH={prefix}/bin  ## on PATH\n",
			[]Annotation{{Key: "#ENV", Value: "PATH={prefix}/bin", Note: "on PATH", Line: 1}}},
		{"an inline comment is stripped", "#DEP: star/2.7  # the aligner\n",
			[]Annotation{{Key: "#DEP", Value: "star/2.7", Line: 1}}},
		{"a shebang is not an annotation", "#!/bin/bash\n", nil},
		{"lower-case prose is not an annotation", "# note: rerun weekly\n", nil},
		{"a bare key without '#' is not an annotation", "TYPE:app\n", nil},
		{"a comment with no colon is not an annotation", "#DESC no colon\n", nil},
		{"a directive is not an annotation", "#SBATCH --time=01:00:00\n", nil},
		{"CRLF", "#DESC: hello\r\n", []Annotation{{Key: "#DESC", Value: "hello", Line: 1}}},
		{"empty", "", nil},
	}
	for _, tt := range tests {
		t.Run(tt.name, func(t *testing.T) {
			got := ScanAnnotations([]byte(tt.text))
			if len(got) != len(tt.want) {
				t.Fatalf("ScanAnnotations = %#v, want %#v", got, tt.want)
			}
			for i := range tt.want {
				if got[i] != tt.want[i] {
					t.Errorf("[%d] = %#v, want %#v", i, got[i], tt.want[i])
				}
			}
		})
	}
}

// A directive carries colons in its value, so it is matched by prefix rather
// than cut on the first colon like an annotation.
func TestDirectivesAreNotAnnotations(t *testing.T) {
	recipe, err := ParseRecipe("x/1.0", strings.NewReader("#SBATCH --time=01:00:00\n#DEP: star/2.7\n"))
	if err != nil {
		t.Fatal(err)
	}
	if len(recipe.Directives) != 1 || recipe.Directives[0] != "#SBATCH --time=01:00:00" {
		t.Errorf("directives = %#v", recipe.Directives)
	}
	if len(recipe.Deps) != 1 || recipe.Deps[0] != "star/2.7" {
		t.Errorf("deps = %#v", recipe.Deps)
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
