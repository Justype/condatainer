package key

import (
	"strings"
	"testing"

	"github.com/Justype/condatainer/internal/artifact/meta"
)

const templateRecipe = `#!/usr/bin/env bash
# Build a STAR index.
#DESC:STAR {star_version} index for GENCODE {gencode_version}
#TARGET:grch38/star/{star_version}/gencode{gencode_version}
#PH:star_version:2.7.11b,2.7.11a
#PH:gencode_version:47-49
#DEP:star/{star_version}
#ENV:STAR_INDEX={prefix}/index   ## pass to --genomeDir
#SBATCH --mem=64G

set -euo pipefail
STAR --runMode genomeGenerate --genomeDir "$CNT_PREFIX/index"
`

// Nothing that only describes a recipe may move its digest — least of all the
// parts tooling rewrites without being asked.
func TestRecipeDigestIgnoresDescription(t *testing.T) {
	base := RecipeDigest([]byte(templateRecipe))
	if !ValidDigest(base) {
		t.Fatalf("digest = %q, not a well-formed digest", base)
	}

	unchanged := map[string]string{
		"a #PH: menu autoupdate grew": strings.Replace(templateRecipe,
			"#PH:star_version:2.7.11b,2.7.11a", "#PH:star_version:2.7.12,2.7.11b,2.7.11a", 1),
		"a reworded #DESC:": strings.Replace(templateRecipe,
			"#DESC:STAR {star_version} index", "#DESC:STAR {star_version} genome index", 1),
		"a changed directive": strings.Replace(templateRecipe, "--mem=64G", "--mem=128G", 1),
		"a reworded comment":  strings.Replace(templateRecipe, "# Build a STAR index.", "# Builds it.", 1),
		"a dropped shebang":   strings.Replace(templateRecipe, "#!/usr/bin/env bash\n", "", 1),
	}
	for what, edited := range unchanged {
		if edited == templateRecipe {
			t.Fatalf("%s: the fixture did not change", what)
		}
		if got := RecipeDigest([]byte(edited)); got != base {
			t.Errorf("%s moved the digest", what)
		}
	}

	moved := map[string]string{
		"a changed command": strings.Replace(templateRecipe, "--runMode", "--runThreadN 8 --runMode", 1),
		"a dropped set -e":  strings.Replace(templateRecipe, "set -euo pipefail\n", "", 1),
	}
	for what, edited := range moved {
		if RecipeDigest([]byte(edited)) == base {
			t.Errorf("%s did not move the digest", what)
		}
	}
}

// A template is stored with its tokens intact and the expansion is never hashed,
// so every variant shares one recipe digest and is told apart by ph= alone. That
// is the whole reason placeholders are a separate group of lines.
func TestTemplateVariantsShareADigestAndDifferByPlaceholders(t *testing.T) {
	key := func(selected map[string]string) (string, string) {
		t.Helper()
		r := Model{
			Kind:         KindIdentity,
			Type:         "data",
			Recipe:       RecipeDigest([]byte(templateRecipe)),
			Placeholders: Placeholders(selected),
		}
		k, err := Key(r)
		if err != nil {
			t.Fatalf("Key: %v", err)
		}
		return r.Recipe, k
	}

	gencode49, key49 := key(map[string]string{"star_version": "2.7.11b", "gencode_version": "49"})
	gencode47, key47 := key(map[string]string{"star_version": "2.7.11b", "gencode_version": "47"})

	if gencode49 != gencode47 {
		t.Error("two variants of one template got different recipe digests")
	}
	if key49 == key47 {
		t.Error("two variants of one template got the same key")
	}

	// Sorted by name, so the record they land in is canonical either way.
	lines := Placeholders(map[string]string{"star_version": "2.7.11b", "gencode_version": "49"})
	if len(lines) != 2 || lines[0].Name != "gencode_version" || lines[1].Name != "star_version" {
		t.Errorf("placeholders are not sorted by name: %v", lines)
	}
	if Placeholders(nil) != nil {
		t.Error("a recipe with no placeholders should contribute no lines")
	}
}

// #ENV: lives in the header, which no digest sees, so these lines are the only
// thing standing between an edited variable and a key that never moves.
func TestEnvLines(t *testing.T) {
	env := []meta.EnvVar{
		{Key: "STAR_INDEX", Value: "{prefix}/index", Note: "pass to --genomeDir"},
		{Key: "GENOME_FASTA", Value: "{prefix}/genome.fa"},
	}
	lines := Env(env)
	if len(lines) != 2 {
		t.Fatalf("lines = %v", lines)
	}
	if lines[0].Key != "GENOME_FASTA" || lines[1].Key != "STAR_INDEX" {
		t.Errorf("lines are not sorted by key: %v", lines)
	}
	// {prefix} stays unsubstituted, which is what keeps the artifact's name out
	// of a record that must not contain it.
	if lines[1].Value != "{prefix}/index" {
		t.Errorf("value = %q, want the token intact", lines[1].Value)
	}
	// The ## note is descriptive; rewording one must not mint a new artifact.
	for _, l := range lines {
		if strings.Contains(l.Value, "genomeDir") {
			t.Errorf("a note reached a record line: %v", l)
		}
	}
	if Env(nil) != nil {
		t.Error("an image with no #ENV: should contribute no lines")
	}
}

// The derivation exists so a build and a later comparison agree. This is that
// agreement, end to end: the lines go into a record and the record hashes.
func TestDerivationFeedsCanonicalPreimage(t *testing.T) {
	r := Model{
		Kind:         KindIdentity,
		Type:         "data",
		Env:          Env([]meta.EnvVar{{Key: "STAR_INDEX", Value: "{prefix}/index", Note: "n"}}),
		Recipe:       RecipeDigest([]byte(templateRecipe)),
		Placeholders: Placeholders(map[string]string{"star_version": "2.7.11b"}),
	}
	data, err := Marshal(r)
	if err != nil {
		t.Fatalf("Marshal: %v", err)
	}
	if !strings.Contains(string(data), "ph=star_version=2.7.11b") {
		t.Errorf("placeholder line missing:\n%s", data)
	}
	if !strings.Contains(string(data), "env=STAR_INDEX={prefix}/index") {
		t.Errorf("env line missing:\n%s", data)
	}
}
