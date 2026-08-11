// Recipe helpers convert inputs shared by the recipe-backed schemes into
// canonical values. Individual scheme files decide whether to include them.
package key

import (
	"sort"

	"github.com/Justype/condatainer/catalog"
	"github.com/Justype/condatainer/internal/artifact/meta"
)

// RecipeDigest returns the digest included by all recipe-backed V1 schemes.
//
// The preimage is the recipe with its whole-line comments removed — what the
// recipe does, with everything that only describes it gone. The header goes with
// them, and that is deliberate: every header item that belongs in a key already
// has its own canonical model field. Hashing the whole file would re-admit
// values deliberately excluded by individual schemes. A #PH: menu that #AUTOUPDATE: grew, a reworded #DESC:, a
// changed scheduler directive: none changes what a build produces, and a key that
// moved for them would mint a new artifact every time a robot edited a recipe.
//
// The input is the recipe as stored at /.cnt/recipe — a template with its tokens
// intact, never the expansion. Which variant was built is carried by Placeholders.
func RecipeDigest(recipe []byte) string {
	return Digest(catalog.StripComments(recipe))
}

// Placeholders converts selected #PH: values into canonical model fields sorted
// by name. All recipe-backed V1 schemes include them because they select the
// concrete variant built from a template.
func Placeholders(selected map[string]string) []PlaceholderValue {
	if len(selected) == 0 {
		return nil
	}
	out := make([]PlaceholderValue, 0, len(selected))
	for name, value := range selected {
		out = append(out, PlaceholderValue{Name: name, Value: value})
	}
	sort.Slice(out, func(i, j int) bool { return out[i].Name < out[j].Name })
	return out
}

// Env converts #ENV: contributions into canonical model fields sorted by key,
// with {prefix} left unsubstituted.
//
// All recipe-backed V1 schemes include them. #ENV: lives in the header, which
// the recipe digest omits. Without these lines editing GENOME_FASTA would move no key at all — and since
// #ENV: is how one recipe finds another's output, that silence would propagate
// into everything built on top.
//
// The ## notes are dropped: they are descriptive, and rewording one must not mint
// a new artifact. Values are kept as well as names, because a value change is itself a
// build-input change. One shared conversion rule is sufficient.
func Env(env []meta.EnvVar) []EnvValue {
	if len(env) == 0 {
		return nil
	}
	out := make([]EnvValue, 0, len(env))
	for _, e := range env {
		out = append(out, EnvValue{Key: e.Key, Value: e.Value})
	}
	sort.Slice(out, func(i, j int) bool { return out[i].Key < out[j].Key })
	return out
}
