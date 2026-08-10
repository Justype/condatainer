// Package key derives the lines an artifact's identity and equivalence records
// are built from. It is the single definition of every rule about what reaches a
// key, so a build and a later comparison cannot disagree about what an artifact's
// keys should be.
package key

import (
	"sort"

	"github.com/Justype/condatainer/catalog"
	"github.com/Justype/condatainer/internal/artifact/meta"
	"github.com/Justype/condatainer/internal/artifact/record"
)

// RecipeDigest returns the digest a recipe build carries in both records.
//
// The preimage is the recipe with its whole-line comments removed — what the
// recipe does, with everything that only describes it gone. The header goes with
// them, and that is deliberate: every header item that belongs in a key already
// has its own record line, and hashing the file too would re-admit the ones left
// out on purpose. A #PH: menu that #AUTOUPDATE: grew, a reworded #DESC:, a
// changed scheduler directive: none changes what a build produces, and a key that
// moved for them would mint a new artifact every time a robot edited a recipe.
//
// The input is the recipe as stored at /.cnt/recipe — a template with its tokens
// intact, never the expansion. Which variant was built is carried by Placeholders.
func RecipeDigest(recipe []byte) string {
	return record.Digest(catalog.StripComments(recipe))
}

// Placeholders converts the selected #PH: values into record lines, sorted by
// name. Both records carry the same ones: a placeholder is a build input under
// either reading, and without them every variant of one template would share a
// single key.
func Placeholders(selected map[string]string) []record.Placeholder {
	if len(selected) == 0 {
		return nil
	}
	out := make([]record.Placeholder, 0, len(selected))
	for name, value := range selected {
		out = append(out, record.Placeholder{Name: name, Value: value})
	}
	sort.Slice(out, func(i, j int) bool { return out[i].Name < out[j].Name })
	return out
}

// Env converts an image's #ENV: contributions into record lines, sorted by key,
// with {prefix} left unsubstituted.
//
// Both records carry them. #ENV: lives in the header, which no digest sees, so
// without these lines editing GENOME_FASTA would move no key at all — and since
// #ENV: is how one recipe finds another's output, that silence would propagate
// into everything built on top.
//
// The ## notes are dropped: they are descriptive, and rewording one must not mint
// a new artifact. Values are kept as well as names, because a value change nearly
// always accompanies a body change that has already moved the digest, and one
// rule for both records is worth more than a distinction that never fires alone.
func Env(env []meta.EnvVar) []record.Env {
	if len(env) == 0 {
		return nil
	}
	out := make([]record.Env, 0, len(env))
	for _, e := range env {
		out = append(out, record.Env{Key: e.Key, Value: e.Value})
	}
	sort.Slice(out, func(i, j int) bool { return out[i].Key < out[j].Key })
	return out
}
