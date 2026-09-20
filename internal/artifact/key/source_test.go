package key

import (
	"strings"
	"testing"

	"github.com/Justype/condatainer/catalog"
	"github.com/Justype/condatainer/internal/artifact/meta"
)

// srcFile builds a fetched-source record whose digest is the hash of text, so a
// test can say "different bytes" without writing hex.
func srcFile(name, text string) meta.SourceFile {
	return meta.SourceFile{Name: name, SHA256: Sum([]byte(text))}
}

func TestSourcesEnterIdentitySortedByName(t *testing.T) {
	gtf, genome := srcFile("gtf", "annotation"), srcFile("genome", "assembly")
	a := Artifact{
		Type:    catalog.TypeData,
		Recipe:  []byte("echo build\n"),
		Fetched: []meta.SourceFile{gtf, genome},
	}

	_, model, err := deriveScriptIdentityV1(a)
	if err != nil {
		t.Fatal(err)
	}
	preimage, err := Marshal(model)
	if err != nil {
		t.Fatal(err)
	}

	want := strings.Join([]string{
		"type=data",
		"recipe=" + RecipeDigest(a.Recipe),
		"src=genome=" + genome.Digest(),
		"src=gtf=" + gtf.Digest(),
		"",
	}, "\n")
	if got := string(preimage); got != want {
		t.Errorf("preimage =\n%s\nwant\n%s", got, want)
	}

	// Declaration order must not move the key.
	a.Fetched = []meta.SourceFile{genome, gtf}
	reordered, _, err := deriveScriptIdentityV1(a)
	if err != nil {
		t.Fatal(err)
	}
	first, _, _ := deriveScriptIdentityV1(Artifact{
		Type: catalog.TypeData, Recipe: a.Recipe,
		Fetched: []meta.SourceFile{gtf, genome},
	})
	if reordered.Ref.SHA256 != first.Ref.SHA256 {
		t.Error("declaration order changed the identity")
	}
}

// A re-cut upstream file is the case this whole scheme exists for: same recipe,
// same placeholders, different bytes.
func TestRecutSourceMovesIdentityButNotEquivalence(t *testing.T) {
	base := Artifact{
		Type:    catalog.TypeData,
		Recipe:  []byte("echo build\n"),
		Fetched: []meta.SourceFile{srcFile("gtf", "release 49")},
	}
	recut := base
	recut.Fetched = []meta.SourceFile{srcFile("gtf", "release 49, corrected")}

	beforeID, _, err := deriveScriptIdentityV1(base)
	if err != nil {
		t.Fatal(err)
	}
	afterID, _, err := deriveScriptIdentityV1(recut)
	if err != nil {
		t.Fatal(err)
	}
	if beforeID.Ref.SHA256 == afterID.Ref.SHA256 {
		t.Error("a re-cut source left the identity unchanged")
	}

	beforeEq, _, err := deriveScriptEquivV1(base)
	if err != nil {
		t.Fatal(err)
	}
	afterEq, _, err := deriveScriptEquivV1(recut)
	if err != nil {
		t.Fatal(err)
	}
	if beforeEq.Ref.SHA256 != afterEq.Ref.SHA256 {
		t.Error("a re-cut source moved equivalence; an installed index would stop substituting")
	}
}

func TestSourceModelRejectsBadInput(t *testing.T) {
	good := SourceValue{Name: "gtf", Digest: Digest([]byte("x"))}
	for _, tc := range []struct {
		name  string
		model Model
		want  string
	}{
		{
			name:  "equivalence carries no sources",
			model: Model{Kind: KindEquiv, Type: catalog.TypeData, Sources: []SourceValue{good}},
			want:  "identity only",
		},
		{
			name: "digest must be sha256",
			model: Model{Kind: KindIdentity, Type: catalog.TypeData,
				Sources: []SourceValue{{Name: "gtf", Digest: "deadbeef"}}},
			want: "not a sha256 digest",
		},
		{
			name: "name must be a token",
			model: Model{Kind: KindIdentity, Type: catalog.TypeData,
				Sources: []SourceValue{{Name: "two words", Digest: good.Digest}}},
			want: "not a usable source name",
		},
		{
			name: "one name once",
			model: Model{Kind: KindIdentity, Type: catalog.TypeData,
				Sources: []SourceValue{good, good}},
			want: "appears more than once",
		},
	} {
		t.Run(tc.name, func(t *testing.T) {
			_, err := Marshal(tc.model)
			if err == nil {
				t.Fatal("expected a rejection")
			}
			if !strings.Contains(err.Error(), tc.want) {
				t.Errorf("error = %v, want it to mention %q", err, tc.want)
			}
		})
	}
}
