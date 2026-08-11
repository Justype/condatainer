package key

import (
	"fmt"
	"strings"
	"testing"

	"github.com/Justype/condatainer/catalog"
	"github.com/Justype/condatainer/internal/artifact/meta"
	"github.com/Justype/condatainer/internal/conda"
)

func TestScriptSchemesHaveCanonicalPreimages(t *testing.T) {
	recipe := []byte("#TARGET:data/{v}\n#PH:v:1,2\n#ENV:B={prefix}/b\n#ENV:A={v}\necho {v}\n")
	starID, starEq := Digest([]byte("star-id")), Digest([]byte("star-eq"))
	samID, samEq := Digest([]byte("sam-id")), Digest([]byte("sam-eq"))
	genomeID, genomeEq := Digest([]byte("genome-id")), Digest([]byte("genome-eq"))
	m := meta.Manifest{
		SchemaVersion: meta.SchemaVersion,
		Name:          "data/1/star/1",
		Type:          catalog.TypeData,
		BuildType:     "script",
		Platform:      meta.NativePlatform(),
		Source: meta.Source{
			Files:        []string{meta.RecipeFileName},
			Placeholders: map[string]string{"v": "1"},
		},
		Dependencies: []meta.Dependency{
			{Name: "star/1", Type: catalog.TypeApp, Identity: starID, Equiv: starEq, Role: meta.RoleApp},
			{Name: "samtools/1", Type: catalog.TypeApp, Identity: samID, Equiv: samEq, Role: meta.RoleHistory},
			{Name: "genome/1", Type: catalog.TypeData, Identity: genomeID, Equiv: genomeEq, Role: meta.RoleData},
		},
	}

	got, err := Generate(m, Sources{meta.RecipeFileName: recipe})
	if err != nil {
		t.Fatal(err)
	}
	wantIdentity := strings.Join([]string{
		"type=data",
		"env=A=1",
		"env=B={prefix}/b",
		"recipe=" + RecipeDigest(recipe),
		"ph=v=1",
		"dep=app samtools/1 " + samID,
		"dep=app star/1 " + starID,
		"dep=data genome/1 " + genomeID,
		"",
	}, "\n")
	wantEquiv := strings.Join([]string{
		"type=data",
		"env=A=1",
		"env=B={prefix}/b",
		"recipe=" + RecipeDigest(recipe),
		"ph=v=1",
		"dep=app star/1",
		"dep=data " + genomeEq,
		"",
	}, "\n")
	if string(got.Identity.Preimage) != wantIdentity {
		t.Errorf("identity preimage:\n%s\nwant:\n%s", got.Identity.Preimage, wantIdentity)
	}
	if string(got.Equiv.Preimage) != wantEquiv {
		t.Errorf("equivalence preimage:\n%s\nwant:\n%s", got.Equiv.Preimage, wantEquiv)
	}
	if got.Identity.Ref.Scheme != string(ScriptIdentityV1) ||
		got.Equiv.Ref.Scheme != string(ScriptEquivV1) {
		t.Errorf("schemes = %s / %s", got.Identity.Ref.Scheme, got.Equiv.Ref.Scheme)
	}
}

func TestDefinitionSchemesHaveCanonicalPreimages(t *testing.T) {
	recipe := []byte("Bootstrap: docker\nFrom: ubuntu:24.04\n")
	upstream := "sha256:" + strings.Repeat("a", 64)
	m := meta.Manifest{
		SchemaVersion: meta.SchemaVersion,
		Name:          "ubuntu24/base",
		Type:          catalog.TypeBase,
		BuildType:     "def",
		Platform:      meta.NativePlatform(),
		Source:        meta.Source{Files: []string{meta.RecipeFileName}},
		Build: meta.Build{From: &meta.From{
			Bootstrap: "docker", Ref: "ubuntu:24.04", Digest: upstream,
		}},
	}
	got, err := Generate(m, Sources{meta.RecipeFileName: recipe})
	if err != nil {
		t.Fatal(err)
	}
	wantIdentity := fmt.Sprintf("type=base\nrecipe=%s\nfrom=%s\n",
		RecipeDigest(recipe), upstream)
	wantEquiv := fmt.Sprintf("type=base\nrecipe=%s\n", RecipeDigest(recipe))
	if string(got.Identity.Preimage) != wantIdentity || string(got.Equiv.Preimage) != wantEquiv {
		t.Errorf("preimages:\n%s\n%s", got.Identity.Preimage, got.Equiv.Preimage)
	}
	if got.Identity.Ref.Scheme != string(DefinitionIdentityV1) ||
		got.Equiv.Ref.Scheme != string(DefinitionEquivV1) {
		t.Errorf("schemes = %s / %s", got.Identity.Ref.Scheme, got.Equiv.Ref.Scheme)
	}
}

func TestCondaSchemesHashStoredExportsDirectly(t *testing.T) {
	explicit := []byte("@EXPLICIT\nhttps://example.test/pkg-1.0-h0.conda\n")
	environment := []byte("channels:\n  - conda-forge\ndependencies:\n  - pkg=1.0\n")
	m := meta.Manifest{
		SchemaVersion: meta.SchemaVersion,
		Name:          "pkg/1.0",
		Type:          catalog.TypeApp,
		BuildType:     "conda",
		Platform:      meta.NativePlatform(),
		Source: meta.Source{Files: []string{
			conda.ExplicitFileName, conda.EnvironmentFileName,
		}},
	}
	got, err := Generate(m, Sources{
		conda.ExplicitFileName: explicit, conda.EnvironmentFileName: environment,
	})
	if err != nil {
		t.Fatal(err)
	}
	if string(got.Identity.Preimage) != string(explicit) || got.Identity.Ref.SHA256 != Sum(explicit) {
		t.Error("conda identity did not hash explicit.txt directly")
	}
	if string(got.Equiv.Preimage) != string(environment) || got.Equiv.Ref.SHA256 != Sum(environment) {
		t.Error("conda equivalence did not hash environment.yml directly")
	}
	if got.Identity.Ref.Scheme != string(CondaExplicitV1) ||
		got.Equiv.Ref.Scheme != string(CondaEnvironmentV1) {
		t.Errorf("schemes = %s / %s", got.Identity.Ref.Scheme, got.Equiv.Ref.Scheme)
	}
}

func TestVerifyRejectsOldAndUnknownSchemes(t *testing.T) {
	recipe := []byte("echo build\n")
	base := meta.Manifest{
		SchemaVersion: meta.SchemaVersion,
		Name:          "tool/1",
		Type:          catalog.TypeApp,
		BuildType:     "script",
		Platform:      meta.NativePlatform(),
		Source:        meta.Source{Files: []string{meta.RecipeFileName}},
	}
	sources := Sources{meta.RecipeFileName: recipe}

	old := base
	old.Keys = meta.Keys{
		Identity: meta.KeyRef{SHA256: strings.Repeat("a", 64)},
		Equiv:    meta.KeyRef{SHA256: strings.Repeat("b", 64)},
	}
	if _, err := Verify(old, sources); err == nil {
		t.Error("old file-backed manifest was accepted")
	}

	generated, err := Generate(base, sources)
	if err != nil {
		t.Fatal(err)
	}
	unknownIdentity := base
	unknownIdentity.Keys = generated.Keys()
	unknownIdentity.Keys.Identity.Scheme = "script-identity-v2"
	if _, err := Verify(unknownIdentity, sources); err == nil {
		t.Error("unknown identity scheme was accepted")
	}

	unknownEquiv := base
	unknownEquiv.Keys = generated.Keys()
	unknownEquiv.Keys.Equiv.Scheme = "script-equiv-v2"
	if _, err := Verify(unknownEquiv, sources); err == nil {
		t.Error("unknown equivalence scheme was accepted")
	}

	mismatchedPair := base
	mismatchedPair.Keys = generated.Keys()
	mismatchedPair.Keys.Equiv.Scheme = string(DefinitionEquivV1)
	if _, err := Verify(mismatchedPair, sources); err == nil {
		t.Error("mismatched known scheme pair was accepted")
	}
}
