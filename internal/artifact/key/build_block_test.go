package key

import (
	"strings"
	"testing"
	"time"

	"github.com/Justype/condatainer/catalog"
	"github.com/Justype/condatainer/internal/artifact/meta"
	"github.com/Justype/condatainer/internal/conda"
)

// Build.Source and Build.Created are diagnostic provenance, and no scheme may
// hash them. Build.Source especially: two labs building one recipe from their
// own mirrors must produce one identity, or nothing deduplicates and a lock pins
// a repository instead of a recipe.
//
// Build.From is deliberately absent from this test — a definition's upstream
// digest *does* reach identity, which is the distinction being guarded.
func TestBuildProvenanceStaysOutOfEveryKey(t *testing.T) {
	recipe := []byte("#DESC:x\necho build\n")
	definition := []byte("Bootstrap: docker\nFrom: ubuntu:24.04\n")
	explicit := []byte("@EXPLICIT\nhttps://example.test/pkg-1.0-h0.conda\n")
	environment := []byte("channels:\n  - conda-forge\ndependencies:\n  - pkg=1.0\n")

	cases := []struct {
		name     string
		manifest meta.Manifest
		sources  Sources
	}{
		{
			name: "script",
			manifest: meta.Manifest{
				SchemaVersion: meta.SchemaVersion,
				Name:          "tool/1",
				Type:          catalog.TypeApp,
				BuildType:     "script",
				Platform:      meta.NativePlatform(),
				Source:        meta.Source{Files: []string{meta.RecipeFileName}},
			},
			sources: Sources{meta.RecipeFileName: recipe},
		},
		{
			name: "script with dependencies",
			manifest: meta.Manifest{
				SchemaVersion: meta.SchemaVersion,
				Name:          "grch38/index/1",
				Type:          catalog.TypeData,
				BuildType:     "script",
				Platform:      meta.NativePlatform(),
				Source:        meta.Source{Files: []string{meta.RecipeFileName}},
				Dependencies: []meta.Dependency{
					{Name: "star/1", Type: catalog.TypeApp, Role: meta.RoleApp,
						Identity: idKey("star id"), Equiv: eqKey("star eq")},
				},
			},
			sources: Sources{meta.RecipeFileName: recipe},
		},
		{
			name: "definition",
			manifest: meta.Manifest{
				SchemaVersion: meta.SchemaVersion,
				Name:          "ubuntu24/base",
				Type:          catalog.TypeBase,
				BuildType:     "def",
				Platform:      meta.NativePlatform(),
				Source:        meta.Source{Files: []string{meta.RecipeFileName}},
				Build: meta.Build{From: &meta.From{
					Bootstrap: "docker", Ref: "ubuntu:24.04",
					Digest: DigestPrefix + strings.Repeat("a", 64),
				}},
			},
			sources: Sources{meta.RecipeFileName: definition},
		},
		{
			name: "conda",
			manifest: meta.Manifest{
				SchemaVersion: meta.SchemaVersion,
				Name:          "pkg/1.0",
				Type:          catalog.TypeApp,
				BuildType:     "conda",
				Platform:      meta.NativePlatform(),
				Source: meta.Source{Files: []string{
					conda.ExplicitFileName, conda.EnvironmentFileName,
				}},
			},
			sources: Sources{
				conda.ExplicitFileName:    explicit,
				conda.EnvironmentFileName: environment,
			},
		},
	}

	for _, tc := range cases {
		t.Run(tc.name, func(t *testing.T) {
			plain, err := Generate(tc.manifest, tc.sources)
			if err != nil {
				t.Fatalf("Generate: %v", err)
			}

			stamped := tc.manifest
			stamped.Build.Source = "https://github.com/Justype/cnt-scripts"
			stamped.Build.Created = time.Date(2026, 7, 21, 9, 30, 0, 0, time.UTC)
			withProvenance, err := Generate(stamped, tc.sources)
			if err != nil {
				t.Fatalf("Generate with provenance: %v", err)
			}

			if plain.Identity.Ref != withProvenance.Identity.Ref {
				t.Errorf("build.source/created moved identity: %+v -> %+v",
					plain.Identity.Ref, withProvenance.Identity.Ref)
			}
			if plain.Equiv.Ref != withProvenance.Equiv.Ref {
				t.Errorf("build.source/created moved equivalence: %+v -> %+v",
					plain.Equiv.Ref, withProvenance.Equiv.Ref)
			}

			// A second build minutes later, from a different mirror, is the same
			// artifact — the case the store depends on to deduplicate at all.
			mirrored := tc.manifest
			mirrored.Build.Source = "https://git.labA.edu/mirror/cnt-scripts"
			mirrored.Build.Created = time.Date(2026, 7, 21, 9, 41, 0, 0, time.UTC)
			other, err := Generate(mirrored, tc.sources)
			if err != nil {
				t.Fatalf("Generate from mirror: %v", err)
			}
			if other.Identity.Ref != withProvenance.Identity.Ref {
				t.Error("one recipe from two collections produced two identities")
			}
		})
	}
}

// A dependency edge carries a complete key, so the scheme is part of the
// preimage on its own. Without this, the change from a bare SHA to a KeyRef
// would be invisible and the two would be interchangeable.
func TestDependencySchemeReachesTheKey(t *testing.T) {
	recipe := []byte("#DESC:x\necho build\n")
	base := meta.Manifest{
		SchemaVersion: meta.SchemaVersion,
		Name:          "grch38/index/1",
		Type:          catalog.TypeData,
		BuildType:     "script",
		Platform:      meta.NativePlatform(),
		Source:        meta.Source{Files: []string{meta.RecipeFileName}},
		Dependencies: []meta.Dependency{
			{Name: "genome/1", Type: catalog.TypeData, Role: meta.RoleData,
				Identity: idKey("genome id"), Equiv: eqKey("genome eq")},
		},
	}
	sources := Sources{meta.RecipeFileName: recipe}

	first, err := Generate(base, sources)
	if err != nil {
		t.Fatal(err)
	}

	// Same digests, different schemes. Nothing else moves.
	reschemed := base
	reschemed.Dependencies = []meta.Dependency{{
		Name: "genome/1", Type: catalog.TypeData, Role: meta.RoleData,
		Identity: meta.KeyRef{Scheme: string(CondaExplicitV1), SHA256: idKey("genome id").SHA256},
		Equiv:    meta.KeyRef{Scheme: string(CondaEnvironmentV1), SHA256: eqKey("genome eq").SHA256},
	}}
	second, err := Generate(reschemed, sources)
	if err != nil {
		t.Fatal(err)
	}

	if first.Identity.Ref == second.Identity.Ref {
		t.Error("a dependency's identity scheme did not reach the identity key")
	}
	if first.Equiv.Ref == second.Equiv.Ref {
		t.Error("a dependency's equivalence scheme did not reach the equivalence key")
	}
}
