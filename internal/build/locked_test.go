package build

import (
	"context"
	"errors"
	"os"
	"path/filepath"
	"strings"
	"testing"

	"github.com/Justype/condatainer/catalog"
	"github.com/Justype/condatainer/internal/artifact/meta"
	"github.com/Justype/condatainer/internal/conda"
)

const lockedRecipe = `#TYPE: app
#DESC: an aligner
#URL: https://example.invalid/star
#ENV:PATH={prefix}/bin
#SBATCH --cpus-per-task=4
echo build
`

func lockedManifest(name string) meta.Manifest {
	return meta.Manifest{
		SchemaVersion: meta.SchemaVersion,
		Name:          name,
		Type:          catalog.TypeApp,
		BuildType:     "script",
		Platform:      meta.NativePlatform(),
		Source:        meta.Source{Files: []string{meta.RecipeFileName}},
		Keys: meta.Keys{
			Identity: meta.KeyRef{Scheme: "script.identity.v1", SHA256: strings.Repeat("a", 64)},
			Equiv:    meta.KeyRef{Scheme: "script.equiv.v1", SHA256: strings.Repeat("b", 64)},
		},
	}
}

// lockedSpec is a minimal valid rebuild of one script artifact into tmp.
func lockedSpec(t *testing.T, recipe string) LockedSpec {
	t.Helper()
	work := t.TempDir()
	return LockedSpec{
		Manifest: lockedManifest("star/2.7.11b"),
		Sources:  map[string][]byte{meta.RecipeFileName: []byte(recipe)},
		Output:   filepath.Join(work, "out.sqf"),
		TmpRoot:  filepath.Join(work, "tmp"),
	}
}

func TestNewLockedObjectBuildsFromTheVendoredRecipe(t *testing.T) {
	spec := lockedSpec(t, lockedRecipe)
	object, err := NewLockedObject(context.Background(), spec)
	if err != nil {
		t.Fatal(err)
	}

	if object.BuildType() != BuildTypeScript {
		t.Errorf("build type = %v, want script", object.BuildType())
	}
	if object.TargetOverlayPath() != spec.Output {
		t.Errorf("target = %q, want the caller's output %q", object.TargetOverlayPath(), spec.Output)
	}
	if !object.locked {
		t.Error("object is not marked locked")
	}

	// The Spec comes from the recipe, so the rebuilt manifest can match.
	got := object.Spec()
	if got.Image.Description != "an aligner" || got.Image.URL != "https://example.invalid/star" {
		t.Errorf("descriptive metadata = %#v", got.Image)
	}
	if len(got.Image.Env) != 1 || got.Image.Env[0].Key != "PATH" {
		t.Errorf("env = %#v, want the recipe's #ENV:", got.Image.Env)
	}
	if got.Image.Prefix != meta.Prefix("star/2.7.11b", catalog.TypeApp) {
		t.Errorf("prefix = %q", got.Image.Prefix)
	}
	if object.ScriptSpecs() == nil || object.ScriptSpecs().Spec.CpusPerTask != 4 {
		t.Errorf("scheduler directives were not read: %#v", object.ScriptSpecs())
	}

	// The runnable copy is written where the build will read it.
	if data, err := os.ReadFile(object.BuildSource()); err != nil {
		t.Fatal(err)
	} else if string(data) != lockedRecipe {
		t.Errorf("runnable recipe = %q", data)
	}
}

// The recipe reaches the image byte for byte, which is what lets its identity be
// regenerated from the checkout.
func TestNewLockedObjectEmbedsTheRecipeVerbatim(t *testing.T) {
	spec := lockedSpec(t, lockedRecipe)
	object, err := NewLockedObject(context.Background(), spec)
	if err != nil {
		t.Fatal(err)
	}
	sources := object.keySources()
	if got := string(sources[meta.RecipeFileName]); got != lockedRecipe {
		t.Fatalf("embedded recipe = %q, want it verbatim", got)
	}
}

// A template runs expanded and embeds the template, exactly as a catalog build
// does. The placeholders come from the manifest, which is the only record of
// which variant this is.
func TestNewLockedObjectExpandsATemplateWithTheRecordedPlaceholders(t *testing.T) {
	template := "#TYPE: data\n#PH: release: 49,50\n#TARGET: gencode/{release}\necho {release}\n"
	spec := lockedSpec(t, template)
	spec.Manifest = lockedManifest("gencode/49")
	spec.Manifest.Type = catalog.TypeData
	spec.Manifest.Source.Placeholders = map[string]string{"release": "49"}
	spec.Manifest.Source.TargetTemplate = "gencode/{release}"
	spec.Sources[meta.RecipeFileName] = []byte(template)

	object, err := NewLockedObject(context.Background(), spec)
	if err != nil {
		t.Fatal(err)
	}
	runnable, err := os.ReadFile(object.BuildSource())
	if err != nil {
		t.Fatal(err)
	}
	if !strings.Contains(string(runnable), "echo 49") {
		t.Errorf("runnable recipe was not expanded: %q", runnable)
	}
	if got := string(object.keySources()[meta.RecipeFileName]); !strings.Contains(got, "echo {release}") {
		t.Errorf("embedded recipe = %q, want the template with its tokens", got)
	}
	if object.Spec().Source.Placeholders["release"] != "49" {
		t.Errorf("placeholders = %#v, want the manifest's", object.Spec().Source.Placeholders)
	}
}

// Dependencies are mounted at the paths the caller supplies. Nothing is looked
// up by name, which is what makes the rebuild reproducible.
func TestNewLockedObjectMountsDependenciesByPath(t *testing.T) {
	work := t.TempDir()
	depPath := filepath.Join(work, "grch38--genome@0123456789ab.sqf")
	if err := os.WriteFile(depPath, []byte("payload"), 0o644); err != nil {
		t.Fatal(err)
	}
	recipe := "#TYPE: data\n#DEP: grch38/genome/gencode49\necho build\n"

	spec := lockedSpec(t, recipe)
	spec.Manifest.Type = catalog.TypeData
	spec.Deps = []LockedDep{{
		Name:     "grch38/genome/gencode49",
		Identity: meta.KeyRef{Scheme: "script.identity.v1", SHA256: strings.Repeat("c", 64)},
		Path:     depPath,
	}}

	object, err := NewLockedObject(context.Background(), spec)
	if err != nil {
		t.Fatal(err)
	}
	deps := object.Dependencies()
	if len(deps) != 1 || deps[0] != depPath {
		t.Fatalf("dependencies = %v, want the supplied path", deps)
	}
	missing, err := object.GetMissingDependencies()
	if err != nil {
		t.Fatal(err)
	}
	if len(missing) != 0 {
		t.Errorf("missing = %v, want none: the file is there", missing)
	}
}

// A recipe declaring dependencies the caller did not supply means the vendored
// records disagree. Refusing here beats a rebuild that mounts nothing and then
// fails on a key mismatch nobody can explain.
func TestNewLockedObjectRefusesADependencyCountMismatch(t *testing.T) {
	spec := lockedSpec(t, "#TYPE: data\n#DEP: grch38/genome/gencode49\necho build\n")
	spec.Manifest.Type = catalog.TypeData

	_, err := NewLockedObject(context.Background(), spec)
	if !errors.Is(err, ErrLockedInvalid) {
		t.Fatalf("err = %v, want ErrLockedInvalid", err)
	}
	if !strings.Contains(err.Error(), "1 dependencies but 0 paths") {
		t.Errorf("err = %v, want it to name both counts", err)
	}
}

func TestNewLockedObjectRefusesAMissingDependencyImage(t *testing.T) {
	spec := lockedSpec(t, "#TYPE: data\n#DEP: grch38/genome/gencode49\necho build\n")
	spec.Manifest.Type = catalog.TypeData
	spec.Deps = []LockedDep{{Name: "grch38/genome/gencode49", Path: filepath.Join(t.TempDir(), "gone.sqf")}}

	if _, err := NewLockedObject(context.Background(), spec); !errors.Is(err, ErrLockedInvalid) {
		t.Fatalf("err = %v, want ErrLockedInvalid", err)
	}
}

func TestNewLockedObjectRefusesARelativeDependencyPath(t *testing.T) {
	spec := lockedSpec(t, "#TYPE: data\n#DEP: grch38/genome/gencode49\necho build\n")
	spec.Manifest.Type = catalog.TypeData
	spec.Deps = []LockedDep{{Name: "grch38/genome/gencode49", Path: "overlays/genome.sqf"}}

	if _, err := NewLockedObject(context.Background(), spec); !errors.Is(err, ErrLockedInvalid) {
		t.Fatalf("err = %v, want ErrLockedInvalid", err)
	}
}

// Nothing could check a rebuild that has no recorded keys to check it against.
func TestNewLockedObjectRefusesAKeylessManifest(t *testing.T) {
	spec := lockedSpec(t, lockedRecipe)
	spec.Manifest.Keys = meta.Keys{}

	err := func() error { _, err := NewLockedObject(context.Background(), spec); return err }()
	if !errors.Is(err, ErrLockedInvalid) {
		t.Fatalf("err = %v, want ErrLockedInvalid", err)
	}
}

func TestNewLockedObjectRefusesAnOccupiedOutput(t *testing.T) {
	spec := lockedSpec(t, lockedRecipe)
	if err := os.WriteFile(spec.Output, []byte("already here"), 0o644); err != nil {
		t.Fatal(err)
	}
	if _, err := NewLockedObject(context.Background(), spec); !errors.Is(err, ErrLockedInvalid) {
		t.Fatalf("err = %v, want ErrLockedInvalid", err)
	}
}

func TestNewLockedObjectRefusesAnUnvendoredRecipe(t *testing.T) {
	spec := lockedSpec(t, lockedRecipe)
	spec.Sources = map[string][]byte{}

	if _, err := NewLockedObject(context.Background(), spec); !errors.Is(err, ErrLockedInvalid) {
		t.Fatalf("err = %v, want ErrLockedInvalid", err)
	}
}

// #INPUT: answers are supplied per invocation because they are never recorded.
// A count that does not match the prompts would feed answers to the wrong reads.
func TestNewLockedObjectChecksInputAnswersAgainstThePrompts(t *testing.T) {
	recipe := "#TYPE: app\n#INPUT: which release?\n#INPUT: which build?\nread -r a\n"
	spec := lockedSpec(t, recipe)
	spec.Manifest.Source.RequiresInput = true

	if _, err := NewLockedObject(context.Background(), spec); !errors.Is(err, ErrLockedInvalid) {
		t.Fatalf("err = %v, want a refusal when no answers were supplied", err)
	}

	spec = lockedSpec(t, recipe)
	spec.Manifest.Source.RequiresInput = true
	spec.Answers = []string{"49", "primary"}
	object, err := NewLockedObject(context.Background(), spec)
	if err != nil {
		t.Fatal(err)
	}
	if got := object.InputAnswers(); len(got) != 2 || got[0] != "49" {
		t.Errorf("answers = %v", got)
	}
}

// A Conda rebuild replays explicit.txt, whose bytes are the recorded identity.
// environment.yml would be a fresh solve against whatever the channels serve now.
func TestNewLockedObjectReplaysTheCondaExplicitExport(t *testing.T) {
	work := t.TempDir()
	spec := LockedSpec{
		Manifest: lockedManifest("myenv/1.0"),
		Sources: map[string][]byte{
			conda.ExplicitFileName:    []byte("@EXPLICIT\nhttps://example.invalid/star-2.7.11b.conda\n"),
			conda.EnvironmentFileName: []byte("name: myenv\n"),
		},
		Output:  filepath.Join(work, "out.sqf"),
		TmpRoot: filepath.Join(work, "tmp"),
	}
	spec.Manifest.BuildType = "conda"
	spec.Manifest.Build.Channels = []string{"conda-forge", "bioconda"}
	spec.Manifest.Source.Files = []string{conda.EnvironmentFileName, conda.ExplicitFileName}

	object, err := NewLockedObject(context.Background(), spec)
	if err != nil {
		t.Fatal(err)
	}
	if object.BuildType() != BuildTypeConda {
		t.Fatalf("build type = %v, want conda", object.BuildType())
	}
	if filepath.Base(object.BuildSource()) != conda.ExplicitFileName {
		t.Errorf("build source = %q, want the explicit export", object.BuildSource())
	}
	source := object.Spec().Source.Conda
	if source == nil {
		t.Fatal("no conda source")
	}
	// Recorded channels, not live config: the manifest is a projection of Spec.
	if len(source.Channels) != 2 || source.Channels[0] != "conda-forge" {
		t.Errorf("channels = %v, want the manifest's", source.Channels)
	}
}

func TestNewLockedObjectRefusesAnUnvendoredCondaExplicitExport(t *testing.T) {
	work := t.TempDir()
	spec := LockedSpec{
		Manifest: lockedManifest("myenv/1.0"),
		Sources:  map[string][]byte{conda.EnvironmentFileName: []byte("name: myenv\n")},
		Output:   filepath.Join(work, "out.sqf"),
		TmpRoot:  filepath.Join(work, "tmp"),
	}
	spec.Manifest.BuildType = "conda"

	if _, err := NewLockedObject(context.Background(), spec); !errors.Is(err, ErrLockedInvalid) {
		t.Fatalf("err = %v, want ErrLockedInvalid", err)
	}
}

func TestNewLockedObjectRefusesAnUnknownBuildType(t *testing.T) {
	spec := lockedSpec(t, lockedRecipe)
	spec.Manifest.BuildType = "cmake"

	if _, err := NewLockedObject(context.Background(), spec); !errors.Is(err, ErrLockedInvalid) {
		t.Fatalf("err = %v, want ErrLockedInvalid", err)
	}
}

// A definition rebuild bootstraps from the digest the lock recorded. Re-resolving
// the reference would follow a tag that has since moved.
func TestNewLockedObjectKeepsTheRecordedUpstreamDigest(t *testing.T) {
	definition := "Bootstrap: docker\nFrom: ubuntu:24.04\n\n%post\n  echo build\n"
	spec := lockedSpec(t, definition)
	spec.Manifest.BuildType = "def"
	spec.Manifest.Type = catalog.TypeOS
	spec.Manifest.Build.From = &meta.From{
		Bootstrap: "docker", Ref: "ubuntu:24.04", Digest: "sha256:" + strings.Repeat("d", 64),
	}
	spec.Sources[meta.RecipeFileName] = []byte(definition)

	object, err := NewLockedObject(context.Background(), spec)
	if err != nil {
		t.Fatal(err)
	}
	if object.BuildType() != BuildTypeDef {
		t.Fatalf("build type = %v, want def", object.BuildType())
	}
	if got := object.Spec().Source.UpstreamDigest(); got != spec.Manifest.Build.From.Digest {
		t.Errorf("upstream digest = %q, want the recorded one", got)
	}

	// resolveUpstream must leave it alone rather than reaching for a registry.
	object.resolveUpstream(context.Background(), []byte(definition))
	if got := object.Spec().Source.UpstreamDigest(); got != spec.Manifest.Build.From.Digest {
		t.Errorf("upstream digest = %q after resolution, want it untouched", got)
	}

	// Re-siting as a .sif build must not leave the recipe behind at the old path.
	if _, err := os.Stat(object.BuildSource()); err != nil {
		t.Errorf("definition source is not where the build will read it: %v", err)
	}
}

// CondaSource selects which vendored export is replayed. environment.yml is the
// second tier restore reaches for when explicit.txt's package URLs have been
// pruned from the channel.
func TestNewLockedObjectReplaysTheNamedCondaExport(t *testing.T) {
	condaSpec := func(t *testing.T, which string) LockedSpec {
		t.Helper()
		work := t.TempDir()
		spec := LockedSpec{
			Manifest: lockedManifest("myenv/1.0"),
			Sources: map[string][]byte{
				conda.ExplicitFileName:    []byte("@EXPLICIT\nhttps://example.invalid/star.conda\n"),
				conda.EnvironmentFileName: []byte("dependencies:\n  - python=3.12.7\n"),
			},
			Output:      filepath.Join(work, "out.sqf"),
			TmpRoot:     filepath.Join(work, "tmp"),
			CondaSource: which,
		}
		spec.Manifest.BuildType = "conda"
		spec.Manifest.Source.Files = []string{conda.EnvironmentFileName, conda.ExplicitFileName}
		return spec
	}

	for _, tc := range []struct{ name, which, want string }{
		{"DefaultsToExplicit", "", conda.ExplicitFileName},
		{"Explicit", conda.ExplicitFileName, conda.ExplicitFileName},
		{"Environment", conda.EnvironmentFileName, conda.EnvironmentFileName},
	} {
		t.Run(tc.name, func(t *testing.T) {
			object, err := NewLockedObject(context.Background(), condaSpec(t, tc.which))
			if err != nil {
				t.Fatal(err)
			}
			if got := filepath.Base(object.BuildSource()); got != tc.want {
				t.Fatalf("build source = %q, want %q", got, tc.want)
			}
		})
	}

	// An export the manifest does not vendor cannot be replayed, and a name that
	// is no Conda export at all is a caller bug rather than a build failure.
	missing := condaSpec(t, conda.EnvironmentFileName)
	delete(missing.Sources, conda.EnvironmentFileName)
	if _, err := NewLockedObject(context.Background(), missing); !errors.Is(err, ErrLockedInvalid) {
		t.Errorf("unvendored export: err = %v, want ErrLockedInvalid", err)
	}
	if _, err := NewLockedObject(context.Background(), condaSpec(t, meta.RecipeFileName)); !errors.Is(err, ErrLockedInvalid) {
		t.Errorf("non-export source: err = %v, want ErrLockedInvalid", err)
	}
}
