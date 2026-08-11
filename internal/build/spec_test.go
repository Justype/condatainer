package build

import (
	"context"
	"errors"
	"os"
	"os/exec"
	"path/filepath"
	"slices"
	"strings"
	"testing"

	"github.com/Justype/condatainer/catalog"
	"github.com/Justype/condatainer/internal/artifact/capsule"
	"github.com/Justype/condatainer/internal/artifact/key"
	"github.com/Justype/condatainer/internal/artifact/meta"
	"github.com/Justype/condatainer/internal/conda"
	"github.com/Justype/condatainer/internal/config"
	"github.com/Justype/condatainer/internal/runtime/container"
)

func TestSourceSpecBuildTypeIsDerived(t *testing.T) {
	tests := []struct {
		name   string
		source SourceSpec
		want   BuildType
	}{
		{"script", SourceSpec{Script: &ScriptSource{}}, BuildTypeScript},
		{"definition", SourceSpec{Definition: &DefinitionSource{}}, BuildTypeDef},
		{"conda", SourceSpec{Conda: &CondaSource{}}, BuildTypeConda},
		{"unset", SourceSpec{}, 0},
	}
	for _, tt := range tests {
		t.Run(tt.name, func(t *testing.T) {
			if got := tt.source.BuildType(); got != tt.want {
				t.Errorf("BuildType() = %v, want %v", got, tt.want)
			}
		})
	}
}

func TestSourceSpecFile(t *testing.T) {
	script := SourceSpec{Script: &ScriptSource{File: SourceFile{Name: "a.sh", Data: []byte("x")}}}
	if f, ok := script.File(); !ok || f.Name != "a.sh" {
		t.Errorf("script file = %+v, %v", f, ok)
	}

	// A single package or a package list has no file to materialize.
	pkg := SourceSpec{Conda: &CondaSource{Package: &CondaPackage{Name: "samtools"}}}
	if _, ok := pkg.File(); ok {
		t.Error("a bare conda package reported a source file")
	}
	list := SourceSpec{Conda: &CondaSource{Packages: []string{"a", "b"}}}
	if _, ok := list.File(); ok {
		t.Error("a conda package list reported a source file")
	}
	yaml := SourceSpec{Conda: &CondaSource{File: &SourceFile{Name: "env.yml"}}}
	if f, ok := yaml.File(); !ok || f.Name != "env.yml" {
		t.Errorf("conda file = %+v, %v", f, ok)
	}
}

func TestSpecManifest(t *testing.T) {
	spec := Spec{
		Image: ImageSpec{
			Name:        "samtools/1.23.1",
			Type:        catalog.TypeApp,
			Description: "SAMtools alignment toolkit",
			URL:         "https://www.htslib.org/",
		},
		Source: SourceSpec{Script: &ScriptSource{}},
	}
	spec.Image.Prefix = "/cnt/samtools/1.23.1"
	spec.Image.Env = []meta.EnvVar{{Key: "SAMTOOLS_DIR", Value: "{prefix}", Note: "install root"}}

	m := spec.Manifest()
	if m.SchemaVersion != meta.SchemaVersion {
		t.Errorf("schema = %d", m.SchemaVersion)
	}
	if m.Name != "samtools/1.23.1" || m.Type != catalog.TypeApp {
		t.Errorf("identity = %q/%q", m.Name, m.Type)
	}
	if m.BuildType != "script" {
		t.Errorf("build_type = %q, want script (derived from the source)", m.BuildType)
	}
	if m.Description != spec.Image.Description || m.URL != spec.Image.URL {
		t.Errorf("descriptive metadata lost: %q / %q", m.Description, m.URL)
	}
	if err := meta.ValidateManifest(m); err != nil {
		t.Errorf("rendered manifest does not validate: %v", err)
	}
	if m.Platform.Arch == "" {
		t.Error("the manifest records no architecture")
	}

	// The mount-time contract comes off the same Spec, and it is the only place
	// the prefix and the environment appear.
	rt := spec.Runtime()
	if err := meta.ValidateRuntime(rt); err != nil {
		t.Errorf("rendered runtime does not validate: %v", err)
	}
	if rt.Prefix != spec.Image.Prefix {
		t.Errorf("runtime prefix = %q, want %q", rt.Prefix, spec.Image.Prefix)
	}
	if len(rt.Env) != 1 || rt.Env[0].Value != "{prefix}" {
		t.Errorf("runtime env = %+v, want {prefix} intact", rt.Env)
	}
	if rt.Description != spec.Image.Description {
		t.Errorf("runtime description = %q", rt.Description)
	}
}

// Both documents are projections of Spec, and Spec has nowhere to put an #INPUT:
// answer. This pins that: answers are execution input and must never be
// embedded, because they routinely carry tokens and licence keys.
func TestSpecMetadataCannotCarryAnswers(t *testing.T) {
	const secret = "s3cret-license-key"

	spec := Spec{
		Image:  ImageSpec{Name: "tool/1.0", Type: catalog.TypeApp, Prefix: "/cnt/tool/1.0"},
		Source: SourceSpec{Script: &ScriptSource{Prompts: []string{"License key:"}}},
	}
	manifest, err := meta.MarshalManifest(spec.Manifest())
	if err != nil {
		t.Fatalf("MarshalManifest: %v", err)
	}
	runtime, err := meta.MarshalRuntime(spec.Runtime())
	if err != nil {
		t.Fatalf("MarshalRuntime: %v", err)
	}
	for _, data := range [][]byte{manifest, runtime} {
		if strings.Contains(string(data), secret) {
			t.Fatalf("an answer reached embedded metadata:\n%s", data)
		}
		// The prompt is a declaration, not an answer, and is also not embedded.
		if strings.Contains(string(data), "License key:") {
			t.Errorf("a prompt declaration reached embedded metadata:\n%s", data)
		}
	}
}

func TestEnvFromRecipeKeepsPrefixToken(t *testing.T) {
	recipe, err := catalog.ParseRecipe("samtools/1.23.1", strings.NewReader(
		"#DESC:SAMtools\n"+
			"#ENV:SAMTOOLS_DIR={prefix}  ## install root\n"+
			"#ENV:PATH_EXTRA={prefix}/bin\n",
	))
	if err != nil {
		t.Fatalf("ParseRecipe: %v", err)
	}

	env := envFromRecipe(recipe.Env)
	if len(env) != 2 {
		t.Fatalf("env = %+v", env)
	}
	// {prefix} must survive into the image: the install prefix is not known when
	// the image is built, only when it is loaded.
	if env[0].Value != "{prefix}" {
		t.Errorf("env[0] = %q, want the token intact", env[0].Value)
	}
	if env[0].Note != "install root" {
		t.Errorf("note = %q", env[0].Note)
	}
	if got := env[1].Resolved(meta.Prefix("samtools/1.23.1", catalog.TypeApp)); got != "/cnt/samtools/1.23.1/bin" {
		t.Errorf("resolved = %q", got)
	}
}

func TestTargetFor(t *testing.T) {
	target := targetFor("/images/samtools--1.23.1.sqf")
	if target.Path != "/images/samtools--1.23.1.sqf" {
		t.Errorf("Path = %q", target.Path)
	}
	if target.Lock != target.Path+".lock" {
		t.Errorf("Lock = %q", target.Lock)
	}
	if target.Prepared != "" {
		t.Errorf("Prepared = %q, want empty until the build claims one", target.Prepared)
	}
}

func TestSourceFileName(t *testing.T) {
	if got := sourceFileName("grch38/genome/gencode", false); got != "cnt--grch38--genome--gencode.sh" {
		t.Errorf("script name = %q", got)
	}
	if got := sourceFileName("ubuntu24/base", true); got != "cnt--ubuntu24--base.def" {
		t.Errorf("def name = %q", got)
	}
}

// Resolution has to produce a Spec that renders valid metadata, from a real
// catalog lookup rather than a hand-built struct. This is what proves the
// descriptive metadata mapping — #DESC: to Description, #URL: to URL, #ENV: to
// the runtime env — actually runs end to end.
func TestResolvedSpecRendersValidManifest(t *testing.T) {
	setTestSource(t, writeRecipe(t, "samtools/1.23.1", strings.Join([]string{
		"#!/usr/bin/env bash",
		"#DESC:SAMtools alignment toolkit",
		"#URL:https://www.htslib.org/",
		"#ENV:SAMTOOLS_DIR={prefix}  ## install root",
		"#ENV:PATH_EXTRA={prefix}/bin",
		"echo build",
		"",
	}, "\n")))

	imagesDir := t.TempDir()
	obj, err := NewBuildObject(context.Background(), "samtools/1.23.1", false, imagesDir, false)
	if err != nil {
		t.Fatalf("NewBuildObject: %v", err)
	}

	spec := obj.Spec()
	if spec.Image.Name != "samtools/1.23.1" {
		t.Errorf("name = %q", spec.Image.Name)
	}
	if spec.Image.Description != "SAMtools alignment toolkit" {
		t.Errorf("description = %q", spec.Image.Description)
	}
	if spec.Image.URL != "https://www.htslib.org/" {
		t.Errorf("url = %q", spec.Image.URL)
	}
	if spec.Source.BuildType() != BuildTypeScript {
		t.Errorf("build type = %v", spec.Source.BuildType())
	}
	if len(spec.Dependencies) != 0 {
		t.Errorf("dependencies = %v, want none: an app is self-contained", spec.Dependencies)
	}
	if spec.Source.Script == nil || len(spec.Source.Script.File.Data) == 0 {
		t.Error("the resolved recipe text was not captured")
	}

	m := obj.Manifest()
	if err := meta.ValidateManifest(m); err != nil {
		t.Fatalf("resolved manifest does not validate: %v", err)
	}
	if m.BuildType != "script" {
		t.Errorf("build_type = %q", m.BuildType)
	}

	rt := obj.Runtime()
	if err := meta.ValidateRuntime(rt); err != nil {
		t.Fatalf("resolved runtime does not validate: %v", err)
	}
	if len(rt.Env) != 2 || rt.Env[0].Value != "{prefix}" {
		t.Errorf("env = %+v, want {prefix} intact", rt.Env)
	}
	if rt.Env[0].Note != "install root" {
		t.Errorf("note = %q", rt.Env[0].Note)
	}
}

// Only a data recipe may declare #DEP:, so that is where dependency capture is
// exercised — and an app that declares one stops its own build.
func TestResolvedSpecDependencies(t *testing.T) {
	t.Run("data carries its deps", func(t *testing.T) {
		setTestSource(t, writeRecipe(t, "grch38/gtf/49", strings.Join([]string{
			"#DESC:GENCODE 49 annotation",
			"#DEP:samtools/1.23.1",
			"echo build",
			"",
		}, "\n")))

		obj, err := NewBuildObject(context.Background(), "grch38/gtf/49", false, t.TempDir(), false)
		if err != nil {
			t.Fatalf("NewBuildObject: %v", err)
		}
		spec := obj.Spec()
		if spec.Image.Type != catalog.TypeData {
			t.Fatalf("type = %q, want data", spec.Image.Type)
		}
		if len(spec.Dependencies) != 1 || spec.Dependencies[0] != "samtools/1.23.1" {
			t.Errorf("dependencies = %v", spec.Dependencies)
		}
	})

	t.Run("an app declaring one is rejected", func(t *testing.T) {
		setTestSource(t, writeRecipe(t, "samtools/1.23.1", "#DEP:zlib/1.3\necho build\n"))

		_, err := NewBuildObject(context.Background(), "samtools/1.23.1", false, t.TempDir(), false)
		if !errors.Is(err, catalog.ErrInvalidRecipe) {
			t.Fatalf("err = %v, want ErrInvalidRecipe", err)
		}
	})
}

// A name no collection provides is conda's. It gets an app spec with no
// description, and its source records which of the three conda modes applies.
func TestResolvedSpecForCondaFallback(t *testing.T) {
	setTestSource(t, t.TempDir())

	obj, err := NewBuildObject(context.Background(), "numpy/2.1.0", false, t.TempDir(), false)
	if err != nil {
		t.Fatalf("NewBuildObject: %v", err)
	}

	spec := obj.Spec()
	if spec.Source.BuildType() != BuildTypeConda {
		t.Fatalf("build type = %v, want conda", spec.Source.BuildType())
	}
	if spec.Source.Conda.Package == nil {
		t.Fatalf("conda source = %+v, want a single package", spec.Source.Conda)
	}
	if got := spec.Source.Conda.Package; got.Name != "numpy" || got.Version != "2.1.0" {
		t.Errorf("package = %+v", got)
	}
	// A conda environment is self-contained software whatever its name shape.
	if spec.Image.Type != catalog.TypeApp {
		t.Errorf("type = %q, want app", spec.Image.Type)
	}
	if err := meta.ValidateManifest(obj.Manifest()); err != nil {
		t.Errorf("conda manifest does not validate: %v", err)
	}
	if err := meta.ValidateRuntime(obj.Runtime()); err != nil {
		t.Errorf("conda runtime does not validate: %v", err)
	}
}

// envMap turns a KEY=VALUE slice into a map for assertions.
func envMap(t *testing.T, settings []string) map[string]string {
	t.Helper()
	out := map[string]string{}
	for _, kv := range settings {
		k, v, ok := strings.Cut(kv, "=")
		if !ok {
			t.Fatalf("malformed env entry %q", kv)
		}
		if _, dup := out[k]; dup {
			t.Errorf("%s is set more than once", k)
		}
		out[k] = v
	}
	return out
}

func TestBuildEnvFromSpec(t *testing.T) {
	spec := Spec{
		Image:  ImageSpec{Name: "samtools/1.23.1", Type: catalog.TypeApp, Prefix: "/cnt/samtools/1.23.1"},
		Source: SourceSpec{Script: &ScriptSource{}},
	}

	env := envMap(t, buildEnv(spec, Options{}))
	if env["CNT_TYPE"] != "app" {
		t.Errorf("CNT_TYPE = %q", env["CNT_TYPE"])
	}
	if env["CNT_PREFIX"] != "/cnt/samtools/1.23.1" {
		t.Errorf("CNT_PREFIX = %q", env["CNT_PREFIX"])
	}
	// Scratch is backend-neutral: a recipe never learns which mode it runs under.
	if env["CNT_TMP"] != ScratchPath || env["TMPDIR"] != ScratchPath {
		t.Errorf("scratch = %q / %q, want %s", env["CNT_TMP"], env["TMPDIR"], ScratchPath)
	}
	if env["IN_CONDATAINER"] != "1" {
		t.Errorf("IN_CONDATAINER = %q", env["IN_CONDATAINER"])
	}
	// The scheduler-normalized names are the documented exception to CNT_.
	if _, ok := env["NCPUS"]; !ok {
		t.Error("NCPUS missing; resource vars are part of the contract")
	}
}

// CNT_PREFIX has to be the same path the image records, or a payload bakes in a
// location that does not exist when the image is loaded.
func TestBuildEnvPrefixMatchesRuntime(t *testing.T) {
	for _, name := range []string{"samtools/1.23.1", "grch38/star/2.7.11b/gencode47"} {
		spec := Spec{
			Image:  ImageSpec{Name: name, Type: catalog.TypeApp, Prefix: meta.Prefix(name, catalog.TypeApp)},
			Source: SourceSpec{Script: &ScriptSource{}},
		}
		env := envMap(t, buildEnv(spec, Options{}))
		if got, want := env["CNT_PREFIX"], spec.Runtime().Prefix; got != want {
			t.Errorf("%s: CNT_PREFIX = %q, recorded prefix = %q", name, got, want)
		}
	}
}

// The environment is a projection of Spec, so an image with no type still gets
// a usable one rather than an empty CNT_TYPE the recipe would have to handle.
func TestBuildEnvDefaultsType(t *testing.T) {
	env := envMap(t, buildEnv(Spec{Image: ImageSpec{Name: "tool/1.0"}}, Options{}))
	if env["CNT_TYPE"] != "app" {
		t.Errorf("CNT_TYPE = %q, want app", env["CNT_TYPE"])
	}
	if env["CNT_PREFIX"] != "/cnt/tool/1.0" {
		t.Errorf("CNT_PREFIX = %q", env["CNT_PREFIX"])
	}
}

// CNT_NAME is the complete name and there is no CNT_VERSION: not every image
// has one version axis, so splitting at the last path segment would invent a
// boundary the name does not have.
func TestBuildEnvNameIsCompleteAndUnsplit(t *testing.T) {
	for _, name := range []string{
		"samtools/1.23.1",
		"grch38/star/2.7.11b/gencode47-101",
		"ubuntu24/base",
	} {
		spec := Spec{Image: ImageSpec{Name: name, Type: catalog.TypeApp}}
		env := envMap(t, buildEnv(spec, Options{}))
		if env["CNT_NAME"] != name {
			t.Errorf("CNT_NAME = %q, want the complete name %q", env["CNT_NAME"], name)
		}
		if v, ok := env["CNT_VERSION"]; ok {
			t.Errorf("CNT_VERSION = %q; it was removed", v)
		}
	}
}

// The stored recipe is the template, tokens intact: it is what a rebuild starts
// from and what both keys hash. The expansion exists only in the workspace, and
// which variant this is survives as manifest.source.placeholders.
func TestResolvedTemplateEmbedsTheTemplate(t *testing.T) {
	setTestSource(t, writeRecipe(t, "grch38/star-gencode", strings.Join([]string{
		"#DESC:STAR {star_version} index for GENCODE {gencode_version}",
		"#TARGET:grch38/star/{star_version}/gencode{gencode_version}",
		"#PH:star_version:2.7.11b,2.7.11a",
		"#PH:gencode_version:47-49",
		"echo building {star_version} against {gencode_version}",
		"",
	}, "\n")))

	obj, err := NewBuildObject(context.Background(),
		"grch38/star/2.7.11b/gencode49", false, t.TempDir(), false)
	if err != nil {
		t.Fatalf("NewBuildObject: %v", err)
	}
	spec := obj.Spec()

	embedded, ok := spec.Source.RecipeFile()
	if !ok {
		t.Fatal("a script build embeds no recipe")
	}
	if !strings.Contains(string(embedded.Data), "echo building {star_version} against {gencode_version}") {
		t.Errorf("the embedded recipe was expanded:\n%s", embedded.Data)
	}

	// The workspace copy is the one that runs, and it is expanded.
	ran, err := os.ReadFile(obj.BuildSource())
	if err != nil {
		t.Fatalf("reading the materialized recipe: %v", err)
	}
	if !strings.Contains(string(ran), "echo building 2.7.11b against 49") {
		t.Errorf("the recipe that runs was not expanded:\n%s", ran)
	}

	// Which variant it is lives in the manifest, since the recipe no longer says.
	m := obj.Manifest()
	if got := m.Source.Placeholders; got["star_version"] != "2.7.11b" || got["gencode_version"] != "49" {
		t.Errorf("placeholders = %v", got)
	}
	if m.Source.TargetTemplate == "" {
		t.Error("the manifest records no target template")
	}
	if len(m.Source.Files) != 1 || m.Source.Files[0] != meta.RecipeFileName {
		t.Errorf("source files = %v, want the embedded recipe named", m.Source.Files)
	}
	if m.Source.RequiresInput {
		t.Error("a recipe with no #INPUT: should not claim it needs one")
	}
}

// A plain recipe has no placeholders and no target, so the source block says
// only what was embedded.
func TestResolvedNonTemplateSourceBlock(t *testing.T) {
	setTestSource(t, writeRecipe(t, "samtools/1.23.1", "#DESC:SAMtools\necho build\n"))

	obj, err := NewBuildObject(context.Background(), "samtools/1.23.1", false, t.TempDir(), false)
	if err != nil {
		t.Fatalf("NewBuildObject: %v", err)
	}
	m := obj.Manifest()
	if len(m.Source.Placeholders) != 0 || m.Source.TargetTemplate != "" {
		t.Errorf("source = %+v, want no template fields", m.Source)
	}
	if m.Source.RequiresInput {
		t.Error("a recipe with no #INPUT: should not claim it needs one")
	}
}

// requires_input records that a rebuild needs a human, and that is the *only*
// thing an #INPUT: may leave behind: the prompt is not embedded and neither is
// the answer.
func TestSourceBlockRecordsOnlyThatInputIsNeeded(t *testing.T) {
	spec := Spec{
		Image: ImageSpec{Name: "vendor/tool/1.0", Type: catalog.TypeApp, Prefix: "/cnt/vendor/tool/1.0"},
		Source: SourceSpec{
			Script:        &ScriptSource{Prompts: []string{"paste the licence key"}},
			RequiresInput: true,
		},
	}
	m := spec.Manifest()
	if !m.Source.RequiresInput {
		t.Error("requires_input was not recorded")
	}
	data, err := meta.MarshalManifest(m)
	if err != nil {
		t.Fatal(err)
	}
	if strings.Contains(string(data), "licence key") {
		t.Errorf("an #INPUT: prompt reached the manifest:\n%s", data)
	}
}

// A Conda build records the channels it was offered, in priority order. They are
// captured into Spec at resolution rather than read from config when the
// manifest renders, so the manifest stays a projection of Spec and nothing else.
func TestCondaSpecCapturesChannels(t *testing.T) {
	prev := config.Global.Build.Channels
	config.Global.Build.Channels = []string{"conda-forge", "bioconda"}
	t.Cleanup(func() { config.Global.Build.Channels = prev })

	setTestSource(t, t.TempDir())
	obj, err := NewBuildObject(context.Background(), "numpy/2.1.0", false, t.TempDir(), false)
	if err != nil {
		t.Fatalf("NewBuildObject: %v", err)
	}
	if got := obj.Spec().Source.Conda.Channels; !slices.Equal(got, []string{"conda-forge", "bioconda"}) {
		t.Fatalf("channels = %v", got)
	}
	if got := obj.Manifest().Build.Channels; !slices.Equal(got, []string{"conda-forge", "bioconda"}) {
		t.Errorf("manifest channels = %v", got)
	}

	// Changing config afterwards must not change what this build records.
	config.Global.Build.Channels = []string{"nvidia"}
	if got := obj.Manifest().Build.Channels; !slices.Equal(got, []string{"conda-forge", "bioconda"}) {
		t.Errorf("manifest channels followed config after resolution: %v", got)
	}
}

// A recipe build has no channels to record, so the block stays out of the JSON
// entirely rather than appearing empty.
func TestNonCondaManifestHasNoBuildBlock(t *testing.T) {
	spec := Spec{
		Image:  ImageSpec{Name: "samtools/1.23.1", Type: catalog.TypeApp, Prefix: "/cnt/samtools/1.23.1"},
		Source: SourceSpec{Script: &ScriptSource{}},
	}
	data, err := meta.MarshalManifest(spec.Manifest())
	if err != nil {
		t.Fatal(err)
	}
	if strings.Contains(string(data), `"build":`) {
		t.Errorf("an empty build block was written:\n%s", data)
	}
}

// The manifest names exactly the files that were staged, so a reader never has
// to probe for what an image carries.
func TestManifestNamesOnlyStagedSources(t *testing.T) {
	b := newPackObject(t, catalog.TypeApp)
	if got := b.Manifest().Source.Files; len(got) != 0 {
		t.Errorf("files = %v before anything was captured", got)
	}

	b.embedSource(SourceFile{Name: "explicit.txt", Data: []byte("@EXPLICIT\n")})
	b.embedSource(SourceFile{Name: "environment.yml", Data: []byte("channels: []\n")})
	if got := b.Manifest().Source.Files; !slices.Equal(got, []string{"explicit.txt", "environment.yml"}) {
		t.Errorf("files = %v", got)
	}

	// Re-capturing replaces rather than duplicating.
	b.embedSource(SourceFile{Name: "explicit.txt", Data: []byte("@EXPLICIT\nhttps://x\n")})
	if got := b.Manifest().Source.Files; len(got) != 2 {
		t.Errorf("files = %v, want the earlier capture replaced", got)
	}

	// A record is derived from the sources, not one of them: manifest.keys names
	// it, and source.files must not.
}

// A recipe build writes both records, and manifest.keys names files whose
// digests reproduce the keys — the property that lets a reader verify without
// knowing which build type produced the image.
func TestRecipeBuildRecordsKeys(t *testing.T) {
	setTestSource(t, writeRecipe(t, "samtools/1.23.1", strings.Join([]string{
		"#DESC:SAMtools",
		"#ENV:SAMTOOLS_DIR={prefix}  ## install root",
		"echo build",
		"",
	}, "\n")))

	obj, err := NewBuildObject(context.Background(), "samtools/1.23.1", false, t.TempDir(), false)
	if err != nil {
		t.Fatalf("NewBuildObject: %v", err)
	}
	if err := obj.deriveKeys(context.Background()); err != nil {
		t.Fatalf("deriveKeys: %v", err)
	}

	m := obj.Manifest()
	if m.Keys.Identity.Scheme != string(key.ScriptIdentityV1) ||
		m.Keys.Equiv.Scheme != string(key.ScriptEquivV1) {
		t.Fatalf("keys = %+v", m.Keys)
	}

	for _, f := range obj.embedded {
		if strings.HasSuffix(f.Name, ".record") {
			t.Errorf("recipe build staged obsolete file %s", f.Name)
		}
	}

	// An app has no dependencies, so identity and equivalence have the same
	// canonical bytes and SHA. Their scheme fields keep the two keys distinct.
	if m.Keys.Identity.SHA256 != m.Keys.Equiv.SHA256 {
		t.Errorf("untagged app keys have different SHA values: %+v", m.Keys)
	}
	derived, err := key.Verify(m, obj.keySources())
	if err != nil {
		t.Fatal(err)
	}
	identity := derived.IdentityModel
	if identity == nil {
		t.Fatal("script scheme produced no semantic identity")
	}
	if len(identity.Env) != 1 || identity.Env[0].Key != "SAMTOOLS_DIR" {
		t.Errorf("env = %v, want the #ENV: contribution", identity.Env)
	}
	if identity.Recipe == "" {
		t.Error("the record carries no recipe digest")
	}
	if len(m.Dependencies) != 0 || m.ProvenanceComplete != nil {
		t.Errorf("an app recorded dependencies: %v / %v", m.Dependencies, m.ProvenanceComplete)
	}
}

// A base is keyed by its definition and the upstream it bootstrapped from,
// resolved before the build — which is what lets a SIF carry its own keys.
func TestBaseIsKeyedByItsDefinition(t *testing.T) {
	upstream := "sha256:" + strings.Repeat("a", 64)
	b := newPackObject(t, catalog.TypeBase)
	b.spec.Source = SourceSpec{Definition: &DefinitionSource{
		File: SourceFile{Name: meta.RecipeFileName, Data: []byte("Bootstrap: docker\nFrom: ubuntu:24.04\n")},
		From: &meta.From{Bootstrap: "docker", Ref: "ubuntu:24.04", Digest: upstream},
	}}
	b.embedSource(b.spec.Source.Definition.File)
	if err := b.deriveKeys(context.Background()); err != nil {
		t.Fatalf("deriveKeys: %v", err)
	}

	m := b.Manifest()
	if m.Keys.Identity.Scheme != string(key.DefinitionIdentityV1) ||
		m.Keys.Equiv.Scheme != string(key.DefinitionEquivV1) {
		t.Fatalf("a base was not keyed: %+v", m.Keys)
	}
	if m.Build.From == nil || m.Build.From.URI() != "docker://ubuntu:24.04" || m.Build.From.Digest != upstream {
		t.Errorf("build.from = %+v", m.Build.From)
	}

	derived, err := key.Verify(m, b.keySources())
	if err != nil {
		t.Fatal(err)
	}
	identity := derived.IdentityModel
	if identity.Type != catalog.TypeBase {
		t.Errorf("type = %q, want base", identity.Type)
	}
	if identity.From != upstream {
		t.Errorf("from = %q, want the resolved upstream", identity.From)
	}

	// The upstream identifies the build without deciding substitution.
	if strings.Contains(string(derived.Equiv.Preimage), "from=") {
		t.Errorf("the upstream digest reached the equivalence preimage:\n%s", derived.Equiv.Preimage)
	}
}

// A Conda app writes no records: everything one could hold is ruled out for it,
// so the manifest names the exports themselves and hashing them gives the keys.
func TestCondaKeysNameTheExports(t *testing.T) {
	b := newPackObject(t, catalog.TypeApp)
	explicit := []byte("@EXPLICIT\nhttps://conda.anaconda.org/conda-forge/linux-64/bzip2-1.0.8-hda65f42_9.conda\n")
	environment := []byte("channels:\n  - conda-forge\ndependencies:\n  - bzip2=1.0.8\n")
	b.embedSource(SourceFile{Name: conda.ExplicitFileName, Data: explicit})
	b.embedSource(SourceFile{Name: conda.EnvironmentFileName, Data: environment})

	if err := b.deriveKeys(context.Background()); err != nil {
		t.Fatalf("deriveKeys: %v", err)
	}
	m := b.Manifest()
	if m.Keys.Identity.Scheme != string(key.CondaExplicitV1) ||
		m.Keys.Identity.SHA256 != key.Sum(explicit) {
		t.Errorf("identity = %+v", m.Keys.Identity)
	}
	if m.Keys.Equiv.Scheme != string(key.CondaEnvironmentV1) ||
		m.Keys.Equiv.SHA256 != key.Sum(environment) {
		t.Errorf("equiv = %+v", m.Keys.Equiv)
	}
	for _, f := range b.embedded {
		if strings.HasSuffix(f.Name, ".record") {
			t.Errorf("a conda build staged %s", f.Name)
		}
	}
}

// A capture that failed leaves no keys rather than a claim about a file the
// image does not carry.
func TestCondaWithoutExportsRecordsNoKeys(t *testing.T) {
	b := newPackObject(t, catalog.TypeApp)
	if err := b.deriveKeys(context.Background()); err != nil {
		t.Fatalf("deriveKeys: %v", err)
	}
	if b.Manifest().Keys != (meta.Keys{}) {
		t.Errorf("keys were claimed with nothing captured: %+v", b.Manifest().Keys)
	}
}

// A data build embeds the closure it was built from, read out of its
// dependencies' own images rather than re-derived from the catalog.
func TestDataBuildComposesACapsule(t *testing.T) {
	if _, err := exec.LookPath("mksquashfs"); err != nil {
		t.Skip("mksquashfs not available")
	}

	// A dependency image carrying scheme-backed sources and a capsule entry of its own.
	imagesDir := t.TempDir()
	depRoot := t.TempDir()
	depMeta := filepath.Join(depRoot, meta.DirName)
	complete := true
	depRecipe := []byte("echo gtf\n")
	depManifest := meta.Manifest{
		SchemaVersion:      meta.SchemaVersion,
		Name:               "grch38/gtf/49",
		Type:               catalog.TypeData,
		BuildType:          "script",
		Platform:           meta.NativePlatform(),
		Source:             meta.Source{Files: []string{meta.RecipeFileName}},
		ProvenanceComplete: &complete,
	}
	derived, err := key.Generate(depManifest, key.Sources{meta.RecipeFileName: depRecipe})
	if err != nil {
		t.Fatal(err)
	}
	depManifest.Keys = derived.Keys()
	if err := meta.StageManifest(depMeta, depManifest); err != nil {
		t.Fatal(err)
	}
	if err := meta.StageRuntime(depMeta, meta.Runtime{
		SchemaVersion: meta.SchemaVersion, Name: "grch38/gtf/49", Type: catalog.TypeData,
		Platform: meta.NativePlatform(), Prefix: "/cnt/grch38/gtf/49",
	}); err != nil {
		t.Fatal(err)
	}
	if err := meta.StageBytes(depMeta, meta.RecipeFileName, depRecipe); err != nil {
		t.Fatal(err)
	}
	inheritedRecipe := []byte("echo genome\n")
	inheritedManifest := meta.Manifest{
		SchemaVersion: meta.SchemaVersion,
		Name:          "grch38/genome", Type: catalog.TypeData, BuildType: "script",
		Platform: meta.NativePlatform(),
		Source:   meta.Source{Files: []string{meta.RecipeFileName}},
	}
	inheritedDerived, err := key.Generate(inheritedManifest, key.Sources{meta.RecipeFileName: inheritedRecipe})
	if err != nil {
		t.Fatal(err)
	}
	inheritedManifest.Keys = inheritedDerived.Keys()
	inherited := filepath.Join(depMeta, capsule.DirName,
		capsule.EntryName(inheritedManifest.Name, inheritedManifest.Keys.Identity.Digest()))
	if err := os.MkdirAll(inherited, 0o755); err != nil {
		t.Fatal(err)
	}
	if err := meta.StageManifest(inherited, inheritedManifest); err != nil {
		t.Fatal(err)
	}
	if err := meta.StageBytes(inherited, meta.RecipeFileName, inheritedRecipe); err != nil {
		t.Fatal(err)
	}
	depImage := filepath.Join(imagesDir, "grch38--gtf--49.sqf")
	if out, err := exec.Command("mksquashfs", depRoot, depImage, "-no-progress", "-noappend", "-quiet", "-no-xattrs").CombinedOutput(); err != nil {
		t.Fatalf("mksquashfs: %v\n%s", err, out)
	}

	// GlobalDataPaths caches resolved search paths, so setting the env var only
	// takes effect on a cold cache — set the resolved list, as the rest of this
	// package's tests do.
	prevPaths := config.GlobalDataPaths
	config.GlobalDataPaths.ImagesDirs = []string{imagesDir}
	container.InvalidateInstalledOverlaysCache()
	t.Cleanup(func() {
		config.GlobalDataPaths = prevPaths
		container.InvalidateInstalledOverlaysCache()
	})

	setTestSource(t, writeRecipe(t, "grch38/star/index", strings.Join([]string{
		"#DESC:index",
		"#DEP:grch38/gtf/49",
		"echo build",
		"",
	}, "\n")))

	obj, err := NewBuildObject(context.Background(), "grch38/star/index", false, t.TempDir(), false)
	if err != nil {
		t.Fatalf("NewBuildObject: %v", err)
	}
	if err := os.RemoveAll(obj.ws.MetaDir); err != nil {
		t.Fatal(err)
	}
	if err := obj.deriveKeys(context.Background()); err != nil {
		t.Fatalf("deriveKeys: %v", err)
	}

	m := obj.Manifest()
	if len(m.Dependencies) != 1 {
		t.Fatalf("dependencies = %+v", m.Dependencies)
	}
	dep := m.Dependencies[0]
	if dep.Name != "grch38/gtf/49" || dep.Type != catalog.TypeData || dep.Role != meta.RoleData {
		t.Errorf("dependency = %+v", dep)
	}
	if dep.Records == meta.Unrecorded {
		t.Errorf("a recorded dependency was read as unrecorded: %+v", dep)
	}
	if m.ProvenanceComplete == nil || !*m.ProvenanceComplete {
		t.Errorf("provenance_complete = %v, want true", m.ProvenanceComplete)
	}

	entries, err := capsule.Entries(filepath.Join(obj.ws.MetaDir, capsule.DirName))
	if err != nil {
		t.Fatalf("Entries: %v", err)
	}
	var dirs []string
	for _, e := range entries {
		dirs = append(dirs, e.Dir)
	}
	// The dependency itself, and the entry it carried, unioned into one flat set.
	if len(dirs) != 2 || !strings.HasPrefix(dirs[1], "grch38--gtf--49@") ||
		dirs[0] != capsule.EntryName(inheritedManifest.Name, inheritedManifest.Keys.Identity.Digest()) {
		t.Errorf("capsule entries = %v", dirs)
	}
}

// #ARCH:noarch is the only way an artifact claims portability. Nothing infers
// it: getting it wrong does not crash, it silently returns wrong answers.
func TestArchReachesTheRecordedPlatform(t *testing.T) {
	tests := []struct {
		name   string
		recipe string
		want   string
	}{
		{"portable", "#DESC:jars\n#ARCH:noarch\necho build\n", meta.ArchNone},
		{"default is strict", "#DESC:compiled\necho build\n", meta.NativeArch()},
		{"explicit native", "#DESC:compiled\n#ARCH:native\necho build\n", meta.NativeArch()},
	}
	for _, tt := range tests {
		t.Run(tt.name, func(t *testing.T) {
			setTestSource(t, writeRecipe(t, "picard/3.1", tt.recipe))
			obj, err := NewBuildObject(context.Background(), "picard/3.1", false, t.TempDir(), false)
			if err != nil {
				t.Fatalf("NewBuildObject: %v", err)
			}
			if got := obj.Runtime().Platform.Arch; got != tt.want {
				t.Errorf("runtime arch = %q, want %q", got, tt.want)
			}
			if got := obj.Manifest().Platform.Arch; got != tt.want {
				t.Errorf("manifest arch = %q, want %q", got, tt.want)
			}
			// The OS label never varies: everything runs in a Linux container.
			if got := obj.Runtime().Platform.OS; got != "linux" {
				t.Errorf("os = %q", got)
			}
		})
	}
}

// Adding #ARCH: later must not invalidate images already built: it is a
// compatibility assertion about the artifact, not a statement about its contents.
func TestArchEntersNoKey(t *testing.T) {
	keyOf := func(t *testing.T, recipe string) string {
		t.Helper()
		setTestSource(t, writeRecipe(t, "picard/3.1", recipe))
		obj, err := NewBuildObject(context.Background(), "picard/3.1", false, t.TempDir(), false)
		if err != nil {
			t.Fatalf("NewBuildObject: %v", err)
		}
		if err := obj.deriveKeys(context.Background()); err != nil {
			t.Fatalf("deriveKeys: %v", err)
		}
		return obj.Manifest().Keys.Identity.SHA256
	}
	before := keyOf(t, "#DESC:jars\necho build\n")
	after := keyOf(t, "#DESC:jars\n#ARCH:noarch\necho build\n")
	if before != after {
		t.Error("adding #ARCH: moved the identity, invalidating every image built before it")
	}
}
