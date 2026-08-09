package build

import (
	"context"
	"strings"
	"testing"

	"github.com/Justype/condatainer/catalog"
	"github.com/Justype/condatainer/internal/image/meta"
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
		Runtime: meta.Runtime{
			Prefix: "/cnt/samtools/1.23.1",
			Env:    []meta.EnvVar{{Key: "SAMTOOLS_DIR", Value: "{prefix}", Note: "install root"}},
		},
	}

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
	if err := meta.Validate(m); err != nil {
		t.Errorf("rendered manifest does not validate: %v", err)
	}
}

// The manifest is a projection of Spec, and Spec has nowhere to put an #INPUT:
// answer. This pins that: answers are execution input and must never be
// embedded, because they routinely carry tokens and licence keys.
func TestSpecManifestCannotCarryAnswers(t *testing.T) {
	const secret = "s3cret-license-key"

	spec := Spec{
		Image:   ImageSpec{Name: "tool/1.0", Type: catalog.TypeApp},
		Source:  SourceSpec{Script: &ScriptSource{Prompts: []string{"License key:"}}},
		Runtime: meta.Runtime{Prefix: "/cnt/tool/1.0"},
	}
	data, err := meta.Marshal(spec.Manifest())
	if err != nil {
		t.Fatalf("Marshal: %v", err)
	}
	if strings.Contains(string(data), secret) {
		t.Fatalf("an answer reached the manifest:\n%s", data)
	}
	// The prompt itself is a declaration, not an answer, and is also not embedded.
	if strings.Contains(string(data), "License key:") {
		t.Errorf("a prompt declaration reached the manifest:\n%s", data)
	}
}

func TestRuntimeFromRecipeKeepsPrefixToken(t *testing.T) {
	recipe, err := catalog.ParseRecipe("samtools/1.23.1", strings.NewReader(
		"#DESC:SAMtools\n"+
			"#ENV:SAMTOOLS_DIR={prefix}  ## install root\n"+
			"#ENV:PATH_EXTRA={prefix}/bin\n",
	))
	if err != nil {
		t.Fatalf("ParseRecipe: %v", err)
	}

	rt := runtimeFromRecipe("samtools/1.23.1", catalog.TypeApp, recipe.Env)
	if rt.Prefix != "/cnt/samtools/1.23.1" {
		t.Errorf("prefix = %q", rt.Prefix)
	}
	if len(rt.Env) != 2 {
		t.Fatalf("env = %+v", rt.Env)
	}
	// {prefix} must survive to the manifest: the install prefix is not known
	// when the image is built, only when it is loaded.
	if rt.Env[0].Value != "{prefix}" {
		t.Errorf("env[0] = %q, want the token intact", rt.Env[0].Value)
	}
	if rt.Env[0].Note != "install root" {
		t.Errorf("note = %q", rt.Env[0].Note)
	}
	if got := rt.Env[1].Resolved(rt.Prefix); got != "/cnt/samtools/1.23.1/bin" {
		t.Errorf("resolved = %q", got)
	}
}

// base and os apply at the container root, so they carry no prefix and a
// {prefix} token in their env has nothing to resolve against.
func TestRuntimeFromRecipeRootTypes(t *testing.T) {
	for _, typ := range []catalog.Type{catalog.TypeBase, catalog.TypeOS} {
		rt := runtimeFromRecipe("ubuntu24/base", typ, nil)
		if rt.Prefix != "" {
			t.Errorf("%s prefix = %q, want empty", typ, rt.Prefix)
		}
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

// Resolution has to produce a Spec that renders a valid manifest, from a real
// catalog lookup rather than a hand-built struct. This is what proves the
// descriptive metadata mapping — #DESC: to Description, #URL: to URL, #ENV: to
// Runtime.Env — actually runs end to end.
func TestResolvedSpecRendersValidManifest(t *testing.T) {
	setTestSource(t, writeRecipe(t, "samtools/1.23.1", strings.Join([]string{
		"#!/usr/bin/env bash",
		"#DESC:SAMtools alignment toolkit",
		"#URL:https://www.htslib.org/",
		"#DEP:zlib/1.3",
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
	if len(spec.Dependencies) != 1 || spec.Dependencies[0] != "zlib/1.3" {
		t.Errorf("dependencies = %v", spec.Dependencies)
	}
	if spec.Source.Script == nil || len(spec.Source.Script.File.Data) == 0 {
		t.Error("the resolved recipe text was not captured")
	}

	m := obj.Manifest()
	if err := meta.Validate(m); err != nil {
		t.Fatalf("resolved manifest does not validate: %v", err)
	}
	if m.BuildType != "script" {
		t.Errorf("build_type = %q", m.BuildType)
	}
	if len(m.Runtime.Env) != 2 || m.Runtime.Env[0].Value != "{prefix}" {
		t.Errorf("env = %+v, want {prefix} intact", m.Runtime.Env)
	}
	if m.Runtime.Env[0].Note != "install root" {
		t.Errorf("note = %q", m.Runtime.Env[0].Note)
	}
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
	if err := meta.Validate(obj.Manifest()); err != nil {
		t.Errorf("conda manifest does not validate: %v", err)
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
		Image:   ImageSpec{Name: "samtools/1.23.1", Type: catalog.TypeApp},
		Source:  SourceSpec{Script: &ScriptSource{}},
		Runtime: meta.Runtime{Prefix: "/cnt/samtools/1.23.1"},
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

// CNT_PREFIX has to be the same path the manifest records, or a payload bakes
// in a location that does not exist when the image is loaded.
func TestBuildEnvPrefixMatchesManifest(t *testing.T) {
	for _, name := range []string{"samtools/1.23.1", "grch38/star/2.7.11b/gencode47"} {
		spec := Spec{
			Image:   ImageSpec{Name: name, Type: catalog.TypeApp},
			Source:  SourceSpec{Script: &ScriptSource{}},
			Runtime: runtimeFromRecipe(name, catalog.TypeApp, nil),
		}
		env := envMap(t, buildEnv(spec, Options{}))
		if got, want := env["CNT_PREFIX"], spec.Manifest().Runtime.Prefix; got != want {
			t.Errorf("%s: CNT_PREFIX = %q, manifest prefix = %q", name, got, want)
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
