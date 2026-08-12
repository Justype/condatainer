package compare

import (
	"os"
	"os/exec"
	"path/filepath"
	"strings"
	"testing"

	"github.com/Justype/condatainer/catalog"
	"github.com/Justype/condatainer/internal/artifact/key"
	"github.com/Justype/condatainer/internal/artifact/meta"
	"github.com/Justype/condatainer/internal/conda"
)

func requireSquashfsTools(t *testing.T) {
	t.Helper()
	for _, bin := range []string{"mksquashfs", "unsquashfs"} {
		if _, err := exec.LookPath(bin); err != nil {
			t.Skipf("%s not available", bin)
		}
	}
}

// build stages an image the way a build backend does, then packs it — so the
// whole read path is exercised against a real archive rather than a directory.
type build struct {
	name     string
	typ      catalog.Type
	arch     string
	format   string
	env      []meta.EnvVar
	recipe   string
	from     string
	ph       map[string]string
	deps     []key.Dep
	explicit string // set for a Conda app instead of a recipe
	environ  string
	tamper   func(dir string)
}

func pack(t *testing.T, b build) string {
	t.Helper()
	root := t.TempDir()
	dir := filepath.Join(root, meta.DirName)

	arch := b.arch
	if arch == "" {
		arch = meta.NativeArch()
	}
	format := b.format
	if format == "" {
		format = "script"
	}
	if err := meta.StageRuntime(dir, meta.Runtime{
		SchemaVersion: meta.SchemaVersion,
		Name:          b.name,
		Type:          b.typ,
		Platform:      meta.Platform{OS: "linux", Arch: arch},
		Prefix:        meta.Prefix(b.name, b.typ),
		Env:           b.env,
	}); err != nil {
		t.Fatal(err)
	}

	manifest := meta.Manifest{
		SchemaVersion: meta.SchemaVersion,
		Name:          b.name,
		Type:          b.typ,
		BuildType:     format,
		Platform:      meta.Platform{OS: "linux", Arch: arch},
	}
	sources := key.Sources{}

	switch {
	case b.explicit != "":
		if err := meta.StageBytes(dir, conda.ExplicitFileName, []byte(b.explicit)); err != nil {
			t.Fatal(err)
		}
		if err := meta.StageBytes(dir, conda.EnvironmentFileName, []byte(b.environ)); err != nil {
			t.Fatal(err)
		}
		manifest.Source.Files = []string{conda.ExplicitFileName, conda.EnvironmentFileName}
		sources[conda.ExplicitFileName] = []byte(b.explicit)
		sources[conda.EnvironmentFileName] = []byte(b.environ)
	case b.recipe != "":
		artifact := key.Artifact{
			Name: b.name, Type: b.typ, Env: b.env,
			Recipe: []byte(b.recipe), Placeholders: b.ph, Deps: b.deps,
			From: b.from,
		}
		if err := meta.StageBytes(dir, meta.RecipeFileName, []byte(b.recipe)); err != nil {
			t.Fatal(err)
		}
		manifest.Source.Files = []string{meta.RecipeFileName}
		manifest.Source.Placeholders = b.ph
		sources[meta.RecipeFileName] = []byte(b.recipe)
		manifest.Dependencies, manifest.ProvenanceComplete = key.Manifest(artifact)
		if b.from != "" {
			manifest.Build.From = &meta.From{Bootstrap: "docker", Ref: "example:latest", Digest: b.from}
		}
	}
	if len(sources) > 0 {
		derived, err := key.Generate(manifest, sources)
		if err != nil {
			t.Fatal(err)
		}
		manifest.Keys = derived.Keys()
	}

	if err := meta.StageManifest(dir, manifest); err != nil {
		t.Fatal(err)
	}
	if b.tamper != nil {
		b.tamper(dir)
	}

	payload := filepath.Join(root, "cnt", b.name)
	if err := os.MkdirAll(payload, 0o755); err != nil {
		t.Fatal(err)
	}
	if err := os.WriteFile(filepath.Join(payload, "data"), []byte("x"), 0o644); err != nil {
		t.Fatal(err)
	}

	out := filepath.Join(t.TempDir(), "image.sqf")
	cmd := exec.Command("mksquashfs", root, out, "-no-progress", "-noappend", "-quiet", "-no-xattrs")
	if output, err := cmd.CombinedOutput(); err != nil {
		t.Fatalf("mksquashfs: %v\n%s", err, output)
	}
	return out
}

func read(t *testing.T, b build) Artifact {
	t.Helper()
	a, err := Read(pack(t, b))
	if err != nil {
		t.Fatalf("Read: %v", err)
	}
	return a
}

const starRecipe = "#TARGET:grch38/star/2.7.11b/gencode{gencode_version}\n#PH:gencode_version:47,49\n#ENV:STAR_INDEX={prefix}/index  ## for --genomeDir\n#DESC:index\nSTAR --runMode genomeGenerate\n"

// depKey is a dependency's complete key. Comparison holds an edge to both
// halves, so a fixture that supplied only a digest would not exercise it.
func depKey(scheme key.Scheme, seed string) meta.KeyRef {
	return meta.KeyRef{Scheme: string(scheme), SHA256: key.Sum([]byte(seed))}
}

func starIndex() build {
	return build{
		name:   "grch38/star/2.7.11b/gencode49",
		typ:    catalog.TypeData,
		env:    []meta.EnvVar{{Key: "STAR_INDEX", Value: "{prefix}/index", Note: "for --genomeDir"}},
		recipe: starRecipe,
		ph:     map[string]string{"gencode_version": "49"},
		deps: []key.Dep{
			{Name: "star/2.7.11b", Type: catalog.TypeApp,
				Identity: depKey(key.ScriptIdentityV1, "star id"), Equiv: depKey(key.ScriptEquivV1, "star eq")},
			{Name: "samtools/1.23.1", Type: catalog.TypeApp,
				Identity: depKey(key.ScriptIdentityV1, "sam id"), Equiv: depKey(key.ScriptEquivV1, "sam eq")},
		},
	}
}

func TestVerdicts(t *testing.T) {
	requireSquashfsTools(t)
	want := read(t, starIndex())

	t.Run("the same build is exact", func(t *testing.T) {
		if got := Compare(want, read(t, starIndex())); got.Verdict != Exact {
			t.Errorf("verdict = %s (%s)", got.Verdict, got.Reason)
		}
	})

	t.Run("a rebuilt history-only dependency is equivalent", func(t *testing.T) {
		b := starIndex()
		b.deps[1].Identity = depKey(key.ScriptIdentityV1, "samtools rebuilt")
		b.deps[1].Equiv = depKey(key.ScriptEquivV1, "samtools rebuilt eq")

		got := Compare(want, read(t, b))
		if got.Verdict != Equivalent {
			t.Fatalf("verdict = %s (%s)", got.Verdict, got.Reason)
		}
		// The diff still names which dependency moved, from the manifest's list —
		// the equivalence model deliberately dropped it.
		if !hasField(got.Diffs, "dep:samtools/1.23.1") {
			t.Errorf("diffs = %v", got.Diffs)
		}
	})

	t.Run("a changed recipe is different", func(t *testing.T) {
		b := starIndex()
		b.recipe = strings.Replace(starRecipe, "genomeGenerate", "genomeGenerate --sjdbOverhang 100", 1)

		got := Compare(want, read(t, b))
		if got.Verdict != Different {
			t.Fatalf("verdict = %s", got.Verdict)
		}
		if !hasField(got.Diffs, "recipe") {
			t.Errorf("diffs = %v, want the recipe named", got.Diffs)
		}
	})

	t.Run("a changed placeholder is different and names itself", func(t *testing.T) {
		b := starIndex()
		b.ph = map[string]string{"gencode_version": "47"}

		got := Compare(want, read(t, b))
		if got.Verdict != Different {
			t.Fatalf("verdict = %s", got.Verdict)
		}
		if !hasDiff(got.Diffs, Diff{Field: "ph:gencode_version", Want: "49", Got: "47"}) {
			t.Errorf("diffs = %v", got.Diffs)
		}
	})

	t.Run("a comment moves nothing", func(t *testing.T) {
		b := starIndex()
		b.recipe = strings.Replace(starRecipe, "#DESC:index", "#DESC:reworded\n# explain it", 1)
		if got := Compare(want, read(t, b)); got.Verdict != Exact {
			t.Errorf("verdict = %s, want exact", got.Verdict)
		}
	})
}

// Gates run before keys, and a gate failure is a verdict with a specific diff
// rather than a fallthrough to equivalence.
func TestGates(t *testing.T) {
	requireSquashfsTools(t)
	want := read(t, starIndex())

	t.Run("a foreign architecture is refused", func(t *testing.T) {
		b := starIndex()
		b.arch = "aarch64"
		got := compareOn(want, read(t, b), "x86_64")
		if got.Verdict != Different || !hasField(got.Diffs, "arch") {
			t.Errorf("verdict = %s, diffs = %v", got.Verdict, got.Diffs)
		}
	})

	t.Run("noarch passes on any host", func(t *testing.T) {
		b := starIndex()
		b.arch = meta.ArchNone
		portable := read(t, b)
		// Compared against itself on a host it was not built on.
		if got := compareOn(portable, portable, "riscv64"); got.Verdict != Exact {
			t.Errorf("verdict = %s (%s)", got.Verdict, got.Reason)
		}
	})

	t.Run("a different type is never interchangeable", func(t *testing.T) {
		b := starIndex()
		b.typ = catalog.TypeApp
		b.deps = nil
		got := Compare(want, read(t, b))
		if got.Verdict != Different || !hasField(got.Diffs, "type") {
			t.Errorf("verdict = %s, diffs = %v", got.Verdict, got.Diffs)
		}
	})

	t.Run("a different name is never compared", func(t *testing.T) {
		b := starIndex()
		b.name = "grch38/star/2.7.11b/gencode47"
		got := Compare(want, read(t, b))
		if got.Verdict != Different || !hasField(got.Diffs, "name") {
			t.Errorf("verdict = %s, diffs = %v", got.Verdict, got.Diffs)
		}
	})

	// The runtime gate is free integrity: runtime.json was never the hashed
	// preimage, so an edited one cannot agree with the regenerated model.
	t.Run("an edited runtime.json is unverifiable", func(t *testing.T) {
		b := starIndex()
		b.tamper = func(dir string) {
			rt := meta.Runtime{
				SchemaVersion: meta.SchemaVersion, Name: b.name, Type: b.typ,
				Platform: meta.NativePlatform(), Prefix: meta.Prefix(b.name, b.typ),
				Env: []meta.EnvVar{{Key: "STAR_INDEX", Value: "/somewhere/else"}},
			}
			if err := meta.StageRuntime(dir, rt); err != nil {
				t.Fatal(err)
			}
		}
		got := Compare(want, read(t, b))
		if got.Verdict != Unverifiable {
			t.Errorf("verdict = %s, want unverifiable", got.Verdict)
		}
	})
}

// A key that does not match its own file is exactly what verification exists to
// catch — the manifest is what the publisher wrote, not proof.
func TestTamperedKeysAreUnverifiable(t *testing.T) {
	requireSquashfsTools(t)
	want := read(t, starIndex())

	b := starIndex()
	b.tamper = func(dir string) {
		if err := os.WriteFile(filepath.Join(dir, meta.RecipeFileName),
			[]byte(strings.Replace(starRecipe, "genomeGenerate", "tampered", 1)), 0o644); err != nil {
			t.Fatal(err)
		}
	}
	got := Compare(want, read(t, b))
	if got.Verdict != Unverifiable {
		t.Errorf("verdict = %s, want unverifiable", got.Verdict)
	}
	if !strings.Contains(got.Reason, "derives identity") {
		t.Errorf("reason = %q, want a derived identity mismatch", got.Reason)
	}
}

// An image whose manifest names a file it does not carry makes no claim, rather
// than being reported as different.
func TestMissingKeyFileIsUnverifiable(t *testing.T) {
	requireSquashfsTools(t)
	want := read(t, starIndex())

	b := starIndex()
	b.tamper = func(dir string) {
		if err := os.Remove(filepath.Join(dir, meta.RecipeFileName)); err != nil {
			t.Fatal(err)
		}
	}
	if got := Compare(want, read(t, b)); got.Verdict != Unverifiable {
		t.Errorf("verdict = %s, want unverifiable", got.Verdict)
	}
}

// Conda splits its keys where a recipe cannot: a build string moves identity and
// not equivalence, a version moves both.
func TestCondaVerdicts(t *testing.T) {
	requireSquashfsTools(t)

	condaApp := func(explicit, environ string) build {
		return build{
			name: "cutadapt/5.0", typ: catalog.TypeApp, format: "conda",
			explicit: explicit, environ: environ,
		}
	}
	const environ = "channels:\n  - conda-forge\ndependencies:\n  - cutadapt=5.0\n"
	base := condaApp("@EXPLICIT\nhttps://conda.anaconda.org/conda-forge/linux-64/cutadapt-5.0-py312h0_0.conda\n", environ)
	want := read(t, base)

	t.Run("a rebuilt package is equivalent", func(t *testing.T) {
		rebuilt := condaApp("@EXPLICIT\nhttps://conda.anaconda.org/conda-forge/linux-64/cutadapt-5.0-py312h0_1.conda\n", environ)
		got := Compare(want, read(t, rebuilt))
		if got.Verdict != Equivalent {
			t.Fatalf("verdict = %s (%s)", got.Verdict, got.Reason)
		}
		if !hasField(got.Diffs, "packages") {
			t.Errorf("diffs = %v", got.Diffs)
		}
	})

	t.Run("a new version is different", func(t *testing.T) {
		newer := condaApp(
			"@EXPLICIT\nhttps://conda.anaconda.org/conda-forge/linux-64/cutadapt-5.1-py312h0_0.conda\n",
			"channels:\n  - conda-forge\ndependencies:\n  - cutadapt=5.1\n")
		got := Compare(want, read(t, newer))
		if got.Verdict != Different {
			t.Fatalf("verdict = %s", got.Verdict)
		}
		if !hasField(got.Diffs, "environment") {
			t.Errorf("diffs = %v", got.Diffs)
		}
	})

	// The digests are over different kinds of file, so a match across formats
	// would be meaningless.
	t.Run("formats are never compared", func(t *testing.T) {
		script := build{name: "cutadapt/5.0", typ: catalog.TypeApp, recipe: "echo build\n"}
		got := Compare(want, read(t, script))
		if got.Verdict != Different || !hasField(got.Diffs, "build_type") {
			t.Errorf("verdict = %s, diffs = %v", got.Verdict, got.Diffs)
		}
	})
}

// An image carrying no keys makes no claim: unverifiable rather than different.
func TestAnUnkeyedImageIsUnverifiable(t *testing.T) {
	requireSquashfsTools(t)
	a := read(t, build{name: "ubuntu24/base", typ: catalog.TypeBase, format: "def"})
	if got := Compare(a, a); got.Verdict != Unverifiable {
		t.Errorf("verdict = %s, want unverifiable", got.Verdict)
	}
}

// A base compares like anything else: rebuilt against a moved upstream it is a
// different recorded build that still substitutes, and the diff says which.
func TestBaseComparesByDefinitionAndUpstream(t *testing.T) {
	requireSquashfsTools(t)
	def := "Bootstrap: docker\nFrom: ubuntu:24.04\n"
	january := build{name: "ubuntu24/base", typ: catalog.TypeBase, format: "def",
		recipe: def, from: "sha256:" + strings.Repeat("a", 64)}
	june := january
	june.from = "sha256:" + strings.Repeat("b", 64)

	want, got := read(t, january), read(t, june)

	if r := Compare(want, want); r.Verdict != Exact {
		t.Errorf("a base did not match itself: %s (%s)", r.Verdict, r.Reason)
	}

	r := Compare(want, got)
	if r.Verdict != Equivalent {
		t.Fatalf("verdict = %s (%s), want equivalent", r.Verdict, r.Reason)
	}
	var named bool
	for _, d := range r.Diffs {
		if d.Field == "from" {
			named = true
		}
	}
	if !named {
		t.Errorf("the moved upstream was not named: %v", r.Diffs)
	}
}

func TestMountAllowed(t *testing.T) {
	native := meta.Runtime{Platform: meta.NativePlatform()}
	if err := MountAllowed(native); err != nil {
		t.Errorf("a native image was refused: %v", err)
	}
	portable := meta.Runtime{Platform: meta.Platform{OS: "linux", Arch: meta.ArchNone}}
	if err := MountAllowed(portable); err != nil {
		t.Errorf("a noarch image was refused: %v", err)
	}
	foreign := meta.Runtime{Platform: meta.Platform{OS: "linux", Arch: "s390x"}}
	if err := MountAllowed(foreign); err == nil {
		t.Error("a foreign-architecture image was allowed to mount")
	}
}

func hasField(diffs []Diff, field string) bool {
	for _, d := range diffs {
		if d.Field == field {
			return true
		}
	}
	return false
}

func hasDiff(diffs []Diff, want Diff) bool {
	for _, d := range diffs {
		if d == want {
			return true
		}
	}
	return false
}
