package key

import (
	"strings"
	"testing"

	"github.com/Justype/condatainer/catalog"
	"github.com/Justype/condatainer/internal/artifact/meta"
)

func digestOf(s string) string { return Digest([]byte(s)) }

// starIndex is the worked example: a STAR index whose name mentions the tool
// that produced it, built against a GTF dataset and two apps.
func starIndex() Artifact {
	return Artifact{
		Name:         "grch38/star/2.7.11b/gencode49-101",
		Type:         catalog.TypeData,
		Env:          []meta.EnvVar{{Key: "STAR_INDEX_DIR", Value: "{prefix}", Note: "for --genomeDir"}},
		Recipe:       []byte("#DESC:x\nSTAR --runMode genomeGenerate\n"),
		Placeholders: map[string]string{"star_version": "2.7.11b", "read_length": "101"},
		Deps: []Dep{
			{Name: "star/2.7.11b", Type: catalog.TypeApp, Identity: digestOf("star id"), Equiv: digestOf("star eq")},
			{Name: "samtools/1.23.1", Type: catalog.TypeApp, Identity: digestOf("sam id"), Equiv: digestOf("sam eq")},
			{Name: "grch38/gtf-gencode/49", Type: catalog.TypeData, Identity: digestOf("gtf id"), Equiv: digestOf("gtf eq")},
		},
	}
}

func lines(t *testing.T, r Model) []string {
	t.Helper()
	data, err := Marshal(r)
	if err != nil {
		t.Fatalf("Marshal: %v", err)
	}
	return strings.Split(strings.TrimSuffix(string(data), "\n"), "\n")
}

func has(got []string, want string) bool {
	for _, line := range got {
		if line == want {
			return true
		}
	}
	return false
}

// deriveModelsForTest exposes the models selected by the individual scheme
// implementations so their contribution rules can be tested field by field.
func deriveModelsForTest(a Artifact) (Model, Model, error) {
	switch a.Type {
	case catalog.TypeOS, catalog.TypeBase:
		_, identity, err := deriveDefinitionIdentityV1(a)
		if err != nil {
			return Model{}, Model{}, err
		}
		_, equiv, err := deriveDefinitionEquivV1(a)
		return identity, equiv, err
	default:
		_, identity, err := deriveScriptIdentityV1(a)
		if err != nil {
			return Model{}, Model{}, err
		}
		_, equiv, err := deriveScriptEquivV1(a)
		return identity, equiv, err
	}
}

// Identity pins every direct dependency; equivalence keeps only what decides
// substitution. That difference is the reason there are two schemes.
func TestSchemeRulesProjectDependencies(t *testing.T) {
	a := starIndex()
	identity, equiv, err := deriveModelsForTest(a)
	if err != nil {
		t.Fatalf("derive schemes: %v", err)
	}

	id := lines(t, identity)
	for _, want := range []string{
		"dep=app star/2.7.11b " + digestOf("star id"),
		"dep=app samtools/1.23.1 " + digestOf("sam id"),
		"dep=data grch38/gtf-gencode/49 " + digestOf("gtf id"),
	} {
		if !has(id, want) {
			t.Errorf("identity missing %q:\n%s", want, strings.Join(id, "\n"))
		}
	}

	eq := lines(t, equiv)
	// STAR is in the name, so its version is the contract — and nothing else
	// about it, not its identity and not its own equivalence.
	if !has(eq, "dep=app star/2.7.11b") {
		t.Errorf("equiv missing the named app:\n%s", strings.Join(eq, "\n"))
	}
	// Samtools is not in the name: mounted, recorded in identity, and silent here.
	for _, line := range eq {
		if strings.Contains(line, "samtools") {
			t.Errorf("a history-only dependency reached equivalence: %q", line)
		}
	}
	// A data dependency contributes its equivalence and no name.
	if !has(eq, "dep=data "+digestOf("gtf eq")) {
		t.Errorf("equiv missing the data dependency:\n%s", strings.Join(eq, "\n"))
	}
	for _, line := range eq {
		if strings.Contains(line, "gtf-gencode") {
			t.Errorf("a data dependency's name reached equivalence: %q", line)
		}
	}
}

// The two models share every field that describes the recipe. What one asks of a
// dependency differs; what the recipe is does not.
func TestSchemeRulesShareTheirSourceLines(t *testing.T) {
	identity, equiv, err := deriveModelsForTest(starIndex())
	if err != nil {
		t.Fatal(err)
	}
	if identity.Recipe != equiv.Recipe {
		t.Error("the two models disagree about the recipe")
	}
	if len(identity.Placeholders) != len(equiv.Placeholders) {
		t.Error("the two models disagree about the placeholders")
	}
	if len(identity.Env) != 1 || identity.Env[0].Key != "STAR_INDEX_DIR" {
		t.Errorf("env = %v", identity.Env)
	}
	// No canonical model may carry the artifact's name or its prefix.
	for _, r := range []Model{identity, equiv} {
		for _, line := range lines(t, r) {
			if strings.Contains(line, "gencode49-101") && !strings.HasPrefix(line, "ph=") {
				t.Errorf("the artifact's name reached a model: %q", line)
			}
			if strings.HasPrefix(line, "name=") || strings.HasPrefix(line, "prefix=") {
				t.Errorf("a forbidden field reached a model: %q", line)
			}
		}
	}
}

// Nothing may fail because a dependency predates this format. An unrecorded one
// is stable and comparable; it simply carries a weaker claim, said out loud.
func TestUnrecordedDependencies(t *testing.T) {
	a := starIndex()
	a.Deps = []Dep{
		{Name: "star/2.7.11b", Type: catalog.TypeApp},           // named, no records
		{Name: "samtools/1.23.1", Type: catalog.TypeApp},        // history, no records
		{Name: "grch38/gtf-gencode/49", Type: catalog.TypeData}, // data, no records
	}
	identity, equiv, err := deriveModelsForTest(a)
	if err != nil {
		t.Fatalf("derive schemes: %v", err)
	}

	id := lines(t, identity)
	for _, want := range []string{
		"dep=app star/2.7.11b unrecorded",
		"dep=app samtools/1.23.1 unrecorded",
		"dep=data grch38/gtf-gencode/49 unrecorded",
	} {
		if !has(id, want) {
			t.Errorf("identity missing %q:\n%s", want, strings.Join(id, "\n"))
		}
	}

	eq := lines(t, equiv)
	// A named app can still state its version honestly.
	if !has(eq, "dep=app star/2.7.11b") {
		t.Errorf("equiv lost a named app:\n%s", strings.Join(eq, "\n"))
	}
	// A data dependency with no digest says which one it was instead.
	if !has(eq, "dep=data unrecorded grch38/gtf-gencode/49") {
		t.Errorf("equiv missing the unrecorded data dependency:\n%s", strings.Join(eq, "\n"))
	}

	// Two builds against the same unrecorded dependency agree.
	again, _, err := deriveModelsForTest(a)
	if err != nil {
		t.Fatal(err)
	}
	first, _ := Key(identity)
	second, _ := Key(again)
	if first != second {
		t.Error("an unrecorded dependency produced unstable records")
	}

	deps, complete := Manifest(a)
	if complete == nil || *complete {
		t.Error("provenance_complete should be false with an unrecorded dependency")
	}
	for _, d := range deps {
		if d.Records != meta.Unrecorded {
			t.Errorf("%s: records = %q", d.Name, d.Records)
		}
	}
}

// What moves each key is the design in one test.
func TestWhatMovesEachKey(t *testing.T) {
	keys := func(a Artifact) (string, string) {
		t.Helper()
		identity, equiv, err := deriveModelsForTest(a)
		if err != nil {
			t.Fatalf("derive schemes: %v", err)
		}
		i, err := Key(identity)
		if err != nil {
			t.Fatal(err)
		}
		e, err := Key(equiv)
		if err != nil {
			t.Fatal(err)
		}
		return i, e
	}
	baseID, baseEq := keys(starIndex())

	t.Run("a history-only dependency moves identity alone", func(t *testing.T) {
		a := starIndex()
		a.Deps[1].Identity = digestOf("samtools rebuilt")
		a.Deps[1].Equiv = digestOf("samtools rebuilt eq")
		id, eq := keys(a)
		if id == baseID {
			t.Error("identity did not move")
		}
		if eq != baseEq {
			t.Error("equivalence moved on a history-only dependency")
		}
	})

	t.Run("rebuilding a named app at the same version moves identity alone", func(t *testing.T) {
		a := starIndex()
		a.Deps[0].Identity = digestOf("star rebuilt")
		a.Deps[0].Equiv = digestOf("star rebuilt eq")
		id, eq := keys(a)
		if id == baseID {
			t.Error("identity did not move")
		}
		if eq != baseEq {
			t.Error("equivalence moved when a named app was rebuilt at one version")
		}
	})

	t.Run("a data dependency's equivalence moves both", func(t *testing.T) {
		a := starIndex()
		a.Deps[2].Identity = digestOf("gtf id 2")
		a.Deps[2].Equiv = digestOf("gtf eq 2")
		id, eq := keys(a)
		if id == baseID || eq == baseEq {
			t.Error("a data dependency change did not move both keys")
		}
	})

	t.Run("renaming a data dependency moves neither", func(t *testing.T) {
		a := starIndex()
		a.Deps[2].Name = "grch38/annotation/49"
		_, eq := keys(a)
		if eq != baseEq {
			t.Error("equivalence followed a data dependency's name")
		}
	})

	t.Run("the recipe body moves both", func(t *testing.T) {
		a := starIndex()
		a.Recipe = []byte("#DESC:x\nSTAR --runMode genomeGenerate --sjdbOverhang 100\n")
		id, eq := keys(a)
		if id == baseID || eq == baseEq {
			t.Error("a recipe change did not move both keys")
		}
	})

	t.Run("a comment moves neither", func(t *testing.T) {
		a := starIndex()
		a.Recipe = []byte("#DESC:reworded\n# explain it\nSTAR --runMode genomeGenerate\n")
		id, eq := keys(a)
		if id != baseID || eq != baseEq {
			t.Error("a comment-only edit moved a key")
		}
	})

	t.Run("an env value moves both", func(t *testing.T) {
		a := starIndex()
		a.Env[0].Value = "{prefix}/index"
		id, eq := keys(a)
		if id == baseID || eq == baseEq {
			t.Error("an #ENV: change did not move both keys")
		}
	})

	t.Run("an env note moves neither", func(t *testing.T) {
		a := starIndex()
		a.Env[0].Note = "reworded"
		id, eq := keys(a)
		if id != baseID || eq != baseEq {
			t.Error("rewording an #ENV: note moved a key")
		}
	})
}

// The name decides which tools count, so it is a gate on the projection rather
// than a field. Component matching is exact and contiguous.
func TestProjectionFollowsTheName(t *testing.T) {
	tests := []struct {
		name     string
		artifact string
		dep      string
		wantRole string
	}{
		{"clean components", "grch38/star/2.7.11b/gencode49", "star/2.7.11b", meta.RoleApp},
		{"glued version misses", "grch38/star2.7.11b/gencode49", "star/2.7.11b", meta.RoleHistory},
		{"unmentioned tool", "grch38/genome/gencode", "samtools/1.23.1", meta.RoleHistory},
		{"an OS name carries its distro", "grch38/pytorch/2.9/emb", "ubuntu24/pytorch/2.9", meta.RoleHistory},
		{"the distro spelled out", "grch38/ubuntu24/pytorch/2.9/emb", "ubuntu24/pytorch/2.9", meta.RoleApp},
	}
	for _, tt := range tests {
		t.Run(tt.name, func(t *testing.T) {
			d := Dep{Name: tt.dep, Type: catalog.TypeApp, Identity: digestOf("x"), Equiv: digestOf("y")}
			if got := Role(tt.artifact, d); got != tt.wantRole {
				t.Errorf("Role = %q, want %q", got, tt.wantRole)
			}
			a := Artifact{Name: tt.artifact, Type: catalog.TypeData, Recipe: []byte("echo\n"), Deps: []Dep{d}}
			_, equiv, err := deriveModelsForTest(a)
			if err != nil {
				t.Fatal(err)
			}
			inEquiv := len(equiv.Deps) > 0
			if want := tt.wantRole == meta.RoleApp; inEquiv != want {
				t.Errorf("in equivalence = %v, want %v", inEquiv, want)
			}
		})
	}
}

// An OS dependency projects exactly like an app; only the type on the line
// differs.
func TestOSDependencyProjectsLikeAnApp(t *testing.T) {
	a := Artifact{
		Name:   "grch38/ubuntu24/pytorch/2.9/embeddings",
		Type:   catalog.TypeData,
		Recipe: []byte("echo\n"),
		Deps: []Dep{
			{Name: "ubuntu24/pytorch/2.9", Type: catalog.TypeOS, Identity: digestOf("os id"), Equiv: digestOf("os eq")},
		},
	}
	identity, equiv, err := deriveModelsForTest(a)
	if err != nil {
		t.Fatal(err)
	}
	if !has(lines(t, identity), "dep=os ubuntu24/pytorch/2.9 "+digestOf("os id")) {
		t.Errorf("identity:\n%s", strings.Join(lines(t, identity), "\n"))
	}
	if !has(lines(t, equiv), "dep=os ubuntu24/pytorch/2.9") {
		t.Errorf("equiv:\n%s", strings.Join(lines(t, equiv), "\n"))
	}
}

func TestSchemeRulesRejectInvalidSubject(t *testing.T) {
	app := Artifact{
		Name: "samtools/1.23.1", Type: catalog.TypeApp, Recipe: []byte("echo\n"),
		Deps: []Dep{{Name: "zlib/1.3", Type: catalog.TypeApp, Identity: digestOf("z")}},
	}
	if _, _, err := deriveModelsForTest(app); err == nil {
		t.Error("an app was given dependencies")
	}
}

// A dependency pinned in one record and not the other would be a half-claim. The
// pair is all-or-nothing, so an image offering one key counts as unrecorded.
func TestHalfRecordedDependencyIsUnrecorded(t *testing.T) {
	a := starIndex()
	a.Deps = []Dep{{Name: "grch38/gtf-gencode/49", Type: catalog.TypeData, Identity: digestOf("only identity")}}

	identity, equiv, err := deriveModelsForTest(a)
	if err != nil {
		t.Fatalf("derive schemes: %v", err)
	}
	if !has(lines(t, identity), "dep=data grch38/gtf-gencode/49 unrecorded") {
		t.Errorf("identity:\n%s", strings.Join(lines(t, identity), "\n"))
	}
	if !has(lines(t, equiv), "dep=data unrecorded grch38/gtf-gencode/49") {
		t.Errorf("equiv:\n%s", strings.Join(lines(t, equiv), "\n"))
	}

	// Whatever was known still reaches the manifest, marked for what it is.
	deps, complete := Manifest(a)
	if len(deps) != 1 || deps[0].Records != meta.Unrecorded {
		t.Errorf("dependencies = %+v", deps)
	}
	if deps[0].Identity != digestOf("only identity") {
		t.Error("the manifest dropped the one key that was known")
	}
	if complete == nil || *complete {
		t.Error("provenance_complete should be false")
	}
}

// A base is identified like any other definition build; what makes it special is
// only that nothing may depend on it.
func TestBaseIsIdentifiedLikeAnOS(t *testing.T) {
	upstream := digestOf("ubuntu:24.04")
	a := Artifact{
		Name:   "ubuntu24/base",
		Type:   catalog.TypeBase,
		Recipe: []byte("Bootstrap: docker\nFrom: ubuntu:24.04\n"),
		From:   upstream,
	}
	identity, equiv, err := deriveModelsForTest(a)
	if err != nil {
		t.Fatalf("a base was refused records: %v", err)
	}
	if !has(lines(t, identity), "type=base") || !has(lines(t, equiv), "type=base") {
		t.Errorf("identity:\n%s", strings.Join(lines(t, identity), "\n"))
	}
	if !has(lines(t, identity), "from="+upstream) {
		t.Errorf("the upstream digest is missing from identity:\n%s", strings.Join(lines(t, identity), "\n"))
	}
}

// An upstream security push makes a different recorded build that still
// substitutes, which is what asking for a tag meant.
func TestUpstreamMovesIdentityAndNotEquivalence(t *testing.T) {
	recipe := []byte("Bootstrap: docker\nFrom: ubuntu:24.04\n")
	january := Artifact{Name: "ubuntu24/os", Type: catalog.TypeOS, Recipe: recipe, From: digestOf("january")}
	june := january
	june.From = digestOf("june")

	janID, janEq, err := deriveModelsForTest(january)
	if err != nil {
		t.Fatal(err)
	}
	junID, junEq, err := deriveModelsForTest(june)
	if err != nil {
		t.Fatal(err)
	}

	idA, err := Key(janID)
	if err != nil {
		t.Fatal(err)
	}
	idB, err := Key(junID)
	if err != nil {
		t.Fatal(err)
	}
	if idA == idB {
		t.Error("two upstreams produced one identity")
	}

	eqA, err := Key(janEq)
	if err != nil {
		t.Fatal(err)
	}
	eqB, err := Key(junEq)
	if err != nil {
		t.Fatal(err)
	}
	if eqA != eqB {
		t.Error("an upstream digest reached the equivalence key")
	}
}

// An unreachable registry is a weaker claim, not a failed build — and distinct
// from a definition with no upstream at all, which writes no line.
func TestUnresolvedUpstreamIsRecordedAsUnrecorded(t *testing.T) {
	recipe := []byte("Bootstrap: docker\nFrom: ubuntu:24.04\n")
	unresolved, _, err := deriveModelsForTest(Artifact{Name: "u/os", Type: catalog.TypeOS, Recipe: recipe, From: meta.Unrecorded})
	if err != nil {
		t.Fatal(err)
	}
	if !has(lines(t, unresolved), "from="+meta.Unrecorded) {
		t.Errorf("identity:\n%s", strings.Join(lines(t, unresolved), "\n"))
	}

	none, _, err := deriveModelsForTest(Artifact{Name: "u/os", Type: catalog.TypeOS, Recipe: recipe})
	if err != nil {
		t.Fatal(err)
	}
	for _, line := range lines(t, none) {
		if strings.HasPrefix(line, "from=") {
			t.Errorf("a definition with no upstream wrote %q", line)
		}
	}
}
