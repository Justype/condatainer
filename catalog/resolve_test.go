package catalog

import (
	"os"
	"path/filepath"
	"slices"
	"strings"
	"testing"
)

// depSource writes a collection with a small dependency graph:
//
//	app/1.0 -> lib/2.0 -> base-tool/1.0
//	other/1.0 -> lib/2.0            (a diamond)
//	app/1.0 -> samtools/1.23.1>=1.10 (no recipe: a fallback's problem)
func depSource(t *testing.T) Catalog {
	t.Helper()
	root := t.TempDir()
	write := func(rel, body string) {
		p := filepath.Join(root, "recipes", filepath.FromSlash(rel))
		if err := os.MkdirAll(filepath.Dir(p), 0o755); err != nil {
			t.Fatal(err)
		}
		if err := os.WriteFile(p, []byte(body), 0o644); err != nil {
			t.Fatal(err)
		}
	}
	write("app/1.0", "#DESC:app\n#DEP:lib/2.0\n#DEP:samtools/1.23.1>=1.10\n")
	write("other/1.0", "#DESC:other\n#DEP:lib/2.0\n")
	write("lib/2.0", "#DESC:lib\n#DEP:base-tool/1.0\n")
	write("lib/1.0", "#DESC:old lib\n")
	write("base-tool/1.0", "#DESC:base tool\n")

	cat, err := Open(t.Context(), []Spec{{Name: "local", Base: root}}, Cache{})
	if err != nil {
		t.Fatal(err)
	}
	return cat
}

func names(nodes []Node) []string {
	out := make([]string, len(nodes))
	for i, n := range nodes {
		out[i] = n.Name()
	}
	return out
}

func TestResolveOrder(t *testing.T) {
	cat := depSource(t)
	plan, err := cat.Resolve(t.Context(), []string{"app/1.0"}, nil)
	if err != nil {
		t.Fatal(err)
	}

	got := names(plan.Order)
	want := []string{"base-tool/1.0", "lib/2.0", "samtools/1.23.1", "app/1.0"}
	if !slices.Equal(got, want) {
		t.Fatalf("Order = %v, want %v (dependencies first)", got, want)
	}
	// Nothing is installed, so everything is missing.
	if len(plan.Missing) != len(plan.Order) {
		t.Errorf("Missing = %v, want all of Order", names(plan.Missing))
	}
}

// A name no source provides is a node with a nil Entry, not a failed walk.
func TestResolveNotProvided(t *testing.T) {
	cat := depSource(t)
	plan, err := cat.Resolve(t.Context(), []string{"app/1.0"}, nil)
	if err != nil {
		t.Fatal(err)
	}

	var conda *Node
	for i, n := range plan.Order {
		if n.Dep.Name == "samtools" {
			conda = &plan.Order[i]
		}
	}
	if conda == nil {
		t.Fatal("the unprovided dep was dropped from the plan")
	}
	if conda.Entry != nil || conda.Source != nil {
		t.Errorf("unprovided node carries an entry: %+v", conda)
	}
	// The dep survives intact, since it is what the caller hands its fallback.
	if conda.Dep.Version != "1.23.1" || conda.Dep.Op != ">=" || conda.Dep.Min != "1.10" {
		t.Errorf("dep = %+v, want the constraint preserved", conda.Dep)
	}
}

func TestResolveDiamond(t *testing.T) {
	cat := depSource(t)
	plan, err := cat.Resolve(t.Context(), []string{"app/1.0", "other/1.0"}, nil)
	if err != nil {
		t.Fatal(err)
	}

	got := names(plan.Order)
	if n := slices.Index(got, "lib/2.0"); n < 0 {
		t.Fatalf("Order = %v", got)
	}
	if c := strings.Count(strings.Join(got, " "), "lib/2.0"); c != 1 {
		t.Errorf("lib/2.0 appears %d times, want deduplicated", c)
	}
	// Every dependency still precedes its dependents.
	for _, pair := range [][2]string{
		{"base-tool/1.0", "lib/2.0"}, {"lib/2.0", "app/1.0"}, {"lib/2.0", "other/1.0"},
	} {
		if slices.Index(got, pair[0]) > slices.Index(got, pair[1]) {
			t.Errorf("%s must come before %s in %v", pair[0], pair[1], got)
		}
	}
}

func TestResolveCycle(t *testing.T) {
	root := t.TempDir()
	write := func(rel, body string) {
		p := filepath.Join(root, "recipes", filepath.FromSlash(rel))
		if err := os.MkdirAll(filepath.Dir(p), 0o755); err != nil {
			t.Fatal(err)
		}
		if err := os.WriteFile(p, []byte(body), 0o644); err != nil {
			t.Fatal(err)
		}
	}
	write("a/1.0", "#DEP:b/1.0\n")
	write("b/1.0", "#DEP:a/1.0\n")

	cat, err := Open(t.Context(), []Spec{{Name: "local", Base: root}}, Cache{})
	if err != nil {
		t.Fatal(err)
	}
	if _, err := cat.Resolve(t.Context(), []string{"a/1.0"}, nil); err == nil ||
		!strings.Contains(err.Error(), "cycle") {
		t.Errorf("err = %v, want a cycle report", err)
	}
}

// Installed satisfies a dep and stops the walk: what is already there needs no
// dependencies rebuilt beneath it.
func TestResolveHave(t *testing.T) {
	cat := depSource(t)
	have := func(name string) []string {
		if name == "lib" {
			return []string{"1.0", "2.0"}
		}
		return nil
	}
	plan, err := cat.Resolve(t.Context(), []string{"app/1.0"}, have)
	if err != nil {
		t.Fatal(err)
	}

	var lib *Node
	for i, n := range plan.Order {
		if n.Dep.Name == "lib" {
			lib = &plan.Order[i]
		}
	}
	if lib == nil || lib.Installed != "2.0" {
		t.Fatalf("lib node = %+v, want Installed 2.0", lib)
	}
	if slices.Contains(names(plan.Order), "base-tool/1.0") {
		t.Error("an installed dep should not pull its own dependencies")
	}
	if slices.Contains(names(plan.Missing), "lib/2.0") {
		t.Error("an installed dep should not be Missing")
	}
}

// With no preferred version, the newest in range wins — installed first.
func TestResolveNewestInRange(t *testing.T) {
	cat := depSource(t)

	plan, err := cat.Resolve(t.Context(), []string{"lib"}, nil)
	if err != nil {
		t.Fatal(err)
	}
	if got := names(plan.Order); !slices.Contains(got, "lib/2.0") {
		t.Errorf("Order = %v, want the newest available", got)
	}

	// An older installed copy is preferred over a newer available one, so a
	// bare dep does not rebuild whenever upstream moves.
	have := func(name string) []string {
		if name == "lib" {
			return []string{"1.0"}
		}
		return nil
	}
	plan, err = cat.Resolve(t.Context(), []string{"lib"}, have)
	if err != nil {
		t.Fatal(err)
	}
	if plan.Order[0].Installed != "1.0" || len(plan.Missing) != 0 {
		t.Errorf("node = %+v, missing = %v, want the installed 1.0 reused",
			plan.Order[0], names(plan.Missing))
	}
}

// A constraint bounds both ends: [min, preferred].
func TestResolveConstraintBounds(t *testing.T) {
	cat := depSource(t)
	have := func(string) []string { return []string{"1.09", "1.30"} }

	plan, err := cat.Resolve(t.Context(), []string{"samtools/1.23.1>=1.10"}, have)
	if err != nil {
		t.Fatal(err)
	}
	// 1.30 exceeds the preferred version and 1.09 is below the minimum, so
	// neither satisfies and the node stays unprovided.
	if plan.Order[0].Installed != "" {
		t.Errorf("Installed = %q, want neither candidate accepted", plan.Order[0].Installed)
	}
}

// A template dep is substituted from the parent's expansion before it has a
// name at all.
func TestResolveTemplateDep(t *testing.T) {
	root := t.TempDir()
	write := func(rel, body string) {
		p := filepath.Join(root, "recipes", filepath.FromSlash(rel))
		if err := os.MkdirAll(filepath.Dir(p), 0o755); err != nil {
			t.Fatal(err)
		}
		if err := os.WriteFile(p, []byte(body), 0o644); err != nil {
			t.Fatal(err)
		}
	}
	write("grch38/gtf-gencode", "#TARGET:grch38/gtf/{gencode_version}\n#PH:gencode_version:47-49\n")
	write("grch38/star-gencode",
		"#TARGET:grch38/star/{star_version}/gencode{gencode_version}\n"+
			"#PH:star_version:2.7.11b\n#PH:gencode_version:47-49\n"+
			"#DEP:grch38/gtf/{gencode_version}\n")

	cat, err := Open(t.Context(), []Spec{{Name: "local", Base: root}}, Cache{})
	if err != nil {
		t.Fatal(err)
	}
	plan, err := cat.Resolve(t.Context(), []string{"grch38/star/2.7.11b/gencode48"}, nil)
	if err != nil {
		t.Fatal(err)
	}
	got := names(plan.Order)
	if !slices.Contains(got, "grch38/gtf/48") {
		t.Errorf("Order = %v, want the dep substituted with the parent's 48", got)
	}
	if slices.Index(got, "grch38/gtf/48") != 0 {
		t.Errorf("Order = %v, want the dep first", got)
	}
}

// A base is the container root a build runs inside, not an input it is composed
// from: it may be a root of the walk, never an edge in it.
func TestResolveRejectsABaseDependency(t *testing.T) {
	root := t.TempDir()
	write := func(rel, body string) {
		p := filepath.Join(root, "recipes", filepath.FromSlash(rel))
		if err := os.MkdirAll(filepath.Dir(p), 0o755); err != nil {
			t.Fatal(err)
		}
		if err := os.WriteFile(p, []byte(body), 0o644); err != nil {
			t.Fatal(err)
		}
	}
	write("ubuntu24/base.def", "#DESC:root\nBootstrap: docker\n")
	write("ubuntu24/r-essential.def", "#DESC:os layer\nBootstrap: docker\n")
	write("grch38/gtf/49", "#DESC:annotation\n#DEP:ubuntu24/base\n")
	write("grch38/vcf/49", "#DESC:variants\n#DEP:ubuntu24/r-essential\n")

	cat, err := Open(t.Context(), []Spec{{Name: "local", Base: root}}, Cache{})
	if err != nil {
		t.Fatal(err)
	}

	if _, err := cat.Resolve(t.Context(), []string{"grch38/gtf/49"}, nil); err == nil {
		t.Error("a base was accepted as a dependency")
	} else if !strings.Contains(err.Error(), "ubuntu24/base") {
		t.Errorf("err = %v, want it to name the base", err)
	}

	// Building the base itself is exactly what a root is for.
	if _, err := cat.Resolve(t.Context(), []string{"ubuntu24/base"}, nil); err != nil {
		t.Errorf("building a base was rejected: %v", err)
	}

	// An OS layer is a legal dependency; only base is refused.
	if _, err := cat.Resolve(t.Context(), []string{"grch38/vcf/49"}, nil); err != nil {
		t.Errorf("an os dependency was rejected: %v", err)
	}
}
