package build

import (
	"os"
	"path/filepath"
	"testing"

	"github.com/Justype/condatainer/internal/config"
)

// withInstalledBase installs a base image and points every image search at it.
// Its bytes are never read: the callers under test resolve a path.
func withInstalledBase(t *testing.T) string {
	t.Helper()
	dir := filepath.Join(t.TempDir(), "images")
	if err := os.MkdirAll(dir, 0o755); err != nil {
		t.Fatal(err)
	}
	base := filepath.Join(dir, "ubuntu24--base.sif")
	if err := os.WriteFile(base, []byte("SIF"), 0o644); err != nil {
		t.Fatal(err)
	}

	prevPaths, prevBase := config.GlobalDataPaths, config.Global.Base
	config.GlobalDataPaths.ImagesDirs = []string{dir}
	config.Global.Base = "ubuntu24"
	t.Cleanup(func() { config.GlobalDataPaths, config.Global.Base = prevPaths, prevBase })
	return base
}

// A definition bootstraps its own root from the definition itself, so a plan of
// only definitions — the base's own build among them — resolves no base. This
// is what lets the base be built when none exists yet.
func TestGraphSkipsBaseForDefinitionOnlyPlans(t *testing.T) {
	// No base configured and none installed, so any resolution attempt fails.
	prevBase := config.Global.Base
	config.Global.Base = ""
	t.Cleanup(func() { config.Global.Base = prevBase })

	def := &BuildObject{spec: Spec{Image: ImageSpec{Name: "ubuntu24/base"}}, buildType: BuildTypeDef}
	bg := &BuildGraph{graph: map[string]*BuildObject{def.NameVersion(): def}}

	if err := bg.resolveBase(t.Context()); err != nil {
		t.Fatalf("a definition-only plan tried to resolve a base: %v", err)
	}
	if def.spec.Base != "" {
		t.Errorf("Base = %q, want empty — a definition has no base to run inside", def.spec.Base)
	}
}

// Script and Conda builds run their install and their packing inside the base,
// so the graph records it on them before any node runs.
func TestGraphRecordsBaseOnDependents(t *testing.T) {
	base := withInstalledBase(t)

	script := &BuildObject{spec: Spec{Image: ImageSpec{Name: "samtools/1.23.1"}}, buildType: BuildTypeScript}
	conda := &BuildObject{spec: Spec{Image: ImageSpec{Name: "numpy/2.1.0"}}, buildType: BuildTypeConda}
	bg := &BuildGraph{
		update: true, // treat both as missing without touching the installed set
		graph: map[string]*BuildObject{
			script.NameVersion(): script,
			conda.NameVersion():  conda,
		},
	}

	if err := bg.resolveBase(t.Context()); err != nil {
		t.Fatalf("resolveBase: %v", err)
	}
	for _, obj := range []*BuildObject{script, conda} {
		if obj.spec.Base != base {
			t.Errorf("%s Base = %q, want %q", obj.NameVersion(), obj.spec.Base, base)
		}
	}
}

// A base already chosen for the build is kept: the graph resolves once, and a
// backend re-checking must not pick a different root mid-build.
func TestResolveBaseKeepsAnAlreadyChosenBase(t *testing.T) {
	withInstalledBase(t)

	b := &BuildObject{spec: Spec{Image: ImageSpec{Name: "samtools/1.23.1"}}, buildType: BuildTypeScript}
	b.spec.Base = "/somewhere/else.sif"

	if err := b.resolveBase(t.Context()); err != nil {
		t.Fatalf("resolveBase: %v", err)
	}
	if b.spec.Base != "/somewhere/else.sif" {
		t.Errorf("Base = %q, want the one already chosen", b.spec.Base)
	}
}
