package build

import (
	"context"
	"path/filepath"
	"strconv"
	"strings"
	"testing"

	ocispec "github.com/opencontainers/image-spec/specs-go/v1"

	"github.com/Justype/condatainer/catalog"
	"github.com/Justype/condatainer/internal/artifact/meta"
	"github.com/Justype/condatainer/internal/config"
	"github.com/Justype/condatainer/internal/registry"
)

func prebuiltObject(t *testing.T, recipe []byte, endpoints ...string) *BuildObject {
	t.Helper()
	b := &BuildObject{
		spec: Spec{
			Image: ImageSpec{Name: "demo/1", Type: catalog.TypeApp},
			Source: SourceSpec{Script: &ScriptSource{File: SourceFile{
				Name: meta.RecipeFileName, Data: recipe,
			}}},
		},
		tgt: targetFor(filepath.Join(t.TempDir(), "demo--1.sqf")),
		catalogSource: &catalog.Source{Desc: catalog.Descriptor{OCI: catalog.OCI{
			Pull: endpoints, Audience: "restricted",
		}}},
	}
	b.embedSource(SourceFile{Name: meta.RecipeFileName, Data: recipe})
	return b
}

func TestTryPrebuiltUsesOrderedEndpointAndEquivalence(t *testing.T) {
	b := prebuiltObject(t, []byte("echo demo\n"), "mirror.invalid/lab", "origin.invalid/lab")
	want, err := b.prebuiltEquivalence(t.Context())
	if err != nil {
		t.Fatal(err)
	}

	oldResolve, oldPull := resolvePrebuilt, pullPrebuilt
	t.Cleanup(func() { resolvePrebuilt, pullPrebuilt = oldResolve, oldPull })
	var resolved []string
	resolvePrebuilt = func(_ context.Context, base, repo, tag string) (ocispec.Descriptor, map[string]string, error) {
		resolved = append(resolved, base+"/"+repo+":"+tag)
		if len(resolved) == 1 {
			return ocispec.Descriptor{}, nil, registry.ErrNotFound
		}
		return ocispec.Descriptor{}, map[string]string{
			registry.AnnTitle:       "demo/1",
			registry.AnnSchema:      strconv.Itoa(meta.SchemaVersion),
			registry.AnnEquivScheme: want.Scheme,
			registry.AnnEquivSHA:    want.SHA256,
		}, nil
	}
	pulled := false
	pullPrebuilt = func(_ context.Context, base, repo string, _ ocispec.Descriptor, _ map[string]string, dest string, _ registry.Kind) error {
		pulled = base == "origin.invalid/lab" && repo == "demo" && dest == b.tgt.Path
		return nil
	}

	got, err := b.tryPrebuilt(t.Context())
	if err != nil {
		t.Fatal(err)
	}
	if got != prebuiltResult(true) || !pulled || len(resolved) != 2 {
		t.Fatalf("result=%v pulled=%v resolved=%v", got, pulled, resolved)
	}
}

func TestTryPrebuiltSkipsEquivalenceMismatch(t *testing.T) {
	b := prebuiltObject(t, []byte("echo current\n"), "registry.invalid/lab")
	oldResolve, oldPull := resolvePrebuilt, pullPrebuilt
	t.Cleanup(func() { resolvePrebuilt, pullPrebuilt = oldResolve, oldPull })
	resolvePrebuilt = func(context.Context, string, string, string) (ocispec.Descriptor, map[string]string, error) {
		return ocispec.Descriptor{}, map[string]string{
			registry.AnnTitle:       "demo/1",
			registry.AnnSchema:      strconv.Itoa(meta.SchemaVersion),
			registry.AnnEquivScheme: "script-equiv-v1",
			registry.AnnEquivSHA:    "different",
		}, nil
	}
	pullPrebuilt = func(context.Context, string, string, ocispec.Descriptor, map[string]string, string, registry.Kind) error {
		t.Fatal("pull called for mismatched candidate")
		return nil
	}

	if pulled, err := b.tryPrebuilt(t.Context()); err != nil || pulled {
		t.Fatalf("result=%v error=%v, want a local build", pulled, err)
	}
}

func TestTryPrebuiltFallsBackWhenEndpointsUnavailable(t *testing.T) {
	b := prebuiltObject(t, []byte("echo demo\n"), "mirror.invalid/lab", "origin.invalid/lab")
	oldResolve := resolvePrebuilt
	t.Cleanup(func() { resolvePrebuilt = oldResolve })
	resolvePrebuilt = func(context.Context, string, string, string) (ocispec.Descriptor, map[string]string, error) {
		return ocispec.Descriptor{}, nil, registry.ErrUnavailable
	}

	if pulled, err := b.tryPrebuilt(t.Context()); err != nil || pulled {
		t.Fatalf("result=%v error=%v", pulled, err)
	}
}

// --no-prebuilt reaches no registry: nothing is resolved and the build carries on.
func TestTryPrebuiltSkippedByNoPrebuilt(t *testing.T) {
	b := prebuiltObject(t, []byte("echo demo\n"), "origin.invalid/lab")

	oldResolve := resolvePrebuilt
	t.Cleanup(func() { resolvePrebuilt = oldResolve })
	resolvePrebuilt = func(context.Context, string, string, string) (ocispec.Descriptor, map[string]string, error) {
		t.Error("a registry was consulted under --no-prebuilt")
		return ocispec.Descriptor{}, nil, registry.ErrNotFound
	}
	prev := config.Global.Build.SkipPrebuilt
	config.Global.Build.SkipPrebuilt = true
	t.Cleanup(func() { config.Global.Build.SkipPrebuilt = prev })

	got, err := b.tryPrebuilt(t.Context())
	if err != nil || got != prebuiltResult(false) {
		t.Fatalf("result=%v err=%v, want a local build", got, err)
	}
}

// annotationsFor is what a registry serves for an artifact built with equiv.
func annotationsFor(equiv meta.KeyRef) map[string]string {
	return map[string]string{
		registry.AnnTitle:       "demo/1",
		registry.AnnSchema:      strconv.Itoa(meta.SchemaVersion),
		registry.AnnEquivScheme: equiv.Scheme,
		registry.AnnEquivSHA:    equiv.SHA256,
	}
}

func TestPlanPrebuiltDecides(t *testing.T) {
	oldResolve := resolvePrebuilt
	t.Cleanup(func() { resolvePrebuilt = oldResolve })

	tests := []struct {
		name      string
		published func(want meta.KeyRef) meta.KeyRef
		known     bool
		want      prebuiltChoice
	}{
		{"matching artifact is pulled", func(w meta.KeyRef) meta.KeyRef { return w }, true, prebuiltPull},
		{"other recipe is built", func(meta.KeyRef) meta.KeyRef { return meta.KeyRef{Scheme: "script-equiv-v1", SHA256: "other"} }, true, prebuiltNone},
		{"pending dependency stays undecided", func(w meta.KeyRef) meta.KeyRef { return w }, false, prebuiltUndecided},
	}
	for _, tt := range tests {
		t.Run(tt.name, func(t *testing.T) {
			b := prebuiltObject(t, []byte("echo demo\n"), "origin.invalid/lab")
			want, err := b.prebuiltEquivalence(t.Context())
			if err != nil {
				t.Fatal(err)
			}
			resolvePrebuilt = func(context.Context, string, string, string) (ocispec.Descriptor, map[string]string, error) {
				return ocispec.Descriptor{}, annotationsFor(tt.published(want)), nil
			}
			if err := b.planPrebuilt(t.Context(), tt.known); err != nil {
				t.Fatal(err)
			}
			if b.prebuilt.choice != tt.want {
				t.Fatalf("choice = %v, want %v", b.prebuilt.choice, tt.want)
			}
		})
	}
}

// A planned node reaches the registry for its pull and no further: the
// candidate was found when planning ran, and a build plan skips the lookup.
func TestTryPrebuiltFollowsThePlan(t *testing.T) {
	b := prebuiltObject(t, []byte("echo demo\n"), "origin.invalid/lab")
	want, err := b.prebuiltEquivalence(t.Context())
	if err != nil {
		t.Fatal(err)
	}
	oldResolve, oldPull := resolvePrebuilt, pullPrebuilt
	t.Cleanup(func() { resolvePrebuilt, pullPrebuilt = oldResolve, oldPull })
	resolves := 0
	resolvePrebuilt = func(context.Context, string, string, string) (ocispec.Descriptor, map[string]string, error) {
		resolves++
		return ocispec.Descriptor{}, annotationsFor(want), nil
	}
	pulls := 0
	pullPrebuilt = func(context.Context, string, string, ocispec.Descriptor, map[string]string, string, registry.Kind) error {
		pulls++
		return nil
	}

	if err := b.planPrebuilt(t.Context(), true); err != nil {
		t.Fatal(err)
	}
	if got, err := b.tryPrebuilt(t.Context()); err != nil || got != prebuiltResult(true) || resolves != 1 || pulls != 1 {
		t.Fatalf("planned pull: result=%v err=%v resolves=%d pulls=%d", got, err, resolves, pulls)
	}

	b.prebuilt = prebuiltPlan{choice: prebuiltNone}
	if got, err := b.tryPrebuilt(t.Context()); err != nil || got || resolves != 1 || pulls != 1 {
		t.Fatalf("planned build: result=%v err=%v resolves=%d pulls=%d", got, err, resolves, pulls)
	}
}

func planNode(name string, typ catalog.Type, choice prebuiltChoice, deps ...string) *BuildObject {
	return &BuildObject{
		spec:     Spec{Image: ImageSpec{Name: name, Type: typ}, Dependencies: deps},
		prebuilt: prebuiltPlan{choice: choice},
	}
}

func TestDependencyPlan(t *testing.T) {
	key := meta.KeyRef{Scheme: "script-equiv-v1", SHA256: "abc"}
	keyed := planNode("data/1", catalog.TypeData, prebuiltUndecided)
	keyed.plannedEquiv = key
	bg := &BuildGraph{graph: map[string]*BuildObject{
		"data/1": keyed,
		"data/2": planNode("data/2", catalog.TypeData, prebuiltUndecided), // equivalence not derived yet
		"tool/1": planNode("tool/1", catalog.TypeApp, prebuiltUndecided),
	}}

	planned, ok := bg.dependencyPlan(planNode("top/1", catalog.TypeData, prebuiltUndecided, "data/1", "tool/1"))
	if !ok || planned["data/1"].Equiv != key || planned["tool/1"].Type != catalog.TypeApp {
		t.Errorf("planned = %+v, ok = %v", planned, ok)
	}
	if _, ok := bg.dependencyPlan(planNode("top/1", catalog.TypeData, prebuiltUndecided, "data/2")); ok {
		t.Error("a data dependency with no derived equivalence counted as known")
	}
	if _, ok := bg.dependencyPlan(planNode("top/1", catalog.TypeData, prebuiltUndecided, "elsewhere/1")); ok {
		t.Error("a dependency neither installed nor planned counted as known")
	}
}

func TestPruneDependencies(t *testing.T) {
	build := func(topChoice prebuiltChoice) (*BuildGraph, []*BuildObject, *BuildObject) {
		tool := planNode("tool/1", catalog.TypeApp, prebuiltUndecided)
		data := planNode("data/1", catalog.TypeData, prebuiltUndecided)
		top := planNode("top/1", catalog.TypeData, topChoice, "tool/1", "data/1")
		graph := map[string]*BuildObject{"tool/1": tool, "data/1": data, "top/1": top}
		return &BuildGraph{graph: graph}, []*BuildObject{tool, data, top}, top
	}

	bg, order, top := build(prebuiltPull)
	if kept := bg.pruneDependencies(order, []string{"top/1"}); len(kept) != 1 || kept[0] != top ||
		len(top.prunedDeps) != 2 || !order[0].pruned || top.pruned {
		t.Errorf("a pulled root kept %d nodes, pruned %v", len(kept), top.prunedDeps)
	}

	bg, order, top = build(prebuiltNone)
	if kept := bg.pruneDependencies(order, []string{"top/1"}); len(kept) != 3 || len(top.prunedDeps) != 0 {
		t.Errorf("a built root kept %d nodes, pruned %v", len(kept), top.prunedDeps)
	}

	// A dependency another root asked for stays, and is not reported as left out.
	bg, order, top = build(prebuiltPull)
	if kept := bg.pruneDependencies(order, []string{"top/1", "tool/1"}); len(kept) != 2 ||
		len(top.prunedDeps) != 1 || top.prunedDeps[0] != "data/1" {
		t.Errorf("kept %d nodes, pruned %v", len(kept), top.prunedDeps)
	}
}

// A dependency not installed yet enters the equivalence by its planned key.
func TestPlannedDependencyEntersEquivalence(t *testing.T) {
	keyOf := func(sha string) meta.KeyRef {
		b := prebuiltObject(t, []byte("echo demo\n"), "origin.invalid/lab")
		b.spec.Image.Type = catalog.TypeData
		b.spec.Dependencies = []string{"data/1"}
		b.plannedDeps = map[string]plannedDep{
			"data/1": {Name: "data/1", Type: catalog.TypeData, Equiv: meta.KeyRef{Scheme: "script-equiv-v1", SHA256: sha}},
		}
		got, err := b.prebuiltEquivalence(t.Context())
		if err != nil {
			t.Fatal(err)
		}
		return got
	}
	if keyOf(strings.Repeat("a", 64)) != keyOf(strings.Repeat("a", 64)) {
		t.Error("the same planned dependency gave two keys")
	}
	if keyOf(strings.Repeat("a", 64)) == keyOf(strings.Repeat("b", 64)) {
		t.Error("a different planned dependency key gave the same equivalence")
	}
}
