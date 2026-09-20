package build

import (
	"context"
	"errors"
	"path/filepath"
	"strconv"
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

func TestTryPrebuiltRejectsEquivalenceMismatch(t *testing.T) {
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

	if _, err := b.tryPrebuilt(t.Context()); !errors.Is(err, registry.ErrInvalidArtifact) {
		t.Fatalf("error = %v, want ErrInvalidArtifact", err)
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
