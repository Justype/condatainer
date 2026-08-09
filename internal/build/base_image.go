package build

import (
	"context"
	"fmt"
	"path/filepath"

	"github.com/Justype/condatainer/catalog"
	"github.com/Justype/condatainer/internal/config"
	"github.com/Justype/condatainer/internal/image/meta"
	"github.com/Justype/condatainer/internal/runtime/apptainer"
	"github.com/Justype/condatainer/internal/utils"
)

// ResolveBase returns a base image usable as a container root, building the
// configured base first when none is installed. Script and Conda builds run
// inside it, so it is an implicit prerequisite, not a recipe #DEP:.
func ResolveBase(ctx context.Context) (string, error) {
	if path := config.FindBaseImage(); path != "" {
		if err := meta.CheckBase(path); err != nil {
			return "", err
		}
		return path, nil
	}
	if err := buildBase(ctx, false); err != nil {
		return "", err
	}
	path := config.FindBaseImage()
	if path == "" {
		return "", fmt.Errorf("base image build succeeded but installed nothing")
	}
	return path, nil
}

// RebuildBase rebuilds the configured base image, replacing an installed copy.
func RebuildBase(ctx context.Context) error {
	return buildBase(ctx, true)
}

// buildBase runs the base through the same BuildObject every other image goes
// through; only its type differs, which is what keeps the SIF.
func buildBase(ctx context.Context, update bool) error {
	obj, err := newBaseObject(ctx, update)
	if err != nil {
		return err
	}
	return obj.Build(ctx, false)
}

// resolveBase records the container root this build runs inside, resolving and
// building the configured base on first use.
func (b *BuildObject) resolveBase(ctx context.Context) error {
	if b.spec.Base != "" {
		return nil
	}
	base, err := ResolveBase(ctx)
	if err != nil {
		return fmt.Errorf("cannot build %s: %w", b.spec.Image.Name, err)
	}
	b.spec.Base = base
	return nil
}

// newBaseObject creates the BuildObject for the base image, resolved through the
// catalog like any other name and differing only in its type. update=true
// rebuilds even when one is already installed.
func newBaseObject(ctx context.Context, update bool) (*BuildObject, error) {
	if err := apptainer.EnsureApptainer(); err != nil {
		return nil, err
	}

	// The base recipe is resolved through the catalog like any other name.
	// Config `base` wins; a source's default_base fills in when it is unset,
	// which is safe here because the catalog is open either way.
	cat, err := config.OpenCatalog(ctx)
	if err != nil {
		return nil, err
	}
	config.EnsureBase(cat)
	nameVersion := config.BaseRecipeNameFrom(cat)
	if nameVersion == "" {
		return nil, fmt.Errorf("no base configured: set `base`, or configure a source declaring default_base")
	}

	// The same function FindBaseImage searches with, so the builder cannot write
	// a name the finder will not look for.
	targetOverlayPath, err := config.GetBaseImageWritePath()
	if err != nil {
		return nil, fmt.Errorf("no writable directory for base image: %w", err)
	}
	if err := utils.MkdirAllShared(filepath.Dir(targetOverlayPath)); err != nil {
		return nil, fmt.Errorf("failed to create images directory: %w", err)
	}
	if abs, err := filepath.Abs(targetOverlayPath); err == nil {
		targetOverlayPath = abs
	}

	tmpDir := resolveTmpDirForDef()
	if absDir, err := filepath.Abs(tmpDir); err == nil {
		tmpDir = absDir
	}
	base := &BuildObject{
		spec:      Spec{Image: ImageSpec{Name: nameVersion, Type: catalog.TypeBase}},
		tgt:       targetFor(targetOverlayPath),
		submitJob: false, // base image is always built locally
		update:    update,
	}

	// Resolve build source using the same mechanism as regular build objects.
	_, isContainer, err := resolveBuildSource(ctx, base, tmpDir)
	if err != nil {
		return nil, fmt.Errorf("failed to resolve base image build source: %w", err)
	}
	if !isContainer {
		return nil, fmt.Errorf("base recipe %s was not found, or is not a definition", nameVersion)
	}

	// The same siting every other definition build gets; a base differs only in
	// keeping the .sif instead of extracting it.
	base.asDefinitionBuild()

	// After capture, which takes the type from the recipe: a base is the
	// container root whatever its recipe says, and has no install prefix.
	base.spec.Image.Type = catalog.TypeBase
	base.spec.Runtime = meta.Runtime{}
	return base, nil
}
