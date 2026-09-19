package build

import (
	"context"
	"fmt"

	"github.com/Justype/condatainer/internal/config"
	"github.com/Justype/condatainer/internal/runtime/apptainer"
)

// ResolveBase returns the container root to build inside, building the
// configured default root first when none is installed. Script and Conda
// builds run inside it, so it is an implicit prerequisite, not a recipe
// #DEP:.
func ResolveBase(ctx context.Context) (string, error) {
	if path := config.FindBaseImage(); path != "" {
		return path, nil
	}
	if err := buildConfiguredRoot(ctx); err != nil {
		return "", err
	}
	path := config.FindBaseImage()
	if path == "" {
		return "", fmt.Errorf("root image build succeeded but installed nothing")
	}
	return path, nil
}

// resolveBase records the container root this build runs inside, resolving
// and building the configured default root on first use.
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

// buildConfiguredRoot builds the configured default root — config
// `default_distro`, or a source's default_distro — through the exact same
// catalog/build path every other name takes: it is an ordinary os recipe,
// not a special build type.
func buildConfiguredRoot(ctx context.Context) error {
	if _, err := apptainer.ForBuild(); err != nil {
		return err
	}

	// Config `default_distro` wins; a source's default_distro fills in when
	// it is unset.
	cat, err := config.OpenCatalog(ctx)
	if err != nil {
		return err
	}
	config.EnsureDefaultDistro(cat)
	nameVersion := config.BaseRecipeNameFrom(cat)
	if nameVersion == "" {
		return fmt.Errorf("no default distro configured: set `default_distro`, or configure a source declaring default_distro")
	}

	imagesDir, err := config.GetWritableImagesDir()
	if err != nil {
		return fmt.Errorf("no writable directory for the root image: %w", err)
	}

	obj, err := NewBuildObject(ctx, nameVersion, false, imagesDir, false)
	if err != nil {
		return fmt.Errorf("failed to resolve root recipe %s: %w", nameVersion, err)
	}
	// A root image is always built locally, never submitted to a scheduler:
	// script/Conda builds waiting on it run locally too, in the same process.
	graph, err := NewBuildGraph(ctx, []*BuildObject{obj}, imagesDir, false, false)
	if err != nil {
		return fmt.Errorf("failed to plan root build %s: %w", nameVersion, err)
	}
	return graph.Run(ctx)
}
