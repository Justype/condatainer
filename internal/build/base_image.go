package build

import (
	"context"
	"fmt"
	"os"
	"path/filepath"
	"strings"

	"github.com/Justype/condatainer/internal/apptainer"
	"github.com/Justype/condatainer/internal/config"
	"github.com/Justype/condatainer/internal/logging"
	"github.com/Justype/condatainer/internal/utils"
)

// BaseImageBuildObject wraps BuildObject but overrides IsInstalled() and Build()
// to search all configured image paths and skip SquashFS extraction.
type BaseImageBuildObject struct {
	*BuildObject
}

// EnsureBaseImage is the fast-path entry point for ensuring the base image exists.
// If the base image is already installed and update=false, returns immediately
// without constructing any build object or fetching remote metadata.
func EnsureBaseImage(ctx context.Context, update bool) error {
	if err := apptainer.EnsureApptainer(); err != nil {
		return err
	}
	if !update && config.FindBaseImage() != "" {
		return nil
	}
	obj, err := NewBaseImageBuildObject(ctx, update)
	if err != nil {
		return err
	}
	return obj.Build(ctx, false)
}

// IsInstalled checks all image search paths for the base image.
func (b *BaseImageBuildObject) IsInstalled() bool {
	return config.FindBaseImage() != ""
}

// Build overrides DefBuildObject.Build for base images.
// Unlike regular def builds, the SIF produced by apptainer is kept directly as the
// final image — no SquashFS extraction step is performed.
func (b *BaseImageBuildObject) Build(ctx context.Context, buildDeps bool) error {
	targetPath, finalPath := buildOverlayPaths(b.BuildObject)
	log := logging.FromContext(ctx)

	if !b.update && b.IsInstalled() {
		return nil
	}

	// Before starting any build work, check that no exec/run is holding the existing base image.
	if skip, err := checkShouldBuild(b.BuildObject); err != nil {
		return err
	} else if skip {
		// IsInstalled already handled the skip case above; checkShouldBuild handles the lock check.
		return nil
	}

	log.Info("building base image", "image", filepath.Base(targetPath), "source", b.buildSource)

	done := watchContext(ctx, "base image build")
	defer close(done)

	log.Info("running apptainer build", "source", b.buildSource)

	// Ensure the tmp directory exists before apptainer tries to write the SIF there.
	if err := utils.EnsureTmpSubdir(b.tmpDir); err != nil {
		return fmt.Errorf("failed to create tmp dir %s: %w", b.tmpDir, err)
	}

	buildOpts := &apptainer.BuildOptions{
		Force:     false,
		NoCleanup: false,
	}

	if err := apptainer.Build(ctx, b.tmpOverlayPath, b.buildSource, buildOpts); err != nil {
		b.Cleanup(true)
		if apptainer.IsBuildCancelled(err) {
			log.Info("build cancelled, base image unchanged", "image", filepath.Base(targetPath))
			return ErrBuildCancelled
		}
		return fmt.Errorf("failed to build SIF from %s: %w", b.buildSource, err)
	}

	// Move the SIF to its final location (no SquashFS extraction needed).
	if err := os.Rename(b.tmpOverlayPath, finalPath); err != nil {
		os.Remove(finalPath) //nolint:errcheck
		b.Cleanup(true)
		return fmt.Errorf("failed to move SIF to %s: %w", finalPath, err)
	}

	if err := utils.MakeExecutable(finalPath); err != nil {
		log.Debug("failed to set permissions", "path", finalPath, "err", err)
	}

	if err := atomicInstall(finalPath, targetPath, b.update); err != nil {
		return err
	}

	log.Info("base image ready", "kind", "success", "path", targetPath)
	b.Cleanup(false)
	return nil
}

// NewBaseImageBuildObject creates a BuildObject for the base image using the same
// path resolution as regular build objects. The base recipe is resolved through
// the catalog like any other module.
//
// When update=false the build is skipped if the image already exists anywhere in
// the configured search paths. When update=true the image is always rebuilt
// (written to .new then atomically renamed).
func NewBaseImageBuildObject(ctx context.Context, update bool) (*BaseImageBuildObject, error) {
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

	imagesDir, err := config.GetWritableImagesDir()
	if err != nil {
		return nil, fmt.Errorf("no writable directory for base image: %w", err)
	}

	if err := utils.MkdirAllShared(imagesDir); err != nil {
		return nil, fmt.Errorf("failed to create images directory: %w", err)
	}

	targetOverlayPath := filepath.Join(imagesDir, strings.ReplaceAll(nameVersion, "/", "--")+".sif")
	if abs, err := filepath.Abs(targetOverlayPath); err == nil {
		targetOverlayPath = abs
	}

	tmpDir := resolveTmpDirForDef()
	if absDir, err := filepath.Abs(tmpDir); err == nil {
		tmpDir = absDir
	}
	tmpOverlayPath, _ := buildTmpPaths(nameVersion, tmpDir, ".img")

	base := &BuildObject{
		nameVersion:       nameVersion,
		submitJob:         false, // base image is always built locally
		tmpDir:            tmpDir,
		tmpOverlayPath:    tmpOverlayPath,
		targetOverlayPath: targetOverlayPath,
		update:            update,
	}

	// Resolve build source using the same mechanism as regular build objects.
	_, isContainer, err := resolveBuildSource(ctx, base, tmpDir)
	if err != nil {
		return nil, fmt.Errorf("failed to resolve base image build source: %w", err)
	}
	if !isContainer {
		return nil, fmt.Errorf("base image build source not found or is not a .def file")
	}

	base.buildType = BuildTypeDef
	return &BaseImageBuildObject{BuildObject: base}, nil
}
