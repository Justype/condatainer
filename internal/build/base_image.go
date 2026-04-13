package build

import (
	"context"
	"fmt"
	"os"
	"path/filepath"
	"strings"

	"github.com/Justype/condatainer/internal/apptainer"
	"github.com/Justype/condatainer/internal/config"
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
	obj, err := NewBaseImageBuildObject(update)
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
	styledImage := utils.StyleName(filepath.Base(targetPath))

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

	utils.PrintMessage("Building base image %s (local build) from %s", styledImage, utils.StylePath(b.buildSource))

	done := watchContext(ctx, "base image build")
	defer close(done)

	// Try to download prebuilt .sif from GitHub releases if build source is remote.
	if b.isRemote {
		if writableDir, err := config.GetWritableImagesDir(); err == nil {
			if filepath.Dir(targetPath) == writableDir {
				downloadPath := buildFinalPath(targetPath, b.update)
				if tryDownloadPrebuiltSif(ctx, b.nameVersion, downloadPath, b.prebuiltLink) {
					if err := atomicInstall(downloadPath, targetPath, b.update); err != nil {
						return err
					}
					b.Cleanup(false)
					return nil
				}
				if b.update {
					os.Remove(downloadPath) //nolint:errcheck
				}
			}
		}
	}

	utils.PrintMessage("Running apptainer build from %s", utils.StylePath(b.buildSource))

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
			utils.PrintMessage("Build cancelled for %s. Base image unchanged.", styledImage)
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

	if err := os.Chmod(finalPath, utils.PermExec); err != nil {
		utils.PrintDebug("Failed to set permissions on %s: %v", finalPath, err)
	}

	if err := atomicInstall(finalPath, targetPath, b.update); err != nil {
		return err
	}

	utils.PrintSuccess("Finished base image %s", utils.StylePath(targetPath))
	b.Cleanup(false)
	return nil
}

// NewBaseImageBuildObject creates a BuildObject for the base image using the same
// path resolution as regular build objects. base_image.def is searched in the
// build-scripts directories (same structure as any other build script) and downloaded
// from GitHub if not found locally. When the def comes from GitHub, BaseImageBuildObject
// also tries to download a prebuilt .sif before falling back to a local build.
//
// When update=false the build is skipped if the image already exists anywhere in
// the configured search paths. When update=true the image is always rebuilt
// (written to .new then atomically renamed).
func NewBaseImageBuildObject(update bool) (*BaseImageBuildObject, error) {
	if err := apptainer.EnsureApptainer(); err != nil {
		return nil, err
	}

	// nameVersion mirrors the .sif filename convention, e.g. "ubuntu24/base_image"
	nameVersion := config.Global.DefaultDistro + "/base_image"

	imagesDir, err := config.GetWritableImagesDir()
	if err != nil {
		return nil, fmt.Errorf("no writable directory for base image: %w", err)
	}

	if err := os.MkdirAll(imagesDir, utils.PermDir); err != nil {
		return nil, fmt.Errorf("failed to create images directory: %w", err)
	}

	targetOverlayPath := filepath.Join(imagesDir, strings.ReplaceAll(nameVersion, "/", "--")+".sif")
	if abs, err := filepath.Abs(targetOverlayPath); err == nil {
		targetOverlayPath = abs
	}

	tmpDir := resolveTmpDirForDef()
	tmpOverlayPath, _ := buildTmpPaths(nameVersion, tmpDir, ".img")

	base := &BuildObject{
		nameVersion:       nameVersion,
		submitJob:         false, // base image is always built locally
		tmpOverlayPath:    tmpOverlayPath,
		targetOverlayPath: targetOverlayPath,
		update:            update,
	}

	// Resolve build source using the same mechanism as regular build objects.
	// Finds base_image.def in local build-scripts dirs, or downloads from GitHub.
	_, isContainer, err := resolveBuildSource(base, tmpDir)
	if err != nil {
		return nil, fmt.Errorf("failed to resolve base image build source: %w", err)
	}
	if !isContainer {
		return nil, fmt.Errorf("base image build source not found or is not a .def file")
	}

	base.buildType = BuildTypeDef
	return &BaseImageBuildObject{BuildObject: base}, nil
}

// tryDownloadPrebuiltSif attempts to download a prebuilt .sif base image from the given prebuiltLink base URL.
func tryDownloadPrebuiltSif(ctx context.Context, nameVersion, destPath, prebuiltLink string) bool {
	return tryDownloadPrebuilt(ctx, nameVersion, destPath, "sif", prebuiltLink, utils.DownloadExecutable)
}
