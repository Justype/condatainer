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

// BaseImageBuildObject wraps DefBuildObject but overrides IsInstalled() to search
// all configured image paths (not just the writable directory).
type BaseImageBuildObject struct {
	*DefBuildObject
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

// Build overrides DefBuildObject.Build to check all configured image search paths
// (not just the writable directory) before deciding whether a build is needed.
func (b *BaseImageBuildObject) Build(ctx context.Context, buildDeps bool) error {
	if !b.update && b.IsInstalled() {
		return nil
	}
	return b.DefBuildObject.Build(ctx, buildDeps)
}

// NewBaseImageBuildObject creates a BuildObject for the base image using the same
// path resolution as regular build objects. base_image.def is searched in the
// build-scripts directories (same structure as any other build script) and downloaded
// from GitHub if not found locally. When the def comes from GitHub, DefBuildObject
// also tries to download a prebuilt overlay before falling back to a local build.
//
// When update=false the build is skipped if the image already exists anywhere in
// the configured search paths. When update=true the image is always rebuilt
// (written to .new then atomically renamed).
func NewBaseImageBuildObject(update bool) (*BaseImageBuildObject, error) {
	if err := apptainer.EnsureApptainer(); err != nil {
		return nil, err
	}

	// nameVersion mirrors the .sqf filename convention, e.g. "ubuntu24/base_image"
	nameVersion := config.Global.DefaultDistro + "/base_image"

	imagesDir, err := config.GetWritableImagesDir()
	if err != nil {
		return nil, fmt.Errorf("no writable directory for base image: %w", err)
	}

	if err := os.MkdirAll(imagesDir, utils.PermDir); err != nil {
		return nil, fmt.Errorf("failed to create images directory: %w", err)
	}

	targetOverlayPath := filepath.Join(imagesDir, strings.ReplaceAll(nameVersion, "/", "--")+".sqf")
	if abs, err := filepath.Abs(targetOverlayPath); err == nil {
		targetOverlayPath = abs
	}

	tmpDir := config.GetWritableTmpDir()
	tmpOverlayPath := getTmpOverlayPath(nameVersion, tmpDir)

	base := &BaseBuildObject{
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

	defObj := &DefBuildObject{BaseBuildObject: base}
	return &BaseImageBuildObject{DefBuildObject: defObj}, nil
}
