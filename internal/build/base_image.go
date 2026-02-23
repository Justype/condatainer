package build

import (
	"context"
	"fmt"
	"os"
	"path/filepath"
	"runtime"
	"strings"
	"sync"

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

// Build overrides DefBuildObject.Build for base images.
// Unlike regular def builds, the SIF produced by apptainer is kept directly as the
// final image — no SquashFS extraction step is performed.
func (b *BaseImageBuildObject) Build(ctx context.Context, buildDeps bool) error {
	targetPath := b.targetOverlayPath
	if absTarget, err := filepath.Abs(targetPath); err == nil {
		targetPath = absTarget
	}

	styledImage := utils.StyleName(filepath.Base(targetPath))

	if !b.update && b.IsInstalled() {
		return nil
	}

	utils.PrintMessage("Building base image %s (local build) from %s", styledImage, utils.StylePath(b.buildSource))

	done := make(chan struct{})
	defer close(done)

	var cleanupOnce sync.Once
	cleanupFunc := func() {
		cleanupOnce.Do(func() {
			utils.PrintMessage("Cleaning up temporary files for %s...", styledImage)
			b.Cleanup(true)
		})
	}

	go func() {
		select {
		case <-ctx.Done():
			utils.PrintWarning("Build cancelled. Interrupting build...")
		case <-done:
			return
		}
	}()

	// Try to download prebuilt .sif from GitHub releases if build source is remote.
	if b.isRemote {
		if writableDir, err := config.GetWritableImagesDir(); err == nil {
			if filepath.Dir(targetPath) == writableDir {
				downloadPath := targetPath
				if b.update {
					downloadPath = targetPath + ".new"
				}
				if tryDownloadPrebuiltSif(b.nameVersion, downloadPath) {
					if b.update {
						os.Remove(targetPath) //nolint:errcheck
						if err := os.Rename(downloadPath, targetPath); err != nil {
							os.Remove(downloadPath) //nolint:errcheck
							return fmt.Errorf("failed to replace base image %s: %w", targetPath, err)
						}
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

	buildOpts := &apptainer.BuildOptions{
		Force:     false,
		NoCleanup: false,
	}

	if err := apptainer.Build(ctx, b.tmpOverlayPath, b.buildSource, buildOpts); err != nil {
		cleanupFunc()
		if apptainer.IsBuildCancelled(err) {
			utils.PrintMessage("Build cancelled for %s. Base image unchanged.", styledImage)
			return ErrBuildCancelled
		}
		return fmt.Errorf("failed to build SIF from %s: %w", b.buildSource, err)
	}

	// Determine final output path; write to .new in update mode for atomic replacement.
	finalPath := targetPath
	if b.update {
		finalPath = targetPath + ".new"
	}

	// Move the SIF to its final location (no SquashFS extraction needed).
	if err := os.Rename(b.tmpOverlayPath, finalPath); err != nil {
		os.Remove(finalPath) //nolint:errcheck
		cleanupFunc()
		return fmt.Errorf("failed to move SIF to %s: %w", finalPath, err)
	}

	if err := os.Chmod(finalPath, utils.PermExec); err != nil {
		utils.PrintDebug("Failed to set permissions on %s: %v", finalPath, err)
	}

	if b.update {
		os.Remove(targetPath) //nolint:errcheck
		if err := os.Rename(finalPath, targetPath); err != nil {
			os.Remove(finalPath) //nolint:errcheck
			return fmt.Errorf("failed to replace base image %s: %w", targetPath, err)
		}
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

// tryDownloadPrebuiltSif attempts to download a prebuilt .sif base image from GitHub releases.
func tryDownloadPrebuiltSif(nameVersion, destPath string) bool {
	arch := runtime.GOARCH

	archMap := map[string]string{
		"amd64": "x86_64",
		"arm64": "aarch64",
	}

	archName, ok := archMap[arch]
	if !ok {
		return false
	}

	normalized := utils.NormalizeNameVersion(nameVersion)
	parts := strings.SplitN(normalized, "/", 2)
	if len(parts) != 2 {
		return false
	}
	sifFilename := parts[1] + "_" + archName + ".sif"
	url := fmt.Sprintf("%s/%s/%s", config.PrebuiltBaseURL, parts[0], sifFilename)

	if !utils.URLExists(url) {
		return false
	}

	utils.PrintMessage("Found pre-built base image %s. Downloading...", utils.StyleName(normalized))

	if err := utils.DownloadExecutable(url, destPath); err != nil {
		utils.PrintWarning("Download failed. Falling back to local build.")
		return false
	}

	utils.PrintSuccess("Pre-built base image %s downloaded.", utils.StyleName(normalized))
	return true
}
