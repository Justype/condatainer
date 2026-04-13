package build

import (
	"context"
	"fmt"
	"os"
	"path/filepath"

	"github.com/Justype/condatainer/internal/apptainer"
	"github.com/Justype/condatainer/internal/config"
	"github.com/Justype/condatainer/internal/utils"
)

// buildDef implements the Apptainer .def build workflow on BuildObject.
// Workflow:
//  1. Check if overlay already exists (skip if yes)
//  2. Try downloading prebuilt overlay if available
//  3. Build SIF from .def file using apptainer build --fakeroot
//  4. Extract SquashFS partition from SIF
//  5. Set permissions
//  6. Cleanup
func (b *BuildObject) buildDef(ctx context.Context) error {
	targetPath, finalPath := buildOverlayPaths(b)
	styledOverlay := utils.StyleName(filepath.Base(targetPath))

	if skip, err := checkShouldBuild(b); skip || err != nil {
		return err
	}

	if err := b.createBuildLock(); err != nil {
		return err
	}
	defer b.removeBuildLock()

	utils.PrintMessage("Building overlay %s (%s build) from %s", styledOverlay, utils.StyleAction(buildModeLabel(b)), utils.StylePath(b.buildSource))

	done := watchContext(ctx, "def build")
	defer close(done)

	// Try to download prebuilt overlay first, fall back to building if not available.
	// Only attempt download if the build script source is remote.
	if b.isRemote {
		if writableDir, err := config.GetWritableImagesDir(); err == nil {
			if filepath.Dir(targetPath) == writableDir {
				downloadPath := buildFinalPath(targetPath, b.update)
				if tryDownloadPrebuiltOverlay(ctx, b.nameVersion, downloadPath, b.prebuiltLink) {
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

	// Ensure the tmp directory exists before apptainer tries to write the SIF there.
	if err := utils.EnsureTmpSubdir(b.tmpDir); err != nil {
		return fmt.Errorf("failed to create tmp dir %s: %w", b.tmpDir, err)
	}

	// For PL template .def files, substitute {key} placeholders before passing to
	// Apptainer. Apptainer reads directives like From: verbatim — no env expansion.
	defSource := b.buildSource
	if len(b.vars) > 0 {
		subPath, err := substituteTemplateFile(b.buildSource, b.vars, b.tmpDir)
		if err != nil {
			return fmt.Errorf("failed to substitute placeholders in .def file: %w", err)
		}
		defer os.Remove(subPath) //nolint:errcheck
		defSource = subPath
	}

	utils.PrintMessage("Running apptainer build from %s", utils.StylePath(b.buildSource))

	buildOpts := &apptainer.BuildOptions{
		Force:     false,
		NoCleanup: false,
	}

	if err := apptainer.Build(ctx, b.tmpOverlayPath, defSource, buildOpts); err != nil {
		b.Cleanup(true)
		if apptainer.IsBuildCancelled(err) {
			utils.PrintMessage("Build cancelled for %s. Overlay unchanged.", styledOverlay)
			return ErrBuildCancelled
		}
		return fmt.Errorf("failed to build SIF from %s: %w", b.buildSource, err)
	}

	utils.PrintMessage("Extracting SquashFS to %s", utils.StylePath(finalPath))

	if err := apptainer.DumpSifToSquashfs(ctx, b.tmpOverlayPath, finalPath); err != nil {
		os.Remove(finalPath) //nolint:errcheck
		b.Cleanup(true)
		if apptainer.IsBuildCancelled(err) {
			utils.PrintMessage("Build cancelled for %s. Overlay unchanged.", styledOverlay)
			return ErrBuildCancelled
		}
		return fmt.Errorf("failed to dump SquashFS from SIF: %w", err)
	}

	if err := os.Chmod(finalPath, utils.PermFile); err != nil {
		utils.PrintDebug("Failed to set permissions on %s: %v", finalPath, err)
	}

	if err := atomicInstall(finalPath, targetPath, b.update); err != nil {
		return err
	}

	utils.PrintSuccess("Finished overlay %s", utils.StylePath(targetPath))
	b.Cleanup(false)
	return nil
}

// tryDownloadPrebuiltOverlay attempts to download a prebuilt overlay from the given prebuiltLink base URL.
func tryDownloadPrebuiltOverlay(ctx context.Context, nameVersion, destPath, prebuiltLink string) bool {
	return tryDownloadPrebuilt(ctx, nameVersion, destPath, "sqf", prebuiltLink, utils.DownloadFile)
}
