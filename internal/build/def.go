package build

import (
	"context"
	"fmt"
	"os"
	"path/filepath"
	"runtime"
	"strings"

	"github.com/Justype/condatainer/internal/apptainer"
	"github.com/Justype/condatainer/internal/config"
	"github.com/Justype/condatainer/internal/utils"
)

// DefBuildObject implements BuildObject for Apptainer definition files
type DefBuildObject struct {
	*BaseBuildObject
}

// newDefBuildObject creates a DefBuildObject from base
func newDefBuildObject(base *BaseBuildObject) (*DefBuildObject, error) {
	return &DefBuildObject{BaseBuildObject: base}, nil
}

// Type returns the build type
func (d *DefBuildObject) Type() BuildType { return BuildTypeDef }

// String returns a string representation
func (d *DefBuildObject) String() string {
	return fmt.Sprintf(`BuildObject:
		name_version: %s
		build_source_type: Apptainer File
		build_source: %s
		dependencies: %v
		script_specs: %v
		tmp_overlay_path: %s
		target_overlay_path: %s
		cnt_dir_path: %s`,
		d.nameVersion,
		d.buildSource,
		d.dependencies,
		d.scriptSpecs,
		d.tmpOverlayPath,
		d.targetOverlayPath,
		d.cntDirPath,
	)
}

// Build implements the Apptainer .def build workflow
// Workflow:
//  1. Check if overlay already exists (skip if yes)
//  2. Try downloading prebuilt overlay if available
//  3. Build SIF from .def file using apptainer build --fakeroot
//  4. Extract SquashFS partition from SIF
//  5. Set permissions
//  6. Cleanup
func (d *DefBuildObject) Build(ctx context.Context, buildDeps bool) error {
	targetPath, finalPath := buildOverlayPaths(d.BaseBuildObject)
	styledOverlay := utils.StyleName(filepath.Base(targetPath))

	if skip, err := checkShouldBuild(d.BaseBuildObject); skip || err != nil {
		return err
	}

	if err := d.createBuildLock(); err != nil {
		return err
	}
	defer d.removeBuildLock()

	utils.PrintMessage("Building overlay %s (%s build) from %s", styledOverlay, utils.StyleAction(buildModeLabel(d.BaseBuildObject)), utils.StylePath(d.buildSource))

	done := watchContext(ctx, "def build")
	defer close(done)

	// Try to download prebuilt overlay first, fall back to building if not available
	// Only attempt download if the build script source is remote
	if d.isRemote {
		if writableDir, err := config.GetWritableImagesDir(); err == nil {
			if filepath.Dir(targetPath) == writableDir {
				downloadPath := buildFinalPath(targetPath, d.update)
				if tryDownloadPrebuiltOverlay(d.nameVersion, downloadPath) {
					if err := atomicInstall(downloadPath, targetPath, d.update); err != nil {
						return err
					}
					d.Cleanup(false)
					return nil
				}
				if d.update {
					os.Remove(downloadPath) //nolint:errcheck // clean up partial download
				}
			}
		}
	}

	utils.PrintMessage("Running apptainer build from %s", utils.StylePath(d.buildSource))

	// Build SIF using apptainer build --fakeroot
	// Note: tmpOverlayPath is actually a SIF file at this stage, not sqf yet
	buildOpts := &apptainer.BuildOptions{
		Force:     false,
		NoCleanup: false,
	}

	if err := apptainer.Build(ctx, d.tmpOverlayPath, d.buildSource, buildOpts); err != nil {
		d.Cleanup(true)
		if apptainer.IsBuildCancelled(err) {
			utils.PrintMessage("Build cancelled for %s. Overlay unchanged.", styledOverlay)
			return ErrBuildCancelled
		}
		return fmt.Errorf("failed to build SIF from %s: %w", d.buildSource, err)
	}

	utils.PrintMessage("Extracting SquashFS to %s", utils.StylePath(finalPath))

	// Extract SquashFS from SIF
	if err := apptainer.DumpSifToSquashfs(ctx, d.tmpOverlayPath, finalPath); err != nil {
		os.Remove(finalPath) //nolint:errcheck
		d.Cleanup(true)
		if apptainer.IsBuildCancelled(err) {
			utils.PrintMessage("Build cancelled for %s. Overlay unchanged.", styledOverlay)
			return ErrBuildCancelled
		}
		return fmt.Errorf("failed to dump SquashFS from SIF: %w", err)
	}

	// Set permissions on output file
	if err := os.Chmod(finalPath, 0o664); err != nil {
		utils.PrintDebug("Failed to set permissions on %s: %v", finalPath, err)
	}

	if err := atomicInstall(finalPath, targetPath, d.update); err != nil {
		return err
	}

	utils.PrintSuccess("Finished overlay %s", utils.StylePath(targetPath))
	d.Cleanup(false)
	return nil
}

// tryDownloadPrebuiltOverlay attempts to download a prebuilt overlay from GitHub releases
func tryDownloadPrebuiltOverlay(nameVersion, destPath string) bool {
	arch := runtime.GOARCH

	// Map Go arch names to the format used in releases
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
	overlayFilename := parts[1] + "_" + archName + ".sqf"
	url := fmt.Sprintf("%s/%s/%s", config.Global.PrebuiltLink, parts[0], overlayFilename)

	// Check if prebuilt exists before attempting download
	if !utils.URLExists(url) {
		return false
	}

	utils.PrintMessage("Found pre-built %s. Downloading...", utils.StyleName(normalized))

	if err := utils.DownloadFile(url, destPath); err != nil {
		utils.PrintWarning("Download failed. Falling back to local build.")
		return false
	}

	utils.PrintSuccess("Pre-built %s downloaded.", utils.StyleName(normalized))
	return true
}
