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
	"github.com/Justype/condatainer/internal/overlay"
	"github.com/Justype/condatainer/internal/utils"
)

// DefBuildObject implements BuildObject for Apptainer definition files
type DefBuildObject struct {
	*BaseBuildObject
}

// newDefBuildObject creates a DefBuildObject from base
func newDefBuildObject(base *BaseBuildObject) (*DefBuildObject, error) {
	return &DefBuildObject{
		BaseBuildObject: base,
	}, nil
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
	targetOverlayPath := d.targetOverlayPath
	if absTarget, err := filepath.Abs(targetOverlayPath); err == nil {
		targetOverlayPath = absTarget
	}

	styledOverlay := utils.StyleName(filepath.Base(targetOverlayPath))

	// Check if already exists (skip only when not in update mode)
	if !d.update {
		if _, err := os.Stat(targetOverlayPath); err == nil {
			utils.PrintMessage("Overlay %s already exists at %s. Skipping creation.",
				styledOverlay, utils.StylePath(targetOverlayPath))
			return nil
		}
	}

	// Before starting any build work, check that no exec/run holds the existing overlay.
	if d.update && utils.FileExists(targetOverlayPath) {
		if lock, err := overlay.AcquireLock(targetOverlayPath, true); err != nil {
			return fmt.Errorf("cannot update %s: %w", d.nameVersion, err)
		} else {
			lock.Close()
		}
	}

	buildMode := "local"
	if d.RequiresScheduler() {
		buildMode = "sbatch"
	}
	utils.PrintMessage("Building overlay %s (%s build) from %s", styledOverlay, utils.StyleAction(buildMode), utils.StylePath(d.buildSource))

	done := make(chan struct{})

	// Ensure cleanup on function exit
	defer close(done)

	// Setup cleanup function with sync.Once to prevent double cleanup
	var cleanupOnce sync.Once
	cleanupFunc := func() {
		cleanupOnce.Do(func() {
			utils.PrintMessage("Cleaning up temporary files for %s...", styledOverlay)
			d.Cleanup(true)
		})
	}

	// Monitor for context cancellation
	go func() {
		select {
		case <-ctx.Done():
			utils.PrintWarning("Build cancelled. Interrupting build...")
			// Do not call cleanupFunc() here to avoid race condition with process termination
			// cleanupFunc() will be called after apptainer.Build returns
		case <-done:
			return
		}
	}()

	// Try to download prebuilt overlay first, fall back to building if not available
	// Only attempt download if the build script source is remote
	if d.isRemote {
		if writableDir, err := config.GetWritableImagesDir(); err == nil {
			if filepath.Dir(targetOverlayPath) == writableDir {
				// In update mode write to .new so the existing overlay is intact until success
				downloadPath := targetOverlayPath
				if d.update {
					downloadPath = targetOverlayPath + ".new"
				}
				if tryDownloadPrebuiltOverlay(d.nameVersion, downloadPath) {
					if d.update {
						os.Remove(targetOverlayPath) //nolint:errcheck
						if err := os.Rename(downloadPath, targetOverlayPath); err != nil {
							os.Remove(downloadPath) //nolint:errcheck
							return fmt.Errorf("failed to replace overlay %s: %w", targetOverlayPath, err)
						}
					}
					// Cleanup downloaded remote build script (if any)
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
		cleanupFunc()
		if apptainer.IsBuildCancelled(err) {
			utils.PrintMessage("Build cancelled for %s. Overlay unchanged.", styledOverlay)
			return ErrBuildCancelled
		}
		return fmt.Errorf("failed to build SIF from %s: %w", d.buildSource, err)
	}

	// Determine final output path: write to .new when updating for atomic replacement
	finalPath := targetOverlayPath
	if d.update {
		finalPath = targetOverlayPath + ".new"
	}

	utils.PrintMessage("Extracting SquashFS to %s", utils.StylePath(finalPath))

	// Extract SquashFS from SIF
	if err := apptainer.DumpSifToSquashfs(ctx, d.tmpOverlayPath, finalPath); err != nil {
		os.Remove(finalPath) //nolint:errcheck
		cleanupFunc()
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

	// Atomic replacement: remove old and rename .new → target
	if d.update {
		os.Remove(targetOverlayPath) //nolint:errcheck
		if err := os.Rename(finalPath, targetOverlayPath); err != nil {
			os.Remove(finalPath) //nolint:errcheck
			return fmt.Errorf("failed to replace overlay %s: %w", targetOverlayPath, err)
		}
	}

	utils.PrintSuccess("Finished overlay %s", utils.StylePath(targetOverlayPath))

	// Cleanup temp files
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
