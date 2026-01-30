package build

import (
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
	return &DefBuildObject{
		BaseBuildObject: base,
	}, nil
}

// Type check implementations
func (d *DefBuildObject) IsConda() bool { return false }
func (d *DefBuildObject) IsDef() bool   { return true }
func (d *DefBuildObject) IsShell() bool { return false }
func (d *DefBuildObject) IsRef() bool   { return false }

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
func (d *DefBuildObject) Build(buildDeps bool) error {
	targetOverlayPath := d.targetOverlayPath
	if absTarget, err := filepath.Abs(targetOverlayPath); err == nil {
		targetOverlayPath = absTarget
	}

	styledOverlay := utils.StyleName(filepath.Base(targetOverlayPath))

	// Check if already exists
	if _, err := os.Stat(targetOverlayPath); err == nil {
		utils.PrintMessage("Overlay %s already exists at %s. Skipping creation.",
			styledOverlay, utils.StylePath(targetOverlayPath))
		// Ensure it's in the whitelist (handles existing .def overlays)
		UpdateDefBuiltWhitelist(d.nameVersion)
		return nil
	}

	buildMode := "local"
	if d.RequiresScheduler() {
		buildMode = "sbatch"
	}
	utils.PrintMessage("Building overlay %s (%s build) from %s", styledOverlay, utils.StyleAction(buildMode), utils.StylePath(d.buildSource))

	// Try to download prebuilt overlay first, fall back to building if not available
	if writableDir, err := config.GetWritableImagesDir(); err == nil {
		if filepath.Dir(targetOverlayPath) == writableDir {
			if tryDownloadPrebuiltOverlay(d.nameVersion, targetOverlayPath) {
				// Mark as .def-built in whitelist
				UpdateDefBuiltWhitelist(d.nameVersion)
				// Cleanup downloaded remote build script (if any)
				d.Cleanup(false)
				return nil
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

	if err := apptainer.Build(d.tmpOverlayPath, d.buildSource, buildOpts); err != nil {
		if apptainer.IsBuildCancelled(err) {
			utils.PrintMessage("Build cancelled for %s. Overlay unchanged.", styledOverlay)
			d.Cleanup(true)
			return ErrBuildCancelled
		}
		d.Cleanup(true)
		return fmt.Errorf("failed to build SIF from %s: %w", d.buildSource, err)
	}

	utils.PrintMessage("Extracting SquashFS to %s", utils.StylePath(targetOverlayPath))

	// Extract SquashFS from SIF
	if err := apptainer.DumpSifToSquashfs(d.tmpOverlayPath, targetOverlayPath); err != nil {
		d.Cleanup(true)
		return fmt.Errorf("failed to dump SquashFS from SIF: %w", err)
	}

	// Set permissions on target overlay
	if err := os.Chmod(targetOverlayPath, 0o664); err != nil {
		utils.PrintDebug("Failed to set permissions on %s: %v", targetOverlayPath, err)
	}

	utils.PrintSuccess("Finished overlay %s", utils.StylePath(targetOverlayPath))

	// Mark this overlay as .def-built in the whitelist
	UpdateDefBuiltWhitelist(d.nameVersion)

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
	overlayFilename := strings.ReplaceAll(normalized, "/", "--") + "_" + archName + ".sqf"
	url := fmt.Sprintf("%s/%s", config.PrebuiltBaseURL, overlayFilename)

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
