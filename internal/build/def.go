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
//  2. Try downloading prebuilt overlay if available (TODO)
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
		return nil
	}

	buildMode := "local"
	if d.RequiresScheduler() {
		buildMode = "sbatch"
	}
	utils.PrintMessage("Building overlay %s (%s build) from %s", styledOverlay, utils.StyleAction(buildMode), utils.StylePath(d.buildSource))

	// Check if prebuilt overlay is available and try to download it
	prebuiltOverlays := map[string]bool{
		"build-essential": true,
	}
	if prebuiltOverlays[d.nameVersion] && filepath.Dir(targetOverlayPath) == config.Global.ImagesDir {
		if tryDownloadPrebuiltOverlay(d.nameVersion, targetOverlayPath) {
			return nil
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
		utils.PrintError("Apptainer build failed for %s: %v", styledOverlay, err)
		d.Cleanup(true)
		return fmt.Errorf("failed to build SIF from %s: %w", d.buildSource, err)
	}

	utils.PrintMessage("Extracting SquashFS to %s", utils.StylePath(targetOverlayPath))

	// Extract SquashFS from SIF
	if err := apptainer.DumpSifToSquashfs(d.tmpOverlayPath, targetOverlayPath); err != nil {
		utils.PrintError("Failed to export SquashFS for %s: %v", styledOverlay, err)
		d.Cleanup(true)
		return fmt.Errorf("failed to dump SquashFS from SIF: %w", err)
	}

	// Set permissions on target overlay
	if err := os.Chmod(targetOverlayPath, 0o664); err != nil {
		utils.PrintDebug("Failed to set permissions on %s: %v", targetOverlayPath, err)
	}

	utils.PrintMessage("Finished overlay %s", styledOverlay)
	utils.PrintDebug("Overlay %s created at %s. Removing temporary overlay...",
		filepath.Base(targetOverlayPath), utils.StylePath(targetOverlayPath))

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
		utils.PrintWarning("Pre-built overlays not available for architecture: %s", arch)
		return false
	}

	normalized := utils.NormalizeNameVersion(nameVersion)
	overlayFilename := strings.ReplaceAll(normalized, "/", "--") + "_" + archName + ".sqf"
	url := fmt.Sprintf("https://github.com/Justype/condatainer/releases/download/v1.0.5/%s", overlayFilename)

	utils.PrintMessage("Attempting to download pre-built overlay for %s...", utils.StyleName(normalized))

	if err := utils.DownloadFile(url, destPath); err != nil {
		utils.PrintWarning("Failed to download pre-built overlay: %v", err)
		return false
	}

	utils.PrintSuccess("Pre-built %s downloaded.", utils.StyleName(normalized))
	return true
}
