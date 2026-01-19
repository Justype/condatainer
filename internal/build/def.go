package build

import (
	"fmt"
	"os"
	"path/filepath"

	"github.com/Justype/condatainer/internal/apptainer"
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
    sbatch_flags: %v
    tmp_overlay_path: %s
    target_overlay_path: %s
    cnt_dir_path: %s`,
		d.nameVersion,
		d.buildSource,
		d.dependencies,
		d.sbatchFlags,
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
	// Check if already exists
	if _, err := os.Stat(d.targetOverlayPath); err == nil {
		utils.PrintDebug("Overlay %s already exists at %s. Skipping creation.",
			filepath.Base(d.targetOverlayPath), utils.StylePath(d.targetOverlayPath))
		return nil
	}

	// TODO: Check if prebuilt overlay is available and try to download it
	// if d.nameVersion in Config.PREBUILT_OVERLAYS {
	//     if tryDownloadPrebuiltOverlay(d.nameVersion) {
	//         return nil
	//     }
	// }

	utils.PrintDebug("Building Apptainer SIF from %s", utils.StylePath(d.buildSource))

	// Build SIF using apptainer build --fakeroot
	// Note: tmpOverlayPath is actually a SIF file at this stage, not sqf yet
	buildOpts := &apptainer.BuildOptions{
		Force:     false,
		NoCleanup: false,
	}

	if err := apptainer.Build(d.tmpOverlayPath, d.buildSource, buildOpts); err != nil {
		d.Cleanup(true)
		return fmt.Errorf("failed to build SIF from %s: %w", d.buildSource, err)
	}

	utils.PrintDebug("Dumping SquashFS partition to %s", utils.StylePath(d.targetOverlayPath))

	// Extract SquashFS from SIF
	if err := apptainer.DumpSifToSquashfs(d.tmpOverlayPath, d.targetOverlayPath); err != nil {
		d.Cleanup(true)
		return fmt.Errorf("failed to dump SquashFS from SIF: %w", err)
	}

	// Set permissions on target overlay
	if err := os.Chmod(d.targetOverlayPath, 0o664); err != nil {
		utils.PrintDebug("Failed to set permissions on %s: %v", d.targetOverlayPath, err)
	}

	utils.PrintDebug("Overlay %s created at %s. Removing temporary overlay...",
		filepath.Base(d.targetOverlayPath), utils.StylePath(d.targetOverlayPath))

	// Cleanup temp files
	d.Cleanup(false)

	return nil
}
