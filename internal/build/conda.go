package build

import (
	"fmt"
	"os"
	"os/exec"
	"path/filepath"
	"strings"

	"github.com/Justype/condatainer/internal/utils"
)

// CondaBuildObject implements BuildObject for conda packages
type CondaBuildObject struct {
	*BaseBuildObject
	packageName    string
	packageVersion string
}

// newCondaBuildObject creates a CondaBuildObject from base
func newCondaBuildObject(base *BaseBuildObject) (*CondaBuildObject, error) {
	parts := strings.Split(base.nameVersion, "/")
	if len(parts) != 2 {
		return nil, fmt.Errorf("conda package must be in format name/version, got: %s", base.nameVersion)
	}

	return &CondaBuildObject{
		BaseBuildObject: base,
		packageName:     parts[0],
		packageVersion:  parts[1],
	}, nil
}

// Type check implementations
func (c *CondaBuildObject) IsConda() bool { return true }
func (c *CondaBuildObject) IsDef() bool   { return false }
func (c *CondaBuildObject) IsShell() bool { return false }
func (c *CondaBuildObject) IsRef() bool   { return false }

// String returns a string representation
func (c *CondaBuildObject) String() string {
	return fmt.Sprintf(`BuildObject:
    name_version: %s
    build_source_type: Conda Package
    build_source: %s
    dependencies: %v
    sbatch_flags: %v
    tmp_overlay_path: %s
    target_overlay_path: %s
    cnt_dir_path: %s`,
		c.nameVersion,
		c.buildSource,
		c.dependencies,
		c.sbatchFlags,
		c.tmpOverlayPath,
		c.targetOverlayPath,
		c.cntDirPath,
	)
}

// Build implements the conda package build workflow
// Workflow:
//  1. Check if overlay already exists (skip if yes)
//  2. Create temporary ext3 overlay
//  3. Run micromamba create inside container with overlay
//  4. Set permissions
//  5. Pack to SquashFS
//  6. Cleanup
func (c *CondaBuildObject) Build(buildDeps bool) error {
	// Check if already exists
	if _, err := os.Stat(c.targetOverlayPath); err == nil {
		utils.PrintDebug("Overlay %s already exists at %s. Skipping creation.",
			filepath.Base(c.targetOverlayPath), utils.StylePath(c.targetOverlayPath))
		return nil
	}

	// Create temporary overlay
	if err := c.CreateTmpOverlay(false); err != nil {
		return fmt.Errorf("temporary overlay for %s already exists. Maybe a build is still running?", c.nameVersion)
	}

	utils.PrintDebug("Building SquashFS overlay at %s with packages: %s=%s",
		utils.StylePath(c.targetOverlayPath), c.packageName, c.packageVersion)

	// Build the apptainer exec command
	// TODO: These paths and settings should come from config
	baseImage := os.Getenv("CONDATAINER_BASE_IMAGE")
	if baseImage == "" {
		baseImage = "base.sif" // Default
	}
	apptainerBin := os.Getenv("SYSTEM_APPTAINER_BIN")
	if apptainerBin == "" {
		apptainerBin = "apptainer" // Default
	}
	condatainerDir := os.Getenv("CONDATAINER_DIR")
	if condatainerDir == "" {
		condatainerDir = filepath.Dir(c.targetOverlayPath)
	}

	relativePath := c.nameVersion

	// Build the bash script to run inside container
	bashScript := fmt.Sprintf(`
mkdir -p $TMPDIR
echo "Creating conda environment %s=%s in overlay..."
micromamba create -r /ext3/tmp -c conda-forge -c bioconda -q -y -p /cnt/%s %s=%s

echo "Setting permissions..."
find /cnt -type f -exec chmod ug+rw,o+r {} \;
find /cnt -type d -exec chmod ug+rwx,o+rx {} \;

echo "Packing overlay to SquashFS..."
mksquashfs /cnt %s -processors %d -keep-as-directory -comp gzip -b 1M
`,
		c.packageName, c.packageVersion,
		relativePath,
		c.packageName, c.packageVersion,
		filepath.Base(c.targetOverlayPath), c.ncpus,
	)

	// Build exec command args
	args := []string{
		"exec",
		"--env", "TMPDIR=/ext3/tmp",
		"--overlay", c.tmpOverlayPath,
		// TODO: Add bind and GPU args from config
		baseImage,
		"/bin/bash", "-c",
		bashScript,
	}

	cmd := exec.Command(apptainerBin, args...)
	utils.PrintDebug("[BUILD] Creating overlay with command: %s %s", apptainerBin, strings.Join(args, " "))

	// Run the build
	output, err := cmd.CombinedOutput()
	if err != nil {
		utils.PrintDebug("Build failed with output: %s", string(output))
		c.Cleanup(true)
		return fmt.Errorf("failed to build conda package %s: %w", c.nameVersion, err)
	}

	// Set permissions on target overlay
	if err := os.Chmod(c.targetOverlayPath, 0o664); err != nil {
		utils.PrintDebug("Failed to set permissions on %s: %v", c.targetOverlayPath, err)
	}

	utils.PrintDebug("Overlay %s created at %s. Removing temporary overlay...",
		filepath.Base(c.targetOverlayPath), utils.StylePath(c.targetOverlayPath))

	// Cleanup temp files
	c.Cleanup(false)

	return nil
}
