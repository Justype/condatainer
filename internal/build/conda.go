package build

import (
	"fmt"
	"os"
	"os/exec"
	"path/filepath"
	"strings"

	"github.com/Justype/condatainer/internal/config"
	"github.com/Justype/condatainer/internal/utils"
)

// CondaBuildObject implements BuildObject for conda packages
// Supports three modes:
// 1. Single package (name/version): packageName and packageVersion are set
// 2. Multiple packages (-n flag): buildSource contains comma-separated package specs
// 3. YAML file (-p flag): buildSource contains path to .yml/.yaml file
type CondaBuildObject struct {
	*BaseBuildObject
	packageName    string // Primary package name (for single package mode)
	packageVersion string // Primary package version (for single package mode)
}

// newCondaBuildObject creates a CondaBuildObject from base
func newCondaBuildObject(base *BaseBuildObject) (*CondaBuildObject, error) {
	parts := strings.Split(base.nameVersion, "/")

	// If buildSource is explicitly set (via -n flag), allow custom names without version
	// Otherwise, validate that conda packages follow name/version format
	var packageName, packageVersion string

	if base.buildSource != "" {
		// Custom buildSource provided (YAML file or comma-separated packages)
		// Use the full nameVersion as packageName, set version to "env"
		if len(parts) == 2 {
			packageName = parts[0]
			packageVersion = parts[1]
		} else {
			// Single name without version (e.g., "nvim_test")
			packageName = base.nameVersion
			packageVersion = "env"
		}
	} else {
		// Standard conda package - must be in name/version format
		if len(parts) != 2 {
			return nil, fmt.Errorf("conda package must be in format name/version, got: %s", base.nameVersion)
		}
		packageName = parts[0]
		packageVersion = parts[1]
	}

	return &CondaBuildObject{
		BaseBuildObject: base,
		packageName:     packageName,
		packageVersion:  packageVersion,
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
		script_specs: %v
		tmp_overlay_path: %s
		target_overlay_path: %s
		cnt_dir_path: %s`,
		c.nameVersion,
		c.buildSource,
		c.dependencies,
		c.scriptSpecs,
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

	// Determine build mode and create appropriate micromamba command
	var installCmd string
	var bindPaths []string

	if strings.HasSuffix(c.buildSource, ".yml") || strings.HasSuffix(c.buildSource, ".yaml") {
		// Mode 3: YAML file (-p prefix -f environment.yml)
		// The conda-meta folder will be automatically created by micromamba
		absFilePath, err := filepath.Abs(c.buildSource)
		if err != nil {
			return fmt.Errorf("failed to get absolute path for %s: %w", c.buildSource, err)
		}
		bindPaths = append(bindPaths, filepath.Dir(absFilePath))
		installCmd = fmt.Sprintf("micromamba create -r /ext3/tmp -c conda-forge -c bioconda -q -y -p /cnt/%s -f %s", c.nameVersion, absFilePath)
		utils.PrintDebug("Building SquashFS overlay at %s from YAML file %s",
			utils.StylePath(c.targetOverlayPath), utils.StylePath(absFilePath))
	} else if c.buildSource != "" {
		// Mode 2: Multiple packages (-n name pkg1 pkg2 ...)
		// buildSource contains comma-separated package specs like "nvim,nodejs"
		packages := strings.Split(c.buildSource, ",")
		// Convert slashes to equals for conda syntax (samtools/1.16 -> samtools=1.16)
		for i, pkg := range packages {
			packages[i] = strings.ReplaceAll(strings.TrimSpace(pkg), "/", "=")
		}
		installCmd = fmt.Sprintf("micromamba create -r /ext3/tmp -c conda-forge -c bioconda -q -y -p /cnt/%s %s",
			c.nameVersion, strings.Join(packages, " "))
		utils.PrintDebug("Building SquashFS overlay at %s with packages: %s",
			utils.StylePath(c.targetOverlayPath), strings.Join(packages, ", "))
	} else {
		// Mode 1: Single package (name/version)
		installCmd = fmt.Sprintf("micromamba create -r /ext3/tmp -c conda-forge -c bioconda -q -y -p /cnt/%s %s=%s",
			c.nameVersion, c.packageName, c.packageVersion)
		utils.PrintDebug("Building SquashFS overlay at %s with package: %s=%s",
			utils.StylePath(c.targetOverlayPath), c.packageName, c.packageVersion)
	}

	// Build the apptainer exec command using config values
	baseImage := config.GetBaseImage()
	apptainerBin := config.Global.ApptainerBin

	utils.PrintDebug("Using base image: %s", baseImage)
	utils.PrintDebug("Using apptainer binary: %s", apptainerBin)

	// Build the bash script to run inside container
	bashScript := fmt.Sprintf(`
trap 'exit 130' INT TERM

mkdir -p $TMPDIR
echo "Creating conda environment in overlay..."
%s

echo "Setting permissions..."
find /cnt -type f -exec chmod ug+rw,o+r {} \;
find /cnt -type d -exec chmod ug+rwx,o+rx {} \;

echo "Packing overlay to SquashFS..."
mksquashfs /cnt %s -processors %d -keep-as-directory %s -b 1M
`,
		installCmd,
		c.targetOverlayPath, c.ncpus, config.Global.Build.CompressArgs,
	)

	// Build exec command args
	args := []string{
		"exec",
		"--env", "TMPDIR=/ext3/tmp",
		"--overlay", c.tmpOverlayPath,
	}

	// Always bind the target overlay directory so mksquashfs can write there
	bindPaths = append(bindPaths, filepath.Dir(c.targetOverlayPath))

	// Add bind paths
	for _, path := range bindPaths {
		args = append(args, "--bind", path)
	}

	args = append(args,
		baseImage,
		"/bin/bash", "-c",
		bashScript,
	)

	cmd := exec.Command(apptainerBin, args...)
	utils.PrintDebug("[BUILD] Creating overlay with command: %s %s", apptainerBin, strings.Join(args, " "))

	// Run the build - show output in real-time
	cmd.Stdout = os.Stdout
	cmd.Stderr = os.Stderr

	// Use signal-aware runner to ensure cleanup runs on Ctrl+C
	if err := runCommandWithSignalHandling(cmd); err != nil {
		c.Cleanup(true)
		if isCancelledByUser(err) {
			return ErrBuildCancelled
		}
		return fmt.Errorf("failed to build conda package %s: %w", c.nameVersion, err)
	}

	// Set permissions on target overlay
	if err := os.Chmod(c.targetOverlayPath, 0o664); err != nil {
		utils.PrintDebug("Failed to set permissions on %s: %v", c.targetOverlayPath, err)
	}

	utils.PrintSuccess("Finished overlay %s", utils.StylePath(c.targetOverlayPath))

	// Cleanup temp files
	c.Cleanup(false)

	return nil
}
