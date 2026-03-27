package build

import (
	"context"
	"fmt"
	"os"
	"path/filepath"
	"strings"

	"github.com/Justype/condatainer/internal/config"
	"github.com/Justype/condatainer/internal/exec"
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
		// Channel-annotated input (e.g. "bioconda::star/2.7.11b"): use full spec for micromamba.
		if base.condaChannelPkg != "" {
			packageName = base.condaChannelPkg
		}
	}

	return &CondaBuildObject{
		BaseBuildObject: base,
		packageName:     packageName,
		packageVersion:  packageVersion,
	}, nil
}

// Type returns the build type
func (c *CondaBuildObject) Type() BuildType { return BuildTypeConda }

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
//  2. Create build workspace (ext3 overlay or host dirs)
//  3. Run micromamba create inside container
//  4. Pack to SquashFS and perform atomic install
//  5. Cleanup
func (c *CondaBuildObject) Build(ctx context.Context, buildDeps bool) error {
	targetPath, finalPath := buildOverlayPaths(c.BaseBuildObject)
	styledOverlay := utils.StyleName(filepath.Base(targetPath))

	if skip, err := checkShouldBuild(c.BaseBuildObject); skip || err != nil {
		return err
	}

	if err := c.createBuildLock(); err != nil {
		return err
	}
	defer c.removeBuildLock()

	utils.PrintMessage("Building overlay %s (%s build)", styledOverlay, utils.StyleAction(buildModeLabel(c.BaseBuildObject)))

	if err := prepareBuildWorkspace(ctx, c.BaseBuildObject, config.Global.Build.UseTmpOverlay); err != nil {
		return err
	}

	if err := c.runCondaBuild(ctx, finalPath); err != nil {
		c.Cleanup(true)
		return err
	}

	if err := os.Chmod(finalPath, utils.PermFile); err != nil {
		utils.PrintDebug("Failed to set permissions on %s: %v", finalPath, err)
	}

	if err := atomicInstall(finalPath, targetPath, c.update); err != nil {
		return err
	}

	utils.PrintSuccess("Finished overlay %s", utils.StylePath(targetPath))
	c.Cleanup(false)
	return nil
}

// buildChannelFlags builds the micromamba -c flags from the configured channels.
func buildChannelFlags() string {
	var parts []string
	for _, ch := range config.Global.Build.Channels {
		parts = append(parts, "-c "+ch)
	}
	return strings.Join(parts, " ")
}

// buildInstallCmd returns the micromamba install command string and any extra bind paths.
// Handles three modes: YAML file, comma-separated packages, and single package.
func (c *CondaBuildObject) buildInstallCmd() (cmd string, extraBindPaths []string, err error) {
	var quietFlag string
	if utils.QuietMode {
		quietFlag = "-q"
	}
	channelFlags := buildChannelFlags()

	if utils.IsYaml(c.buildSource) {
		// Mode 3: YAML file (-p prefix -f environment.yml)
		absFilePath, err := filepath.Abs(c.buildSource)
		if err != nil {
			return "", nil, fmt.Errorf("failed to get absolute path for %s: %w", c.buildSource, err)
		}
		extraBindPaths = []string{filepath.Dir(absFilePath)}
		cmd = fmt.Sprintf("micromamba create -r /ext3/tmp %s -y %s -p /cnt/%s -f %s",
			channelFlags, quietFlag, c.nameVersion, absFilePath)
		utils.PrintMessage("Populating overlay %s via %s",
			utils.StyleName(filepath.Base(c.targetOverlayPath)), "environment.yml")
	} else if c.buildSource != "" {
		// Mode 2: Multiple packages (-n name pkg1 pkg2 ...)
		packages := strings.Split(c.buildSource, ",")
		for i, pkg := range packages {
			packages[i] = strings.ReplaceAll(strings.TrimSpace(pkg), "/", "=")
		}
		cmd = fmt.Sprintf("micromamba create -r /ext3/tmp %s -y %s -p /cnt/%s %s",
			channelFlags, quietFlag, c.nameVersion, strings.Join(packages, " "))
		utils.PrintMessage("Populating overlay %s via %s",
			utils.StyleName(filepath.Base(c.targetOverlayPath)),
			fmt.Sprintf("micromamba (%s)", strings.Join(packages, ", ")))
	} else {
		// Mode 1: Single package (name/version)
		cmd = fmt.Sprintf("micromamba create -r /ext3/tmp %s -y %s -p /cnt/%s %s=%s",
			channelFlags, quietFlag, c.nameVersion, c.packageName, c.packageVersion)
		utils.PrintMessage("Populating overlay %s via %s",
			utils.StyleName(filepath.Base(c.targetOverlayPath)),
			fmt.Sprintf("micromamba (%s=%s)", c.packageName, c.packageVersion))
	}

	return cmd, extraBindPaths, nil
}

// buildCondaExecOpts constructs exec.Options for the conda build container run.
func (c *CondaBuildObject) buildCondaExecOpts(finalPath string) (exec.Options, error) {
	installCmd, extraBindPaths, err := c.buildInstallCmd()
	if err != nil {
		return exec.Options{}, err
	}

	var echoPrefix string
	if utils.QuietMode {
		echoPrefix = ": #"
	} else {
		echoPrefix = "echo"
	}

	bashScript := fmt.Sprintf(`
trap 'exit 130' INT TERM
set -e

mkdir -p $TMPDIR
%[1]s "Creating conda environment in overlay..."
%[2]s

if [ -z "$(ls -A /cnt 2>/dev/null)" ]; then
    %[1]s "Conda environment is empty, nothing to pack."
    exit 1
fi

%[1]s "Packing overlay to SquashFS..."
mksquashfs /cnt %[3]s -processors %[4]d -b %[5]s -keep-as-directory -all-root %[6]s
`,
		echoPrefix, installCmd,
		finalPath, c.effectiveNcpus(), config.Global.Build.BlockSize, config.Global.Build.CompressArgs,
	)

	// Always bind the target overlay directory so mksquashfs can write there.
	bindPaths := append(extraBindPaths, filepath.Dir(c.targetOverlayPath))

	var opts exec.Options
	if !config.Global.Build.UseTmpOverlay {
		buildTmpDir := getBuildTmpDir(c.BaseBuildObject)
		bindPaths = append(bindPaths,
			buildTmpDir+":/ext3/tmp",
			c.cntDirPath+":/cnt",
		)
		opts = exec.Options{
			BaseImage:      config.GetBaseImage(),
			ApptainerBin:   config.Global.ApptainerBin,
			Overlays:       []string{},
			BindPaths:      bindPaths,
			EnvSettings:    []string{"TMPDIR=/ext3/tmp"},
			Command:        []string{"/bin/bash", "-c", bashScript},
			HidePrompt:     true,
			WritableImg:    false,
			ApptainerFlags: []string{"--writable-tmpfs"},
			PassThruStdin:  true,
		}
	} else {
		opts = exec.Options{
			BaseImage:     config.GetBaseImage(),
			ApptainerBin:  config.Global.ApptainerBin,
			Overlays:      []string{c.tmpOverlayPath},
			BindPaths:     bindPaths,
			EnvSettings:   []string{"TMPDIR=/ext3/tmp"},
			Command:       []string{"/bin/bash", "-c", bashScript},
			HidePrompt:    true,
			WritableImg:   true,
			PassThruStdin: true,
		}
	}

	return opts, nil
}

// runCondaBuild sets up a context watcher, builds exec opts, and runs the conda build.
// The caller is responsible for calling Cleanup(true) if an error is returned.
func (c *CondaBuildObject) runCondaBuild(ctx context.Context, finalPath string) error {
	opts, err := c.buildCondaExecOpts(finalPath)
	if err != nil {
		return err
	}

	utils.PrintDebug("[BUILD] Creating overlay for %s with options: overlays=%v, bindPaths=%v",
		c.nameVersion, opts.Overlays, opts.BindPaths)

	done := watchContext(ctx, "conda build")
	defer close(done)

	if err := exec.Run(ctx, opts); err != nil {
		if isCancelledByUser(err) {
			return ErrBuildCancelled
		}
		return fmt.Errorf("failed to build conda package %s: %w", c.nameVersion, err)
	}

	return nil
}
