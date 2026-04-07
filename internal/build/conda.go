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

// buildConda implements the conda package build workflow on BuildObject.
// Workflow:
//  1. Check if overlay already exists (skip if yes)
//  2. Create build workspace (ext3 overlay or host dirs)
//  3. Run micromamba create inside container
//  4. Pack to SquashFS and perform atomic install
//  5. Cleanup
func (b *BuildObject) buildConda(ctx context.Context) error {
	targetPath, finalPath := buildOverlayPaths(b)
	styledOverlay := utils.StyleName(filepath.Base(targetPath))

	if skip, err := checkShouldBuild(b); skip || err != nil {
		return err
	}

	if err := b.createBuildLock(); err != nil {
		return err
	}
	defer b.removeBuildLock()

	utils.PrintMessage("Building overlay %s (%s build)", styledOverlay, utils.StyleAction(buildModeLabel(b)))

	if err := prepareBuildWorkspace(ctx, b, config.Global.Build.UseTmpOverlay); err != nil {
		return err
	}

	if err := b.runCondaBuild(ctx, finalPath); err != nil {
		b.Cleanup(true)
		return err
	}

	if err := os.Chmod(finalPath, utils.PermFile); err != nil {
		utils.PrintDebug("Failed to set permissions on %s: %v", finalPath, err)
	}

	if err := atomicInstall(finalPath, targetPath, b.update); err != nil {
		return err
	}

	b.saveCondaEnvFile(targetPath)

	utils.PrintSuccess("Finished overlay %s", utils.StylePath(targetPath))
	b.Cleanup(false)
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
func (b *BuildObject) buildInstallCmd() (cmd string, extraBindPaths []string, err error) {
	var quietFlag string
	if utils.QuietMode {
		quietFlag = "-q"
	}
	channelFlags := buildChannelFlags()

	if utils.IsYaml(b.buildSource) {
		// Mode 3: YAML file (-p prefix -f environment.yml)
		absFilePath, err := filepath.Abs(b.buildSource)
		if err != nil {
			return "", nil, fmt.Errorf("failed to get absolute path for %s: %w", b.buildSource, err)
		}
		extraBindPaths = []string{filepath.Dir(absFilePath)}
		cmd = fmt.Sprintf("micromamba create -r /ext3/tmp %s -y %s -p /cnt/%s -f %s",
			channelFlags, quietFlag, b.nameVersion, absFilePath)
		utils.PrintMessage("Populating overlay %s via %s",
			utils.StyleName(filepath.Base(b.targetOverlayPath)), "environment.yml")
	} else if b.buildSource != "" {
		// Mode 2: Multiple packages (-n name pkg1 pkg2 ...)
		packages := strings.Split(b.buildSource, ",")
		for i, pkg := range packages {
			packages[i] = strings.ReplaceAll(strings.TrimSpace(pkg), "/", "=")
		}
		cmd = fmt.Sprintf("micromamba create -r /ext3/tmp %s -y %s -p /cnt/%s %s",
			channelFlags, quietFlag, b.nameVersion, strings.Join(packages, " "))
		utils.PrintMessage("Populating overlay %s via %s",
			utils.StyleName(filepath.Base(b.targetOverlayPath)),
			fmt.Sprintf("micromamba (%s)", strings.Join(packages, ", ")))
	} else {
		// Mode 1: Single package (name/version)
		cmd = fmt.Sprintf("micromamba create -r /ext3/tmp %s -y %s -p /cnt/%s %s=%s",
			channelFlags, quietFlag, b.nameVersion, b.packageName, b.packageVersion)
		utils.PrintMessage("Populating overlay %s via %s",
			utils.StyleName(filepath.Base(b.targetOverlayPath)),
			fmt.Sprintf("micromamba (%s=%s)", b.packageName, b.packageVersion))
	}

	return cmd, extraBindPaths, nil
}

// buildCondaExecOpts constructs exec.Options for the conda build container run.
func (b *BuildObject) buildCondaExecOpts(finalPath string) (exec.Options, error) {
	installCmd, extraBindPaths, err := b.buildInstallCmd()
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
		finalPath, b.effectiveNcpus(), config.Global.Build.BlockSize, config.Global.Build.CompressArgs,
	)

	// Always bind the target overlay directory so mksquashfs can write there.
	bindPaths := append(extraBindPaths, filepath.Dir(b.targetOverlayPath))

	var opts exec.Options
	if !config.Global.Build.UseTmpOverlay {
		buildTmpDir := getBuildTmpDir(b)
		bindPaths = append(bindPaths,
			buildTmpDir+":/ext3/tmp",
			b.cntDirPath+":/cnt",
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
			Overlays:      []string{b.tmpOverlayPath},
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
func (b *BuildObject) runCondaBuild(ctx context.Context, finalPath string) error {
	opts, err := b.buildCondaExecOpts(finalPath)
	if err != nil {
		return err
	}

	utils.PrintDebug("[BUILD] Creating overlay for %s with options: overlays=%v, bindPaths=%v",
		b.nameVersion, opts.Overlays, opts.BindPaths)

	done := watchContext(ctx, "conda build")
	defer close(done)

	if err := exec.Run(ctx, opts); err != nil {
		if isCancelledByUser(err) {
			return ErrBuildCancelled
		}
		return fmt.Errorf("failed to build conda package %s: %w", b.nameVersion, err)
	}

	return nil
}
