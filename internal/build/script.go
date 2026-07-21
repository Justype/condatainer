package build

import (
	"context"
	"fmt"
	"os"
	"path/filepath"
	"strings"

	"log/slog"

	"github.com/Justype/condatainer/internal/config"
	"github.com/Justype/condatainer/internal/container"
	execpkg "github.com/Justype/condatainer/internal/exec"
	"github.com/Justype/condatainer/internal/logging"
	"github.com/Justype/condatainer/internal/scheduler"
	"github.com/Justype/condatainer/internal/utils"
)

// getAllBaseDirs returns all configured base directories that actually exist
func getAllBaseDirs() []string {
	var dirs []string

	addIfExists := func(dir string) {
		if dir == "" {
			return
		}
		if _, err := os.Stat(dir); err == nil {
			dirs = append(dirs, dir)
		}
	}

	for _, entry := range config.GetExtraImageDirs() {
		path, _ := config.ParseDirEntry(entry)
		addIfExists(path)
	}
	for _, path := range config.GetExtraBuildDirs() {
		addIfExists(path)
	}
	addIfExists(config.GetExtraRootDir())
	addIfExists(config.GetRootDir())
	addIfExists(config.GetScratchDataDir())
	addIfExists(config.GetUserDataDir())

	return dirs
}

// buildScript implements the shell script build workflow on BuildObject.
// Workflow:
//  1. Check if overlay already exists (skip if yes)
//  2. Create temporary ext3 overlay or host dirs
//  3. Check and build missing dependencies if buildDeps=true
//  4. Run shell script inside container with overlays
//  5. Create SquashFS from the output directory
//  6. Extract and save ENV variables if present
//  7. Set permissions, atomic install, and cleanup
func (b *BuildObject) buildScript(ctx context.Context, buildDeps bool) error {
	targetPath, finalPath := buildOverlayPaths(b)
	log := logging.FromContext(ctx)

	if skip, err := checkShouldBuild(b); skip || err != nil {
		return err
	}

	if err := b.createBuildLock(); err != nil {
		return err
	}
	defer b.removeBuildLock()

	log.Info("building overlay", "overlay", filepath.Base(targetPath), "mode", buildModeLabel(b))

	if err := prepareBuildWorkspace(ctx, b, config.Global.Build.UseTmpOverlay); err != nil {
		return err
	}

	log.Info("populating overlay", "overlay", filepath.Base(targetPath), "source", b.buildSource)

	if err := b.buildDependencies(ctx, buildDeps); err != nil {
		return err
	}

	if err := b.runBuildScript(ctx); err != nil {
		b.Cleanup(true)
		return err
	}

	if err := b.packOutput(ctx, finalPath); err != nil {
		return err
	}

	b.saveEnvFile(targetPath)

	utils.ShareWithParentGroup(finalPath)

	if err := atomicInstall(finalPath, targetPath, b.update); err != nil {
		return err
	}

	log.Info("overlay ready", "kind", "success", "path", targetPath)
	b.Cleanup(false)
	return nil
}

// buildDependencies checks for missing dependencies and optionally builds them.
func (b *BuildObject) buildDependencies(ctx context.Context, buildDeps bool) error {
	missingDeps, err := b.GetMissingDependencies()
	if err != nil {
		return fmt.Errorf("failed to check dependencies: %w", err)
	}
	if len(missingDeps) == 0 {
		return nil
	}

	depList := strings.Join(missingDeps, ", ")

	if !buildDeps {
		logging.FromContext(ctx).Error("missing dependencies", "overlay", filepath.Base(b.targetOverlayPath), "deps", depList)
		return fmt.Errorf("missing dependencies for %s: %s. Please install them first", b.nameVersion, depList)
	}

	logging.FromContext(ctx).Info("building missing dependencies", "overlay", filepath.Base(b.targetOverlayPath), "deps", depList)

	writableImagesDir, err := config.GetWritableImagesDir()
	if err != nil {
		return fmt.Errorf("no writable images directory found: %w", err)
	}

	for _, dep := range missingDeps {
		// Strip any version constraint (e.g. "samtools/1.21>=1.16" → "samtools/1.21")
		// before passing to NewBuildObject, which expects a plain name/version.
		preferredNV, _, _ := utils.SplitDepConstraint(dep)
		depObj, err := NewBuildObject(ctx, preferredNV, false, writableImagesDir, config.GetWritableTmpDir(), false)
		if err != nil {
			return fmt.Errorf("failed to create build object for dependency %s: %w", preferredNV, err)
		}
		if err := depObj.Build(ctx, false); err != nil {
			return fmt.Errorf("failed to build dependency %s: %w", preferredNV, err)
		}
	}

	logging.FromContext(ctx).Info("all dependencies built", "kind", "success", "overlay", filepath.Base(b.targetOverlayPath))
	return nil
}

// buildExecOpts constructs the exec.Options for running the build script inside the container.
func (b *BuildObject) buildExecOpts() (execpkg.Options, execpkg.IO, error) {
	nameVersionParts := strings.Split(b.nameVersion, "/")
	name := nameVersionParts[0]
	version := "env"
	if len(nameVersionParts) >= 2 {
		version = nameVersionParts[1]
	}

	effRS := buildEffectiveResourceSpec(b.scriptSpecs)
	buildRS := &scheduler.ResourceSpec{
		Nodes:        1,
		TasksPerNode: 1,
		CpusPerTask:  b.effectiveNcpus(),
		MemPerCpuMB:  effRS.MemPerCpuMB,
		MemPerNodeMB: effRS.MemPerNodeMB,
	}

	envSettings := append(
		scheduler.ResourceEnvVars(buildRS),
		fmt.Sprintf("app_name=%s", name),
		fmt.Sprintf("version=%s", version),
		fmt.Sprintf("app_name_version=%s", b.nameVersion),
		"tmp_dir=/ext3/tmp",
		"TMPDIR=/ext3/tmp",
		"IN_CONDATAINER=1",
	)

	isRef := b.buildType == BuildTypeRef

	if !config.Global.Build.UseTmpOverlay {
		if isRef {
			targetDir := filepath.Join(b.cntDirPath, b.nameVersion)
			envSettings = append(envSettings, fmt.Sprintf("target_dir=%s", targetDir))
		} else {
			envSettings = append(envSettings, fmt.Sprintf("target_dir=/cnt/%s", b.nameVersion))
		}
	} else {
		if isRef {
			if err := utils.MkdirAllShared(b.cntDirPath); err != nil {
				return execpkg.Options{}, execpkg.IO{}, fmt.Errorf("failed to create cnt_dir %s: %w", b.cntDirPath, err)
			}
			targetDir := filepath.Join(b.cntDirPath, b.nameVersion)
			envSettings = append(envSettings, fmt.Sprintf("target_dir=%s", targetDir))
		} else {
			envSettings = append(envSettings, fmt.Sprintf("target_dir=/cnt/%s", b.nameVersion))
		}
	}

	bashScript := fmt.Sprintf(`
trap 'exit 130' INT TERM

mkdir -p $TMPDIR
bash %s
if [ $? -ne 0 ]; then
    echo "Build script %s failed."
    exit 1
fi
`, b.buildSource, b.buildSource)

	overlayArgs, err := GetOverlayArgsFromDependencies(b.dependencies)
	if err != nil {
		slog.Default().Warn("failed to get overlay args from dependencies", "err", err)
	}

	var overlays []string
	if config.Global.Build.UseTmpOverlay {
		overlays = []string{b.tmpOverlayPath}
	}
	for i := 0; i < len(overlayArgs); i += 2 {
		if overlayArgs[i] == "--overlay" && i+1 < len(overlayArgs) {
			overlays = append(overlays, overlayArgs[i+1])
		}
	}

	bindDirs := container.DeduplicateBindPaths(getAllBaseDirs())
	if !config.Global.Build.UseTmpOverlay {
		buildTmpDir := getBuildTmpDir(b)
		if isRef {
			bindDirs = append(bindDirs, buildTmpDir+":/ext3/tmp", b.cntDirPath+":"+b.cntDirPath)
		} else {
			bindDirs = append(bindDirs, buildTmpDir+":/ext3/tmp", b.cntDirPath+":/cnt")
		}
	}

	opts := execpkg.Options{
		BaseImage:    config.GetBaseImage(),
		ApptainerBin: config.Global.ApptainerBin,
		Overlays:     overlays,
		BindPaths:    bindDirs,
		EnvSettings:  envSettings,
		Command:      []string{"/bin/bash", "-c", bashScript},
		HidePrompt:   true,
		WritableImg:  config.Global.Build.UseTmpOverlay,
	}
	if !config.Global.Build.UseTmpOverlay {
		opts.ApptainerFlags = []string{"--writable-tmpfs"}
	}

	var ioStreams execpkg.IO
	if len(b.interactiveInputs) > 0 {
		inputStr := strings.Join(b.interactiveInputs, "\n") + "\n"
		opts.PassThruStdin = true
		ioStreams.Stdin = strings.NewReader(inputStr)
	} else {
		opts.PassThruStdin = true
	}

	return opts, ioStreams, nil
}

// runBuildScript sets up a context watcher and runs the build script.
func (b *BuildObject) runBuildScript(ctx context.Context) error {
	if len(b.vars) > 0 {
		origSource := b.buildSource
		subPath, err := substituteTemplateFile(origSource, b.vars, b.tmpDir)
		if err != nil {
			return fmt.Errorf("failed to substitute placeholders in build script: %w", err)
		}
		// Keep origSource alive so saveEnvFile can read #ENV: from it after this
		// function returns. Cleanup() will remove it via isRemote. Only the
		// substituted tmp script (subPath) needs to be removed here.
		defer func() {
			os.Remove(subPath) //nolint:errcheck
			b.buildSource = origSource
		}()
		b.buildSource = subPath
	}

	opts, ioStreams, err := b.buildExecOpts()
	if err != nil {
		return err
	}

	logging.FromContext(ctx).Debug("running build script",
		"name", b.nameVersion, "overlays", opts.Overlays, "bindPaths", opts.BindPaths, "passThruStdin", opts.PassThruStdin)

	done := watchContext(ctx, "build script")
	defer close(done)

	if ioStreams.IsZero() {
		ioStreams = execpkg.IOFromContext(ctx)
	}
	if err := execpkg.Run(ctx, opts, ioStreams); err != nil {
		if isCancelledByUser(err) {
			return ErrBuildCancelled
		}
		return fmt.Errorf("build script %s failed: %w", b.buildSource, err)
	}

	return nil
}

// packOutput determines the squashfs source directory and creates the SquashFS file.
func (b *BuildObject) packOutput(ctx context.Context, finalPath string) error {
	isRef := b.buildType == BuildTypeRef
	var targetDir string
	if isRef {
		targetDir = filepath.Join(b.cntDirPath, b.nameVersion)
	}

	log := logging.FromContext(ctx)

	if !config.Global.Build.UseTmpOverlay {
		checkDir := b.cntDirPath
		if isRef && targetDir != "" {
			checkDir = targetDir
		}
		if entries, err := os.ReadDir(checkDir); err != nil || len(entries) == 0 {
			b.Cleanup(true)
			return fmt.Errorf("build script did not create any files in %s", checkDir)
		}
		log.Info("creating SquashFS", "source", b.cntDirPath, "overlay", filepath.Base(b.targetOverlayPath))
		if err := createSquashfs(ctx, b, isRef, b.cntDirPath, finalPath); err != nil {
			b.Cleanup(true)
			return err
		}
	} else if isRef {
		if entries, err := os.ReadDir(targetDir); err != nil || len(entries) == 0 {
			b.Cleanup(true)
			return fmt.Errorf("overlay build script did not create any files in %s", targetDir)
		}
		log.Info("creating SquashFS", "source", b.cntDirPath, "overlay", filepath.Base(b.targetOverlayPath))
		if err := createSquashfs(ctx, b, isRef, b.cntDirPath, finalPath); err != nil {
			b.Cleanup(true)
			return err
		}
	} else {
		log.Info("preparing SquashFS from /cnt", "overlay", filepath.Base(b.targetOverlayPath))
		if err := createSquashfs(ctx, b, isRef, "/cnt", finalPath); err != nil {
			b.Cleanup(true)
			return err
		}
	}

	return nil
}

// saveEnvFile extracts #ENV: declarations and #WHATIS: from the build script and writes the .env file.
func (b *BuildObject) saveEnvFile(targetPath string) {
	envDict, err := GetEnvDictFromBuildScript(b.buildSource)
	if err != nil {
		slog.Default().Warn("failed to extract ENV from build script", "err", err)
		return
	}
	whatis := utils.GetWhatIsFromScript(b.buildSource)
	if len(b.vars) > 0 {
		whatis = utils.InterpolateVars(whatis, b.vars)
		for key, entry := range envDict {
			entry.Value = utils.InterpolateVars(entry.Value, b.vars)
			entry.Note = utils.InterpolateVars(entry.Note, b.vars)
			envDict[key] = entry
		}
	}
	if err := SaveEnvFile(targetPath, envDict, b.nameVersion, whatis); err != nil {
		slog.Default().Warn("failed to save ENV file", "err", err)
	}
}

// saveCondaEnvFile fetches a package description from anaconda.org and writes the .env file.
// Skipped for multi-package / YAML builds where packageName is a user-chosen env name, not a real package.
func (b *BuildObject) saveCondaEnvFile(targetPath string) {
	if b.packageName == "" || b.buildSource != "" {
		return
	}
	whatis := utils.FetchCondaSummary(b.packageName, config.Global.Build.Channels)
	if err := SaveEnvFile(targetPath, map[string]EnvEntry{}, b.nameVersion, whatis); err != nil {
		slog.Default().Warn("failed to save ENV file", "err", err)
	}
}
