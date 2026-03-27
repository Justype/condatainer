package build

import (
	"context"
	"fmt"
	"os"
	"path/filepath"
	"strings"

	"github.com/Justype/condatainer/internal/config"
	"github.com/Justype/condatainer/internal/container"
	execpkg "github.com/Justype/condatainer/internal/exec"
	"github.com/Justype/condatainer/internal/scheduler"
	"github.com/Justype/condatainer/internal/utils"
)

// ScriptBuildObject implements BuildObject for shell scripts
type ScriptBuildObject struct {
	*BaseBuildObject
	isRef bool // Only shell scripts can be ref overlays (for large reference datasets)
}

// newScriptBuildObject creates a ScriptBuildObject from base
func newScriptBuildObject(base *BaseBuildObject, isRef bool) (*ScriptBuildObject, error) {
	return &ScriptBuildObject{
		BaseBuildObject: base,
		isRef:           isRef,
	}, nil
}

// Type returns the build type
func (s *ScriptBuildObject) Type() BuildType {
	if s.isRef {
		return BuildTypeRef
	}
	return BuildTypeShell
}

// getAllBaseDirs returns all configured base directories that actually exist
func getAllBaseDirs() []string {
	var dirs []string

	// Helper to add directory if it exists
	addIfExists := func(dir string) {
		if dir == "" {
			return
		}
		if _, err := os.Stat(dir); err == nil {
			dirs = append(dirs, dir)
		}
	}

	// Explicit extra image dirs (strip :ro/:rw marker)
	for _, entry := range config.GetExtraImageDirs() {
		path, _ := config.ParseDirEntry(entry)
		addIfExists(path)
	}

	// Explicit extra build-scripts dirs (plain paths)
	for _, path := range config.GetExtraBuildDirs() {
		addIfExists(path)
	}

	// Extra root dir
	addIfExists(config.GetExtraRootDir())

	// Root dir
	addIfExists(config.GetRootDir())

	// Scratch dir
	addIfExists(config.GetScratchDataDir())

	// User dir
	addIfExists(config.GetUserDataDir())

	return dirs
}

// String returns a string representation
func (s *ScriptBuildObject) String() string {
	return fmt.Sprintf(`BuildObject:
		name_version: %s
		build_source_type: Shell Script
		build_source: %s
		dependencies: %v
		script_specs: %v
		tmp_overlay_path: %s
		target_overlay_path: %s
		cnt_dir_path: %s
		is_ref: %v`,
		s.nameVersion,
		s.buildSource,
		s.dependencies,
		s.scriptSpecs,
		s.tmpOverlayPath,
		s.targetOverlayPath,
		s.cntDirPath,
		s.isRef,
	)
}

// Build implements the shell script build workflow
// Workflow:
//  1. Check if overlay already exists (skip if yes)
//  2. Create temporary ext3 overlay or host dirs
//  3. Check and build missing dependencies if buildDeps=true
//  4. Run shell script inside container with overlays
//  5. Create SquashFS from the output directory
//  6. Extract and save ENV variables if present
//  7. Set permissions, atomic install, and cleanup
func (s *ScriptBuildObject) Build(ctx context.Context, buildDeps bool) error {
	targetPath, finalPath := buildOverlayPaths(s.BaseBuildObject)
	styledOverlay := utils.StyleName(filepath.Base(targetPath))

	if skip, err := checkShouldBuild(s.BaseBuildObject); skip || err != nil {
		return err
	}

	if err := s.createBuildLock(); err != nil {
		return err
	}
	defer s.removeBuildLock()

	utils.PrintMessage("Building overlay %s (%s build)", styledOverlay, utils.StyleAction(buildModeLabel(s.BaseBuildObject)))

	if err := prepareBuildWorkspace(ctx, s.BaseBuildObject, config.Global.Build.UseTmpOverlay); err != nil {
		return err
	}

	utils.PrintMessage("Populating overlay %s via %s", styledOverlay, utils.StyleAction(s.buildSource))

	if err := s.buildDependencies(ctx, buildDeps); err != nil {
		return err
	}

	if err := s.runBuildScript(ctx); err != nil {
		s.Cleanup(true)
		return err
	}

	if err := s.packOutput(ctx, finalPath); err != nil {
		// packOutput calls Cleanup(true) before returning an error.
		return err
	}

	s.saveEnvFile(targetPath)

	if err := os.Chmod(finalPath, utils.PermFile); err != nil {
		utils.PrintDebug("Failed to set permissions on %s: %v", finalPath, err)
	}

	if err := atomicInstall(finalPath, targetPath, s.update); err != nil {
		return err
	}

	utils.PrintSuccess("Finished overlay %s", utils.StylePath(targetPath))
	s.Cleanup(false)
	return nil
}

// buildDependencies checks for missing dependencies and optionally builds them.
func (s *ScriptBuildObject) buildDependencies(ctx context.Context, buildDeps bool) error {
	missingDeps, err := s.GetMissingDependencies()
	if err != nil {
		return fmt.Errorf("failed to check dependencies: %w", err)
	}
	if len(missingDeps) == 0 {
		return nil
	}

	styledOverlay := utils.StyleName(filepath.Base(s.targetOverlayPath))
	depList := strings.Join(missingDeps, ", ")

	if !buildDeps {
		utils.PrintError("Missing dependencies for %s: %s", styledOverlay, utils.StyleNumber(depList))
		return fmt.Errorf("missing dependencies for %s: %s. Please install them first", s.nameVersion, depList)
	}

	utils.PrintMessage("Building missing dependencies for %s: %s", styledOverlay, utils.StyleNumber(depList))

	writableImagesDir, err := config.GetWritableImagesDir()
	if err != nil {
		return fmt.Errorf("no writable images directory found: %w", err)
	}

	for _, dep := range missingDeps {
		depObj, err := NewBuildObject(ctx, dep, false, writableImagesDir, config.GetWritableTmpDir(), false)
		if err != nil {
			return fmt.Errorf("failed to create build object for dependency %s: %w", dep, err)
		}
		if err := depObj.Build(ctx, false); err != nil {
			return fmt.Errorf("failed to build dependency %s: %w", dep, err)
		}
	}

	utils.PrintSuccess("All dependencies for %s built successfully.", styledOverlay)
	return nil
}

// buildExecOpts constructs the exec.Options for running the build script inside the container.
func (s *ScriptBuildObject) buildExecOpts() (execpkg.Options, error) {
	// Parse name and version from nameVersion
	nameVersionParts := strings.Split(s.nameVersion, "/")
	name := nameVersionParts[0]
	version := "env"
	if len(nameVersionParts) >= 2 {
		version = nameVersionParts[1]
	}

	// Build resource spec: builds are always single-node/single-task.
	// effectiveNcpus() derives CPU count from scriptSpecs (folds TasksPerNode) or defaults.
	// Pass both mem fields so GetMemPerTaskMB() on buildRS correctly derives per-task memory.
	effRS := buildEffectiveResourceSpec(s.scriptSpecs)
	buildRS := &scheduler.ResourceSpec{
		Nodes:        1,
		TasksPerNode: 1,
		CpusPerTask:  s.effectiveNcpus(),
		MemPerCpuMB:  effRS.MemPerCpuMB,
		MemPerNodeMB: effRS.MemPerNodeMB,
	}

	// Build environment settings (KEY=VALUE format for exec.Options.EnvSettings)
	envSettings := append(
		scheduler.ResourceEnvVars(buildRS),
		fmt.Sprintf("app_name=%s", name),
		fmt.Sprintf("version=%s", version),
		fmt.Sprintf("app_name_version=%s", s.nameVersion),
		"tmp_dir=/ext3/tmp",
		"TMPDIR=/ext3/tmp",
		"IN_CONDATAINER=1",
	)

	// Set target_dir based on build mode and overlay type.
	// In ext3-mode ref builds: create cnt_dir so the script can write there.
	if !config.Global.Build.UseTmpOverlay {
		if s.isRef {
			targetDir := filepath.Join(s.cntDirPath, s.nameVersion)
			envSettings = append(envSettings, fmt.Sprintf("target_dir=%s", targetDir))
		} else {
			envSettings = append(envSettings, fmt.Sprintf("target_dir=/cnt/%s", s.nameVersion))
		}
	} else {
		if s.isRef {
			if err := os.MkdirAll(s.cntDirPath, utils.PermDir); err != nil {
				return execpkg.Options{}, fmt.Errorf("failed to create cnt_dir %s: %w", s.cntDirPath, err)
			}
			targetDir := filepath.Join(s.cntDirPath, s.nameVersion)
			envSettings = append(envSettings, fmt.Sprintf("target_dir=%s", targetDir))
		} else {
			envSettings = append(envSettings, fmt.Sprintf("target_dir=/cnt/%s", s.nameVersion))
		}
	}

	// Build bash command to execute script
	bashScript := fmt.Sprintf(`
trap 'exit 130' INT TERM

mkdir -p $TMPDIR
bash %s
if [ $? -ne 0 ]; then
    echo "Build script %s failed."
    exit 1
fi
`, s.buildSource, s.buildSource)

	// Collect overlay dependencies
	overlayArgs, err := GetOverlayArgsFromDependencies(s.dependencies)
	if err != nil {
		utils.PrintWarning("Failed to get overlay args from dependencies: %v", err)
	}

	// Build overlay list (dependency overlays; tmp overlay added first in ext3 mode)
	var overlays []string
	if config.Global.Build.UseTmpOverlay {
		overlays = []string{s.tmpOverlayPath}
	}
	for i := 0; i < len(overlayArgs); i += 2 {
		if overlayArgs[i] == "--overlay" && i+1 < len(overlayArgs) {
			overlays = append(overlays, overlayArgs[i+1])
		}
	}

	// Collect bind directories (for accessing overlays and build scripts)
	bindDirs := container.DeduplicateBindPaths(getAllBaseDirs())
	if !config.Global.Build.UseTmpOverlay {
		buildTmpDir := getBuildTmpDir(s.BaseBuildObject)
		if s.isRef {
			// Ref builds: self-bind cntDirPath so the container can write to the absolute host path.
			// /cnt is NOT bound, keeping it free for dependency overlays.
			bindDirs = append(bindDirs, buildTmpDir+":/ext3/tmp", s.cntDirPath+":"+s.cntDirPath)
		} else {
			// App builds: bind cntDirPath as /cnt (no dep overlay conflict for app builds)
			bindDirs = append(bindDirs, buildTmpDir+":/ext3/tmp", s.cntDirPath+":/cnt")
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

	// Enable stdin pass-through for interactive inputs or real-time user interaction
	if len(s.interactiveInputs) > 0 {
		inputStr := strings.Join(s.interactiveInputs, "\n") + "\n"
		opts.PassThruStdin = true
		opts.Stdin = strings.NewReader(inputStr)
	} else {
		opts.PassThruStdin = true
	}

	return opts, nil
}

// runBuildScript sets up a context watcher and runs the build script.
// The caller is responsible for calling Cleanup(true) if an error is returned.
func (s *ScriptBuildObject) runBuildScript(ctx context.Context) error {
	opts, err := s.buildExecOpts()
	if err != nil {
		return err
	}

	utils.PrintDebug("[BUILD] Running build script for %s with options: overlays=%v, bindPaths=%v, passThruStdin=%v",
		s.nameVersion, opts.Overlays, opts.BindPaths, opts.PassThruStdin)

	done := watchContext(ctx, "build script")
	defer close(done)

	if err := execpkg.Run(ctx, opts); err != nil {
		if isCancelledByUser(err) {
			return ErrBuildCancelled
		}
		return fmt.Errorf("build script %s failed: %w", s.buildSource, err)
	}

	return nil
}

// packOutput determines the squashfs source directory and creates the SquashFS file.
// Calls Cleanup(true) and returns an error if packing fails.
func (s *ScriptBuildObject) packOutput(ctx context.Context, finalPath string) error {
	// For ref builds, targetDir is the host path where the script writes files.
	var targetDir string
	if s.isRef {
		targetDir = filepath.Join(s.cntDirPath, s.nameVersion)
	}

	styledOverlay := utils.StyleName(filepath.Base(s.targetOverlayPath))

	if !config.Global.Build.UseTmpOverlay {
		// Dir mode: all build types pack from cntDirPath (host dir bound as /cnt)
		checkDir := s.cntDirPath
		if s.isRef && targetDir != "" {
			checkDir = targetDir
		}
		if entries, err := os.ReadDir(checkDir); err != nil || len(entries) == 0 {
			s.Cleanup(true)
			return fmt.Errorf("build script did not create any files in %s", checkDir)
		}
		utils.PrintMessage("Creating SquashFS from %s for overlay %s", utils.StylePath(s.cntDirPath), styledOverlay)
		if err := createSquashfs(ctx, s.BaseBuildObject, s.isRef, s.cntDirPath, finalPath); err != nil {
			s.Cleanup(true)
			return err
		}
	} else if s.isRef {
		// Ext3 mode, ref overlays: verify files exist, then pack from cnt_dir
		if entries, err := os.ReadDir(targetDir); err != nil || len(entries) == 0 {
			s.Cleanup(true)
			return fmt.Errorf("overlay build script did not create any files in %s", targetDir)
		}
		utils.PrintMessage("Creating SquashFS from %s for overlay %s", utils.StylePath(s.cntDirPath), styledOverlay)
		if err := createSquashfs(ctx, s.BaseBuildObject, s.isRef, s.cntDirPath, finalPath); err != nil {
			s.Cleanup(true)
			return err
		}
	} else {
		// Ext3 mode, app overlays: pack from /cnt inside container via overlay
		utils.PrintMessage("Preparing SquashFS from /cnt for overlay %s", styledOverlay)
		if err := createSquashfs(ctx, s.BaseBuildObject, s.isRef, "/cnt", finalPath); err != nil {
			s.Cleanup(true)
			return err
		}
	}

	return nil
}

// saveEnvFile extracts #ENV: declarations from the build script and writes the .env file.
// Errors are logged as warnings only — env file failures don't abort the build.
func (s *ScriptBuildObject) saveEnvFile(targetPath string) {
	envDict, err := GetEnvDictFromBuildScript(s.buildSource)
	if err != nil {
		utils.PrintWarning("Failed to extract ENV from build script: %v", err)
		return
	}
	if len(envDict) == 0 {
		return
	}
	if err := SaveEnvFile(targetPath, envDict, s.nameVersion); err != nil {
		utils.PrintWarning("Failed to save ENV file: %v", err)
	}
}
