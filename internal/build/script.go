package build

import (
	"context"
	"errors"
	"fmt"
	"os"
	"path/filepath"
	"strings"
	"sync"

	"github.com/Justype/condatainer/internal/config"
	"github.com/Justype/condatainer/internal/container"
	execpkg "github.com/Justype/condatainer/internal/exec"
	"github.com/Justype/condatainer/internal/overlay"
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

	// Extra base dirs
	for _, baseDir := range config.GetExtraBaseDirs() {
		addIfExists(baseDir)
	}

	// Portable dir
	addIfExists(config.GetPortableDataDir())

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
//  2. Create temporary ext3 overlay
//  3. Check and build missing dependencies if buildDeps=true
//  4. Run shell script inside container with overlays
//  5. For ref overlays: verify files exist, create SquashFS from cnt_dir
//  6. For app overlays: set permissions, create SquashFS from /cnt
//  7. Extract and save ENV variables if present
//  8. Set permissions and cleanup
func (s *ScriptBuildObject) Build(ctx context.Context, buildDeps bool) error {
	// Ensure target overlay path is absolute
	targetOverlayPath := s.targetOverlayPath
	if absTarget, err := filepath.Abs(targetOverlayPath); err == nil {
		targetOverlayPath = absTarget
	}

	styledOverlay := utils.StyleName(filepath.Base(targetOverlayPath))

	// Check if already exists (skip only when not in update mode)
	if !s.update {
		if _, err := os.Stat(targetOverlayPath); err == nil {
			utils.PrintMessage("Overlay %s already exists at %s. Skipping creation.",
				styledOverlay, utils.StylePath(targetOverlayPath))
			return nil
		}
	}

	// Before starting any build work, check that no exec/run holds the existing overlay.
	if s.update && utils.FileExists(targetOverlayPath) {
		if lock, err := overlay.AcquireLock(targetOverlayPath, true); err != nil {
			return fmt.Errorf("cannot update %s: %w", s.nameVersion, err)
		} else {
			lock.Close()
		}
	}

	// Acquire build lock to prevent concurrent builds of the same package.
	if err := s.createBuildLock(); err != nil {
		return err
	}
	defer s.removeBuildLock()

	buildMode := "local"
	if s.RequiresScheduler() {
		buildMode = "sbatch"
	}
	utils.PrintMessage("Building overlay %s (%s build)", styledOverlay, utils.StyleAction(buildMode))

	// Create build workspace (host dirs or ext3 overlay depending on config)
	if !config.Global.Build.UseTmpOverlay {
		utils.PrintMessage("Creating build directories for %s", styledOverlay)
		if err := s.CreateBuildDirs(ctx, false); err != nil {
			if errors.Is(err, ErrTmpOverlayExists) {
				utils.PrintWarning("Stale build directory found for %s. Cleaning up...", s.nameVersion)
				s.Cleanup(true)
				if err := s.CreateBuildDirs(ctx, false); err != nil {
					return fmt.Errorf("failed to create build dirs: %w", err)
				}
			} else {
				return fmt.Errorf("failed to create build dirs: %w", err)
			}
		}
	} else {
		utils.PrintMessage("Creating temporary overlay %s", utils.StylePath(s.tmpOverlayPath))
		if err := s.CreateTmpOverlay(ctx, false); err != nil {
			if errors.Is(err, ErrTmpOverlayExists) {
				utils.PrintWarning("Stale temporary overlay found for %s. Cleaning up...", s.nameVersion)
				s.Cleanup(true)
				if err := s.CreateTmpOverlay(ctx, false); err != nil {
					return fmt.Errorf("failed to create temporary overlay: %w", err)
				}
			} else {
				return fmt.Errorf("failed to create temporary overlay: %w", err)
			}
		}
	}

	utils.PrintMessage("Populating overlay %s via %s", styledOverlay, utils.StyleCommand(s.buildSource))

	// Check dependencies
	missingDeps, err := s.GetMissingDependencies()
	if err != nil {
		return fmt.Errorf("failed to check dependencies: %w", err)
	}

	if len(missingDeps) > 0 {
		depList := strings.Join(missingDeps, ", ")
		if buildDeps {
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
		} else {
			utils.PrintError("Missing dependencies for %s: %s", styledOverlay, utils.StyleNumber(depList))
			return fmt.Errorf("missing dependencies for %s: %s. Please install them first", s.nameVersion, depList)
		}
	}

	// Parse name and version
	nameVersionParts := strings.Split(s.nameVersion, "/")
	name := nameVersionParts[0]
	version := "env"
	if len(nameVersionParts) >= 2 {
		version = nameVersionParts[1]
	}

	// Build resource spec: builds are always single-node/single-task.
	// effectiveNcpus() derives CPU count from scriptSpecs (folds TasksPerNode) or defaults.
	// Pass both mem fields so GetMemPerTaskMB() on buildRS correctly derives per-task memory:
	//   MemPerCpuMB × effectiveNcpus()  (when --mem-per-cpu was specified)
	//   MemPerNodeMB / 1                (when --mem was specified; TasksPerNode=1)
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

	// Set target_dir. In ext3-mode, ref overlays use the absolute host path (bound into
	// container via getAllBaseDirs). In dir-mode, ref builds also use the absolute host path
	// to avoid the /cnt bind shadowing dependency overlay content at /cnt. App builds in
	// dir-mode use the container-side /cnt/<nameVersion> path (no dep overlay conflict).
	relativePath := s.nameVersion
	var targetDir string // host path used for post-build checks (ref builds only)
	if !config.Global.Build.UseTmpOverlay {
		if s.isRef {
			// Ref builds: absolute host path so /cnt is not bound and dep overlays are visible
			targetDir = filepath.Join(s.cntDirPath, s.nameVersion)
			envSettings = append(envSettings, fmt.Sprintf("target_dir=%s", targetDir))
		} else {
			// App builds: container-side /cnt path (cntDirPath bound as /cnt)
			envSettings = append(envSettings, fmt.Sprintf("target_dir=/cnt/%s", relativePath))
		}
	} else {
		if s.isRef {
			// Create cnt_dir for ref overlays (large reference datasets)
			if err := os.MkdirAll(s.cntDirPath, 0o775); err != nil {
				return fmt.Errorf("failed to create cnt_dir %s: %w", s.cntDirPath, err)
			}
			targetDir = filepath.Join(s.cntDirPath, s.nameVersion)
			envSettings = append(envSettings, fmt.Sprintf("target_dir=%s", targetDir))
		} else {
			envSettings = append(envSettings, fmt.Sprintf("target_dir=/cnt/%s", relativePath))
		}
	}

	// Add PATH from dependencies
	pathEnv, err := GetPathEnvFromDependencies(s.dependencies)
	if err != nil {
		utils.PrintWarning("Failed to get PATH from dependencies: %v", err)
	} else if pathEnv != "" {
		// Prepend dependency paths to system paths
		pathEnv = pathEnv + ":/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin"
		envSettings = append(envSettings, fmt.Sprintf("PATH=%s", pathEnv))
	}

	// Add overlay env configs from dependencies (already in KEY=VALUE format)
	envConfigs, err := GetOverlayEnvConfigsFromDependencies(s.dependencies)
	if err != nil {
		utils.PrintWarning("Failed to get env configs from dependencies: %v", err)
	} else {
		envSettings = append(envSettings, envConfigs...)
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

	// Convert overlay args to overlay paths (dependency overlays only; tmp overlay added below if ext3 mode)
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
			// Ref builds: self-bind cntDirPath so the container can write to the absolute host
			// path. /cnt is NOT bound, keeping it free for dependency overlays.
			bindDirs = append(bindDirs, buildTmpDir+":/ext3/tmp", s.cntDirPath+":"+s.cntDirPath)
		} else {
			// App builds: bind cntDirPath as /cnt (no dep overlay conflict for app builds)
			bindDirs = append(bindDirs, buildTmpDir+":/ext3/tmp", s.cntDirPath+":/cnt")
		}
	}

	// Create exec options for running build script
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

	// Enable stdin pass-through if there are interactive inputs
	// This allows build scripts to read from stdin for user input (download links, etc.)
	if len(s.interactiveInputs) > 0 {
		// Pre-collected inputs: provide them as stdin
		inputStr := strings.Join(s.interactiveInputs, "\n") + "\n"
		opts.PassThruStdin = true
		opts.Stdin = strings.NewReader(inputStr)
	} else {
		// Allow container to read from user's stdin for real-time interaction
		opts.PassThruStdin = true
	}

	utils.PrintDebug("[BUILD] Running build script for %s with options: overlays=%v, bindPaths=%v, passThruStdin=%v",
		s.nameVersion, opts.Overlays, opts.BindPaths, opts.PassThruStdin)

	// Set up signal handling and cleanup for build script phase
	var cleanupOnce sync.Once
	cleanupFunc := func() {
		cleanupOnce.Do(func() {
			utils.PrintMessage("Cleaning up temporary files for %s...", utils.StyleName(s.nameVersion))
			s.Cleanup(true)
		})
	}

	done := make(chan struct{})
	defer close(done)

	// Monitor signals in background during build script execution
	go func() {
		select {
		case <-ctx.Done():
			utils.PrintWarning("Build cancelled. Interrupting build...")
			// Do not call cleanupFunc() here to avoid race condition with process termination
			// cleanupFunc() will be called after execpkg.Run returns
		case <-done:
			return
		}
	}()

	// Run build script
	if err := execpkg.Run(ctx, opts); err != nil {
		cleanupFunc()
		if isCancelledByUser(err) {
			return ErrBuildCancelled
		}
		return fmt.Errorf("build script %s failed: %w", s.buildSource, err)
	}

	// Build succeeded - stop the signal handler before createSquashfs
	// This prevents the cleanup goroutine from removing tmpOverlayPath
	// while createSquashfs has it mounted, which would cause a Bus error

	// Determine final output path: write to .new when updating for atomic replacement
	finalPath := targetOverlayPath
	if s.update {
		finalPath = targetOverlayPath + ".new"
	}

	// Create SquashFS
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
		utils.PrintMessage("Setting permissions for %s", utils.StylePath(s.cntDirPath))
		if err := utils.ShareToUGORecursive(s.cntDirPath); err != nil {
			utils.PrintWarning("Failed to set permissions for %s: %v", utils.StylePath(s.cntDirPath), err)
		}
		utils.PrintMessage("Creating SquashFS from %s for overlay %s", utils.StylePath(s.cntDirPath), styledOverlay)
		if err := s.createSquashfs(ctx, s.cntDirPath, finalPath); err != nil {
			s.Cleanup(true)
			return err
		}
	} else if s.isRef {
		// Ext3 mode, ref overlays: check files exist, then pack from cnt_dir
		if entries, err := os.ReadDir(targetDir); err != nil || len(entries) == 0 {
			s.Cleanup(true)
			return fmt.Errorf("overlay build script did not create any files in %s", targetDir)
		}

		// Set permissions recursively
		utils.PrintMessage("Setting permissions recursively for %s", utils.StylePath(s.cntDirPath))
		if err := utils.ShareToUGORecursive(s.cntDirPath); err != nil {
			utils.PrintWarning("Failed to set permissions for %s: %v", utils.StylePath(s.cntDirPath), err)
		}

		utils.PrintMessage("Creating SquashFS from %s for overlay %s", utils.StylePath(s.cntDirPath), styledOverlay)
		if err := s.createSquashfs(ctx, s.cntDirPath, finalPath); err != nil {
			s.Cleanup(true)
			return err
		}
	} else {
		// Ext3 mode, app overlays: set permissions and pack from /cnt (inside container via overlay)
		utils.PrintMessage("Preparing SquashFS from /cnt for overlay %s", styledOverlay)
		if err := s.createSquashfs(ctx, "/cnt", finalPath); err != nil {
			s.Cleanup(true)
			return err
		}
	}

	// Extract and save ENV variables from build script (always written directly to targetOverlayPath.env)
	envDict, err := GetEnvDictFromBuildScript(s.buildSource)
	if err != nil {
		utils.PrintWarning("Failed to extract ENV from build script: %v", err)
	} else if len(envDict) > 0 {
		if err := SaveEnvFile(targetOverlayPath, envDict, relativePath); err != nil {
			utils.PrintWarning("Failed to save ENV file: %v", err)
		}
	}

	// Set permissions on output file
	if err := os.Chmod(finalPath, 0o664); err != nil {
		utils.PrintDebug("Failed to set permissions on %s: %v", finalPath, err)
	}

	// Atomic replacement: remove old and rename .new → target
	if s.update {
		os.Remove(targetOverlayPath) //nolint:errcheck
		if err := os.Rename(finalPath, targetOverlayPath); err != nil {
			os.Remove(finalPath) //nolint:errcheck
			return fmt.Errorf("failed to replace overlay %s: %w", targetOverlayPath, err)
		}
	}

	utils.PrintSuccess("Finished overlay %s", utils.StylePath(targetOverlayPath))

	// Cleanup
	s.Cleanup(false)

	return nil
}

// createSquashfs creates a SquashFS file from the given source directory
func (s *ScriptBuildObject) createSquashfs(ctx context.Context, sourceDir, targetOverlayPath string) error {
	targetPath := targetOverlayPath
	if absTarget, err := filepath.Abs(targetPath); err == nil {
		targetPath = absTarget
	}

	// Set up signal handling for SquashFS creation phase
	var cleanupOnce sync.Once
	cleanupFunc := func() {
		cleanupOnce.Do(func() {
			utils.PrintMessage("Cleaning up temporary files for %s...", utils.StyleName(s.nameVersion))
			s.Cleanup(true)
			// Remove partial SquashFS file if it exists
			if _, err := os.Stat(targetPath); err == nil {
				os.Remove(targetPath)
			}
		})
	}

	done := make(chan struct{})
	defer close(done)

	// Monitor signals in background during SquashFS creation
	go func() {
		select {
		case <-ctx.Done():
			utils.PrintWarning("Build cancelled during SquashFS creation. Cleaning up...")
			// Do not call cleanupFunc() here to avoid race condition with process termination
			// cleanupFunc() will be called after execpkg.Run returns
		case <-done:
			return
		}
	}()

	// Build mksquashfs command and exec options
	var bashScript string
	var packOverlays []string
	var packBindDirs []string

	if !config.Global.Build.UseTmpOverlay {
		// Dir mode: sourceDir is a host path (cntDirPath); permissions already set by caller.
		// Bind sourceDir and output dir at their host paths so mksquashfs can access them.
		blockSize := config.Global.Build.BlockSize
		if s.isRef {
			blockSize = config.Global.Build.DataBlockSize
		}
		bashScript = fmt.Sprintf(`
trap 'exit 130' INT TERM
echo "Packing overlay to SquashFS..."
mksquashfs %s %s -processors %d -b %s -keep-as-directory %s
`, sourceDir, targetPath, s.effectiveNcpus(), blockSize, config.Global.Build.CompressArgs)
		packOverlays = []string{}
		packBindDirs = append(container.DeduplicateBindPaths(getAllBaseDirs()), sourceDir, filepath.Dir(targetPath))
	} else if sourceDir == "/cnt" {
		// Ext3 mode, app overlays: set permissions inside container then pack from /cnt
		bashScript = fmt.Sprintf(`
trap 'exit 130' INT TERM
echo "Setting permissions..."
find /cnt -type f -exec chmod ug+rw,o+r {} \;
find /cnt -type d -exec chmod ug+rwx,o+rx {} \;

echo "Packing overlay to SquashFS..."
mksquashfs /cnt %s -processors %d -b %s -keep-as-directory %s
`, targetPath, s.effectiveNcpus(), config.Global.Build.BlockSize, config.Global.Build.CompressArgs)
		packOverlays = []string{s.tmpOverlayPath}
		packBindDirs = container.DeduplicateBindPaths(getAllBaseDirs())
	} else {
		// Ext3 mode, ref overlays: fix permissions on host, pack from host path
		if err := utils.FixPermissionsDefault(sourceDir); err != nil {
			utils.PrintWarning("Failed to fix permissions on %s: %v", sourceDir, err)
		}
		bashScript = fmt.Sprintf(`
trap 'exit 130' INT TERM
echo "Packing overlay to SquashFS..."
mksquashfs %s %s -processors %d -b %s -keep-as-directory %s
`, sourceDir, targetPath, s.effectiveNcpus(), config.Global.Build.DataBlockSize, config.Global.Build.CompressArgs)
		packOverlays = []string{}
		packBindDirs = container.DeduplicateBindPaths(getAllBaseDirs())
	}

	// Create exec options for creating SquashFS
	opts := execpkg.Options{
		BaseImage:    config.GetBaseImage(),
		ApptainerBin: config.Global.ApptainerBin,
		Overlays:     packOverlays,
		BindPaths:    packBindDirs,
		Command:      []string{"/bin/bash", "-c", bashScript},
		HidePrompt:   true,
		WritableImg:  config.Global.Build.UseTmpOverlay,
	}
	if !config.Global.Build.UseTmpOverlay {
		opts.ApptainerFlags = []string{"--writable-tmpfs"}
	}

	utils.PrintDebug("[BUILD] Creating SquashFS for %s with options: overlays=%v, bindPaths=%v",
		s.nameVersion, opts.Overlays, opts.BindPaths)
	utils.PrintMessage("Packing %s into %s", utils.StylePath(sourceDir), utils.StylePath(targetPath))

	// Run the SquashFS creation using exec.Run
	if err := execpkg.Run(ctx, opts); err != nil {
		cleanupFunc()
		if isCancelledByUser(err) {
			return ErrBuildCancelled
		}
		return fmt.Errorf("failed to create SquashFS: %w", err)
	}

	return nil
}
