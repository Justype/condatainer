package build

import (
	"errors"
	"fmt"
	"os"
	"os/exec"
	"path/filepath"
	"strings"

	"github.com/Justype/condatainer/internal/config"
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

// Type check implementations
func (s *ScriptBuildObject) IsConda() bool { return false }
func (s *ScriptBuildObject) IsDef() bool   { return false }
func (s *ScriptBuildObject) IsShell() bool { return true }
func (s *ScriptBuildObject) IsRef() bool   { return s.isRef }

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
func (s *ScriptBuildObject) Build(buildDeps bool) error {
	// Ensure target overlay path is absolute
	targetOverlayPath := s.targetOverlayPath
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
	if s.RequiresScheduler() {
		buildMode = "sbatch"
	}
	utils.PrintMessage("Building overlay %s (%s build)", styledOverlay, utils.StyleAction(buildMode))

	// Create temporary overlay
	utils.PrintMessage("Creating temporary overlay %s", utils.StylePath(s.tmpOverlayPath))
	if err := s.CreateTmpOverlay(false); err != nil {
		if errors.Is(err, ErrTmpOverlayExists) {
			utils.PrintError("Temporary overlay %s already exists. Cancel build.", utils.StylePath(s.tmpOverlayPath))
		}
		return fmt.Errorf("temporary overlay for %s already exists. Maybe a build is still running? %w", s.nameVersion, err)
	}

	utils.PrintMessage("Populating overlay %s via %s", styledOverlay, utils.StyleCommand(s.buildSource))

	// Check dependencies
	missingDeps, err := s.GetMissingDependencies()
	if err != nil {
		return fmt.Errorf("failed to check dependencies: %w", err)
	}

	if len(missingDeps) > 0 {
		depList := strings.Join(missingDeps, ", ")
		utils.PrintError("Missing dependencies for %s: %s", styledOverlay, utils.StyleNumber(depList))
		if buildDeps {
			// TODO: Build dependencies recursively
			return fmt.Errorf("missing dependencies: %s (recursive build not yet implemented)", depList)
		}
		return fmt.Errorf("missing dependencies for %s: %s. Please install them first", s.nameVersion, depList)
	}

	// Parse name and version
	nameVersionParts := strings.Split(s.nameVersion, "/")
	name := nameVersionParts[0]
	version := "env"
	if len(nameVersionParts) >= 2 {
		version = nameVersionParts[1]
	}

	// Build environment settings
	envSettings := []string{
		"--env", fmt.Sprintf("app_name=%s", name),
		"--env", fmt.Sprintf("version=%s", version),
		"--env", fmt.Sprintf("app_name_version=%s", s.nameVersion),
		"--env", "tmp_dir=/ext3/tmp",
		"--env", "TMPDIR=/ext3/tmp",
		"--env", fmt.Sprintf("SLURM_CPUS_PER_TASK=%d", s.ncpus),
		"--env", "IN_CONDATAINER=1",
	}

	// Set target_dir based on ref type
	targetDir := ""
	relativePath := s.nameVersion
	if s.isRef {
		// Create cnt_dir for ref overlays (large reference datasets)
		if err := os.MkdirAll(s.cntDirPath, 0o775); err != nil {
			return fmt.Errorf("failed to create cnt_dir %s: %w", s.cntDirPath, err)
		}
		targetDir = filepath.Join(s.cntDirPath, s.nameVersion)
		envSettings = append(envSettings, "--env", fmt.Sprintf("target_dir=%s", targetDir))
	} else {
		envSettings = append(envSettings, "--env", fmt.Sprintf("target_dir=/cnt/%s", relativePath))
	}

	// TODO: Add PATH and overlay env configs from dependencies
	// envSettings = append(envSettings, "--env", fmt.Sprintf("PATH=%s", getPathEnv(s.dependencies)))

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

	// Build exec command args using config values
	baseImage := config.Global.BaseImage
	apptainerBin := config.Global.ApptainerBin
	condatainerDir := config.Global.BaseDir

	args := []string{"exec"}
	args = append(args, envSettings...)
	// TODO: Add overlay args from dependencies
	args = append(args, "--overlay", s.tmpOverlayPath)
	if condatainerDir != "" {
		args = append(args, "--bind", condatainerDir)
	}
	// TODO: Add GPU args
	args = append(args, baseImage, "/bin/bash", "-c", bashScript)

	cmd := exec.Command(apptainerBin, args...)
	utils.PrintDebug("[BUILD] Running build script with command: %s %s", apptainerBin, strings.Join(args, " "))
	utils.PrintMessage("Executing build script %s for overlay %s", utils.StyleCommand(s.buildSource), styledOverlay)

	// Provide interactive inputs if any
	if len(s.interactiveInputs) > 0 {
		inputStr := strings.Join(s.interactiveInputs, "\n") + "\n"
		cmd.Stdin = strings.NewReader(inputStr)
	}

	// Run the build script - show output in real-time
	cmd.Stdout = os.Stdout
	cmd.Stderr = os.Stderr

	// Use signal-aware runner to ensure cleanup runs on Ctrl+C
	if err := runCommandWithSignalHandling(cmd); err != nil {
		s.Cleanup(true)
		if isCancelledByUser(err) {
			return ErrBuildCancelled
		}
		utils.PrintError("Build script %s failed for overlay %s: %v", utils.StyleCommand(s.buildSource), styledOverlay, err)
		return fmt.Errorf("build script %s failed: %w", s.buildSource, err)
	}

	// Create SquashFS
	if s.isRef {
		// For ref overlays: check files exist, then pack from cnt_dir
		if entries, err := os.ReadDir(targetDir); err != nil || len(entries) == 0 {
			utils.PrintError("No files produced for overlay %s", styledOverlay)
			s.Cleanup(true)
			return fmt.Errorf("overlay build script did not create any files in %s", targetDir)
		}

		// TODO: Set permissions recursively
		// Utils.share_to_ugo_recursive(s.cntDirPath)

		utils.PrintMessage("Creating SquashFS from %s for overlay %s", utils.StylePath(s.cntDirPath), styledOverlay)
		if err := s.createSquashfs(s.cntDirPath, targetOverlayPath); err != nil {
			utils.PrintError("Failed to create SquashFS for overlay %s: %v", styledOverlay, err)
			s.Cleanup(true)
			return err
		}
	} else {
		// For app overlays: set permissions and pack from /cnt
		utils.PrintMessage("Preparing SquashFS from /cnt for overlay %s", styledOverlay)
		if err := s.createSquashfs("/cnt", targetOverlayPath); err != nil {
			utils.PrintError("Failed to create SquashFS for overlay %s: %v", styledOverlay, err)
			s.Cleanup(true)
			return err
		}
	}

	// TODO: Extract and save ENV variables from build script
	// envDict := getEnvDictFromBuildScript(s.buildSource)
	// if len(envDict) > 0 {
	//     saveEnvFile(s.targetOverlayPath + ".env", envDict, relativePath)
	// }

	// Set permissions
	if err := os.Chmod(targetOverlayPath, 0o664); err != nil {
		utils.PrintDebug("Failed to set permissions on %s: %v", targetOverlayPath, err)
	}

	utils.PrintMessage("Finished overlay %s", styledOverlay)
	utils.PrintDebug("Overlay created at %s. Cleaning up...", utils.StylePath(targetOverlayPath))

	// Cleanup
	s.Cleanup(false)

	return nil
}

// createSquashfs creates a SquashFS file from the given source directory
func (s *ScriptBuildObject) createSquashfs(sourceDir, targetOverlayPath string) error {
	// Get config values
	baseImage := config.Global.BaseImage
	apptainerBin := config.Global.ApptainerBin
	condatainerDir := config.Global.BaseDir

	targetPath := targetOverlayPath
	if absTarget, err := filepath.Abs(targetPath); err == nil {
		targetPath = absTarget
	}

	// Build mksquashfs command
	bashScript := ""
	if sourceDir == "/cnt" {
		// For app overlays: set permissions first
		bashScript = fmt.Sprintf(`
trap 'exit 130' INT TERM
echo "Setting permissions..."
find /cnt -type f -exec chmod ug+rw,o+r {} \;
find /cnt -type d -exec chmod ug+rwx,o+rx {} \;
mksquashfs /cnt %s -processors %d -keep-as-directory %s -b 1M
`, targetPath, s.ncpus, config.Global.Build.CompressArgs)
	} else {
		// For ref overlays: just pack the directory
		bashScript = fmt.Sprintf(`
trap 'exit 130' INT TERM
mksquashfs %s %s -processors %d -keep-as-directory %s -b 1M
`, sourceDir, targetPath, s.ncpus, config.Global.Build.CompressArgs)
	}

	args := []string{
		"exec",
		"--overlay", s.tmpOverlayPath,
	}
	if condatainerDir != "" {
		args = append(args, "--bind", condatainerDir)
	}
	args = append(args, baseImage, "/bin/bash", "-c", bashScript)

	cmd := exec.Command(apptainerBin, args...)
	utils.PrintDebug("[BUILD] Creating SquashFS with command: %s %s", apptainerBin, strings.Join(args, " "))
	utils.PrintMessage("Packing %s into %s", utils.StylePath(sourceDir), utils.StylePath(targetPath))

	// Show output in real-time
	cmd.Stdout = os.Stdout
	cmd.Stderr = os.Stderr

	// Use signal-aware runner to ensure cleanup runs on Ctrl+C
	if err := runCommandWithSignalHandling(cmd); err != nil {
		if isCancelledByUser(err) {
			return ErrBuildCancelled
		}
		return fmt.Errorf("failed to create SquashFS: %w", err)
	}

	return nil
}
