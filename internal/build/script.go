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
    sbatch_flags: %v
    tmp_overlay_path: %s
    target_overlay_path: %s
    cnt_dir_path: %s
    is_ref: %v`,
		s.nameVersion,
		s.buildSource,
		s.dependencies,
		s.sbatchFlags,
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
	// Check if already exists
	if _, err := os.Stat(s.targetOverlayPath); err == nil {
		utils.PrintDebug("Overlay %s already exists at %s. Skipping creation.",
			filepath.Base(s.targetOverlayPath), utils.StylePath(s.targetOverlayPath))
		return nil
	}

	// Create temporary overlay
	if err := s.CreateTmpOverlay(false); err != nil {
		return fmt.Errorf("temporary overlay for %s already exists. Maybe a build is still running?", s.nameVersion)
	}

	utils.PrintDebug("Building SquashFS overlay at %s using build source %s",
		utils.StylePath(s.targetOverlayPath), utils.StylePath(s.buildSource))

	// Check dependencies
	missingDeps, err := s.GetMissingDependencies()
	if err != nil {
		return fmt.Errorf("failed to check dependencies: %w", err)
	}

	if len(missingDeps) > 0 {
		if buildDeps {
			utils.PrintDebug("Building missing dependencies for %s: %s",
				s.nameVersion, strings.Join(missingDeps, ", "))
			// TODO: Build dependencies recursively
			// For now, just return error
			return fmt.Errorf("missing dependencies: %s (recursive build not yet implemented)", strings.Join(missingDeps, ", "))
		} else {
			return fmt.Errorf("missing dependencies for %s: %s. Please install them first",
				s.nameVersion, strings.Join(missingDeps, ", "))
		}
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
		"--env", "IN_CONDATINER=1",
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

	// Provide interactive inputs if any
	if len(s.interactiveInputs) > 0 {
		inputStr := strings.Join(s.interactiveInputs, "\n") + "\n"
		cmd.Stdin = strings.NewReader(inputStr)
	}

	// Run the build script - show output in real-time
	cmd.Stdout = os.Stdout
	cmd.Stderr = os.Stderr

	if err := cmd.Run(); err != nil {
		s.Cleanup(true)
		return fmt.Errorf("build script %s failed: %w", s.buildSource, err)
	}

	// Create SquashFS
	if s.isRef {
		// For ref overlays: check files exist, then pack from cnt_dir
		if entries, err := os.ReadDir(targetDir); err != nil || len(entries) == 0 {
			s.Cleanup(true)
			return fmt.Errorf("overlay build script did not create any files in %s", targetDir)
		}

		// TODO: Set permissions recursively
		// Utils.share_to_ugo_recursive(s.cntDirPath)

		utils.PrintDebug("Creating SquashFS file from %s", utils.StylePath(s.cntDirPath))
		if err := s.createSquashfs(s.cntDirPath); err != nil {
			s.Cleanup(true)
			return err
		}
	} else {
		// For app overlays: set permissions and pack from /cnt
		utils.PrintDebug("Creating SquashFS file from /cnt")
		if err := s.createSquashfs("/cnt"); err != nil {
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
	if err := os.Chmod(s.targetOverlayPath, 0o664); err != nil {
		utils.PrintDebug("Failed to set permissions on %s: %v", s.targetOverlayPath, err)
	}

	utils.PrintDebug("Overlay created at %s. Cleaning up...", utils.StylePath(s.targetOverlayPath))

	// Cleanup
	s.Cleanup(false)

	return nil
}

// createSquashfs creates a SquashFS file from the given source directory
func (s *ScriptBuildObject) createSquashfs(sourceDir string) error {
	// Get config values
	baseImage := config.Global.BaseImage
	apptainerBin := config.Global.ApptainerBin
	condatainerDir := config.Global.BaseDir

	// Build mksquashfs command
	bashScript := ""
	if sourceDir == "/cnt" {
		// For app overlays: set permissions first
		bashScript = fmt.Sprintf(`
echo "Setting permissions..."
find /cnt -type f -exec chmod ug+rw,o+r {} \;
find /cnt -type d -exec chmod ug+rwx,o+rx {} \;
mksquashfs /cnt %s -processors %d -keep-as-directory -comp gzip -b 1M
`, filepath.Base(s.targetOverlayPath), s.ncpus)
	} else {
		// For ref overlays: just pack the directory
		bashScript = fmt.Sprintf(`
mksquashfs %s %s -processors %d -keep-as-directory -comp gzip -b 1M
`, sourceDir, filepath.Base(s.targetOverlayPath), s.ncpus)
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

	// Show output in real-time
	cmd.Stdout = os.Stdout
	cmd.Stderr = os.Stderr

	if err := cmd.Run(); err != nil {
		return fmt.Errorf("failed to create SquashFS: %w", err)
	}

	return nil
}
