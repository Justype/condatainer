package cmd

import (
	"context"
	"errors"
	"fmt"
	"os"
	"path/filepath"
	"strings"
	"time"

	"github.com/Justype/condatainer/internal/apptainer"
	"github.com/Justype/condatainer/internal/build"
	"github.com/Justype/condatainer/internal/config"
	execpkg "github.com/Justype/condatainer/internal/exec"
	"github.com/Justype/condatainer/internal/overlay"
	"github.com/Justype/condatainer/internal/scheduler"
	"github.com/Justype/condatainer/internal/utils"
	"github.com/spf13/cobra"
)

var (
	runWritableImg bool
	runBaseImage   string
	runAutoInstall bool
	runEnvSettings []string
	runBindPaths   []string
	runFakeroot    bool
)

var runCmd = &cobra.Command{
	Use:   "run [script]",
	Short: "Run a script and auto-solve the dependencies by #DEP tags",
	Long: `Execute a script with automatic dependency resolution.

The script can contain special comment tags:
  #DEP: package/version  - Declares a dependency
  #CNT args              - Additional arguments to pass to condatainer

Supported #CNT arguments:
  -w, --writable         Make .img overlays writable
  -b, --base-image PATH  Use custom base image
  --env KEY=VALUE        Set environment variable
  --bind HOST:CONTAINER  Bind mount path
  -f, --fakeroot         Run with fakeroot privileges

Dependencies will be automatically loaded, and missing ones can be auto-installed
with the --auto-install flag.`,
	Example: `  condatainer run script.sh              # Run with dependency check
  condatainer run script.sh -a           # Auto-install missing deps
  condatainer run script.sh -w           # Make .img overlays writable
  condatainer run script.sh -b base.sif  # Use custom base image

  # In script.sh, you can use #CNT to set options:
  # #CNT -w
  # #CNT --env MYVAR=value
  # #CNT --bind /data:/mnt/data`,
	Args:         cobra.ExactArgs(1),
	SilenceUsage: true, // Runtime errors should not show usage
	RunE:         runScript,
}

func init() {
	rootCmd.AddCommand(runCmd)
	runCmd.Flags().BoolVarP(&runWritableImg, "writable", "w", false, "Make .img overlays writable (default: read-only)")
	runCmd.Flags().BoolVar(&runWritableImg, "writable-img", false, "Alias for --writable")
	runCmd.Flags().StringVarP(&runBaseImage, "base-image", "b", "", "Base image to use instead of default")
	runCmd.Flags().BoolVarP(&runAutoInstall, "auto-install", "a", false, "Automatically install missing dependencies")
	runCmd.Flags().BoolP("install", "i", false, "Alias for --auto-install")
}

func runScript(cmd *cobra.Command, args []string) error {
	scriptPath := args[0]

	// Check if -i flag was used (alias for -a)
	if installFlag, _ := cmd.Flags().GetBool("install"); installFlag {
		runAutoInstall = true
	}

	// Ensure base image exists
	if err := apptainer.EnsureBaseImage(cmd.Context()); err != nil {
		return err
	}

	// Validate script exists
	if !utils.FileExists(scriptPath) {
		utils.PrintError("Script file %s not found", utils.StylePath(scriptPath))
		return nil
	}

	originScriptPath := scriptPath
	scriptPath, err := filepath.Abs(scriptPath)
	if err != nil {
		utils.PrintError("Failed to get absolute path: %v", err)
		return nil
	}

	// Parse additional arguments from script (#CNT tags)
	scriptArgs, err := parseArgsInScript(scriptPath)
	if err != nil {
		return err
	}

	// Apply #CNT arguments
	if len(scriptArgs) > 0 {
		utils.PrintDebug("[RUN] Additional script arguments found: %v", scriptArgs)
		applyScriptArgs(scriptArgs)
	}

	// Get dependencies from script
	deps, err := utils.GetDependenciesFromScript(scriptPath)
	if err != nil {
		utils.PrintError("Failed to parse dependencies: %v", err)
		return nil
	}

	// Collect build job IDs from auto-install step (populated later if auto-install runs)
	var buildJobIDs []string

	// Check for missing dependencies
	installedOverlays, err := getInstalledOverlaysMap()
	if err != nil {
		return err
	}

	missingDeps := []string{}
	for _, dep := range deps {
		// Check if it's an external overlay file by extension
		isExternalOverlay := utils.IsOverlay(dep)

		if isExternalOverlay {
			// External overlay file - check if file exists
			if !utils.FileExists(dep) {
				missingDeps = append(missingDeps, dep)
			}
		} else {
			// Package name - check if installed
			normalized := utils.NormalizeNameVersion(dep)
			if _, ok := installedOverlays[normalized]; !ok {
				missingDeps = append(missingDeps, dep)
			}
		}
	}

	// Handle missing dependencies
	if len(missingDeps) > 0 {
		// Always show missing dependencies first
		utils.PrintMessage("Missing dependencies:")
		for _, md := range missingDeps {
			utils.PrintMessage("  - %s", utils.StyleWarning(md))
		}

		if !runAutoInstall {
			externalFiles := []string{}
			for _, md := range missingDeps {
				if utils.IsOverlay(md) {
					externalFiles = append(externalFiles, md)
				}
			}

			if len(externalFiles) > 0 {
				fileList := strings.Join(externalFiles, ", ")
				utils.PrintError("External overlay(s) %s not found. Cannot use -a to install.", utils.StyleName(fileList))
				return nil
			}

			utils.PrintHint("Please run %s to install missing dependencies.", utils.StyleAction("condatainer check -a "+originScriptPath))
			return nil
		}

		// Auto-install missing dependencies
		utils.PrintMessage("Attempting to auto-install missing dependencies...")

		// Check for external sqf/overlay files - these can't be auto-installed
		externalFiles := []string{}
		for _, md := range missingDeps {
			if utils.IsOverlay(md) {
				externalFiles = append(externalFiles, md)
			}
		}

		if len(externalFiles) > 0 {
			fileList := strings.Join(externalFiles, ", ")
			utils.PrintError("External overlay(s) %s not found - cannot proceed with auto-installation", utils.StyleName(fileList))
			return nil
		}

		// Create build objects
		imagesDir, err := config.GetWritableImagesDir()
		if err != nil {
			utils.PrintError("No writable images directory found: %v", err)
			return nil
		}
		buildObjects := make([]build.BuildObject, 0, len(missingDeps))
		for _, pkg := range missingDeps {
			bo, err := build.NewBuildObject(pkg, false, imagesDir, config.GetWritableTmpDir())
			if err != nil {
				utils.PrintError("Failed to create build object for %s: %v", utils.StyleName(pkg), err)
				return nil
			}
			buildObjects = append(buildObjects, bo)
		}

		// Build graph and execute
		graph, err := build.NewBuildGraph(buildObjects, imagesDir, config.GetWritableTmpDir(), config.Global.SubmitJob)
		if err != nil {
			utils.PrintError("Failed to create build graph: %v", err)
			return nil
		}

		if err := graph.Run(cmd.Context()); err != nil {
			if errors.Is(err, build.ErrBuildCancelled) ||
				errors.Is(cmd.Context().Err(), context.Canceled) ||
				strings.Contains(err.Error(), "signal: killed") ||
				strings.Contains(err.Error(), "context canceled") {
				utils.PrintWarning("Auto-installation cancelled.")
				return nil
			}
			utils.PrintError("Some overlays failed to install: %v", err)
			return nil
		}

		utils.PrintSuccess("All selected overlays installed/submitted.")

		// Collect scheduler build job IDs (if any) so run job can depend on them
		for _, id := range graph.GetJobIDs() {
			buildJobIDs = append(buildJobIDs, id)
		}
		if len(buildJobIDs) > 0 {
			utils.PrintMessage("Build job(s) submitted: %s", strings.Join(buildJobIDs, ", "))
		}

		// Re-fetch installed overlays after installation
		installedOverlays, err = getInstalledOverlaysMap()
		if err != nil {
			return err
		}
	}

	// Resolve overlay paths and check availability early
	overlays := make([]string, len(deps))
	for i, dep := range deps {
		if utils.IsOverlay(dep) {
			overlays[i] = dep
		} else {
			normalized := utils.NormalizeNameVersion(dep)
			if path, ok := installedOverlays[normalized]; ok {
				overlays[i] = path
			} else {
				overlays[i] = dep
			}
		}
	}

	// Before submitting or running locally, check if .img overlays are available.
	// This prevents submitting jobs that will immediately fail due to locked images.
	for _, ol := range overlays {
		if utils.IsImg(ol) && utils.FileExists(ol) {
			// If writable requested, check for exclusive lock availability.
			if err := overlay.CheckAvailable(ol, runWritableImg); err != nil {
				return err
			}
		}
	}

	// Check if script has scheduler specs and scheduler is available
	scriptSpecs, err := scheduler.ReadScriptSpecsFromPath(scriptPath)
	if err != nil {
		utils.PrintWarning("Failed to parse scheduler specs: %v", err)
	}

	// If script has scheduler specs, try to submit as a job
	if scriptSpecs != nil && len(scriptSpecs.RawFlags) > 0 {
		// Check if job submission is enabled in config
		if !config.Global.SubmitJob {
			utils.PrintNote("Job submission is disabled. Running locally.")
		} else if scheduler.IsInsideJob() {
			// Already inside a scheduler job, run locally (silently)
		} else {
			// Try to detect and use scheduler
			sched, err := scheduler.DetectScheduler()
			if err != nil {
				utils.PrintNote("Script has scheduler specs but scheduler detection failed: %v. Running locally.", err)
			} else if !sched.IsAvailable() {
				utils.PrintNote("Script has scheduler specs but no scheduler is available. Running locally.")
			} else {
				// Scheduler available - submit job
				return submitRunJob(sched, scriptPath, originScriptPath, scriptSpecs, buildJobIDs)
			}
		}
	} else {
		// No scheduler specs found - explicitly note local execution
		utils.PrintNote("No scheduler specs found in script. Running locally.")
	}

	// Prepare execution command
	// Disable module commands and run the script
	executionScript := `module() { :; }
ml() { :; }
export -f module ml
`

	// Check if script is executable
	fileInfo, err := os.Stat(scriptPath)
	if err != nil {
		utils.PrintError("Failed to stat script: %v", err)
		return nil
	}

	if fileInfo.Mode()&0111 != 0 {
		// Script is executable
		executionScript += scriptPath
	} else {
		// Run with bash
		executionScript += "bash " + scriptPath
	}

	// Execute with overlays
	options := execpkg.Options{
		Overlays:     overlays,
		Command:      []string{"/bin/bash", "-c", executionScript},
		WritableImg:  runWritableImg,
		EnvSettings:  runEnvSettings,
		BindPaths:    runBindPaths,
		Fakeroot:     runFakeroot,
		BaseImage:    runBaseImage, // Empty string triggers GetBaseImage() in ensureDefaults()
		ApptainerBin: config.Global.ApptainerBin,
		HidePrompt:   true,
	}

	if err := execpkg.Run(cmd.Context(), options); err != nil {
		// Propagate exit code from container command
		if appErr, ok := err.(*apptainer.ApptainerError); ok {
			if code := appErr.ExitCode(); code >= 0 {
				os.Exit(code)
			}
		}
		return err
	}
	return nil
}

// applyScriptArgs parses #CNT arguments and applies them to run options
func applyScriptArgs(scriptArgs []string) {
	for _, argLine := range scriptArgs {
		parts := strings.Fields(argLine)
		for i := 0; i < len(parts); i++ {
			arg := parts[i]
			switch arg {
			case "-w", "--writable", "--writable-img":
				runWritableImg = true
			case "-f", "--fakeroot":
				runFakeroot = true
			case "-b", "--base-image":
				if i+1 < len(parts) {
					i++
					runBaseImage = parts[i]
				}
			case "--env":
				if i+1 < len(parts) {
					i++
					runEnvSettings = append(runEnvSettings, parts[i])
				}
			case "--bind":
				if i+1 < len(parts) {
					i++
					runBindPaths = append(runBindPaths, parts[i])
				}
			default:
				// Handle --env=VALUE and --bind=VALUE formats
				if strings.HasPrefix(arg, "--env=") {
					runEnvSettings = append(runEnvSettings, strings.TrimPrefix(arg, "--env="))
				} else if strings.HasPrefix(arg, "--bind=") {
					runBindPaths = append(runBindPaths, strings.TrimPrefix(arg, "--bind="))
				} else if strings.HasPrefix(arg, "-b=") {
					runBaseImage = strings.TrimPrefix(arg, "-b=")
				} else if strings.HasPrefix(arg, "--base-image=") {
					runBaseImage = strings.TrimPrefix(arg, "--base-image=")
				}
			}
		}
	}
}

// parseArgsInScript extracts arguments from #CNT comments in the script
func parseArgsInScript(scriptPath string) ([]string, error) {
	data, err := os.ReadFile(scriptPath)
	if err != nil {
		return nil, fmt.Errorf("failed to read script: %w", err)
	}

	args := []string{}
	lines := strings.Split(string(data), "\n")

	for _, line := range lines {
		line = strings.TrimSpace(line)
		if strings.HasPrefix(line, "#CNT") {
			argLine := strings.TrimSpace(strings.TrimPrefix(line, "#CNT"))
			if argLine != "" {
				args = append(args, argLine)
			}
		}
	}

	return args, nil
}

// submitRunJob creates and submits a scheduler job to run the script
func submitRunJob(sched scheduler.Scheduler, scriptPath, originScriptPath string, specs *scheduler.ScriptSpecs, depIDs []string) error {
	info := sched.GetInfo()
	if len(depIDs) > 0 {
		utils.PrintMessage("Submitting %s job to run %s (will wait for build job(s): %s)", info.Type, utils.StylePath(originScriptPath), strings.Join(depIDs, ", "))
	} else {
		utils.PrintMessage("Submitting %s job to run %s", info.Type, utils.StylePath(originScriptPath))
	}

	// Determine log directory - use spec's Stdout path if set, otherwise global log path
	var logsDir string
	if specs.Stdout != "" {
		// Ensure absolute path
		logPath := specs.Stdout
		if !filepath.IsAbs(logPath) {
			if abs, err := filepath.Abs(logPath); err == nil {
				logPath = abs
			}
		}
		logsDir = filepath.Dir(logPath)
	} else {
		// Use global log path
		logsDir = config.Global.LogsDir
		if logsDir == "" {
			logsDir = filepath.Join(os.Getenv("HOME"), "logs")
		}
	}

	// Ensure logs directory exists
	if err := os.MkdirAll(logsDir, 0775); err != nil {
		return fmt.Errorf("failed to create logs directory %s: %w", logsDir, err)
	}

	// Build the condatainer run command (without -a to avoid re-installing in job)
	runCommand := fmt.Sprintf("condatainer run %s", scriptPath)

	// Generate names: short name for job, timestamped name for files
	baseName := filepath.Base(originScriptPath)
	ext := filepath.Ext(baseName)
	nameWithoutExt := strings.TrimSuffix(baseName, ext)
	timestamp := time.Now().Format("20060102_150405")
	fileBaseName := fmt.Sprintf("%s_%s", nameWithoutExt, timestamp)

	// Set job name if not already set in script (short name for squeue)
	if specs.JobName == "" {
		specs.JobName = nameWithoutExt
	}

	// Create job specification (Name used for file naming)
	absScriptPath, err := filepath.Abs(originScriptPath)
	if err != nil {
		absScriptPath = originScriptPath
	}
	jobSpec := &scheduler.JobSpec{
		Name:      fileBaseName,
		Command:   runCommand,
		Specs:     specs,
		DepJobIDs: depIDs,
		Metadata: map[string]string{
			"Script": absScriptPath,
		},
	}

	// Create the job script in the same directory as logs
	jobScriptPath, err := sched.CreateScriptWithSpec(jobSpec, logsDir)
	if err != nil {
		return fmt.Errorf("failed to create job script: %w", err)
	}

	// Submit the job (attach dependency job IDs so run waits for builds)
	jobID, err := sched.Submit(jobScriptPath, depIDs)
	if err != nil {
		return fmt.Errorf("failed to submit job: %w", err)
	}

	utils.PrintSuccess("Submitted %s job %s for %s", info.Type, utils.StyleNumber(jobID), utils.StylePath(originScriptPath))

	// Show script and log file locations
	// utils.PrintHint("Job script: %s", utils.StylePath(jobScriptPath))
	logFile := filepath.Join(logsDir, fmt.Sprintf("%s.log", fileBaseName))
	utils.PrintHint("Log file: %s", utils.StylePath(logFile))

	return nil
}
