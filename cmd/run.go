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

// errRunAborted signals a handled stop (message already printed); caller returns nil.
var errRunAborted = errors.New("run aborted")

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
	if err := apptainer.EnsureBaseImage(cmd.Context(), false, false); err != nil {
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

	// 1. Fast path: already inside a scheduler job — skip spec reading, run locally
	if scheduler.IsInsideJob() {
		if err := processEmbeddedArgs(scriptPath); err != nil {
			return err
		}
		overlays, _, err := checkDepsAndAutoInstall(cmd.Context(), scriptPath, originScriptPath)
		if err != nil {
			if errors.Is(err, errRunAborted) {
				return nil
			}
			return err
		}
		return runLocally(cmd.Context(), scriptPath, overlays)
	}

	// 2. Read specs → resolve the bash script to use for #CNT/#DEP
	contentScript, scriptSpecs := resolveScriptAndSpecs(scriptPath)

	// 3. Embedded args + dependency check/install
	if err := processEmbeddedArgs(contentScript); err != nil {
		return err
	}
	overlays, buildJobIDs, err := checkDepsAndAutoInstall(cmd.Context(), contentScript, originScriptPath)
	if err != nil {
		if errors.Is(err, errRunAborted) {
			return nil
		}
		return err
	}

	// 4. Submit or run locally
	if config.Global.SubmitJob && scheduler.HasSchedulerSpecs(scriptSpecs) {
		if err := validateAndConvertSpecs(scriptSpecs); err != nil {
			return nil
		}
		sched, schedErr := scheduler.DetectScheduler()
		if schedErr != nil {
			utils.PrintNote("Script has scheduler specs but scheduler detection failed: %v. Running locally.", schedErr)
		} else if !sched.IsAvailable() {
			utils.PrintNote("Script has scheduler specs but no scheduler is available. Running locally.")
		} else {
			return submitRunJob(sched, originScriptPath, contentScript, scriptSpecs, buildJobIDs)
		}
	} else if !scheduler.HasSchedulerSpecs(scriptSpecs) {
		utils.PrintNote("No scheduler specs found in script. Running locally.")
	}
	return runLocally(cmd.Context(), contentScript, overlays)
}

// resolveScriptAndSpecs tries to read scheduler specs and resolves the content script path.
// For HTCondor .sub files: specs.ScriptPath is the bash executable — returned as contentScript.
// For SLURM/PBS/LSF: contentScript == scriptPath.
// On failure: prints a warning and returns (scriptPath, nil).
func resolveScriptAndSpecs(scriptPath string) (string, *scheduler.ScriptSpecs) {
	specs, err := scheduler.ReadScriptSpecsFromPath(scriptPath)
	if err != nil {
		utils.PrintWarning("Failed to parse scheduler specs: %v", err)
		return scriptPath, nil
	}
	contentScript := scriptPath
	if specs != nil && specs.ScriptPath != "" && specs.ScriptPath != scriptPath {
		p := specs.ScriptPath
		if !filepath.IsAbs(p) {
			// Resolve relative executable paths against the directory of the submit file
			p = filepath.Join(filepath.Dir(scriptPath), p)
		}
		contentScript = p
	}
	return contentScript, specs
}

// processEmbeddedArgs reads #CNT tags from the script and applies them to run options.
func processEmbeddedArgs(scriptPath string) error {
	scriptArgs, err := parseArgsInScript(scriptPath)
	if err != nil {
		return err
	}
	if len(scriptArgs) > 0 {
		utils.PrintDebug("[RUN] Additional script arguments found: %v", scriptArgs)
		applyScriptArgs(scriptArgs)
	}
	return nil
}

// checkDepsAndAutoInstall resolves #DEP dependencies, checks installed overlays, and optionally
// auto-installs missing ones. Returns overlay paths and build job IDs.
// Returns errRunAborted (message already printed) for handled stop conditions.
func checkDepsAndAutoInstall(ctx context.Context, contentScript, originScriptPath string) (overlays []string, buildJobIDs []string, err error) {
	deps, err := utils.GetDependenciesFromScript(contentScript)
	if err != nil {
		utils.PrintError("Failed to parse dependencies: %v", err)
		return nil, nil, errRunAborted
	}

	installedOverlays, err := getInstalledOverlaysMap()
	if err != nil {
		return nil, nil, err
	}

	missingDeps := []string{}
	for _, dep := range deps {
		if utils.IsOverlay(dep) {
			if !utils.FileExists(dep) {
				missingDeps = append(missingDeps, dep)
			}
		} else {
			normalized := utils.NormalizeNameVersion(dep)
			if _, ok := installedOverlays[normalized]; !ok {
				missingDeps = append(missingDeps, dep)
			}
		}
	}

	if len(missingDeps) > 0 {
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
				return nil, nil, errRunAborted
			}
			utils.PrintHint("Please run %s to install missing dependencies.", utils.StyleAction("condatainer check -a "+originScriptPath))
			return nil, nil, errRunAborted
		}

		// Auto-install missing dependencies
		utils.PrintMessage("Attempting to auto-install missing dependencies...")

		externalFiles := []string{}
		for _, md := range missingDeps {
			if utils.IsOverlay(md) {
				externalFiles = append(externalFiles, md)
			}
		}
		if len(externalFiles) > 0 {
			fileList := strings.Join(externalFiles, ", ")
			utils.PrintError("External overlay(s) %s not found - cannot proceed with auto-installation", utils.StyleName(fileList))
			return nil, nil, errRunAborted
		}

		imagesDir, err := config.GetWritableImagesDir()
		if err != nil {
			utils.PrintError("No writable images directory found: %v", err)
			return nil, nil, errRunAborted
		}
		buildObjects := make([]build.BuildObject, 0, len(missingDeps))
		for _, pkg := range missingDeps {
			bo, err := build.NewBuildObject(pkg, false, imagesDir, config.GetWritableTmpDir())
			if err != nil {
				utils.PrintError("Failed to create build object for %s: %v", utils.StyleName(pkg), err)
				return nil, nil, errRunAborted
			}
			buildObjects = append(buildObjects, bo)
		}

		graph, err := build.NewBuildGraph(buildObjects, imagesDir, config.GetWritableTmpDir(), config.Global.SubmitJob)
		if err != nil {
			utils.PrintError("Failed to create build graph: %v", err)
			return nil, nil, errRunAborted
		}

		if err := graph.Run(ctx); err != nil {
			if errors.Is(err, build.ErrBuildCancelled) ||
				errors.Is(ctx.Err(), context.Canceled) ||
				strings.Contains(err.Error(), "signal: killed") ||
				strings.Contains(err.Error(), "context canceled") {
				utils.PrintWarning("Auto-installation cancelled.")
				return nil, nil, errRunAborted
			}
			utils.PrintError("Some overlays failed to install: %v", err)
			return nil, nil, errRunAborted
		}

		utils.PrintSuccess("All selected overlays installed/submitted.")
		for _, id := range graph.GetJobIDs() {
			buildJobIDs = append(buildJobIDs, id)
		}
		if len(buildJobIDs) > 0 {
			utils.PrintMessage("Build job(s) submitted: %s", strings.Join(buildJobIDs, ", "))
		}

		installedOverlays, err = getInstalledOverlaysMap()
		if err != nil {
			return nil, nil, err
		}
	}

	// Resolve overlay paths
	overlays = make([]string, len(deps))
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

	// Check .img availability before submitting or running locally
	for _, ol := range overlays {
		if utils.IsImg(ol) && utils.FileExists(ol) {
			if err := overlay.CheckAvailable(ol, runWritableImg); err != nil {
				return nil, nil, err
			}
		}
	}

	return overlays, buildJobIDs, nil
}

// runLocally executes the script directly via apptainer with the given overlays.
func runLocally(ctx context.Context, contentScript string, overlays []string) error {
	// Disable module commands and run the script
	executionScript := `module() { :; }
ml() { :; }
export -f module ml
`
	fileInfo, err := os.Stat(contentScript)
	if err != nil {
		utils.PrintError("Failed to stat script: %v", err)
		return nil
	}

	if fileInfo.Mode()&0111 != 0 {
		// Script is executable
		executionScript += contentScript
	} else {
		// Run with bash
		executionScript += "bash " + contentScript
	}

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

	if err := execpkg.Run(ctx, options); err != nil {
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

// validateAndConvertSpecs validates job specs and converts GPU/CPU if needed using the scheduler package.
func validateAndConvertSpecs(specs *scheduler.ScriptSpecs) error {
	validationErr, cpuAdjusted, cpuMsg, gpuConverted, gpuMsg := scheduler.ValidateAndConvertSpecs(specs)

	if validationErr != nil {
		// Validation failed
		utils.PrintError("Job specs validation failed: %v", validationErr)

		// If it's a GPU validation error with suggestions, print them
		if gpuErr, ok := validationErr.(*scheduler.GpuValidationError); ok && len(gpuErr.Suggestions) > 0 {
			utils.PrintMessage("Available GPU options:")
			for _, suggestion := range gpuErr.Suggestions {
				utils.PrintMessage("  - %s", suggestion)
			}
		}

		return validationErr
	}

	if cpuAdjusted {
		utils.PrintWarning("%s", cpuMsg)
	}

	if gpuConverted {
		utils.PrintWarning("%s", gpuMsg)
	}

	return nil
}

// submitRunJob creates and submits a scheduler job to run the script.
// contentScript is the bash script containing #DEP/#CNT directives — for HTCondor this is the
// executable referenced in the .sub file, for other schedulers it equals scriptPath.
func submitRunJob(sched scheduler.Scheduler, originScriptPath, contentScript string, specs *scheduler.ScriptSpecs, depIDs []string) error {
	info := sched.GetInfo()
	if len(depIDs) > 0 {
		utils.PrintMessage("Submitting %s job to run %s (will wait for build job(s): %s)", info.Type, utils.StylePath(originScriptPath), strings.Join(depIDs, ", "))
	} else {
		utils.PrintMessage("Submitting %s job to run %s", info.Type, utils.StylePath(originScriptPath))
	}

	// Determine log directory - use spec's Stdout path if set, otherwise global log path
	var logsDir string
	if specs.Control.Stdout != "" {
		// Ensure absolute path
		logPath := specs.Control.Stdout
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
	if err := os.MkdirAll(logsDir, utils.PermDir); err != nil {
		return fmt.Errorf("failed to create logs directory %s: %w", logsDir, err)
	}

	// Build the condatainer run command pointing at the bash script.
	// When the user provided a .sub file, contentScript is the executable .sh — the submitted
	// job must re-invoke condatainer run on the bash script, not the submit file.
	runCommand := fmt.Sprintf("condatainer run %s", contentScript)

	// Generate names: short name for job, timestamped name for files
	baseName := filepath.Base(originScriptPath)
	ext := filepath.Ext(baseName)
	nameWithoutExt := strings.TrimSuffix(baseName, ext)
	timestamp := time.Now().Format("20060102_150405")
	fileBaseName := fmt.Sprintf("%s_%s", nameWithoutExt, timestamp)

	// Set job name if not already set in script (short name for squeue)
	if specs.Control.JobName == "" {
		specs.Control.JobName = nameWithoutExt
	}

	jobSpec := &scheduler.JobSpec{
		Name:      fileBaseName,
		Command:   runCommand,
		Specs:     specs,
		DepJobIDs: depIDs,
		Metadata:  map[string]string{},
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

	stdoutPath := jobSpec.Specs.Control.Stdout
	stderrPath := jobSpec.Specs.Control.Stderr
	if stdoutPath != "" {
		if stderrPath == "" || stdoutPath == stderrPath {
			utils.PrintHint("Stdout & Stderr => %s", utils.StylePath(stdoutPath))
		} else {
			utils.PrintHint("Stdout => %s", utils.StylePath(stdoutPath))
			utils.PrintHint("Stderr => %s", utils.StylePath(stderrPath))
		}
	}

	return nil
}
