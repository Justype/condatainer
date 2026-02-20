package cmd

import (
	"context"
	"errors"
	"fmt"
	"os"
	osexec "os/exec"
	"path/filepath"
	"sort"
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
	runWritableImg     bool
	runBaseImage       string
	runAutoInstall     bool
	runEnvSettings     []string
	runBindPaths       []string
	runFakeroot        bool
	runParseModuleLoad bool
)

// errRunAborted signals a handled stop (message already printed); caller returns nil.
var errRunAborted = errors.New("run aborted")

var runCmd = &cobra.Command{
	Use:   "run [script]",
	Short: "Run a script and auto-solve the dependencies by #DEP tags",
	Long: `Execute a script with automatic dependency resolution.

The script can contain special comment tags:
  #DEP: package/version   - Declares a dependency
  #DEP: /path/overlay.img - Declares an overlay dependency (.sqf and .img)
  #CNT args               - Additional arguments to pass to condatainer

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
  condatainer run script.sh -b base.sif  # Use custom base image`,
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
	runCmd.Flags().BoolVar(&runParseModuleLoad, "module", false, "Also parse 'module load' / 'ml' lines as dependencies")
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

	// 1. Read specs → resolve the content script (HTCondor: .sub → .sh; others: identity)
	contentScript, scriptSpecs := resolveScriptAndSpecs(scriptPath)

	// 2. Embedded #CNT args + dependency check/install
	if err := processEmbeddedArgs(contentScript); err != nil {
		return err
	}
	if runWritableImg && getNtasks(scriptSpecs) > 1 {
		utils.PrintError("--writable cannot be used with multi-task jobs (ntasks=%d)", getNtasks(scriptSpecs))
		return nil
	}
	overlays, buildJobIDs, err := checkDepsAndAutoInstall(cmd.Context(), contentScript, originScriptPath)
	if err != nil {
		if errors.Is(err, errRunAborted) {
			return nil
		}
		return err
	}

	// 3. In-job or in-container: always run locally (no nested submission)
	if scheduler.IsInsideJob() || config.IsInsideContainer() {
		return runLocally(cmd.Context(), contentScript, overlays, scriptSpecs)
	}

	// 4. Submit if scheduler specs present and scheduler available
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
	return runLocally(cmd.Context(), contentScript, overlays, scriptSpecs)
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
	deps, err := utils.GetDependenciesFromScript(contentScript, config.Global.ParseModuleLoad || runParseModuleLoad)
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

// effectiveResourceSpec resolves the resource spec using the priority chain:
//
//	JobResources (scheduler env, actual allocation) > specs.Spec (directives) > defaults
func effectiveResourceSpec(specs *scheduler.ScriptSpecs) *scheduler.ResourceSpec {
	var jobRes *scheduler.ResourceSpec
	if sched := scheduler.ActiveScheduler(); sched != nil {
		jobRes = sched.GetJobResources()
	}
	return scheduler.ResolveResourceSpec(jobRes, specs)
}

// resourceEnvSettings derives NCPUS, MEM, MEM_MB, MEM_GB from scheduler specs
// and returns them as KEY=VALUE strings suitable for EnvSettings.
// In passthrough mode (Spec == nil) it falls back to live job resources.
// Applies priority chain: JobResources > ScriptSpec > Defaults.
func resourceEnvSettings(specs *scheduler.ScriptSpecs) []string {
	if specs == nil || specs.Spec == nil {
		return liveJobResourceEnvSettings()
	}
	return scheduler.ResourceEnvVars(effectiveResourceSpec(specs))
}

// liveJobResourceEnvSettings returns resource env vars from the active scheduler.
// Returns nil when not in a job or cannot retrieve resources.
func liveJobResourceEnvSettings() []string {
	sched := scheduler.ActiveScheduler()
	if sched == nil {
		return nil
	}
	jobRes := sched.GetJobResources()
	if jobRes == nil {
		return nil
	}
	// Only expose when the scheduler actually provided at least one resource value.
	if jobRes.Nodes == 0 && jobRes.TasksPerNode == 0 && jobRes.CpusPerTask == 0 && jobRes.MemPerNodeMB == 0 {
		return nil
	}
	return scheduler.ResourceEnvVars(jobRes)
}

// runLocally executes the script directly via apptainer with the given overlays.
func runLocally(ctx context.Context, contentScript string, overlays []string, specs *scheduler.ScriptSpecs) error {
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

	// Auto-bind the script's directory so it's accessible inside the container.
	// (Scheduler will copy the script to a temp location)
	scriptDir := filepath.Dir(contentScript)
	bindPaths := append([]string{scriptDir}, runBindPaths...)

	options := execpkg.Options{
		Overlays:     overlays,
		Command:      []string{"/bin/bash", "-c", executionScript},
		WritableImg:  runWritableImg,
		EnvSettings:  append(resourceEnvSettings(specs), runEnvSettings...),
		BindPaths:    bindPaths,
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

// validateAndConvertSpecs validates job specs against cluster limits.
// Prints a descriptive error if any resource exceeds the partition limit.
// For GPU errors, available alternatives are printed as suggestions.
func validateAndConvertSpecs(specs *scheduler.ScriptSpecs) error {
	if err := scheduler.ValidateAndConvertSpecs(specs); err != nil {
		utils.PrintError("Job specs validation failed: %v", err)

		// If it's a GPU validation error with suggestions, print them
		if gpuErr, ok := err.(*scheduler.GpuValidationError); ok && len(gpuErr.Suggestions) > 0 {
			utils.PrintMessage("Available GPU options:")
			for _, suggestion := range gpuErr.Suggestions {
				utils.PrintMessage("  - %s", suggestion)
			}
		}

		return err
	}

	return nil
}

// getNtasks returns the total number of MPI tasks using the priority chain:
// JobResources > ScriptSpec > Defaults. (Skip in passthrough mode)
func getNtasks(specs *scheduler.ScriptSpecs) int {
	if specs == nil || specs.Spec == nil {
		return 1
	}
	rs := effectiveResourceSpec(specs)
	nodes := rs.Nodes
	if nodes <= 0 {
		nodes = 1
	}
	tpn := rs.TasksPerNode
	if tpn <= 0 {
		tpn = 1
	}
	return nodes * tpn
}

// detectMpi probes the current environment for mpiexec, first directly and then
// via the module system.  Returns (available, moduleName, mpiexecPath): moduleName is
// non-empty when a module must be loaded; mpiexecPath is the full path to mpiexec used
// to derive the MPI installation root for bind-mounting into the container.
func detectMpi() (bool, string, string) {
	// Step 1: mpiexec already in PATH?
	if path, err := osexec.LookPath("mpiexec"); err == nil {
		return true, "", path
	}

	// Step 2: query the module system for an openmpi module.
	availCmd := osexec.Command("bash", "-lc", "module avail -t openmpi 2>&1")
	availOut, err := availCmd.Output()
	if err != nil || len(availOut) == 0 {
		return false, "", ""
	}

	var candidates []string
	for _, line := range strings.Split(string(availOut), "\n") {
		line = strings.TrimSpace(line)
		if !strings.HasPrefix(strings.ToLower(line), "openmpi") {
			continue
		}
		// Strip attribute markers like (default), (L), etc.
		if idx := strings.Index(line, "("); idx != -1 {
			line = strings.TrimSpace(line[:idx])
		}
		if line != "" {
			candidates = append(candidates, line)
		}
	}
	if len(candidates) == 0 {
		return false, "", ""
	}

	// Pick the highest version (lexicographic sort is sufficient for openmpi/X.Y.Z).
	sort.Strings(candidates)
	mod := candidates[len(candidates)-1]

	// Step 3: load the module and capture the mpiexec path.
	pathCmd := osexec.Command("bash", "-lc",
		fmt.Sprintf("module purge && module load %s && which mpiexec", mod))
	pathOut, err := pathCmd.Output()
	if err != nil || len(pathOut) == 0 {
		return false, "", ""
	}
	return true, mod, strings.TrimSpace(string(pathOut))
}

// buildMpiRunCommand returns the shell command to embed in the job script.
// When ntasks > 1 it attempts to detect mpiexec (directly or via a module) and
// wraps the condatainer run invocation with mpiexec.
// The user is responsible for installing the same MPI version inside the container.
func buildMpiRunCommand(contentScript string, specs *scheduler.ScriptSpecs) string {
	if getNtasks(specs) <= 1 {
		return fmt.Sprintf("condatainer run %s", contentScript)
	}
	ok, mod, mpiexecPath := detectMpi()
	if !ok {
		utils.PrintWarning("mpiexec not found; running %s without MPI wrapper", utils.StylePath(contentScript))
		return fmt.Sprintf("condatainer run %s", contentScript)
	}
	if mod != "" {
		utils.PrintNote("Detected MPI module: %s", utils.StyleName(mod))
		return fmt.Sprintf("module purge && module load %s && mpiexec condatainer run %s",
			mod, contentScript)
	}
	utils.PrintNote("Detected mpiexec: %s", utils.StylePath(mpiexecPath))
	return fmt.Sprintf("mpiexec condatainer run %s", contentScript)
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
	runCommand := buildMpiRunCommand(contentScript, specs)

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
