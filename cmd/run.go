package cmd

import (
	"context"
	"errors"
	"fmt"
	"os"
	osexec "os/exec"
	"path/filepath"
	"regexp"
	"strconv"
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
	runStdout          string
	runStderr          string
	runAfterOK         string
	runCPU             int
	runMem             string
	runTime            string
	runGPU             string
	runDryRun          bool
	runArray           string
	runArrayLimit      int
	runName            string
)

// errRunAborted signals a handled stop (message already printed); caller returns nil.
var errRunAborted = errors.New("run aborted")

// jobIDRe matches a single scheduler job ID.
// Slurm/LSF: pure digits. PBS: digits followed by .hostname (e.g. 12345.school.edu). HTCondor: ClusterID.ProcessID (e.g. 12345.0).
var jobIDRe = regexp.MustCompile(`^\d+(\.\S+)?$`)

var runCmd = &cobra.Command{
	Use:   "run [flags] <script> [script_args...]",
	Short: "Run a script and auto-solve the dependencies by #DEP tags",
	Long: `Execute a script with automatic dependency resolution.

The script can contain special comment tags:
  #DEP: package/version   - Declares a dependency
  #DEP: /path/overlay.img - Declares an overlay dependency (.sqf and .img)
  #CNT args               - Additional arguments to pass to condatainer

Supported #CNT arguments (also available as CLI flags):
  -w, --writable         Make .img overlays writable
  -b, --base-image PATH  Use custom base image
  --env KEY=VALUE        Set environment variable (repeatable)
  --bind HOST:CONTAINER  Bind mount path (repeatable)
  -f, --fakeroot         Run with fakeroot privileges

Dependencies will be automatically loaded.`,
	Example: `  condatainer run script.sh                        # Run with dependency check
  condatainer run script.sh arg1 arg2              # Pass arguments to the script
  condatainer run -o log/tool1_s1.out run_tool1.sh sample1  # Override stdout
  condatainer run -g a100:2 gpu_job.sh             # Override GPU spec
  condatainer run -c 8 -m 16G -t 2h job.sh        # Override resources

  # Dependency chaining example
  JOB=$(condatainer run run_align.sh sample1)
  condatainer run --afterok "$JOB" run_quant.sh sample1`,
	Args:         cobra.MinimumNArgs(1),
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
	runCmd.Flags().StringVarP(&runStdout, "output", "o", "", "Override job stdout path (creates parent dir if needed)")
	runCmd.Flags().StringVarP(&runStderr, "error", "e", "", "Override job stderr path")
	runCmd.Flags().StringVar(&runAfterOK, "afterok", "", "Depend on job IDs (colon-separated, e.g. 123:456:789)")
	runCmd.Flags().StringArrayVar(&runEnvSettings, "env", nil, "Set environment variable KEY=VALUE (repeatable)")
	runCmd.Flags().StringArrayVar(&runBindPaths, "bind", nil, "Bind mount HOST:CONTAINER (repeatable)")
	runCmd.Flags().BoolVarP(&runFakeroot, "fakeroot", "f", false, "Run with fakeroot privileges")
	runCmd.Flags().IntVarP(&runCPU, "cpu", "c", 0, "Override CPUs per task (e.g. 4)")
	runCmd.Flags().StringVarP(&runMem, "mem", "m", "", "Override memory per node (e.g. 4G, 8192M)")
	runCmd.Flags().StringVarP(&runTime, "time", "t", "", "Override walltime (e.g. 4d12h, 2h30m, 01:30:00)")
	runCmd.Flags().StringVarP(&runGPU, "gpu", "g", "", "Override GPUs per node (e.g. 1, a100:2, a100)")
	runCmd.Flags().BoolVar(&runDryRun, "dry-run", false, "Preview what would happen without executing")
	runCmd.Flags().StringVarP(&runName, "name", "n", "", "Override job name")
	runCmd.Flags().StringVar(&runArray, "array", "", "Input file for array job (one entry per line)")
	runCmd.Flags().IntVar(&runArrayLimit, "array-limit", 0, "Max concurrent array subjobs (0 = unlimited)")
	runCmd.Flags().SetInterspersed(false) // Stop flag parsing after script name; remaining args are passed to the script
}

func runScript(cmd *cobra.Command, args []string) error {
	if cmd.Flags().Changed("afterok") {
		if runAfterOK == "" {
			ExitWithError("--afterok is empty; the upstream job may have failed to submit or run locally")
		}
		for id := range strings.SplitSeq(runAfterOK, ":") {
			if id = strings.TrimSpace(id); !jobIDRe.MatchString(id) {
				ExitWithError("--afterok %q is not a valid job ID", id)
			}
		}
	}

	// Build array spec early so we can validate before doing expensive work
	arraySpec, err := buildArraySpec()
	if err != nil {
		ExitWithError("%v", err)
	}

	scriptPath := args[0]
	scriptArgs := args[1:] // arguments for the script itself

	// Check if -i flag was used (alias for -a)
	if installFlag, _ := cmd.Flags().GetBool("install"); installFlag {
		runAutoInstall = true
	}

	// Ensure base image exists
	if err := ensureBaseImage(cmd.Context()); err != nil {
		return err
	}

	// Validate script exists
	if !utils.FileExists(scriptPath) {
		ExitWithError("Script file %s not found", utils.StylePath(scriptPath))
	}

	originScriptPath := scriptPath
	scriptPath, err = filepath.Abs(scriptPath)
	if err != nil {
		ExitWithError("Failed to get absolute path: %v", err)
	}

	// 1. Read specs → resolve the content script (HTCondor: .sub → .sh; others: identity)
	contentScript, scriptSpecs := resolveScriptAndSpecs(scriptPath)

	// CLI -o/-e/-n override script directives (highest priority)
	if scriptSpecs != nil {
		if runName != "" {
			scriptSpecs.Control.JobName = runName
		}
		if runStdout != "" {
			scriptSpecs.Control.Stdout = runStdout
		}
		if runStderr != "" {
			scriptSpecs.Control.Stderr = runStderr
		}
	}

	// CLI -c/-m/-t/-g resource overrides (highest priority, only when spec is parsed)
	if scriptSpecs != nil && scriptSpecs.Spec == nil && scriptSpecs.HasDirectives {
		// Script is in passthrough mode — resource directives could not be fully parsed.
		// Overrides cannot be applied; warn if the user specified any.
		if runCPU > 0 || runMem != "" || runTime != "" || runGPU != "" {
			utils.PrintWarning("resource overrides (-c/-m/-t/-g) have no effect in passthrough mode; edit the script directives directly")
		}
	} else if scriptSpecs != nil && scriptSpecs.Spec != nil {
		override := &scheduler.ResourceSpec{}
		if runCPU > 0 {
			override.CpusPerTask = runCPU
		}
		if runMem != "" {
			mb, err := utils.ParseMemoryMB(runMem)
			if err != nil {
				ExitWithError("--mem %q: %v", runMem, err)
			}
			override.MemPerNodeMB = mb
		}
		if runTime != "" {
			d, err := utils.ParseWalltime(runTime)
			if err != nil {
				ExitWithError("--time %q: %v", runTime, err)
			}
			override.Time = d
		}
		if runGPU != "" {
			gpu, err := parseGpuFlag(runGPU)
			if err != nil {
				ExitWithError("--gpu %q: %v", runGPU, err)
			}
			override.Gpu = gpu
		}
		scriptSpecs.Spec.Override(override)
	}

	// Conflict: --array flag + native scheduler array directive in the script
	if arraySpec != nil {
		if found := detectNativeArrayDirective(scriptSpecs); found != "" {
			ExitWithError("script contains native array directive %q; remove it and use --array flag instead", found)
		}
	}

	// 2. Embedded #CNT args + dependency check/install
	if err := processEmbeddedArgs(contentScript); err != nil {
		return err
	}
	// Dry run: print summary and exit without executing
	if runDryRun {
		printDryRunSummary(contentScript, originScriptPath, scriptSpecs, scriptArgs, arraySpec)
		return nil
	}

	// Blank lines in the array input file are only warned about in dry-run; error here
	if arraySpec != nil && len(arraySpec.BlankLines) > 0 {
		ExitWithError("--array: blank lines at %v in %s; please remove them", arraySpec.BlankLines, arraySpec.InputFile)
	}

	// Array jobs require scheduler submission; cannot run locally
	if arraySpec != nil && (scheduler.IsInsideJob() || config.IsInsideContainer()) {
		ExitWithError("--array requires scheduler submission; cannot run array jobs locally")
	}
	if runWritableImg && getNtasks(scriptSpecs) > 1 {
		ExitWithError("--writable cannot be used with multi-task jobs (ntasks=%d)", getNtasks(scriptSpecs))
	}
	overlays, depJobIDs, err := checkDepsAndAutoInstall(cmd.Context(), contentScript, originScriptPath)
	if err != nil {
		if errors.Is(err, errRunAborted) {
			os.Exit(ExitCodeError)
		}
		return err
	}

	// 3. In-job or in-container: always run locally (no nested submission)
	if scheduler.IsInsideJob() || config.IsInsideContainer() {
		return runLocally(cmd.Context(), contentScript, overlays, scriptSpecs, scriptArgs)
	}

	// 4. Submit if scheduler specs present and scheduler available
	if config.Global.SubmitJob && scheduler.HasSchedulerSpecs(scriptSpecs) {
		if err := validateAndConvertSpecs(scriptSpecs); err != nil {
			os.Exit(ExitCodeError)
		}
		sched, schedErr := scheduler.DetectScheduler()
		if schedErr != nil {
			utils.PrintNote("Script has scheduler specs but scheduler detection failed: %v. Running locally.", schedErr)
		} else if !sched.IsAvailable() {
			utils.PrintNote("Script has scheduler specs but no scheduler is available. Running locally.")
		} else {
			for _, id := range strings.Split(runAfterOK, ":") {
				if id = strings.TrimSpace(id); id != "" {
					depJobIDs = append(depJobIDs, id)
				}
			}
			return submitRunJob(sched, originScriptPath, contentScript, scriptSpecs, depJobIDs, scriptArgs, arraySpec)
		}
	} else if !scheduler.HasSchedulerSpecs(scriptSpecs) {
		utils.PrintNote("No scheduler specs found in script. Running locally.")
	}
	if runStdout != "" || runStderr != "" || runMem != "" || runTime != "" || runGPU != "" {
		utils.PrintNote("-o/-e/-m/-t/-g are only used for submitted jobs and will be ignored when running locally.")
	}
	return runLocally(cmd.Context(), contentScript, overlays, scriptSpecs, scriptArgs)
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
func checkDepsAndAutoInstall(ctx context.Context, contentScript, originScriptPath string) (overlays []string, depJobIDs []string, err error) {
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
			bo, err := build.NewBuildObject(ctx, pkg, false, imagesDir, config.GetWritableTmpDir(), false)
			if err != nil {
				utils.PrintError("Failed to create build object for %s: %v", utils.StyleName(pkg), err)
				return nil, nil, errRunAborted
			}
			buildObjects = append(buildObjects, bo)
		}

		graph, err := build.NewBuildGraph(ctx, buildObjects, imagesDir, config.GetWritableTmpDir(), config.Global.SubmitJob, false)
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
			depJobIDs = append(depJobIDs, id)
		}
		if len(depJobIDs) > 0 {
			utils.PrintMessage("Build job(s) submitted: %s", strings.Join(depJobIDs, ", "))
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

	return overlays, depJobIDs, nil
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
func runLocally(ctx context.Context, contentScript string, overlays []string, specs *scheduler.ScriptSpecs, scriptArgs []string) error {
	// Disable module commands and run the script
	executionScript := `module() { :; }
ml() { :; }
export -f module ml
`
	fileInfo, err := os.Stat(contentScript)
	if err != nil {
		ExitWithError("Failed to stat script: %v", err)
	}

	if fileInfo.Mode()&0111 != 0 {
		// Script is executable: $0 = contentScript, $@ = scriptArgs
		executionScript += `"$0" "$@"`
	} else {
		// Run with bash: $0 = contentScript, $@ = scriptArgs
		executionScript += `bash "$0" "$@"`
	}

	// Auto-bind the script's directory so it's accessible inside the container.
	// (Scheduler will copy the script to a temp location)
	scriptDir := filepath.Dir(contentScript)
	bindPaths := append([]string{scriptDir}, runBindPaths...)

	options := execpkg.Options{
		Overlays:     overlays,
		Command:      append([]string{"/bin/bash", "-c", executionScript, contentScript}, scriptArgs...),
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

// printDryRunSummary prints what condatainer run would do without executing.
func printDryRunSummary(contentScript, originScript string, specs *scheduler.ScriptSpecs, scriptArgs []string, arraySpec *scheduler.ArraySpec) {
	fmt.Printf("%s %s\n", utils.StyleTitle("Dry run:"), originScript)

	// Dependencies
	baseImg := runBaseImage
	if baseImg == "" {
		baseImg = config.GetBaseImage()
	}

	deps, err := utils.GetDependenciesFromScript(contentScript, config.Global.ParseModuleLoad || runParseModuleLoad)
	if err != nil {
		fmt.Printf("%s Could not parse dependencies: %v\n", utils.StyleError("[ERR]"), err)
	} else {
		installedOverlays, _ := getInstalledOverlaysMap()
		check := utils.StyleSuccess("✓")
		cross := utils.StyleError("✗")

		// First pass: collect results and compute max dep name width for alignment
		type depEntry struct {
			dep, path string
			ok        bool
		}
		entries := make([]depEntry, 0, len(deps))
		for _, dep := range deps {
			var entry depEntry
			entry.dep = dep
			if utils.IsOverlay(dep) {
				entry.ok = utils.FileExists(dep)
				entry.path = dep
			} else {
				normalized := utils.NormalizeNameVersion(dep)
				entry.path, entry.ok = installedOverlays[normalized]
			}
			entries = append(entries, entry)
		}

		found, missing := 0, 0
		if utils.FileExists(baseImg) {
			found++
		} else {
			missing++
		}
		for _, e := range entries {
			if e.ok {
				found++
			} else {
				missing++
			}
		}
		fmt.Printf("%s (%d found, %d missing):\n", utils.StyleTitle("Dependencies:"), found, missing)

		// Base image
		if utils.FileExists(baseImg) {
			fmt.Printf("  Base:      %s %s\n", check, utils.StylePath(baseImg))
		} else {
			fmt.Printf("  Base:      %s %s %s\n", cross, utils.StylePath(baseImg), utils.StyleWarning("(not found)"))
		}

		// External overlays (direct path): show under "External:" header; installed: group by folder
		var externals []string
		dirOrder := []string{}
		dirDeps := map[string][]string{}
		for _, e := range entries {
			if !e.ok {
				continue
			}
			if utils.IsOverlay(e.dep) {
				externals = append(externals, e.dep)
			} else {
				dir := filepath.Dir(e.path)
				if _, exists := dirDeps[dir]; !exists {
					dirOrder = append(dirOrder, dir)
				}
				dirDeps[dir] = append(dirDeps[dir], e.dep)
			}
		}
		if len(externals) > 0 {
			fmt.Printf("  %s\n", "External:")
			for _, dep := range externals {
				fmt.Printf("    %s %s\n", check, utils.StylePath(dep))
			}
		}
		for _, dir := range dirOrder {
			fmt.Printf("  %s\n", utils.StylePath(dir+"/"))
			for _, dep := range dirDeps[dir] {
				fmt.Printf("    %s %s\n", check, dep)
			}
		}
		// Missing deps listed after found
		for _, e := range entries {
			if !e.ok {
				suffix := utils.StyleWarning("(not installed)")
				if utils.IsOverlay(e.dep) {
					suffix = utils.StyleWarning("(not found)")
				}
				fmt.Printf("  %s %s  %s\n", cross, e.dep, suffix)
			}
		}
	}

	// Scheduler specs
	if specs != nil && specs.Spec != nil {
		rs := effectiveResourceSpec(specs)
		fmt.Printf("%s\n", utils.StyleTitle("Resource specs:"))
		if rs.CpusPerTask > 0 {
			fmt.Printf("  CPUs:      %d\n", rs.CpusPerTask)
		}
		if rs.MemPerNodeMB > 0 {
			if rs.MemPerNodeMB >= 1024 && rs.MemPerNodeMB%1024 == 0 {
				fmt.Printf("  Memory:    %d GB\n", rs.MemPerNodeMB/1024)
			} else {
				fmt.Printf("  Memory:    %d MB\n", rs.MemPerNodeMB)
			}
		}
		if rs.Time > 0 {
			total := int64(rs.Time.Seconds())
			fmt.Printf("  Time:      %02d:%02d:%02d\n", total/3600, (total%3600)/60, total%60)
		}
		if rs.Gpu != nil {
			if rs.Gpu.Type != "" {
				fmt.Printf("  GPU:       %s x%d\n", rs.Gpu.Type, rs.Gpu.Count)
			} else {
				fmt.Printf("  GPU:       %d\n", rs.Gpu.Count)
			}
		}
		ntasks := getNtasks(specs)
		if ntasks > 1 {
			if rs.Nodes > 0 {
				fmt.Printf("  Nodes:     %d\n", rs.Nodes)
			}
			if rs.TasksPerNode > 0 {
				fmt.Printf("  Tasks/node: %d\n", rs.TasksPerNode)
			}
			fmt.Printf("  MPI tasks: %d\n", ntasks)
			mpiexecPath, ok := detectMpi()
			if !ok {
				fmt.Printf("  mpiexec:   %s\n", utils.StyleError("not found"))
			} else {
				fmt.Printf("  mpiexec:   %s\n", utils.StylePath(mpiexecPath))
			}
		}

		fmt.Printf("%s\n", utils.StyleTitle("Job control:"))
		scriptBase := strings.TrimSuffix(filepath.Base(originScript), filepath.Ext(originScript))
		if specs.Control.JobName != "" {
			fmt.Printf("  Name:      %s\n", specs.Control.JobName)
		} else {
			fmt.Printf("  Name:      %s (default)\n", scriptBase)
		}
		if specs.Control.Partition != "" {
			fmt.Printf("  Partition: %s\n", specs.Control.Partition)
		}
		logsDir := config.Global.LogsDir
		if logsDir == "" {
			logsDir = filepath.Join(os.Getenv("HOME"), "logs")
		}
		defaultOut := filepath.Join(logsDir, scriptBase+"_<timestamp>.out")
		if specs.Control.Stdout != "" {
			fmt.Printf("  Stdout:    %s\n", specs.Control.Stdout)
		} else {
			fmt.Printf("  Stdout:    %s (default)\n", defaultOut)
		}
		if specs.Control.Stderr != "" {
			fmt.Printf("  Stderr:    %s\n", specs.Control.Stderr)
		} else if specs.Control.Stdout != "" {
			fmt.Printf("  Stderr:    %s (default)\n", specs.Control.Stdout)
		} else {
			fmt.Printf("  Stderr:    %s (default)\n", defaultOut)
		}

		if arraySpec != nil {
			arrayLine := fmt.Sprintf("  Array:     %d subjobs from %s", arraySpec.Count, utils.StylePath(arraySpec.InputFile))
			if arraySpec.Limit > 0 {
				arrayLine += fmt.Sprintf(" (max %d concurrent)", arraySpec.Limit)
			}
			fmt.Println(arrayLine)
			if len(arraySpec.BlankLines) > 0 {
				fmt.Printf("  %s blank lines at %v — remove before submitting\n",
					utils.StyleError("[ERR]"), arraySpec.BlankLines)
			}
		}
		if runAfterOK != "" {
			fmt.Printf("  AfterOK:   %s\n", runAfterOK)
		}

		// Email Notifications
		var events []string
		if specs.Control.EmailOnBegin {
			events = append(events, "BEGIN")
		}
		if specs.Control.EmailOnEnd {
			events = append(events, "END")
		}
		if specs.Control.EmailOnFail {
			events = append(events, "FAIL")
		}
		if len(events) > 0 {
			mailStr := strings.Join(events, ",")
			if specs.Control.MailUser != "" {
				mailStr += " to " + specs.Control.MailUser
			}
			fmt.Printf("  Mail:      %s\n", mailStr)
		}
	}

	// Container args
	hasContainerArgs := len(runEnvSettings) > 0 || len(runBindPaths) > 0 || runFakeroot || runWritableImg
	if hasContainerArgs {
		fmt.Printf("%s\n", utils.StyleTitle("Container args:"))
		for _, env := range runEnvSettings {
			fmt.Printf("  Env:       %s\n", env)
		}
		for _, bind := range runBindPaths {
			fmt.Printf("  Bind:      %s\n", bind)
		}
		if runFakeroot {
			fmt.Printf("  Fakeroot:  yes\n")
		}
		if runWritableImg {
			fmt.Printf("  Writable:  yes\n")
		}
	}

	// Script args (array args prepended, then fixed CLI args)
	hasScriptArgs := arraySpec != nil || len(scriptArgs) > 0
	if hasScriptArgs {
		fmt.Printf("%s\n", utils.StyleTitle("Script args:"))
		if arraySpec != nil {
			fileBase := filepath.Base(arraySpec.InputFile)
			maxLen := 0
			for _, token := range arraySpec.SampleArgs {
				if len(token) > maxLen {
					maxLen = len(token)
				}
			}
			for i, token := range arraySpec.SampleArgs {
				fmt.Printf("  %-*s  (arg%d from %s)\n", maxLen, token, i+1, fileBase)
			}
		}
		for i := 0; i < len(scriptArgs); i++ {
			arg := scriptArgs[i]
			fmt.Printf("  %s\n", arg)
		}
	}

	// Action
	if scheduler.IsInsideJob() {
		fmt.Printf("%s Would run locally (in job)\n", utils.StyleTitle("Action:"))
	} else if config.IsInsideContainer() {
		fmt.Printf("%s Would run locally (in container)\n", utils.StyleTitle("Action:"))
	} else if config.Global.SubmitJob && scheduler.HasSchedulerSpecs(specs) {
		sched, err := scheduler.DetectScheduler()
		if err != nil || !sched.IsAvailable() {
			fmt.Printf("%s Would run locally (scheduler not available)\n", utils.StyleTitle("Action:"))
		} else {
			fmt.Printf("%s Would submit to %s\n", utils.StyleTitle("Action:"), sched.GetInfo().Type)
		}
	} else {
		fmt.Printf("%s Would run locally\n", utils.StyleTitle("Action:"))
	}
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
				// Handle --flag=VALUE formats
				if v, ok := strings.CutPrefix(arg, "--env="); ok {
					runEnvSettings = append(runEnvSettings, v)
				} else if v, ok := strings.CutPrefix(arg, "--bind="); ok {
					runBindPaths = append(runBindPaths, v)
				} else if v, ok := strings.CutPrefix(arg, "-b="); ok {
					runBaseImage = v
				} else if v, ok := strings.CutPrefix(arg, "--base-image="); ok {
					runBaseImage = v
				}
			}
		}
	}
}

// buildArraySpec reads and validates the array input file, returning an ArraySpec.
// All non-empty lines must shell-split into the same number of tokens.
// Returns nil when --array was not provided.
func buildArraySpec() (*scheduler.ArraySpec, error) {
	if runArray == "" {
		return nil, nil
	}
	absFile, err := filepath.Abs(runArray)
	if err != nil {
		return nil, fmt.Errorf("--array: cannot resolve path %q: %w", runArray, err)
	}
	data, err := os.ReadFile(absFile)
	if err != nil {
		return nil, fmt.Errorf("--array: cannot read input file %q: %w", absFile, err)
	}
	var sampleArgs []string
	expectedArgCount := -1
	count := 0
	var blankLines []int
	for lineNum, line := range strings.Split(strings.TrimRight(string(data), "\r\n"), "\n") {
		line = strings.TrimRight(line, "\r") // strip \r from Windows line endings
		if strings.TrimSpace(line) == "" {
			blankLines = append(blankLines, lineNum+1)
			continue
		}
		count++
		tokens := shellSplitLine(line)
		if expectedArgCount == -1 {
			expectedArgCount = len(tokens)
			sampleArgs = tokens
		} else if len(tokens) != expectedArgCount {
			return nil, fmt.Errorf("--array: line %d has %d arg(s) but line 1 has %d; all lines must have the same number of arguments",
				lineNum+1, len(tokens), expectedArgCount)
		}
	}
	if count == 0 {
		return nil, fmt.Errorf("--array: input file %q is empty", absFile)
	}
	return &scheduler.ArraySpec{
		InputFile:  absFile,
		Count:      count,
		Limit:      runArrayLimit,
		ArgCount:   expectedArgCount,
		SampleArgs: sampleArgs,
		BlankLines: blankLines,
	}, nil
}

// parseArgsInScript extracts arguments from #CNT comments in the script
// parseGpuFlag parses --gpu values: "N", "TYPE:N", or "TYPE" (count defaults to 1).
func parseGpuFlag(s string) (*scheduler.GpuSpec, error) {
	if before, after, found := strings.Cut(s, ":"); found {
		// TYPE:N
		count, err := strconv.Atoi(after)
		if err != nil || count <= 0 {
			return nil, fmt.Errorf("invalid GPU count %q (expected TYPE:N, e.g. a100:2)", after)
		}
		return &scheduler.GpuSpec{Type: strings.ToLower(before), Count: count}, nil
	}
	// bare N or bare TYPE
	if count, err := strconv.Atoi(s); err == nil {
		if count <= 0 {
			return nil, fmt.Errorf("GPU count must be > 0")
		}
		return &scheduler.GpuSpec{Count: count}, nil
	}
	// TYPE only — default count 1
	if s == "" {
		return nil, fmt.Errorf("empty GPU spec")
	}
	return &scheduler.GpuSpec{Type: strings.ToLower(s), Count: 1}, nil
}

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

// detectMpi checks whether mpiexec is available in the current PATH.
// Returns the full path and true when found; empty string and false otherwise.
func detectMpi() (string, bool) {
	path, err := osexec.LookPath("mpiexec")
	if err != nil {
		return "", false
	}
	return path, true
}

// buildMpiRunCommand returns the shell command to embed in the job script.
// When ntasks > 1 it attempts to detect mpiexec (directly or via a module) and
// wraps the condatainer run invocation with mpiexec.
// Returns an error when ntasks > 1 but mpiexec cannot be found.
// The user is responsible for installing the same MPI version inside the container.
func buildMpiRunCommand(contentScript string, scriptArgs []string, specs *scheduler.ScriptSpecs, prependArrayArgs bool) (string, error) {
	runCmd := fmt.Sprintf("condatainer run %s", contentScript)
	if prependArrayArgs {
		runCmd += " $ARRAY_ARGS"
	}
	if len(scriptArgs) > 0 {
		quoted := make([]string, len(scriptArgs))
		for i, a := range scriptArgs {
			quoted[i] = shellQuote(a)
		}
		runCmd += " " + strings.Join(quoted, " ")
	}
	if getNtasks(specs) <= 1 {
		return runCmd, nil
	}
	mpiexecPath, ok := detectMpi()
	if !ok {
		return "", fmt.Errorf("mpiexec not found; load the appropriate MPI module before submitting (ntasks=%d)", getNtasks(specs))
	}
	utils.PrintNote("Detected mpiexec: %s", utils.StylePath(mpiexecPath))
	return fmt.Sprintf("mpiexec %s", runCmd), nil
}

// shellQuote returns a single-quoted shell-safe version of s.
// Embedded single quotes are escaped as '\'.
func shellQuote(s string) string {
	return "'" + strings.ReplaceAll(s, "'", `'\''`) + "'"
}

// shellSplitLine splits a shell-like line into tokens.
// Single- and double-quoted strings are kept as one token (quotes stripped).
// Unquoted whitespace is the delimiter.
func shellSplitLine(line string) []string {
	var tokens []string
	var cur strings.Builder
	inSingle, inDouble := false, false
	for _, ch := range line {
		switch {
		case ch == '\'' && !inDouble:
			inSingle = !inSingle
		case ch == '"' && !inSingle:
			inDouble = !inDouble
		case (ch == ' ' || ch == '\t') && !inSingle && !inDouble:
			if cur.Len() > 0 {
				tokens = append(tokens, cur.String())
				cur.Reset()
			}
		default:
			cur.WriteRune(ch)
		}
	}
	if cur.Len() > 0 {
		tokens = append(tokens, cur.String())
	}
	return tokens
}

// detectNativeArrayDirective returns a non-empty string describing a native scheduler
// array directive found in specs, or "" if none. Used to detect conflicts when the
// --array CLI flag is also in use.
func detectNativeArrayDirective(specs *scheduler.ScriptSpecs) string {
	if specs == nil {
		return ""
	}
	for _, f := range specs.RemainingFlags {
		// SLURM: --array=N-M  --array N-M  -a N-M  -a=N-M
		if strings.HasPrefix(f, "--array") || strings.HasPrefix(f, "-a=") ||
			f == "-a" || strings.HasPrefix(f, "-a ") {
			return f
		}
		// PBS: -J N-M (PBS job-name flag is -N, not -J)
		if strings.HasPrefix(f, "-J ") || strings.HasPrefix(f, "-J\t") ||
			strings.HasPrefix(f, "-J=") {
			return f
		}
	}
	// LSF: -J name[N-M] — array range embedded in job name with brackets
	if strings.Contains(specs.Control.JobName, "[") {
		return fmt.Sprintf("-J %s", specs.Control.JobName)
	}
	return ""
}

// submitRunJob creates and submits a scheduler job to run the script.
// contentScript is the bash script containing #DEP/#CNT directives — for HTCondor this is the
// executable referenced in the .sub file, for other schedulers it equals scriptPath.
func submitRunJob(sched scheduler.Scheduler, originScriptPath, contentScript string, specs *scheduler.ScriptSpecs, depIDs []string, scriptArgs []string, arraySpec *scheduler.ArraySpec) error {
	info := sched.GetInfo()

	// Capture separate-output intent before CreateScriptWithSpec overrides Stdout/Stderr to /dev/null
	arraySeparateOutput := arraySpec != nil && specs.Control.Stderr != "" && specs.Control.Stderr != specs.Control.Stdout

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
	runCommand, err := buildMpiRunCommand(contentScript, scriptArgs, specs, arraySpec != nil)
	if err != nil {
		utils.PrintError("%v", err)
		return err
	}

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
		Array:     arraySpec,
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

	fmt.Fprintln(os.Stdout, jobID)
	if len(depIDs) > 0 {
		utils.PrintSuccess("Submitted %s job %s for %s (after: %s)", info.Type, utils.StyleNumber(jobID), utils.StylePath(originScriptPath), strings.Join(depIDs, ", "))
	} else {
		utils.PrintSuccess("Submitted %s job %s for %s", info.Type, utils.StyleNumber(jobID), utils.StylePath(originScriptPath))
	}

	if arraySpec != nil {
		// Array jobs redirect output per-subjob inside the script; show glob pattern
		safeName := strings.ReplaceAll(fileBaseName, "/", "--")
		if arraySeparateOutput {
			utils.PrintMessage("Per-subjob stdout/err => %s", utils.StylePath(filepath.Join(logsDir, safeName+"_*.out/*.err")))
		} else {
			utils.PrintMessage("Per-subjob stdout&err => %s", utils.StylePath(filepath.Join(logsDir, safeName+"_*.log")))
		}
	} else {
		stdoutPath := jobSpec.Specs.Control.Stdout
		stderrPath := jobSpec.Specs.Control.Stderr
		if stdoutPath != "" {
			if stderrPath == "" || stdoutPath == stderrPath {
				utils.PrintMessage("Stdout & Stderr => %s", utils.StylePath(stdoutPath))
			} else {
				utils.PrintMessage("Stdout => %s", utils.StylePath(stdoutPath))
				utils.PrintMessage("Stderr => %s", utils.StylePath(stderrPath))
			}
		}
	}

	return nil
}
