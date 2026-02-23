package build

import (
	"bufio"
	"context"
	"errors"

	"fmt"
	"os"
	"os/exec"
	"path/filepath"
	"strings"
	"time"

	"github.com/Justype/condatainer/internal/config"
	"github.com/Justype/condatainer/internal/overlay"
	"github.com/Justype/condatainer/internal/scheduler"
	"github.com/Justype/condatainer/internal/utils"
)

var ErrTmpOverlayExists = errors.New("temporary overlay already exists")
var ErrBuildCancelled = errors.New("build cancelled by user")

// BuildType represents the type of build target
type BuildType struct {
	IsConda bool
	IsDef   bool
	IsShell bool
	IsRef   bool
}

// Predefined build types
var (
	BuildTypeConda = BuildType{IsConda: true}
	BuildTypeDef   = BuildType{IsDef: true}
	BuildTypeShell = BuildType{IsShell: true}
	BuildTypeRef   = BuildType{IsShell: true, IsRef: true}
)

// String implements fmt.Stringer for BuildType.
// IsRef is checked before IsShell so ref overlays are reported as "ref".
func (bt BuildType) String() string {
	if bt.IsRef {
		return "ref"
	}
	if bt.IsConda {
		return "conda"
	}
	if bt.IsDef {
		return "def"
	}
	if bt.IsShell {
		return "shell"
	}
	return "unknown"
}

// isCancelledByUser checks if the error is due to user cancellation (Ctrl+C)
// Exit code 130 = 128 + SIGINT(2), checks for "signal: killed/interrupt" or context errors
func isCancelledByUser(err error) bool {
	if errors.Is(err, context.Canceled) || errors.Is(err, context.DeadlineExceeded) {
		return true
	}
	errMsg := err.Error()
	if strings.Contains(errMsg, "signal: killed") || strings.Contains(errMsg, "signal: interrupt") {
		return true
	}
	var exitErr *exec.ExitError
	if errors.As(err, &exitErr) {
		// 130 (SIGINT) or -1 (signal killed)
		return exitErr.ExitCode() == 130 || exitErr.ExitCode() == -1
	}
	return false
}

// ScriptSpecs mirrors the scheduler module's job spec metadata.
type ScriptSpecs = scheduler.ScriptSpecs

// buildDefaults holds resource defaults for build operations.
// Set from config at CLI startup via SetBuildDefaults.
var buildDefaults = scheduler.ResourceSpec{
	Nodes:        1,
	TasksPerNode: 1,
	CpusPerTask:  4, // conservative default for local builds
	MemPerNodeMB: 8192,
	Time:         2 * time.Hour,
}

// SetBuildDefaults sets the resource defaults used for build job submissions.
func SetBuildDefaults(d scheduler.ResourceSpec) { buildDefaults = d }

// GetBuildDefaults returns the current build resource defaults.
func GetBuildDefaults() scheduler.ResourceSpec { return buildDefaults }

// buildEffectiveResourceSpec resolves resources for build using the priority chain:
//
//	buildDefaults → scriptSpecs.Spec (when HasDirectives=true) → scheduler job resources
func buildEffectiveResourceSpec(specs *scheduler.ScriptSpecs) *scheduler.ResourceSpec {
	var jobRes *scheduler.ResourceSpec
	if sched := scheduler.ActiveScheduler(); sched != nil {
		jobRes = sched.GetJobResources()
	}
	return scheduler.ResolveResourceSpecFrom(buildDefaults, jobRes, specs)
}

// BuildObject represents a build target with its metadata and dependencies
type BuildObject interface {
	// Core properties
	NameVersion() string
	BuildSource() string
	Dependencies() []string
	IsInstalled() bool

	// Type information
	Type() BuildType

	// Path management
	TmpOverlayPath() string
	TargetOverlayPath() string
	CntDirPath() string

	// Scheduler metadata
	ScriptSpecs() *ScriptSpecs
	RequiresScheduler() bool

	// Build operations
	Build(ctx context.Context, buildDeps bool) error
	GetMissingDependencies() ([]string, error)
	CreateTmpOverlay(ctx context.Context, force bool) error
	Cleanup(failed bool) error

	// Update mode
	Update() bool

	// String representation
	String() string
}

// BaseBuildObject provides common functionality for all build types
type BaseBuildObject struct {
	nameVersion       string
	buildSource       string
	dependencies      []string
	tmpOverlayPath    string
	targetOverlayPath string
	cntDirPath        string
	submitJob         bool // Whether to submit to scheduler (from config at construction time)
	isRemote          bool // Whether build source was downloaded
	update            bool // If true, rebuild even if overlay already exists (atomic .new swap)
	scriptSpecs       *scheduler.ScriptSpecs

	// Interactive inputs for shell scripts
	interactiveInputs []string
}

// Common interface implementations for BaseBuildObject

func (b *BaseBuildObject) NameVersion() string       { return b.nameVersion }
func (b *BaseBuildObject) BuildSource() string       { return b.buildSource }
func (b *BaseBuildObject) Dependencies() []string    { return b.dependencies }
func (b *BaseBuildObject) TmpOverlayPath() string    { return b.tmpOverlayPath }
func (b *BaseBuildObject) TargetOverlayPath() string { return b.targetOverlayPath }
func (b *BaseBuildObject) CntDirPath() string        { return b.cntDirPath }
func (b *BaseBuildObject) ScriptSpecs() *ScriptSpecs { return b.scriptSpecs }
func (b *BaseBuildObject) Update() bool              { return b.update }
func (b *BaseBuildObject) RequiresScheduler() bool {
	return b.submitJob && scheduler.HasSchedulerSpecs(b.scriptSpecs)
}

// effectiveNcpus returns the effective total CPUs (CpusPerTask × TasksPerNode) for this build.
func (b *BaseBuildObject) effectiveNcpus() int {
	rs := buildEffectiveResourceSpec(b.scriptSpecs)
	cpus := rs.CpusPerTask
	if rs.TasksPerNode > 1 {
		cpus *= rs.TasksPerNode
	}
	return cpus
}

// effectiveMemMB returns the effective memory per node in MB for this build.
func (b *BaseBuildObject) effectiveMemMB() int64 {
	return buildEffectiveResourceSpec(b.scriptSpecs).MemPerNodeMB
}

func (b *BaseBuildObject) IsInstalled() bool {
	// Check if target overlay exists
	_, err := os.Stat(b.targetOverlayPath)
	return err == nil
}

func (b *BaseBuildObject) GetMissingDependencies() ([]string, error) {
	missing := []string{}

	// Build a set of installed overlays from all search paths
	installed := make(map[string]bool)

	for _, imagesDir := range config.GetImageSearchPaths() {
		if !utils.DirExists(imagesDir) {
			continue
		}

		entries, err := os.ReadDir(imagesDir)
		if err != nil {
			utils.PrintWarning("Failed to read directory %s: %v", imagesDir, err)
			continue
		}

		for _, entry := range entries {
			if entry.IsDir() {
				continue
			}
			if !utils.IsOverlay(entry.Name()) {
				continue
			}

			// Convert filename to name/version format
			nameVersion := strings.TrimSuffix(entry.Name(), filepath.Ext(entry.Name()))
			// Convert samtools--1.21.sqf to samtools/1.21
			normalized := strings.ReplaceAll(nameVersion, "--", "/")
			installed[normalized] = true
		}
	}

	// Check which dependencies are missing
	for _, dep := range b.dependencies {
		if !installed[dep] {
			missing = append(missing, dep)
		}
	}

	return missing, nil
}

func (b *BaseBuildObject) CreateTmpOverlay(ctx context.Context, force bool) error {
	// Check if tmp overlay already exists
	if _, err := os.Stat(b.tmpOverlayPath); err == nil {
		if !force {
			return fmt.Errorf("%w: %s", ErrTmpOverlayExists, b.tmpOverlayPath)
		}
		// Remove existing if force=true
		if err := os.Remove(b.tmpOverlayPath); err != nil {
			return fmt.Errorf("failed to remove existing tmp overlay: %w", err)
		}
	}

	// Ensure parent directory for tmp overlay exists (mkdir -p)
	parentDir := filepath.Dir(b.tmpOverlayPath)
	if parentDir != "" {
		if err := os.MkdirAll(parentDir, 0o775); err != nil {
			return fmt.Errorf("failed to create tmp overlay parent dir %s: %w", parentDir, err)
		}
	}

	// Create actual ext3 overlay using overlay package
	utils.PrintDebug("Creating temporary overlay at %s", utils.StylePath(b.tmpOverlayPath))

	// Use overlay package to create ext3 overlay
	// For build overlays, we use a temporary size from config (default 20GB)
	// Use "conda" profile for small files, sparse=true for faster creation
	// quiet=true suppresses detailed specs output since users don't need to see those for temp overlays
	if err := overlay.CreateForCurrentUser(ctx, b.tmpOverlayPath, config.Global.Build.TmpSizeMB, "default", true, config.Global.Build.OverlayType, true); err != nil {
		return fmt.Errorf("failed to create temporary overlay: %w", err)
	}

	return nil
}

func (b *BaseBuildObject) Cleanup(failed bool) error {
	// Remove remote build source if downloaded
	if b.isRemote && b.buildSource != "" {
		if err := os.Remove(b.buildSource); err != nil && !os.IsNotExist(err) {
			utils.PrintWarning("Failed to remove remote build source %s: %v", b.buildSource, err)
		} else {
			utils.PrintDebug("Removed remote build source %s", b.buildSource)
		}
	}

	// Remove tmp overlay
	if b.tmpOverlayPath != "" {
		if err := os.Remove(b.tmpOverlayPath); err != nil && !os.IsNotExist(err) {
			utils.PrintWarning("Failed to remove tmp overlay %s: %v", b.tmpOverlayPath, err)
		} else if utils.FileExists(b.tmpOverlayPath) {
			utils.PrintDebug("Removed tmp overlay %s", b.tmpOverlayPath)
		}
	}

	// Remove cnt directory
	if b.cntDirPath != "" {
		cntBaseDir := filepath.Dir(b.cntDirPath)
		if err := os.RemoveAll(cntBaseDir); err != nil && !os.IsNotExist(err) {
			utils.PrintWarning("Failed to remove cnt dir %s: %v", cntBaseDir, err)
		} else if utils.FileExists(cntBaseDir) {
			utils.PrintDebug("Removed cnt dir %s", cntBaseDir)
		}
	}

	// If failed, also remove target overlay
	if failed && b.targetOverlayPath != "" {
		if err := os.Remove(b.targetOverlayPath); err != nil && !os.IsNotExist(err) {
			utils.PrintWarning("Failed to remove target overlay %s: %v", b.targetOverlayPath, err)
		} else if utils.FileExists(b.targetOverlayPath) {
			utils.PrintDebug("Removed target overlay %s", b.targetOverlayPath)
		}
	}

	return nil
}

// parseScriptMetadata extracts dependencies, sbatch flags, and interactive prompts from shell scripts
func (b *BaseBuildObject) parseScriptMetadata() error {
	if b.buildSource == "" {
		return nil
	}

	// Always parse dependencies (needed for dependency graph)
	deps, err := utils.GetDependenciesFromScript(b.buildSource, config.Global.ParseModuleLoad)
	if err != nil {
		return fmt.Errorf("failed to parse dependencies: %w", err)
	}
	b.dependencies = deps

	// Parse interactive prompts from build script and collect user inputs if needed
	prompts, err := utils.GetInteractivePromptsFromScript(b.buildSource)
	if err != nil {
		return fmt.Errorf("failed to parse interactive prompts: %w", err)
	}
	b.interactiveInputs = []string{}
	if len(prompts) > 0 {
		// If --yes flag is set, automatically provide empty responses
		if utils.ShouldAnswerYes() {
			// Provide empty responses for all prompts with --yes
			for range prompts {
				b.interactiveInputs = append(b.interactiveInputs, "")
			}
		} else {
			// If interactive prompts exist, we must be in an interactive shell
			if !utils.IsInteractiveShell() {
				return fmt.Errorf("build script for %s requires interactive input, but no TTY is available", b.nameVersion)
			}

			reader := bufio.NewReader(os.Stdin)
			for _, prompt := range prompts {
				// Replace escaped "\\n" sequences with actual newlines and print the full message lines
				msg := strings.ReplaceAll(prompt, `\\n`, "\n")
				// Also handle single-backslash "\\n" sequences commonly used in scripts
				msg = strings.ReplaceAll(msg, "\\n", "\n")
				for _, line := range strings.Split(msg, "\n") {
					utils.PrintNote("%s", line)
				}
				// Prompt inline without adding a newline
				fmt.Print("Enter here: ")
				input, _ := reader.ReadString('\n')
				input = strings.TrimRight(input, "\r\n")
				if strings.ContainsAny(input, "\r\n") {
					utils.PrintWarning("Multiline input detected. Only the first line will be used.")
					if idx := strings.IndexAny(input, "\r\n"); idx != -1 {
						input = input[:idx]
					}
				}
				b.interactiveInputs = append(b.interactiveInputs, input)
			}
		}
	}

	// Always parse scheduler specs to get ncpus from script (for both local and scheduler builds)
	specs, err := scheduler.ReadScriptSpecsFromPath(b.buildSource)
	if err != nil {
		return err
	}

	// Always store scriptSpecs — effectiveNcpus()/effectiveMemMB() derive values from it.
	b.scriptSpecs = specs

	// Passthrough mode: scheduler directives found but resource parsing failed (unsupported flags).
	// Build cannot proceed without a normalized resource spec.
	if scheduler.IsPassthrough(specs) {
		return fmt.Errorf("build script %s contains unsupported scheduler directives (passthrough mode); remove or fix the unsupported directives", b.buildSource)
	}

	// Resolve using the priority chain: buildDefaults → script → job resources.
	specs.Spec = buildEffectiveResourceSpec(specs)
	return nil
}

// Helper functions

// getCntDirPath returns the container directory path for a name/version
// Format: <tmpDir>/build_<nameVersion>/cnt
func getCntDirPath(nameVersion, tmpDir string) string {
	buildDirName := "build_" + strings.ReplaceAll(nameVersion, "/", "_")
	return filepath.Join(tmpDir, buildDirName, "cnt")
}

// getTmpOverlayPath returns the temporary overlay path
// Format: <tmpDir>/<nameVersion>.sqf (with / replaced by --)
func getTmpOverlayPath(nameVersion, tmpDir string) string {
	filename := strings.ReplaceAll(nameVersion, "/", "--") + ".img"
	return filepath.Join(tmpDir, filename)
}

// NewBuildObject creates a BuildObject from a name/version string
// Format: "name/version" for conda/shell, "name" for def, "prefix/name/version" for ref
// All overlays are stored in imagesDir regardless of type
func NewBuildObject(nameVersion string, external bool, imagesDir, tmpDir string, update bool) (BuildObject, error) {
	normalized := utils.NormalizeNameVersion(nameVersion)
	slashCount := strings.Count(normalized, "/")

	// Determine build type based on slash count
	// 0 slashes: system (def is not defined here)
	// 1 slash: conda or shell (e.g., "numpy/1.24")
	// 2+ slashes: ref shell script (e.g., "genomes/hg38/full")
	isRef := slashCount > 1

	// Make tmpDir absolute
	if absDir, err := filepath.Abs(tmpDir); err == nil {
		tmpDir = absDir
	}

	// Create base object with resolved absolute path (including symlinks)
	targetOverlay := filepath.Join(imagesDir, strings.ReplaceAll(normalized, "/", "--")+".sqf")
	if abs, err := filepath.Abs(targetOverlay); err == nil {
		targetOverlay = abs
	}
	if real, err := filepath.EvalSymlinks(filepath.Dir(targetOverlay)); err == nil {
		targetOverlay = filepath.Join(real, filepath.Base(targetOverlay))
	}

	cntDirPath := getCntDirPath(normalized, tmpDir)
	tmpOverlayPath := getTmpOverlayPath(normalized, tmpDir)

	utils.PrintDebug("[BUILD OBJECT] Creating %s: nameVersion=%s, targetOverlay=%s, tmpOverlay=%s, cntDir=%s", nameVersion, normalized, targetOverlay, tmpOverlayPath, cntDirPath)

	base := &BaseBuildObject{
		nameVersion:       normalized,
		submitJob:         config.Global.SubmitJob,
		cntDirPath:        cntDirPath,
		tmpOverlayPath:    tmpOverlayPath,
		targetOverlayPath: targetOverlay,
		update:            update,
	}

	if external {
		// External builds don't need to resolve build source
		return createConcreteType(base, isRef, tmpDir)
	}

	// Check if already installed or temporary overlay exists.
	// Skip this optimisation in update mode so the correct concrete type is resolved
	// (the build will proceed regardless of install status).
	if !update && (base.IsInstalled() || utils.FileExists(base.tmpOverlayPath)) {
		return newScriptBuildObject(base, isRef)
	}

	// Resolve build source and determine concrete type
	return createConcreteType(base, isRef, tmpDir)
}

// NewCondaObjectWithSource creates a CondaBuildObject with custom buildSource
// This is used for the -n flag to create a single sqf with multiple packages or YAML
// buildSource can be:
//   - Path to YAML file (e.g., "/path/to/environment.yml")
//   - Comma-separated package list (e.g., "nvim,nodejs,samtools/1.16")
func NewCondaObjectWithSource(nameVersion, buildSource string, imagesDir, tmpDir string, update bool) (BuildObject, error) {
	normalized := utils.NormalizeNameVersion(nameVersion)

	// Make tmpDir absolute
	if absDir, err := filepath.Abs(tmpDir); err == nil {
		tmpDir = absDir
	}

	targetOverlay := filepath.Join(imagesDir, strings.ReplaceAll(normalized, "/", "--")+".sqf")
	if abs, err := filepath.Abs(targetOverlay); err == nil {
		targetOverlay = abs
	}

	cntDirPath := getCntDirPath(normalized, tmpDir)
	tmpOverlayPath := getTmpOverlayPath(normalized, tmpDir)

	utils.PrintDebug("[BUILD OBJECT] Creating conda %s: buildSource=%s, targetOverlay=%s, tmpOverlay=%s, cntDir=%s", nameVersion, buildSource, targetOverlay, tmpOverlayPath, cntDirPath)

	base := &BaseBuildObject{
		nameVersion:       normalized,
		buildSource:       buildSource,
		submitJob:         config.Global.SubmitJob,
		cntDirPath:        cntDirPath,
		tmpOverlayPath:    tmpOverlayPath,
		targetOverlayPath: targetOverlay,
		update:            update,
	}

	return newCondaBuildObject(base)
}

// FromExternalSource creates a BuildObject from an external build script or def file
// All overlays are stored in imagesDir regardless of type
func FromExternalSource(targetPrefix, source string, isApptainer bool, imagesDir, tmpDir string) (BuildObject, error) {
	nameVersion := filepath.Base(targetPrefix)
	nameVersion = utils.NormalizeNameVersion(nameVersion)

	// Make tmpDir absolute
	if absDir, err := filepath.Abs(tmpDir); err == nil {
		tmpDir = absDir
	}

	// Determine build type from source file extension
	isDef := isApptainer || strings.HasSuffix(source, ".def")
	isShell := strings.HasSuffix(source, ".sh") || strings.HasSuffix(source, ".bash")

	cntDirPath := getCntDirPath(nameVersion, tmpDir)
	tmpOverlayPath := getTmpOverlayPath(nameVersion, tmpDir)

	utils.PrintDebug("[BUILD OBJECT] Creating external %s: source=%s, targetPrefix=%s, tmpOverlay=%s, cntDir=%s", nameVersion, source, targetPrefix, tmpOverlayPath, cntDirPath)

	base := &BaseBuildObject{
		nameVersion:       nameVersion,
		buildSource:       source,
		submitJob:         config.Global.SubmitJob,
		cntDirPath:        cntDirPath,
		tmpOverlayPath:    tmpOverlayPath,
		targetOverlayPath: targetPrefix + ".sqf",
	}

	// Parse script metadata if it's a shell script
	if isShell {
		if err := base.parseScriptMetadata(); err != nil {
			return nil, err
		}
	}

	var obj BuildObject
	var err error

	if isDef {
		obj, err = newDefBuildObject(base)
	} else if isShell {
		obj, err = newScriptBuildObject(base, false) // External sources are not ref overlays
	} else {
		return nil, fmt.Errorf("unknown source type for %s", source)
	}

	if err != nil {
		return nil, err
	}

	utils.PrintDebug("[EXTERNAL] Created BuildObject from external path: %s", obj.String())
	return obj, nil
}

// createConcreteType creates the appropriate concrete BuildObject type
// It determines whether to create a conda, def, or script build based on:
// - isRef flag (from slash count > 1)
// - build source resolution: .def -> DefBuildObject, shell script -> ScriptBuildObject, not found -> conda
func createConcreteType(base *BaseBuildObject, isRef bool, tmpDir string) (BuildObject, error) {
	// Resolve build source - this determines the actual type based on file extension
	isConda, isContainer, err := resolveBuildSource(base, tmpDir)
	if err != nil {
		return nil, err
	}

	if isConda {
		return newCondaBuildObject(base)
	}

	if isContainer {
		// It's a .def file
		return newDefBuildObject(base)
	}

	// It's a shell script (no extension or .sh/.bash)
	// Parse script metadata
	if err := base.parseScriptMetadata(); err != nil {
		return nil, err
	}

	return newScriptBuildObject(base, isRef)
}

// resolveBuildSource finds the build script and determines the build type
// Returns (isConda, isContainer, error):
// - isConda=true: no build script found, use conda
// - isContainer=true: found .def file
// - both false: found shell script (no extension)
// First checks local and remote build scripts, falls back to conda if not found
func resolveBuildSource(base *BaseBuildObject, tmpDir string) (isConda bool, isContainer bool, err error) {
	// Check if already has a build source
	if base.buildSource != "" {
		// Determine type from extension
		isContainer = strings.HasSuffix(base.buildSource, ".def")
		return false, isContainer, nil
	}

	// Look for build script (local first, then remote)
	info, found := FindBuildScript(base.nameVersion)
	if !found {
		// No build script found, use conda
		utils.PrintDebug("No build script found for %s, using conda", utils.StyleName(base.nameVersion))
		return true, false, nil
	}

	// If remote, download to tmp directory
	if info.IsRemote {
		localPath, err := DownloadRemoteScript(info, tmpDir)
		if err != nil {
			return false, false, fmt.Errorf("failed to download remote build script: %w", err)
		}
		base.buildSource = localPath
		base.isRemote = true
		utils.PrintDebug("Downloaded remote build script to %s", utils.StylePath(localPath))
	} else {
		base.buildSource = info.Path
		utils.PrintDebug("Using local build script at %s", utils.StylePath(info.Path))
	}

	return false, info.IsContainer, nil
}
