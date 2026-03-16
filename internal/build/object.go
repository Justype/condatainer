package build

import (
	"context"
	"errors"
	"fmt"
	"os"
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

// ScriptSpecs mirrors the scheduler module's job spec metadata.
type ScriptSpecs = scheduler.ScriptSpecs

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
	LockPath() string

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
	tmpDir            string // base tmp directory (dynamic; scheduler/TMPDIR-based for conda/app)
	dependencies      []string
	tmpOverlayPath    string
	targetOverlayPath string
	cntDirPath        string
	submitJob         bool   // Whether to submit to scheduler (from config at construction time)
	isRemote          bool   // Whether build source was downloaded
	prebuiltLink      string // per-source prebuilt base URL (from metadata/prebuilt_link); empty = no prebuilt
	update            bool   // If true, rebuild even if overlay already exists (atomic .new swap)
	scriptSpecs       *scheduler.ScriptSpecs

	// Interactive inputs for shell scripts
	interactiveInputs []string
}

// Common interface implementations for BaseBuildObject

func (b *BaseBuildObject) NameVersion() string       { return b.nameVersion }
func (b *BaseBuildObject) BuildSource() string       { return b.buildSource }
func (b *BaseBuildObject) Dependencies() []string    { return b.dependencies }
func (b *BaseBuildObject) TmpDir() string            { return b.tmpDir }
func (b *BaseBuildObject) TmpOverlayPath() string    { return b.tmpOverlayPath }
func (b *BaseBuildObject) TargetOverlayPath() string { return b.targetOverlayPath }
func (b *BaseBuildObject) CntDirPath() string        { return b.cntDirPath }
func (b *BaseBuildObject) ScriptSpecs() *ScriptSpecs { return b.scriptSpecs }
func (b *BaseBuildObject) Update() bool              { return b.update }
func (b *BaseBuildObject) RequiresScheduler() bool {
	return b.submitJob && scheduler.HasSchedulerSpecs(b.scriptSpecs)
}

// BuildLockInfo holds metadata stored inside a build lock file.
type BuildLockInfo struct {
	Type      string `json:"type"`       // "local", "slurm", "pbs", "lsf", or "htcondor"
	JobID     string `json:"job_id"`     // scheduler job ID (empty until submit returns)
	Node      string `json:"node"`       // short hostname where lock was created
	PID       int    `json:"pid"`        // OS PID for local builds; 0 for scheduler
	CreatedAt string `json:"created_at"` // RFC3339 timestamp
}

// buildLockPath returns the lock file path for this build (stable, in imagesDir).
func (b *BaseBuildObject) buildLockPath() string {
	return b.targetOverlayPath + ".lock"
}

// LockPath satisfies the BuildObject interface — returns the lock file path.
func (b *BaseBuildObject) LockPath() string { return b.buildLockPath() }

// writeBuildLock atomically creates the lock file and writes JSON metadata.
// Returns an os.ErrExist-wrapped error if the lock already exists.
func (b *BaseBuildObject) writeBuildLock(info BuildLockInfo) error {
	return acquireBuildLockFile(b.buildLockPath(), info)
}

// readBuildLock reads and parses the lock file JSON.
// An empty or corrupt file (old empty-lock format) returns a zero-value struct.
func (b *BaseBuildObject) readBuildLock() (BuildLockInfo, error) {
	return readBuildLockFile(b.buildLockPath())
}

// updateBuildLock overwrites the lock file contents.
func (b *BaseBuildObject) updateBuildLock(info BuildLockInfo) error {
	return overwriteBuildLockFile(b.buildLockPath(), info)
}

// createBuildLock creates the lock file for a local build.
// If a scheduler job lock already exists for this job, it adopts that lock
// (updating it with the runtime node and PID) instead of failing.
func (b *BaseBuildObject) createBuildLock() error {
	info := BuildLockInfo{
		Type:      "local",
		Node:      shortHostname(),
		PID:       os.Getpid(),
		CreatedAt: time.Now().Format(time.RFC3339),
	}
	if err := b.writeBuildLock(info); err != nil {
		if !os.IsExist(err) {
			return fmt.Errorf("failed to create build lock: %w", err)
		}
		// Lock exists — check if it belongs to our scheduler job.
		existing, readErr := b.readBuildLock()
		if readErr == nil && existing.Type != "local" && existing.Type != "" {
			if myJobID := scheduler.CurrentJobID(); myJobID != "" && existing.JobID == myJobID {
				// Adopt the lock: update with runtime node + PID.
				existing.Node = shortHostname()
				existing.PID = os.Getpid()
				return b.updateBuildLock(existing)
			}
		}
		return fmt.Errorf("build already in progress: lock file exists at %s", b.buildLockPath())
	}
	return nil
}

// removeBuildLock removes the lock file on build completion or failure.
func (b *BaseBuildObject) removeBuildLock() {
	os.Remove(b.buildLockPath()) //nolint:errcheck
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

func (b *BaseBuildObject) IsInstalled() bool {
	_, err := os.Stat(b.targetOverlayPath)
	return err == nil
}

func (b *BaseBuildObject) GetMissingDependencies() ([]string, error) {
	installed := getInstalledOverlays()
	var missing []string
	for _, dep := range b.dependencies {
		if !installed[dep] {
			missing = append(missing, dep)
		}
	}
	return missing, nil
}

func (b *BaseBuildObject) CreateTmpOverlay(ctx context.Context, force bool) error {
	// Check both ext3-mode artifact (.img) and dir-mode artifact (buildDir) for cross-mode stale detection
	buildDir := filepath.Dir(b.cntDirPath)
	stale := utils.FileExists(b.tmpOverlayPath) || (b.cntDirPath != "" && utils.DirExists(buildDir))
	if stale {
		if !force {
			return fmt.Errorf("%w: %s", ErrTmpOverlayExists, b.tmpOverlayPath)
		}
		if err := os.Remove(b.tmpOverlayPath); err != nil && !os.IsNotExist(err) {
			return fmt.Errorf("failed to remove existing tmp overlay: %w", err)
		}
		os.RemoveAll(buildDir) //nolint:errcheck
	}

	// Ensure parent directory for tmp overlay exists (mkdir -p)
	parentDir := filepath.Dir(b.tmpOverlayPath)
	if parentDir != "" {
		if err := os.MkdirAll(parentDir, utils.PermDir); err != nil {
			return fmt.Errorf("failed to create tmp overlay parent dir %s: %w", parentDir, err)
		}
	}

	utils.PrintDebug("Creating temporary overlay at %s", utils.StylePath(b.tmpOverlayPath))

	// Use overlay package to create ext3 overlay.
	// For build overlays, we use a temporary size from config (default 20GB).
	// Use "default" profile, sparse=true for faster creation.
	// quiet=true suppresses detailed specs output since users don't need to see those for temp overlays.
	if err := overlay.CreateWithOptions(ctx, &overlay.CreateOptions{
		Path: b.tmpOverlayPath, SizeMB: config.Global.Build.TmpSizeMB,
		UID: os.Getuid(), GID: os.Getgid(),
		Profile: overlay.ProfileDefault, Sparse: true, FilesystemType: "ext3", Quiet: true,
	}); err != nil {
		return fmt.Errorf("failed to create temporary overlay: %w", err)
	}

	return nil
}

// CreateBuildDirs creates host directories for dir-mode builds (use_tmp_overlay=false).
// Layout: <buildDir>/cnt/ (bound as /cnt) and <buildDir>/tmp/ (bound as /ext3/tmp).
// Checks both dir-mode (buildDir) and ext3-mode (.img) artifacts for cross-mode stale detection.
func (b *BaseBuildObject) CreateBuildDirs(ctx context.Context, force bool) error {
	buildDir := filepath.Dir(b.cntDirPath)
	// Check both dir-mode artifact (buildDir) and ext3-mode artifact (.img)
	stale := utils.DirExists(buildDir) || utils.FileExists(b.tmpOverlayPath)
	if stale {
		if !force {
			return fmt.Errorf("%w: %s", ErrTmpOverlayExists, buildDir)
		}
		os.RemoveAll(buildDir)      //nolint:errcheck — clean dir-mode artifact
		os.Remove(b.tmpOverlayPath) //nolint:errcheck — clean ext3-mode artifact (no-op if "")
	}

	if err := os.MkdirAll(b.cntDirPath, utils.PermDir); err != nil {
		return fmt.Errorf("failed to create build cnt dir %s: %w", b.cntDirPath, err)
	}
	if err := os.MkdirAll(filepath.Join(buildDir, "tmp"), utils.PermDir); err != nil {
		return fmt.Errorf("failed to create build tmp dir: %w", err)
	}
	utils.PrintMessage("Build dir: %s", utils.StylePath(buildDir))
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
		utils.PrintMessage("Cleaning up build directory %s...", cntBaseDir)
		if err := os.RemoveAll(cntBaseDir); err != nil && !os.IsNotExist(err) {
			utils.PrintWarning("Failed to remove cnt dir %s: %v", cntBaseDir, err)
		} else if utils.FileExists(cntBaseDir) {
			utils.PrintDebug("Removed cnt dir %s", cntBaseDir)
		}

		// Remove tmpDir itself if it is now empty (no other builds using it)
		if b.tmpDir != "" && b.tmpDir != cntBaseDir {
			utils.RemoveDirIfEmpty(b.tmpDir)
		}
	}

	// If failed, remove target overlay (or its .new counterpart in update mode).
	if failed && b.targetOverlayPath != "" {
		if b.update {
			// In update mode the build writes to targetOverlayPath+".new"; preserve the
			// existing target so a failed update doesn't destroy the installed overlay.
			newPath := b.targetOverlayPath + ".new"
			if err := os.Remove(newPath); err != nil && !os.IsNotExist(err) {
				utils.PrintWarning("Failed to remove partial new overlay %s: %v", newPath, err)
			}
		} else {
			if err := os.Remove(b.targetOverlayPath); err != nil && !os.IsNotExist(err) {
				utils.PrintWarning("Failed to remove target overlay %s: %v", b.targetOverlayPath, err)
			}
		}
	}

	return nil
}

// parseScriptMetadata extracts dependencies, scheduler specs, and interactive prompts from shell scripts.
func (b *BaseBuildObject) parseScriptMetadata(ctx context.Context) error {
	if b.buildSource == "" {
		return nil
	}
	if err := b.parseDependencies(); err != nil {
		return err
	}
	if err := b.collectInteractiveInputs(ctx); err != nil {
		return err
	}
	return b.resolveResourceSpec()
}

// parseDependencies reads #DEP: lines from the build script and sets b.dependencies.
func (b *BaseBuildObject) parseDependencies() error {
	deps, err := utils.GetDependenciesFromScript(b.buildSource, config.Global.ParseModuleLoad)
	if err != nil {
		return fmt.Errorf("failed to parse dependencies: %w", err)
	}
	b.dependencies = deps
	return nil
}

// collectInteractiveInputs handles #INTERACTIVE: prompts — checks TTY, supports --yes shortcut,
// and reads user input from stdin. Sets b.interactiveInputs.
func (b *BaseBuildObject) collectInteractiveInputs(ctx context.Context) error {
	prompts, err := utils.GetInteractivePromptsFromScript(b.buildSource)
	if err != nil {
		return fmt.Errorf("failed to parse interactive prompts: %w", err)
	}
	b.interactiveInputs = []string{}
	if len(prompts) == 0 {
		return nil
	}

	// If --yes flag is set, automatically provide empty responses
	if utils.ShouldAnswerYes() {
		for range prompts {
			b.interactiveInputs = append(b.interactiveInputs, "")
		}
		return nil
	}

	// Interactive prompts require a TTY
	if !utils.IsInteractiveShell() {
		return fmt.Errorf("build script for %s requires interactive input, but no TTY is available", b.nameVersion)
	}

	for _, prompt := range prompts {
		msg := strings.ReplaceAll(prompt, `\\n`, "\n")
		msg = strings.ReplaceAll(msg, "\\n", "\n")
		for _, line := range strings.Split(msg, "\n") {
			utils.PrintNote("%s", line)
		}
		fmt.Print("Enter here: ")
		input, err := utils.ReadLineContext(ctx)
		if err != nil {
			return err
		}
		b.interactiveInputs = append(b.interactiveInputs, input)
	}
	return nil
}

// resolveResourceSpec parses scheduler directives from the build script and sets b.scriptSpecs.
// Applies the priority chain: buildDefaults → script directives → current job resources.
func (b *BaseBuildObject) resolveResourceSpec() error {
	specs, err := scheduler.ReadScriptSpecsFromPath(b.buildSource)
	if err != nil {
		return err
	}
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

// NewBuildObject creates a BuildObject from a name/version string
// Format: "name/version" for conda/shell, "name" for def, "prefix/name/version" for ref
// All overlays are stored in imagesDir regardless of type
func NewBuildObject(ctx context.Context, nameVersion string, external bool, imagesDir, tmpDir string, update bool) (BuildObject, error) {
	normalized := utils.NormalizeNameVersion(nameVersion)
	slashCount := strings.Count(normalized, "/")

	// Determine build type based on slash count
	// 0 slashes: system (def is not defined here)
	// 1 slash: conda or shell (e.g., "numpy/1.24")
	// 2+ slashes: ref shell script (e.g., "genomes/hg38/full")
	isRef := slashCount > 1

	// Use fast local storage for conda/app builds; keep stable path for ref/data.
	// Def builds will override this in createConcreteType via resolveTmpDirForDef.
	if isRef {
		tmpDir = resolveTmpDirForRef()
	} else {
		tmpDir = resolveTmpDirForConda()
	}

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
		tmpDir:            tmpDir,
		cntDirPath:        cntDirPath,
		tmpOverlayPath:    tmpOverlayPath,
		targetOverlayPath: targetOverlay,
		update:            update,
	}

	if external {
		// External builds don't need to resolve build source
		return createConcreteType(ctx, base, isRef, tmpDir)
	}

	// Check if already installed or a build is currently in progress (lock file exists).
	// Skip this optimisation in update mode so the correct concrete type is resolved.
	if !update {
		if base.IsInstalled() {
			return newScriptBuildObject(base, isRef)
		}
		if utils.FileExists(base.buildLockPath()) {
			info, readErr := base.readBuildLock()
			if readErr != nil {
				// Corrupt or old empty lock → treat as stale.
				utils.PrintWarning("Corrupt build lock for %s; removing.", utils.StyleName(normalized))
				base.removeBuildLock()
			} else if info.JobID != "" && info.JobID == scheduler.CurrentJobID() {
				// Our own scheduler lock — proceed; createBuildLock() will adopt it
			} else if stale, jobStatus, _ := isBuildLockStale(info); stale {
				detail := info.JobID
				if detail == "" {
					detail = fmt.Sprintf("pid=%d", info.PID)
				}
				utils.PrintWarning("Stale build lock for %s (%s); removing.", utils.StyleName(normalized), detail)
				base.removeBuildLock()
			} else {
				statusHint := ""
				switch jobStatus {
				case scheduler.JobStatusPending:
					statusHint = fmt.Sprintf(" (job %s is pending in queue)", info.JobID)
				case scheduler.JobStatusRunning:
					if info.JobID != "" {
						statusHint = fmt.Sprintf(" (job %s is running)", info.JobID)
					} else if info.PID != 0 {
						statusHint = fmt.Sprintf(" (pid %d is running)", info.PID)
					}
				}
				return nil, fmt.Errorf(
					"build lock found for %s%s.\nLock file: %s",
					utils.StyleName(normalized), statusHint, utils.StylePath(base.buildLockPath()),
				)
			}
		}
	}

	// Resolve build source and determine concrete type
	return createConcreteType(ctx, base, isRef, tmpDir)
}

// NewCondaObjectWithSource creates a CondaBuildObject with custom buildSource
// This is used for the -n flag to create a single sqf with multiple packages or YAML
// buildSource can be:
//   - Path to YAML file (e.g., "/path/to/environment.yml")
//   - Comma-separated package list (e.g., "nvim,nodejs,samtools/1.16")
func NewCondaObjectWithSource(nameVersion, buildSource string, imagesDir, tmpDir string, update bool) (BuildObject, error) {
	normalized := utils.NormalizeNameVersion(nameVersion)

	// Conda builds always use fast local storage
	tmpDir = resolveTmpDirForConda()

	// Make tmpDir absolute
	if absDir, err := filepath.Abs(tmpDir); err == nil {
		tmpDir = absDir
	}

	targetOverlay := filepath.Join(imagesDir, strings.ReplaceAll(normalized, "/", "--")+".sqf")
	if abs, err := filepath.Abs(targetOverlay); err == nil {
		targetOverlay = abs
	}

	tmpOverlayPath, cntDirPath := buildTmpPaths(normalized, tmpDir, ".img")

	utils.PrintDebug("[BUILD OBJECT] Creating conda %s: buildSource=%s, targetOverlay=%s, tmpOverlay=%s, cntDir=%s", nameVersion, buildSource, targetOverlay, tmpOverlayPath, cntDirPath)

	base := &BaseBuildObject{
		nameVersion:       normalized,
		buildSource:       buildSource,
		submitJob:         config.Global.SubmitJob,
		tmpDir:            tmpDir,
		cntDirPath:        cntDirPath,
		tmpOverlayPath:    tmpOverlayPath,
		targetOverlayPath: targetOverlay,
		update:            update,
	}

	return newCondaBuildObject(base)
}

// FromExternalSource creates a BuildObject from an external build script or def file
// All overlays are stored in imagesDir regardless of type
func FromExternalSource(ctx context.Context, targetPrefix, source string, isApptainer bool, imagesDir string) (BuildObject, error) {
	nameVersion := filepath.Base(targetPrefix)
	nameVersion = utils.NormalizeNameVersion(nameVersion)

	// Determine build type from source file extension
	isDef := isApptainer || strings.HasSuffix(source, ".def")
	isShell := strings.HasSuffix(source, ".sh") || strings.HasSuffix(source, ".bash")
	externalType := "app"
	if isShell {
		parsedType, err := utils.GetExternalBuildTypeFromScript(source)
		if err != nil {
			return nil, fmt.Errorf("failed to parse external build type: %w", err)
		}
		externalType = parsedType
	}

	// Build artifacts go next to the target (user controls target location)
	targetDir := resolveTmpDirForExternal(filepath.Dir(targetPrefix), externalType)
	if absDir, err := filepath.Abs(targetDir); err == nil {
		targetDir = absDir
	}
	// Def builds produce a SIF; shell builds use ext3 (.img) or dir-mode (no overlay).
	var ext string
	if isDef {
		ext = ".sif"
	} else if config.Global.Build.UseTmpOverlay {
		ext = ".img"
	}
	tmpOverlayPath, cntDirPath := buildTmpPaths(nameVersion, targetDir, ext)

	utils.PrintDebug("[BUILD OBJECT] Creating external %s: source=%s, targetPrefix=%s, tmpOverlay=%s, cntDir=%s", nameVersion, source, targetPrefix, tmpOverlayPath, cntDirPath)

	base := &BaseBuildObject{
		nameVersion:       nameVersion,
		buildSource:       source,
		submitJob:         config.Global.SubmitJob,
		tmpDir:            targetDir,
		cntDirPath:        cntDirPath,
		tmpOverlayPath:    tmpOverlayPath,
		targetOverlayPath: targetPrefix + ".sqf",
	}

	// Parse script metadata if it's a shell script
	if isShell {
		if err := base.parseScriptMetadata(ctx); err != nil {
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
func createConcreteType(ctx context.Context, base *BaseBuildObject, isRef bool, tmpDir string) (BuildObject, error) {
	// Resolve build source - this determines the actual type based on file extension
	isConda, isContainer, err := resolveBuildSource(base, tmpDir)
	if err != nil {
		return nil, err
	}

	if isConda {
		return newCondaBuildObject(base)
	}

	if isContainer {
		// Internal def builds use the stable writable tmp dir (not fast local scratch).
		defTmpDir := resolveTmpDirForDef()
		if absDir, err := filepath.Abs(defTmpDir); err == nil {
			defTmpDir = absDir
		}
		base.tmpDir = defTmpDir
		base.tmpOverlayPath, base.cntDirPath = buildTmpPaths(base.nameVersion, defTmpDir, ".sif")
		return newDefBuildObject(base)
	}

	// It's a shell script (no extension or .sh/.bash)
	// Parse script metadata
	if err := base.parseScriptMetadata(ctx); err != nil {
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

	// If remote, skip download when target already exists and we're not updating.
	// This avoids unnecessary network requests for installed overlays (e.g. dep
	// graph traversal where bg.update is propagated to already-installed nodes).
	if info.IsRemote && !base.update && base.IsInstalled() {
		utils.PrintDebug("Skipping remote download for %s: target already exists", utils.StyleName(base.nameVersion))
		return false, info.IsContainer, nil
	}

	// If remote, download to tmp directory
	if info.IsRemote {
		dlDir := tmpDir
		if info.IsContainer {
			// Def builds use GetWritableTmpDir; download there so buildSource
			// and tmpOverlayPath (SIF) are in the same directory.
			dlDir = resolveTmpDirForDef()
		}
		localPath, err := DownloadRemoteScript(info, dlDir)
		if err != nil {
			return false, false, fmt.Errorf("failed to download remote build script: %w", err)
		}
		base.buildSource = localPath
		base.isRemote = true
		base.prebuiltLink = info.PrebuiltLink
		utils.PrintDebug("Downloaded remote build script to %s", utils.StylePath(localPath))
	} else {
		base.buildSource = info.Path
		utils.PrintDebug("Using local build script at %s", utils.StylePath(info.Path))
	}

	return false, info.IsContainer, nil
}
