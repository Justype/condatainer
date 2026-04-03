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

// BuildType represents the type of build target.
type BuildType int

// Predefined build types.
const (
	BuildTypeConda BuildType = iota + 1
	BuildTypeDef
	BuildTypeShell
	BuildTypeRef
)

// String implements fmt.Stringer for BuildType.
func (bt BuildType) String() string {
	switch bt {
	case BuildTypeConda:
		return "conda"
	case BuildTypeDef:
		return "def"
	case BuildTypeShell:
		return "shell"
	case BuildTypeRef:
		return "ref"
	}
	return "unknown"
}

// ScriptSpecs mirrors the scheduler module's job spec metadata.
type ScriptSpecs = scheduler.ScriptSpecs

// BuildObject holds all state and implements all build operations.
type BuildObject struct {
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
	condaChannelPkg   string // channel-annotated package spec, e.g. "bioconda::star"; set when input uses "::" notation

	// Interactive inputs for shell scripts
	interactiveInputs []string

	// Placeholder variable values for PL (template) scripts, e.g. {"star_version": "2.7.11b"}.
	// Injected as env vars at build time and interpolated into dependency names.
	vars map[string]string

	// Build type and conda-specific fields
	buildType      BuildType
	packageName    string // conda: primary package name
	packageVersion string // conda: primary package version
}

// Common interface implementations for BuildObject

func (b *BuildObject) NameVersion() string         { return b.nameVersion }
func (b *BuildObject) BuildSource() string         { return b.buildSource }
func (b *BuildObject) Dependencies() []string      { return b.dependencies }
func (b *BuildObject) TmpDir() string              { return b.tmpDir }
func (b *BuildObject) TmpOverlayPath() string      { return b.tmpOverlayPath }
func (b *BuildObject) TargetOverlayPath() string   { return b.targetOverlayPath }
func (b *BuildObject) CntDirPath() string          { return b.cntDirPath }
func (b *BuildObject) ScriptSpecs() *ScriptSpecs   { return b.scriptSpecs }
func (b *BuildObject) Update() bool                { return b.update }
func (b *BuildObject) Type() BuildType             { return b.buildType }
func (b *BuildObject) InteractiveInputs() []string { return b.interactiveInputs }

func (b *BuildObject) String() string {
	return fmt.Sprintf(`BuildObject:
		name_version: %s
		build_source_type: %s
		build_source: %s
		dependencies: %v
		script_specs: %v
		tmp_overlay_path: %s
		target_overlay_path: %s
		cnt_dir_path: %s`,
		b.nameVersion, b.buildType, b.buildSource,
		b.dependencies, b.scriptSpecs,
		b.tmpOverlayPath, b.targetOverlayPath, b.cntDirPath,
	)
}

// Build dispatches to the appropriate build implementation based on buildType.
func (b *BuildObject) Build(ctx context.Context, buildDeps bool) error {
	switch b.buildType {
	case BuildTypeConda:
		return b.buildConda(ctx)
	case BuildTypeDef:
		return b.buildDef(ctx)
	default: // BuildTypeShell, BuildTypeRef
		return b.buildScript(ctx, buildDeps)
	}
}

// setupCondaFields parses packageName and packageVersion from nameVersion/buildSource.
// Must be called before setting buildType = BuildTypeConda.
func (b *BuildObject) setupCondaFields() error {
	parts := strings.Split(b.nameVersion, "/")
	if b.buildSource != "" {
		// Custom buildSource (YAML or comma-separated packages): name without version is OK.
		if len(parts) == 2 {
			b.packageName = parts[0]
			b.packageVersion = parts[1]
		} else {
			b.packageName = b.nameVersion
			b.packageVersion = "env"
		}
	} else {
		// Standard conda package: must be name/version.
		if len(parts) != 2 {
			if b.condaChannelPkg != "" {
				return fmt.Errorf("channel-annotated package requires a version (e.g. %s=1.0)", b.condaChannelPkg)
			}
			return fmt.Errorf("conda package must be in format name/version, got: %s", b.nameVersion)
		}
		b.packageName = parts[0]
		b.packageVersion = parts[1]
		if b.condaChannelPkg != "" {
			b.packageName = b.condaChannelPkg
		}
	}
	return nil
}
func (b *BuildObject) RequiresScheduler() bool {
	return b.submitJob && (config.Global.Build.AlwaysSubmit || scheduler.HasSchedulerSpecs(b.scriptSpecs))
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
func (b *BuildObject) buildLockPath() string {
	return b.targetOverlayPath + ".lock"
}

// LockPath returns the lock file path.
func (b *BuildObject) LockPath() string { return b.buildLockPath() }

// writeBuildLock atomically creates the lock file and writes JSON metadata.
// Returns an os.ErrExist-wrapped error if the lock already exists.
func (b *BuildObject) writeBuildLock(info BuildLockInfo) error {
	return acquireBuildLockFile(b.buildLockPath(), info)
}

// readBuildLock reads and parses the lock file JSON.
// An empty or corrupt file (old empty-lock format) returns a zero-value struct.
func (b *BuildObject) readBuildLock() (BuildLockInfo, error) {
	return readBuildLockFile(b.buildLockPath())
}

// updateBuildLock overwrites the lock file contents.
func (b *BuildObject) updateBuildLock(info BuildLockInfo) error {
	return overwriteBuildLockFile(b.buildLockPath(), info)
}

// createBuildLock creates the lock file for a local build.
// If a scheduler job lock already exists for this job, it adopts that lock
// (updating it with the runtime node and PID) instead of failing.
func (b *BuildObject) createBuildLock() error {
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
func (b *BuildObject) removeBuildLock() {
	os.Remove(b.buildLockPath()) //nolint:errcheck
}

// effectiveNcpus returns the effective total CPUs (CpusPerTask × TasksPerNode) for this build.
func (b *BuildObject) effectiveNcpus() int {
	rs := buildEffectiveResourceSpec(b.scriptSpecs)
	cpus := rs.CpusPerTask
	if rs.TasksPerNode > 1 {
		cpus *= rs.TasksPerNode
	}
	return cpus
}

func (b *BuildObject) IsInstalled() bool {
	_, err := os.Stat(b.targetOverlayPath)
	return err == nil
}

func (b *BuildObject) GetMissingDependencies() ([]string, error) {
	installed := getInstalledOverlays()
	var missing []string
	for _, dep := range b.dependencies {
		nameVersion, op, minVersion := utils.SplitDepConstraint(dep)
		if op == "" {
			// Exact match (existing behaviour).
			if !installed[nameVersion] {
				missing = append(missing, dep)
			}
			continue
		}
		// Constraint present: accept any installed version of the same package
		// that satisfies op+minVersion and does not exceed the preferred version.
		name := nameVersion
		preferredVer := ""
		if idx := strings.LastIndex(nameVersion, "/"); idx >= 0 {
			name = nameVersion[:idx]
			preferredVer = nameVersion[idx+1:]
		}
		prefix := name + "/"
		satisfied := false
		for key := range installed {
			if after, ok := strings.CutPrefix(key, prefix); ok {
				installedVer := after
				if utils.DepSatisfiedByVersion(installedVer, op, minVersion, preferredVer) {
					satisfied = true
					break
				}
			}
		}
		if !satisfied {
			missing = append(missing, nameVersion) // build the preferred version
		}
	}
	return missing, nil
}

func (b *BuildObject) CreateTmpOverlay(ctx context.Context, force bool) error {
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
func (b *BuildObject) CreateBuildDirs(ctx context.Context, force bool) error {
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

	// Create cnt-$USER leaf first with appropriate permissions (0700 under /tmp, 0775 elsewhere).
	tmpBase := b.tmpDir
	if tmpBase == "" {
		tmpBase = filepath.Dir(buildDir)
	}
	if err := utils.EnsureTmpSubdir(tmpBase); err != nil {
		return fmt.Errorf("failed to create tmp dir %s: %w", tmpBase, err)
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

func (b *BuildObject) Cleanup(failed bool) error {
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
func (b *BuildObject) parseScriptMetadata(ctx context.Context) error {
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
// Skips parsing if dependencies were already populated from remote metadata.
func (b *BuildObject) parseDependencies() error {
	if b.dependencies != nil {
		// Dependencies were pre-populated from metadata; apply var interpolation if needed.
		if len(b.vars) > 0 {
			for i, dep := range b.dependencies {
				b.dependencies[i] = utils.InterpolateVars(dep, b.vars)
			}
		}
		return nil
	}
	deps, err := utils.GetDependenciesFromScript(b.buildSource, config.Global.ParseModuleLoad)
	if err != nil {
		return fmt.Errorf("failed to parse dependencies: %w", err)
	}
	if len(b.vars) > 0 {
		for i, dep := range deps {
			deps[i] = utils.InterpolateVars(dep, b.vars)
		}
	}
	b.dependencies = deps
	return nil
}

// collectInteractiveInputs handles #INTERACTIVE: prompts — checks TTY, supports --yes shortcut,
// and reads user input from stdin. Sets b.interactiveInputs.
func (b *BuildObject) collectInteractiveInputs(ctx context.Context) error {
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

	// Interactive prompts require a TTY or piped stdin (e.g. scheduler job with embedded heredoc)
	if !utils.IsInteractiveShell() && !utils.IsStdinPiped() {
		return fmt.Errorf("build script for %s requires interactive input, but no TTY is available", b.nameVersion)
	}

	for _, prompt := range prompts {
		if len(b.vars) > 0 {
			prompt = utils.InterpolateVars(prompt, b.vars)
		}
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
func (b *BuildObject) resolveResourceSpec() error {
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
func NewBuildObject(ctx context.Context, nameVersion string, external bool, imagesDir, tmpDir string, update bool) (*BuildObject, error) {
	normalized := utils.NormalizeNameVersion(nameVersion)

	// Handle channel annotation (e.g. "bioconda::star/2.7.11b"):
	// strip the channel prefix for path/naming; keep it for the micromamba spec.
	var condaChannelPkg string
	if colonIdx := strings.Index(normalized, "::"); colonIdx != -1 {
		channel := normalized[:colonIdx]
		rest := normalized[colonIdx+2:] // "star/2.7.11b" or "star"
		pkgName := rest
		if before, _, ok := strings.Cut(rest, "/"); ok {
			pkgName = before
		}
		condaChannelPkg = channel + "::" + pkgName
		normalized = rest // strip channel prefix for sqf naming and env path
	}

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

	base := &BuildObject{
		nameVersion:       normalized,
		submitJob:         config.Global.SubmitJob,
		tmpDir:            tmpDir,
		cntDirPath:        cntDirPath,
		tmpOverlayPath:    tmpOverlayPath,
		targetOverlayPath: targetOverlay,
		update:            update,
		condaChannelPkg:   condaChannelPkg,
	}

	if external {
		// External builds don't need to resolve build source
		return createConcreteType(ctx, base, isRef, tmpDir)
	}

	// Check if already installed or a build is currently in progress (lock file exists).
	// Skip this optimisation in update mode so the correct concrete type is resolved.
	if !update {
		if base.IsInstalled() {
			if isRef {
				base.buildType = BuildTypeRef
			} else {
				base.buildType = BuildTypeShell
			}
			return base, nil
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
func NewCondaObjectWithSource(nameVersion, buildSource string, imagesDir, tmpDir string, update bool) (*BuildObject, error) {
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

	base := &BuildObject{
		nameVersion:       normalized,
		buildSource:       buildSource,
		submitJob:         config.Global.SubmitJob,
		tmpDir:            tmpDir,
		cntDirPath:        cntDirPath,
		tmpOverlayPath:    tmpOverlayPath,
		targetOverlayPath: targetOverlay,
		update:            update,
	}

	if err := base.setupCondaFields(); err != nil {
		return nil, err
	}
	base.buildType = BuildTypeConda
	return base, nil
}

// FromExternalSource creates a BuildObject from an external build script or def file
// All overlays are stored in imagesDir regardless of type
func FromExternalSource(ctx context.Context, targetPrefix, source string, isApptainer bool, imagesDir string) (*BuildObject, error) {
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

	base := &BuildObject{
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

	if isDef {
		base.buildType = BuildTypeDef
	} else if isShell {
		base.buildType = BuildTypeShell // External sources are not ref overlays
	} else {
		return nil, fmt.Errorf("unknown source type for %s", source)
	}

	utils.PrintDebug("[EXTERNAL] Created BuildObject from external path: %s", base.String())
	return base, nil
}

// createConcreteType creates the appropriate concrete BuildObject type
// It determines whether to create a conda, def, or script build based on:
// - isRef flag (from slash count > 1)
// - build source resolution: .def -> BuildTypeDef, shell script -> BuildTypeShell/Ref, not found -> conda
func createConcreteType(ctx context.Context, base *BuildObject, isRef bool, tmpDir string) (*BuildObject, error) {
	// Resolve build source - this determines the actual type based on file extension
	isConda, isContainer, err := resolveBuildSource(base, tmpDir)
	if err != nil {
		return nil, err
	}

	if isConda {
		if err := base.setupCondaFields(); err != nil {
			return nil, err
		}
		base.buildType = BuildTypeConda
		return base, nil
	}

	if isContainer {
		// Internal def builds use the stable writable tmp dir (not fast local scratch).
		defTmpDir := resolveTmpDirForDef()
		if absDir, err := filepath.Abs(defTmpDir); err == nil {
			defTmpDir = absDir
		}
		base.tmpDir = defTmpDir
		base.tmpOverlayPath, base.cntDirPath = buildTmpPaths(base.nameVersion, defTmpDir, ".sif")
		base.buildType = BuildTypeDef
		return base, nil
	}

	// It's a shell script (no extension or .sh/.bash)
	if err := base.parseScriptMetadata(ctx); err != nil {
		return nil, err
	}
	if isRef {
		base.buildType = BuildTypeRef
	} else {
		base.buildType = BuildTypeShell
	}
	return base, nil
}

// resolveBuildSource finds the build script and determines the build type
// Returns (isConda, isContainer, error):
// - isConda=true: no build script found, use conda
// - isContainer=true: found .def file
// - both false: found shell script (no extension)
// First checks local and remote build scripts, falls back to conda if not found
func resolveBuildSource(base *BuildObject, tmpDir string) (isConda bool, isContainer bool, err error) {
	// Channel-annotated packages (e.g. "bioconda::star") always go through conda.
	if base.condaChannelPkg != "" {
		return true, false, nil
	}

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

	// Pre-populate dependencies from metadata when available (avoids re-parsing the script file).
	// nil means "not set" (local script or metadata without deps field); []string{} means "no deps".
	if info.Deps != nil {
		base.dependencies = info.Deps
	}

	// Capture placeholder variable values for PL scripts.
	if vars := info.CurrentVars(); len(vars) > 0 {
		base.vars = vars
	}

	return false, info.IsContainer, nil
}

// substituteTemplateFile reads srcPath, replaces all {key} tokens using vars,
// writes the result to a temp file in tmpDir, and returns the temp path.
// The caller is responsible for removing the temp file when done.
func substituteTemplateFile(srcPath string, vars map[string]string, tmpDir string) (string, error) {
	data, err := os.ReadFile(srcPath)
	if err != nil {
		return "", fmt.Errorf("failed to read template file: %w", err)
	}
	content := utils.InterpolateVars(string(data), vars)
	tmpPath := filepath.Join(tmpDir, "pl-"+filepath.Base(srcPath))
	if err := os.WriteFile(tmpPath, []byte(content), utils.PermFile); err != nil {
		return "", fmt.Errorf("failed to write substituted file: %w", err)
	}
	return tmpPath, nil
}
