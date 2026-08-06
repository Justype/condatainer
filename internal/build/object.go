package build

import (
	"context"
	"errors"
	"fmt"
	"os"
	"path/filepath"
	"strings"
	"time"

	"log/slog"

	"github.com/Justype/condatainer/catalog"
	"github.com/Justype/condatainer/internal/config"
	"github.com/Justype/condatainer/internal/logging"
	"github.com/Justype/condatainer/internal/overlay"
	"github.com/Justype/condatainer/internal/scheduler"
	"github.com/Justype/condatainer/internal/utils"
)

var ErrTmpOverlayExists = errors.New("temporary overlay already exists")
var ErrBuildCancelled = errors.New("build cancelled by user")

// BuildType is how a target is built — the shape of its source, nothing more.
// What the payload *is* (base, os, app, data) is its catalog.Kind.
type BuildType int

// Predefined build types.
const (
	BuildTypeConda  BuildType = iota + 1 // micromamba
	BuildTypeDef                         // apptainer definition file
	BuildTypeScript                      // recipe run as a shell script
)

// String implements fmt.Stringer for BuildType.
func (bt BuildType) String() string {
	switch bt {
	case BuildTypeConda:
		return "conda"
	case BuildTypeDef:
		return "def"
	case BuildTypeScript:
		return "script"
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
	submitJob         bool // Whether to submit to scheduler (from config at construction time)
	tempSource        bool // Whether buildSource is a temp file this object wrote
	update            bool // If true, rebuild even if overlay already exists (atomic .new swap)
	scriptSpecs       *scheduler.ScriptSpecs
	condaChannelPkg   string // channel-annotated package spec, e.g. "bioconda::star"; set when input uses "::" notation

	// What the payload is. Decides the tmp root and the pack block size, and is
	// exported to the build as CNT_KIND. Authoritative from the recipe's entry;
	// derived from the name shape when no source provides one.
	kind catalog.Kind

	// #INPUT: prompts and the answers collected for them, same order.
	inputs            []string
	interactiveInputs []string

	// Placeholder values for a template recipe, e.g. {"star_version": "2.7.11b"}.
	// Handed to the catalog, which expands the recipe before it is written out.
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
	default: // BuildTypeScript
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
		parsed, err := catalog.ParseDep(dep)
		if err != nil {
			continue
		}
		if parsed.Op == "" {
			// Exact match (existing behaviour).
			if !installed[parsed.NameVersion()] {
				missing = append(missing, dep)
			}
			continue
		}
		// Constraint present: accept any installed version of the same package
		// that the dep admits — at or above the minimum, never above the
		// preferred version.
		prefix := parsed.Name + "/"
		satisfied := false
		for key := range installed {
			if version, ok := strings.CutPrefix(key, prefix); ok && parsed.Satisfies(version) {
				satisfied = true
				break
			}
		}
		if !satisfied {
			missing = append(missing, parsed.NameVersion()) // build the preferred version
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
		if err := utils.MkdirAllShared(parentDir); err != nil {
			return fmt.Errorf("failed to create tmp overlay parent dir %s: %w", parentDir, err)
		}
	}

	logging.FromContext(ctx).Debug("creating temporary overlay", "path", b.tmpOverlayPath)

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
	if err := utils.MkdirAllShared(b.cntDirPath); err != nil {
		return fmt.Errorf("failed to create build cnt dir %s: %w", b.cntDirPath, err)
	}
	buildTmpDir := filepath.Join(buildDir, "tmp")
	if err := utils.MkdirAllShared(buildTmpDir); err != nil {
		return fmt.Errorf("failed to create build tmp dir: %w", err)
	}
	logging.FromContext(ctx).Info("build dir created", "path", buildDir)
	return nil
}

// retargetWorkspace re-sites the build workspace when the recipe's kind implies
// a different tmp root than the name shape did. A no-op unless the recipe
// declared #TYPE:, which is the only way the two disagree.
func (b *BuildObject) retargetWorkspace() {
	root := tmpRootForKind(b.kind)
	if abs, err := filepath.Abs(root); err == nil {
		root = abs
	}
	if root == b.tmpDir {
		return
	}
	b.tmpDir = root
	b.tmpOverlayPath, b.cntDirPath = buildTmpPaths(b.nameVersion, root, ".img")
}

// Cleanup removes the build workspace (materialized recipe, tmp overlay, build dir), plus
// the partial target overlay on failure. Announces the work when there is something
// to remove; a no-op cleanup stays silent.
func (b *BuildObject) Cleanup(failed bool) error {
	log := slog.Default()

	willClean := (b.tempSource && b.buildSource != "") || b.tmpOverlayPath != "" || b.cntDirPath != ""
	if willClean {
		log.Info("cleaning up temporary files")
	}

	// Remove the materialized recipe
	if b.tempSource && b.buildSource != "" {
		if err := os.Remove(b.buildSource); err != nil && !os.IsNotExist(err) {
			log.Warn("failed to remove materialized recipe", "path", b.buildSource, "err", err)
		} else {
			log.Debug("removed materialized recipe", "path", b.buildSource)
		}
	}

	// Remove tmp overlay
	if b.tmpOverlayPath != "" {
		if err := os.Remove(b.tmpOverlayPath); err != nil && !os.IsNotExist(err) {
			log.Warn("failed to remove tmp overlay", "path", b.tmpOverlayPath, "err", err)
		}
	}

	// Remove cnt directory
	if b.cntDirPath != "" {
		cntBaseDir := filepath.Dir(b.cntDirPath)
		log.Debug("cleaning up build directory", "path", cntBaseDir)
		if err := os.RemoveAll(cntBaseDir); err != nil && !os.IsNotExist(err) {
			log.Warn("failed to remove cnt dir", "path", cntBaseDir, "err", err)
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
				log.Warn("failed to remove partial new overlay", "path", newPath, "err", err)
			}
		} else {
			if err := os.Remove(b.targetOverlayPath); err != nil && !os.IsNotExist(err) {
				log.Warn("failed to remove target overlay", "path", b.targetOverlayPath, "err", err)
			}
		}
	}

	if willClean {
		log.Info("temporary files cleaned", "kind", "success")
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
	if err := b.parseInputs(); err != nil {
		return err
	}
	if err := b.collectInteractiveInputs(ctx); err != nil {
		return err
	}
	return b.resolveResourceSpec()
}

// parseDependencies reads #DEP: lines from the build script and sets b.dependencies.
// Skips parsing when the catalog already materialized the recipe, whose deps
// arrive expanded.
func (b *BuildObject) parseDependencies() error {
	if b.dependencies != nil {
		return nil
	}
	deps, err := utils.GetDependenciesFromScript(b.buildSource, config.Global.ParseModuleLoad)
	if err != nil {
		return fmt.Errorf("failed to parse dependencies: %w", err)
	}
	b.dependencies = deps
	return nil
}

// parseInputs reads #INPUT: lines from the build source and sets b.inputs.
// Skips parsing when the catalog already materialized the recipe.
func (b *BuildObject) parseInputs() error {
	if b.inputs != nil || b.tempSource {
		return nil
	}
	file, err := os.Open(b.buildSource)
	if err != nil {
		return fmt.Errorf("failed to open %s: %w", b.buildSource, err)
	}
	defer file.Close()
	recipe, err := catalog.ParseRecipe(b.buildSource, file)
	if err != nil {
		return fmt.Errorf("failed to parse inputs: %w", err)
	}
	b.inputs = recipe.Inputs
	return nil
}

// collectInteractiveInputs asks for every #INPUT: the recipe declares — checks
// TTY, supports the --yes shortcut, and reads answers from stdin.
//
// Answers stay in declaration order and are fed to the recipe on stdin, which
// is the only channel that carries them literally: apptainer shell-evaluates
// env values, so a $ or a backtick in a pasted URL would be mangled or run.
func (b *BuildObject) collectInteractiveInputs(ctx context.Context) error {
	b.interactiveInputs = []string{}
	if len(b.inputs) == 0 {
		return nil
	}

	// If --yes flag is set, automatically provide empty responses
	if utils.ShouldAnswerYes() {
		for range b.inputs {
			b.interactiveInputs = append(b.interactiveInputs, "")
		}
		return nil
	}

	// Prompts require a TTY or piped stdin (e.g. scheduler job with embedded heredoc)
	if !utils.IsInteractiveShell() && !utils.IsStdinPiped() {
		return fmt.Errorf("recipe for %s requires input, but no TTY is available", b.nameVersion)
	}

	log := logging.FromContext(ctx)
	for _, prompt := range b.inputs {
		msg := strings.ReplaceAll(prompt, `\\n`, "\n")
		msg = strings.ReplaceAll(msg, "\\n", "\n")
		for _, line := range strings.Split(msg, "\n") {
			log.Info(line, "kind", "note")
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
	normalized := catalog.Normalize(nameVersion)

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

	// Provisional kind from the name shape. The recipe's entry is authoritative
	// and may override it via #TYPE:, but the workspace has to be sited before
	// anything is looked up — createConcreteType re-points it if the kind moves.
	kind := catalog.DeriveKind(normalized, "", false, "")

	// Use fast local storage for app builds; keep a stable path for data.
	// Def builds will override this in createConcreteType via resolveTmpDirForDef.
	tmpDir = tmpRootForKind(kind)

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

	logging.FromContext(ctx).Debug("creating build object",
		"input", nameVersion, "nameVersion", normalized,
		"targetOverlay", targetOverlay, "tmpOverlay", tmpOverlayPath, "cntDir", cntDirPath)

	base := &BuildObject{
		nameVersion:       normalized,
		submitJob:         config.Global.SubmitJob,
		tmpDir:            tmpDir,
		cntDirPath:        cntDirPath,
		tmpOverlayPath:    tmpOverlayPath,
		targetOverlayPath: targetOverlay,
		update:            update,
		condaChannelPkg:   condaChannelPkg,
		kind:              kind,
	}

	if external {
		// External builds don't need to resolve build source
		return createConcreteType(ctx, base, tmpDir)
	}

	// Check if already installed or a build is currently in progress (lock file exists).
	// Skip this optimisation in update mode so the correct concrete type is resolved.
	if !update {
		if base.IsInstalled() {
			base.buildType = BuildTypeScript
			return base, nil
		}
		if utils.FileExists(base.buildLockPath()) {
			info, readErr := base.readBuildLock()
			if readErr != nil {
				// Corrupt or old empty lock → treat as stale.
				logging.FromContext(ctx).Warn("corrupt build lock, removing", "name", normalized)
				base.removeBuildLock()
			} else if info.JobID != "" && info.JobID == scheduler.CurrentJobID() {
				// Our own scheduler lock — proceed; createBuildLock() will adopt it
			} else if stale, jobStatus, _ := isBuildLockStale(info); stale {
				detail := info.JobID
				if detail == "" {
					detail = fmt.Sprintf("pid=%d", info.PID)
				}
				logging.FromContext(ctx).Warn("stale build lock, removing", "name", normalized, "detail", detail)
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
					normalized, statusHint, base.buildLockPath(),
				)
			}
		}
	}

	// Resolve build source and determine concrete type
	return createConcreteType(ctx, base, tmpDir)
}

// NewCondaObjectWithSource creates a CondaBuildObject with custom buildSource
// This is used for the -n flag to create a single sqf with multiple packages or YAML
// buildSource can be:
//   - Path to YAML file (e.g., "/path/to/environment.yml")
//   - Comma-separated package list (e.g., "nvim,nodejs,samtools/1.16")
func NewCondaObjectWithSource(nameVersion, buildSource string, imagesDir, tmpDir string, update bool) (*BuildObject, error) {
	normalized := catalog.Normalize(nameVersion)

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

	slog.Default().Debug("creating conda build object",
		"nameVersion", nameVersion, "buildSource", buildSource,
		"targetOverlay", targetOverlay, "tmpOverlay", tmpOverlayPath, "cntDir", cntDirPath)

	base := &BuildObject{
		nameVersion:       normalized,
		buildSource:       buildSource,
		submitJob:         config.Global.SubmitJob,
		tmpDir:            tmpDir,
		cntDirPath:        cntDirPath,
		tmpOverlayPath:    tmpOverlayPath,
		targetOverlayPath: targetOverlay,
		update:            update,
		kind:              catalog.KindApp,
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
	nameVersion = catalog.Normalize(nameVersion)

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

	logging.FromContext(ctx).Debug("creating external build object",
		"nameVersion", nameVersion, "source", source,
		"targetPrefix", targetPrefix, "tmpOverlay", tmpOverlayPath, "cntDir", cntDirPath)

	base := &BuildObject{
		nameVersion:       nameVersion,
		buildSource:       source,
		submitJob:         config.Global.SubmitJob,
		tmpDir:            targetDir,
		cntDirPath:        cntDirPath,
		tmpOverlayPath:    tmpOverlayPath,
		targetOverlayPath: targetPrefix + ".sqf",
		kind:              catalog.DeriveKind(nameVersion, "", isDef, externalType),
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
		base.buildType = BuildTypeScript
	} else {
		return nil, fmt.Errorf("unknown source type for %s", source)
	}

	logging.FromContext(ctx).Debug("created external build object", "obj", base.String())
	return base, nil
}

// createConcreteType creates the appropriate concrete BuildObject type
// It determines whether to create a conda, def, or script build based on:
// - build source resolution: .def -> BuildTypeDef, shell script -> BuildTypeScript, not found -> conda
func createConcreteType(ctx context.Context, base *BuildObject, tmpDir string) (*BuildObject, error) {
	// Resolve build source - this determines the actual type based on file extension
	isConda, isContainer, err := resolveBuildSource(ctx, base, tmpDir)
	if err != nil {
		return nil, err
	}

	if isConda {
		if err := base.setupCondaFields(); err != nil {
			return nil, err
		}
		// A conda environment is self-contained software, whatever its name shape.
		base.kind = catalog.KindApp
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

	// The entry may have declared a kind the name shape does not imply, which
	// moves where the build works.
	base.retargetWorkspace()

	// It's a shell script (no extension or .sh/.bash)
	if err := base.parseScriptMetadata(ctx); err != nil {
		return nil, err
	}
	base.buildType = BuildTypeScript
	return base, nil
}

// resolveBuildSource resolves the module through the catalog and materializes
// its recipe. Returns (isConda, isContainer, error).
//
// A name no source provides is conda's — the normal path for a package with no
// recipe anywhere, not a failure. Everything else is written to a temp file:
// the catalog hands back recipe text already expanded, so a template needs no
// separate substitution pass and a remote source no separate download.
func resolveBuildSource(ctx context.Context, base *BuildObject, tmpDir string) (isConda bool, isContainer bool, err error) {
	// Channel-annotated packages (e.g. "bioconda::star") always go through conda.
	if base.condaChannelPkg != "" {
		return true, false, nil
	}

	// An explicit build source (a path or docker:// URI) bypasses resolution.
	if base.buildSource != "" {
		return false, strings.HasSuffix(base.buildSource, ".def"), nil
	}

	cat, err := config.OpenCatalog(ctx)
	if err != nil {
		return false, false, err
	}
	match, found, err := cat.Lookup(ctx, base.nameVersion)
	if err != nil {
		return false, false, err
	}
	// After the lookup, which is what populates Err: an unreachable source is
	// skipped, and the next one — or conda — answers in its place.
	config.WarnUnreachableSources(ctx, cat)
	if !found {
		slog.Default().Debug("no recipe found, using conda", "name", base.nameVersion)
		return true, false, nil
	}
	isContainer = strings.HasSuffix(match.Entry.Path, ".def")
	base.vars = match.Vars
	base.kind = match.Entry.Kind

	// Nothing to materialize for an overlay that already exists and is not
	// being updated — dependency walks reach installed nodes routinely.
	if !base.update && base.IsInstalled() {
		slog.Default().Debug("target already exists, skipping recipe fetch", "name", base.nameVersion)
		return false, isContainer, nil
	}

	recipe, err := cat.Open(ctx, base.nameVersion, base.vars)
	if err != nil {
		return false, false, fmt.Errorf("failed to read recipe for %s: %w", base.nameVersion, err)
	}

	dir := tmpDir
	if isContainer {
		// Def builds keep buildSource and the SIF in one directory.
		dir = resolveTmpDirForDef()
	}
	path, err := writeRecipeFile(recipe, dir, isContainer)
	if err != nil {
		return false, false, err
	}
	base.buildSource = path
	base.tempSource = true
	base.dependencies = recipe.Deps
	base.inputs = recipe.Inputs
	slog.Default().Debug("materialized recipe", "path", path, "source", match.Source.Name)

	return false, isContainer, nil
}

// writeRecipeFile writes an expanded recipe to a temp file the build can run.
func writeRecipeFile(recipe *catalog.Recipe, dir string, isContainer bool) (string, error) {
	if err := utils.MkdirAllShared(dir); err != nil {
		return "", fmt.Errorf("failed to create tmp directory: %w", err)
	}
	name := "cnt--" + strings.ReplaceAll(recipe.Name, "/", "--")
	if isContainer {
		name += ".def"
	} else {
		name += ".sh"
	}
	path := filepath.Join(dir, name)

	file, err := utils.CreateFileWritable(path)
	if err != nil {
		return "", fmt.Errorf("failed to create %s: %w", path, err)
	}
	defer file.Close()
	if _, err := file.Write(recipe.Text); err != nil {
		return "", fmt.Errorf("failed to write %s: %w", path, err)
	}
	utils.ShareWithParentGroup(path)
	return path, nil
}
