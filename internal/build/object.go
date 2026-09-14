package build

import (
	"bytes"
	"context"
	"errors"
	"fmt"
	"os"
	"path/filepath"
	"slices"
	"strings"
	"time"

	"log/slog"

	"github.com/Justype/condatainer/catalog"
	"github.com/Justype/condatainer/internal/artifact/meta"
	"github.com/Justype/condatainer/internal/config"
	"github.com/Justype/condatainer/internal/image/ext3"
	"github.com/Justype/condatainer/internal/image/producer"
	"github.com/Justype/condatainer/internal/logging"
	"github.com/Justype/condatainer/internal/scheduler"
	"github.com/Justype/condatainer/internal/utils"
)

var ErrTmpOverlayExists = errors.New("temporary overlay already exists")
var ErrBuildCancelled = errors.New("build cancelled by user")

// BuildType is how a target is built — the shape of its source, nothing more.
// What the payload *is* (os, app, data, env) is its catalog.Type.
//
// It is meta.BuildType, the type the manifest records. The unset value is "",
// which SourceSpec.BuildType returns until a source is resolved. There is no
// constant for meta.BuildTypeSnapshot: a freeze packs an overlay, it does not
// build one.
type BuildType = meta.BuildType

// Predefined build types.
const (
	BuildTypeConda  = meta.BuildTypeConda  // micromamba
	BuildTypeDef    = meta.BuildTypeDef    // apptainer definition file
	BuildTypeScript = meta.BuildTypeScript // recipe run as a shell script
)

// ScriptSpecs mirrors the scheduler module's job spec metadata.
type ScriptSpecs = scheduler.ScriptSpecs

// BuildObject holds all state and implements all build operations.
type BuildObject struct {
	buildSource string
	submitJob   bool // Whether to submit to scheduler (from config at construction time)
	tempSource  bool // Whether buildSource is a temp file this object wrote
	update      bool // If true, rebuild even if overlay already exists (atomic .new swap)
	// storeOverflow files the finished build under its identity instead of the
	// bare name, so a second build of a name is kept beside the first rather
	// than skipped. Set per object: it is what the user asked to build, never
	// its dependencies, which stay ordinary and are skipped when installed.
	storeOverflow bool
	// locked marks a rebuild described by a project lock rather than the
	// catalog: dependencies are mounted at supplied paths, the upstream digest
	// is the recorded one, and the output belongs to the caller. See locked.go.
	locked bool
	// lockedKeys are the keys the lock recorded. Only their schemes are used —
	// a rebuild derives with those rather than with the current pair, so a
	// scheme that shipped since does not make the result incomparable.
	lockedKeys      meta.Keys
	scriptSpecs     *scheduler.ScriptSpecs
	condaChannelPkg string // channel-annotated package spec, e.g. "bioconda::star"; set when input uses "::" notation

	// The recipe's #INPUT: questions and the answers, same order. Answers reach
	// the recipe on stdin and go no further — they are not in Spec, so they
	// cannot enter the manifest.
	inputPrompts []string
	inputAnswers []string

	// Placeholder values for a template recipe, e.g. {"star_version": "2.7.11b"}.
	// Handed to the catalog, which expands the recipe before it is written out.
	vars map[string]string
	// catalogSource is the exact collection whose recipe was selected. Its
	// descriptor supplies automatic prebuilt endpoints; local and Conda builds
	// leave it nil.
	catalogSource *catalog.Source

	// Build type and conda-specific fields
	buildType      BuildType
	packageName    string // conda: primary package name
	packageVersion string // conda: primary package version

	// What this build produces, and the only source for the name and type: the
	// manifest and the build environment both read them here.
	spec Spec

	// Where the work happens. Derived as one set by workspaceFor, so re-siting
	// it cannot move some paths and leave others behind. Removed after the build.
	ws Workspace

	// Where the image lands. Prepared is filled in once the build lock is held,
	// since it derives from the lock owner.
	tgt Target

	// embedded are the rebuild sources staged into /.cnt verbatim beside the two
	// metadata documents: a recipe at resolution or a Conda build's exports once
	// the environment exists.
	embedded []embeddedFile

	// What the build learned once its sources and dependencies were known.
	keys         meta.Keys
	dependencies []meta.Dependency
	// depImagePaths is the image each recorded edge's keys were read from, so
	// the capsule is composed from the same file rather than from a second
	// name resolution that could land elsewhere.
	depImagePaths      map[string]string
	provenanceComplete *bool
	buildTools         meta.BuildTools

	// foreignRoot is set only by FromForeignRoot: a .sif or sandbox directory
	// condatainer did not build, packed by buildForeign instead of buildDef.
	foreignRoot *foreignRoot
}

// embeddedFile is one rebuild source staged into /.cnt and named by
// manifest.source.files.
type embeddedFile struct {
	SourceFile
	isSource bool
}

// embedSource records a file the image was built from.
func (b *BuildObject) embedSource(file SourceFile) {
	b.embed(embeddedFile{SourceFile: file, isSource: true})
}

// embed stages a file into /.cnt, replacing any earlier one of the same name so
// a re-resolved build does not stage two.
func (b *BuildObject) embed(file embeddedFile) {
	for i, have := range b.embedded {
		if have.Name == file.Name {
			b.embedded[i] = file
			return
		}
	}
	b.embedded = append(b.embedded, file)
}

// Spec returns the resolved description of the image this build produces.
func (b *BuildObject) Spec() Spec { return b.spec }

// Manifest returns what this build's image records about itself, including the
// sources it embeds — which are known only once they have been captured.
func (b *BuildObject) Manifest() meta.Manifest {
	m := b.spec.Manifest()
	for _, file := range b.embedded {
		if file.isSource {
			m.Source.Files = append(m.Source.Files, file.Name)
		}
	}
	m.Keys = b.keys
	m.Dependencies = b.dependencies
	m.ProvenanceComplete = b.provenanceComplete
	m.Build.Tools = b.buildTools
	return m
}

// Runtime returns the mount-time contract this build embeds in its image.
func (b *BuildObject) Runtime() meta.Runtime { return b.spec.Runtime() }

// Common interface implementations for BuildObject

func (b *BuildObject) NameVersion() string       { return b.spec.Image.Name }
func (b *BuildObject) BuildSource() string       { return b.buildSource }
func (b *BuildObject) Dependencies() []string    { return b.spec.Dependencies }
func (b *BuildObject) TmpDir() string            { return b.ws.Root }
func (b *BuildObject) TmpOverlayPath() string    { return b.ws.Overlay }
func (b *BuildObject) TargetOverlayPath() string { return b.tgt.Path }
func (b *BuildObject) CntDirPath() string        { return b.ws.CntDir }
func (b *BuildObject) ScriptSpecs() *ScriptSpecs { return b.scriptSpecs }
func (b *BuildObject) Update() bool              { return b.update }

// SetStoreOverflow files this build under its identity rather than the bare
// name. See the storeOverflow field.
func (b *BuildObject) SetStoreOverflow(v bool) { b.storeOverflow = v }

// StoreOverflow reports whether this build is filed by identity.
func (b *BuildObject) StoreOverflow() bool    { return b.storeOverflow }
func (b *BuildObject) BuildType() BuildType   { return b.buildType }
func (b *BuildObject) InputAnswers() []string { return b.inputAnswers }

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
		b.spec.Image.Name, b.buildType, b.buildSource,
		b.spec.Dependencies, b.scriptSpecs,
		b.ws.Overlay, b.tgt.Path, b.ws.CntDir,
	)
}

// Build dispatches to the appropriate build implementation based on buildType.
func (b *BuildObject) Build(ctx context.Context, buildDeps bool) error {
	switch b.buildType {
	case BuildTypeConda:
		return b.buildConda(ctx)
	case BuildTypeDef:
		if b.foreignRoot != nil {
			return b.buildForeign(ctx)
		}
		return b.buildDef(ctx)
	default: // BuildTypeScript
		return b.buildScript(ctx, buildDeps)
	}
}

// setupCondaFields parses packageName and packageVersion from nameVersion/buildSource.
// Must be called before setting buildType = BuildTypeConda.
func (b *BuildObject) setupCondaFields() error {
	parts := strings.Split(b.spec.Image.Name, "/")
	if b.buildSource != "" {
		// Custom buildSource (YAML or comma-separated packages): name without version is OK.
		if len(parts) == 2 {
			b.packageName = parts[0]
			b.packageVersion = parts[1]
		} else {
			b.packageName = b.spec.Image.Name
			b.packageVersion = "env"
		}
	} else {
		// Standard conda package: must be name/version.
		if len(parts) != 2 {
			if b.condaChannelPkg != "" {
				return fmt.Errorf("channel-annotated package requires a version (e.g. %s=1.0)", b.condaChannelPkg)
			}
			return fmt.Errorf("conda package must be in format name/version, got: %s", b.spec.Image.Name)
		}
		b.packageName = parts[0]
		b.packageVersion = parts[1]
		if b.condaChannelPkg != "" {
			b.packageName = b.condaChannelPkg
		}
	}
	return nil
}

// setCondaSpec fills the Spec for a Conda build. A Conda environment carries its
// own libraries, so it is always an app contributing only its prefix — there is
// no recipe to take #ENV: or a description from.
func (b *BuildObject) setCondaSpec() {
	b.spec.Image.Type = catalog.TypeApp
	b.spec.Image.Prefix = meta.Prefix(b.spec.Image.Name, catalog.TypeApp)
	b.spec.Image.Env = nil

	src := &CondaSource{Channels: slices.Clone(config.Global.Build.Channels)}
	switch {
	case b.buildSource != "" && utils.IsCondaFile(b.buildSource):
		data, err := os.ReadFile(b.buildSource)
		if err != nil {
			// The file is read again at execution; a failure here only costs
			// the captured copy, so resolution does not fail on it.
			slog.Default().Debug("could not capture conda input file", "path", b.buildSource, "err", err)
		}
		src.File = &SourceFile{Name: filepath.Base(b.buildSource), Data: data}
	case b.buildSource != "":
		src.Packages = strings.Split(b.buildSource, ",")
	default:
		src.Package = &CondaPackage{
			Name:           b.packageName,
			Version:        b.packageVersion,
			ChannelPackage: b.condaChannelPkg,
		}
	}
	b.spec.Source = SourceSpec{Conda: src}
}

func (b *BuildObject) RequiresScheduler() bool {
	return b.submitJob && (config.Global.Build.AlwaysSubmit || scheduler.HasSchedulerSpecs(b.scriptSpecs))
}

// BuildLockInfo holds metadata stored inside a build lock file.
type BuildLockInfo = producer.Info

// LockPath returns the lock file path.
func (b *BuildObject) LockPath() string { return b.tgt.Lock }

// writeBuildLock atomically creates the lock file and writes JSON metadata.
// Returns an os.ErrExist-wrapped error if the lock already exists.
func (b *BuildObject) writeBuildLock(info BuildLockInfo) error {
	return acquireBuildLockFile(b.tgt.Lock, info)
}

// readBuildLock reads and parses the lock file JSON.
// An empty or corrupt file (old empty-lock format) returns a zero-value struct.
func (b *BuildObject) readBuildLock() (BuildLockInfo, error) {
	return readBuildLockFile(b.tgt.Lock)
}

// updateBuildLock overwrites the lock file contents.
func (b *BuildObject) updateBuildLock(info BuildLockInfo) error {
	return overwriteBuildLockFile(b.tgt.Lock, info)
}

// clearStaleLock removes a build lock whose owner is gone, and the partial output
// it left behind. A live owner is an error instead. Runs in update mode too,
// which is how a user retries after a killed build.
func (b *BuildObject) clearStaleLock(ctx context.Context) error {
	if !utils.FileExists(b.tgt.Lock) {
		return nil
	}
	name := b.spec.Image.Name

	info, readErr := b.readBuildLock()
	if readErr != nil {
		// Corrupt or old empty lock → treat as stale.
		logging.FromContext(ctx).Warn("corrupt build lock, removing", "name", name)
		b.removeBuildLock()
		return nil
	}
	if info.JobID != "" && info.JobID == scheduler.CurrentJobID() {
		return nil // our own scheduler lock; createBuildLock() adopts it
	}

	stale, jobStatus, _ := isBuildLockStale(info)
	if stale {
		detail := info.JobID
		if detail == "" {
			detail = fmt.Sprintf("pid=%d", info.PID)
		}
		logging.FromContext(ctx).Warn("stale build lock, removing", "name", name, "detail", detail)
		b.removeBuildLock()
		b.removeOrphanedOutput(ctx, info)
		b.removeOwnerWorkspace(info)
		return nil
	}

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
	return fmt.Errorf("build lock found for %s%s.\nLock file: %s", name, statusHint, b.tgt.Lock)
}

// createBuildLock creates the lock file for a local build.
// If a scheduler job lock already exists for this job, it adopts that lock
// (updating it with the runtime node and PID) instead of failing.
func (b *BuildObject) createBuildLock() error {
	info := BuildLockInfo{
		Runner:    "local",
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
		if readErr == nil && existing.Runner != "local" && existing.Runner != "" {
			if myJobID := scheduler.CurrentJobID(); myJobID != "" && existing.JobID == myJobID {
				// Adopt the lock: update with runtime node + PID.
				existing.Node = shortHostname()
				existing.PID = os.Getpid()
				if err := b.updateBuildLock(existing); err != nil {
					return err
				}
				b.setPreparedPath(existing)
				if err := b.adoptWorkspace(existing); err != nil {
					b.removeBuildLock()
					return err
				}
				return nil
			}
		}
		return fmt.Errorf("build already in progress: lock file exists at %s", b.tgt.Lock)
	}
	b.setPreparedPath(info)
	if err := b.adoptWorkspace(info); err != nil {
		b.removeBuildLock()
		return err
	}
	return nil
}

// adoptWorkspace re-sites all temporary build state under the identity in the
// lock. Local builds normally already use this path; scheduler jobs move their
// process-private resolved recipe to the adopted scheduler-job workspace.
func (b *BuildObject) adoptWorkspace(info BuildLockInfo) error {
	next := workspaceForOwner(b.spec.Image.Name, b.ws.BaseRoot, b.ws.imageExt, b.ws.isDef, info)
	if next.Root == b.ws.Root {
		return nil
	}
	old := b.ws
	if utils.DirExists(old.Root) {
		if err := utils.MkdirAllShared(filepath.Dir(next.Root)); err != nil {
			return err
		}
		if err := os.Rename(old.Root, next.Root); err != nil {
			return fmt.Errorf("failed to adopt build workspace %s: %w", next.Root, err)
		}
	}
	b.ws = next
	if b.tempSource && b.buildSource == old.Source {
		b.buildSource = next.Source
	}
	return nil
}

func (b *BuildObject) removeOwnerWorkspace(info BuildLockInfo) {
	ws := workspaceForOwner(b.spec.Image.Name, b.ws.BaseRoot, b.ws.imageExt, b.ws.isDef, info)
	os.RemoveAll(ws.Root) //nolint:errcheck
	utils.RemoveDirIfEmpty(filepath.Dir(ws.Root))
}

// setPreparedPath records where this build writes its output, derived from the
// lock it actually holds. An adopted scheduler lock names the job, so the
// submitter and the job that runs it agree on the path.
func (b *BuildObject) setPreparedPath(info BuildLockInfo) {
	b.tgt.Prepared = preparedPathFor(b.tgt.Path, info)
}

// PreparedPath is where this build writes its output before installing it.
// Empty until the build lock is held.
func (b *BuildObject) PreparedPath() string { return b.tgt.Prepared }

// removeBuildLock removes the lock file on build completion or failure.
func (b *BuildObject) removeBuildLock() {
	os.Remove(b.tgt.Lock) //nolint:errcheck
}

// removeOrphanedOutput removes the partial output a dead lock owner left behind.
// The path is recomputed from the lock, which is why prepared paths derive from
// the owner: a build killed mid-write never cleans up after itself.
func (b *BuildObject) removeOrphanedOutput(ctx context.Context, info BuildLockInfo) {
	orphan := preparedPathFor(b.tgt.Path, info)
	if err := os.Remove(orphan); err == nil {
		logging.FromContext(ctx).Warn("removed partial output from stale build", "path", orphan)
	}
}

// effectiveNcpus returns the effective total CPUs (CpusPerTask × TasksPerNode) for this build.
func (b *BuildObject) effectiveNcpus() int {
	rs := EffectiveResourceSpec(b.scriptSpecs)
	cpus := rs.CpusPerTask
	if rs.TasksPerNode > 1 {
		cpus *= rs.TasksPerNode
	}
	return cpus
}

// IsInstalled reports whether this image is already built. The configured
// default root (config.BaseRecipeName, still conventionally named ".../base")
// is searched across every image path, not just its target, so one from a
// shared install is not rebuilt into the user's own directory.
func (b *BuildObject) IsInstalled() bool {
	if b.spec.Image.Name == config.BaseRecipeName() {
		return config.FindBaseImage() != ""
	}
	_, err := os.Stat(b.tgt.Path)
	return err == nil
}

func (b *BuildObject) GetMissingDependencies() ([]string, error) {
	installed := getInstalledOverlays()
	var missing []string
	for _, dep := range b.spec.Dependencies {
		// An overlay path is satisfied by the file being there. The installed map
		// is keyed by name, so a path would never match it and would look
		// permanently missing.
		if utils.IsOverlay(dep) || utils.IsSif(dep) {
			if !utils.FileExists(dep) {
				missing = append(missing, dep)
			}
			continue
		}
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
	buildDir := b.ws.BuildDir
	stale := utils.FileExists(b.ws.Overlay) || (b.ws.CntDir != "" && utils.DirExists(buildDir))
	if stale {
		if !force {
			return fmt.Errorf("%w: %s", ErrTmpOverlayExists, b.ws.Overlay)
		}
		if err := os.Remove(b.ws.Overlay); err != nil && !os.IsNotExist(err) {
			return fmt.Errorf("failed to remove existing tmp overlay: %w", err)
		}
		os.RemoveAll(buildDir) //nolint:errcheck
	}

	// Ensure parent directory for tmp overlay exists (mkdir -p)
	parentDir := filepath.Dir(b.ws.Overlay)
	if parentDir != "" {
		if err := utils.MkdirAllShared(parentDir); err != nil {
			return fmt.Errorf("failed to create tmp overlay parent dir %s: %w", parentDir, err)
		}
	}

	logging.FromContext(ctx).Debug("creating temporary overlay", "path", b.ws.Overlay)

	// Sparse for speed, quiet because a scratch image's specs are noise.
	if err := ext3.CreateWithOptions(ctx, &ext3.CreateOptions{
		Path: b.ws.Overlay, SizeMB: config.Global.Build.AppTmpOverlaySizeMB,
		UID: os.Getuid(), GID: os.Getgid(),
		Profile: ext3.ProfileDefault, Sparse: true, FilesystemType: "ext3", Quiet: true,
	}); err != nil {
		return fmt.Errorf("failed to create temporary overlay: %w", err)
	}

	return nil
}

// CreateBuildDirs creates host directories for dir-mode builds (app_tmp_overlay=false).
// Layout: <buildDir>/cnt/ (bound as /cnt) and <buildDir>/tmp/ (bound as ScratchPath).
// Checks both dir-mode (buildDir) and ext3-mode (.img) artifacts for cross-mode stale detection.
func (b *BuildObject) CreateBuildDirs(ctx context.Context, force bool) error {
	buildDir := b.ws.BuildDir
	// Check both dir-mode artifact (buildDir) and ext3-mode artifact (.img)
	stale := utils.DirExists(buildDir) || utils.FileExists(b.ws.Overlay)
	if stale {
		if !force {
			return fmt.Errorf("%w: %s", ErrTmpOverlayExists, buildDir)
		}
		os.RemoveAll(buildDir)  //nolint:errcheck — clean dir-mode artifact
		os.Remove(b.ws.Overlay) //nolint:errcheck — clean ext3-mode artifact (no-op if "")
	}

	if err := ensureWorkspaceRoot(b); err != nil {
		return err
	}
	if err := utils.MkdirAllShared(b.ws.CntDir); err != nil {
		return fmt.Errorf("failed to create build cnt dir %s: %w", b.ws.CntDir, err)
	}
	buildTmpDir := filepath.Join(buildDir, "tmp")
	if err := utils.MkdirAllShared(buildTmpDir); err != nil {
		return fmt.Errorf("failed to create build tmp dir: %w", err)
	}
	logging.FromContext(ctx).Info("build dir created", "path", buildDir)
	return nil
}

// retargetWorkspace re-sites the build workspace when the recipe's type implies
// a different tmp root than the name shape did. A no-op unless the recipe
// declared #TYPE:, which is the only way the two disagree.
func (b *BuildObject) retargetWorkspace() {
	typ := b.spec.Image.Type
	root := tmpRootForType(typ)
	if abs, err := filepath.Abs(root); err == nil {
		root = abs
	}
	// Re-derived from the type, not carried over: a guess of app corrected to
	// data must lose its ext3 image, not just move it to another root.
	b.ws = workspaceFor(b.spec.Image.Name, root, appExt3ScratchExt(typ), false)
}

// Cleanup removes the build workspace (materialized recipe, tmp overlay, build dir), plus
// the partial target overlay on failure. Announces the work when there is something
// to remove; a no-op cleanup stays silent.
func (b *BuildObject) Cleanup(failed bool) error {
	log := slog.Default()

	willClean := (b.tempSource && b.buildSource != "") || b.ws.Overlay != "" || b.ws.CntDir != ""
	if willClean {
		log.Info("cleaning up temporary files")
	}

	// Every generated input and intermediate belongs to this producer root.
	if b.ws.Root != "" {
		log.Debug("cleaning up build workspace", "path", b.ws.Root)
		if err := os.RemoveAll(b.ws.Root); err != nil && !os.IsNotExist(err) {
			log.Warn("failed to remove build workspace", "path", b.ws.Root, "err", err)
		}
		utils.RemoveDirIfEmpty(filepath.Dir(b.ws.Root))
		utils.RemoveDirIfEmpty(b.ws.BaseRoot)
	}

	// On failure remove only this build's own output. The installed image is
	// never touched: a build writes to preparedPath and installs by rename, so
	// whatever is at targetOverlayPath is a complete image someone may be using.
	if failed && b.tgt.Prepared != "" {
		if err := os.Remove(b.tgt.Prepared); err != nil && !os.IsNotExist(err) {
			log.Warn("failed to remove partial build output", "path", b.tgt.Prepared, "err", err)
		}
	}

	if willClean {
		log.Info("temporary files cleaned")
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
	if err := b.parseInputPrompts(); err != nil {
		return err
	}
	if err := b.collectInputAnswers(ctx); err != nil {
		return err
	}
	return b.resolveResourceSpec()
}

// parseDependencies reads #DEP: lines from the build script into the Spec.
// Skips parsing when the catalog already materialized the recipe, whose deps
// arrive expanded.
func (b *BuildObject) parseDependencies() error {
	if b.spec.Dependencies != nil {
		return nil
	}
	deps, err := utils.GetDependenciesFromScript(b.buildSource)
	if err != nil {
		return fmt.Errorf("failed to parse dependencies: %w", err)
	}
	b.spec.Dependencies = deps
	return nil
}

// parseInputPrompts reads #INPUT: lines from the build source and sets b.inputPrompts.
// Skips parsing when the catalog already materialized the recipe.
func (b *BuildObject) parseInputPrompts() error {
	if b.inputPrompts != nil || b.tempSource {
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
	b.inputPrompts = recipe.Inputs
	return nil
}

// collectInputAnswers asks for every #INPUT: the recipe declares, in order.
// Answers reach the recipe on stdin, not an env var: apptainer shell-evaluates
// env values, so a $ or backtick in a pasted URL would be mangled or run.
func (b *BuildObject) collectInputAnswers(ctx context.Context) error {
	b.inputAnswers = []string{}
	if len(b.inputPrompts) == 0 {
		return nil
	}

	// If --yes flag is set, automatically provide empty responses
	if utils.ShouldAnswerYes() {
		for range b.inputPrompts {
			b.inputAnswers = append(b.inputAnswers, "")
		}
		return nil
	}

	// Prompts require a TTY or piped stdin (e.g. scheduler job with embedded heredoc)
	if !utils.IsInteractiveShell() && !utils.IsStdinPiped() {
		return fmt.Errorf("recipe for %s requires input, but no TTY is available", b.spec.Image.Name)
	}

	log := logging.FromContext(ctx)
	for _, prompt := range b.inputPrompts {
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
		b.inputAnswers = append(b.inputAnswers, input)
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
		return fmt.Errorf("build script %s has unsupported scheduler directives; remove or fix them", b.buildSource)
	}

	// Resolve using the priority chain: buildDefaults → script → job resources.
	specs.Spec = EffectiveResourceSpec(specs)
	return nil
}

// appExt3ScratchExt returns the ext3 scratch-image extension for a build: ".img"
// for an app under build.app_tmp_overlay, "" for every other type whatever the config
// says. See the README's Workspace Strategy.
func appExt3ScratchExt(typ catalog.Type) string {
	if typ != catalog.TypeApp || !config.Global.Build.AppTmpOverlay {
		return ""
	}
	return ".img"
}

// NewBuildObject creates a BuildObject from a name/version string
// Format: "name/version" for conda/shell, "name" for def, "prefix/name/version" for ref
// All overlays are stored in imagesDir regardless of type
func NewBuildObject(ctx context.Context, nameVersion string, external bool, imagesDir string, update bool) (*BuildObject, error) {
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

	// Provisional type from the name shape. The recipe's entry is authoritative
	// and may override it via #TYPE:, but the workspace has to be sited before
	// anything is looked up — createConcreteType re-points it if the type moves.
	typ := catalog.DeriveType(normalized, "", false, "")

	// Use fast local storage for app builds; keep a stable path for data.
	// A definition build re-sites this in asDefinitionBuild once its type is known.
	tmpDir := tmpRootForType(typ)

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

	ws := workspaceFor(normalized, tmpDir, appExt3ScratchExt(typ), false)

	logging.FromContext(ctx).Debug("creating build object",
		"input", nameVersion, "nameVersion", normalized,
		"targetOverlay", targetOverlay, "tmpOverlay", ws.Overlay, "cntDir", ws.CntDir)

	base := &BuildObject{
		spec:            Spec{Image: ImageSpec{Name: normalized, Type: typ}},
		ws:              ws,
		tgt:             targetFor(targetOverlay),
		submitJob:       config.Global.SubmitJob,
		update:          update,
		condaChannelPkg: condaChannelPkg,
	}

	if external {
		// External builds don't need to resolve build source
		return createConcreteType(ctx, base, tmpDir)
	}

	// Already installed: skip resolving the concrete type. Not in update mode,
	// where the type has to be resolved to rebuild it.
	if !update && base.IsInstalled() {
		base.buildType = BuildTypeScript
		return base, nil
	}

	if err := base.clearStaleLock(ctx); err != nil {
		return nil, err
	}

	// Resolve build source and determine concrete type
	return createConcreteType(ctx, base, tmpDir)
}

// NewCondaObjectWithSource creates a Conda BuildObject packing one image from a
// custom buildSource: a YAML or spec file path, or a comma-separated package
// list. This is the -n flag.
func NewCondaObjectWithSource(nameVersion, buildSource, imagesDir string, update bool) (*BuildObject, error) {
	normalized := catalog.Normalize(nameVersion)

	// A conda environment is an app, so it always gets fast local scratch.
	tmpDir := tmpRootForType(catalog.TypeApp)

	// Make tmpDir absolute
	if absDir, err := filepath.Abs(tmpDir); err == nil {
		tmpDir = absDir
	}

	targetOverlay := filepath.Join(imagesDir, strings.ReplaceAll(normalized, "/", "--")+".sqf")
	if abs, err := filepath.Abs(targetOverlay); err == nil {
		targetOverlay = abs
	}

	ws := workspaceFor(normalized, tmpDir, appExt3ScratchExt(catalog.TypeApp), false)

	slog.Default().Debug("creating conda build object",
		"nameVersion", nameVersion, "buildSource", buildSource,
		"targetOverlay", targetOverlay, "tmpOverlay", ws.Overlay, "cntDir", ws.CntDir)

	base := &BuildObject{
		spec:        Spec{Image: ImageSpec{Name: normalized, Type: catalog.TypeApp}},
		ws:          ws,
		tgt:         targetFor(targetOverlay),
		buildSource: buildSource,
		submitJob:   config.Global.SubmitJob,
		update:      update,
	}

	if err := base.setupCondaFields(); err != nil {
		return nil, err
	}
	base.buildType = BuildTypeConda
	base.setCondaSpec()
	return base, nil
}

// captureLocalSourceSpec fills the Spec from a local script or definition the
// user supplied rather than the catalog, parsed with the same reader so a recipe
// means the same either way. A parse failure only costs descriptive metadata.
func (b *BuildObject) captureLocalSourceSpec(isDef bool) {
	// Drop anything a previous source contributed; the recipe below re-supplies it.
	b.spec.Image.Description = ""
	b.spec.Image.URL = ""
	b.spec.Image.License = ""
	b.spec.Image.Redistribute = nil
	b.spec.Image.Prefix = meta.Prefix(b.spec.Image.Name, b.spec.Image.Type)
	b.spec.Image.Env = nil

	data, err := os.ReadFile(b.buildSource)
	if err != nil {
		slog.Default().Debug("could not read build source for metadata", "path", b.buildSource, "err", err)
		return
	}
	file := SourceFile{Name: filepath.Base(b.buildSource), Data: data}
	if isDef {
		b.spec.Source = SourceSpec{Definition: &DefinitionSource{File: file}}
	} else {
		b.spec.Source = SourceSpec{Script: &ScriptSource{File: file, Prompts: b.inputPrompts}}
	}
	b.embedSource(SourceFile{Name: meta.RecipeFileName, Data: data})

	recipe, err := catalog.ParseRecipe(b.buildSource, bytes.NewReader(data))
	if err != nil {
		slog.Default().Debug("could not parse build source for metadata", "path", b.buildSource, "err", err)
		return
	}
	b.spec.Image.Description = recipe.Description
	b.spec.Image.URL = recipe.URL
	b.spec.Image.License = recipe.License
	b.spec.Image.Redistribute = recipe.Redistributable()
	b.spec.Image.Prefix = meta.Prefix(b.spec.Image.Name, b.spec.Image.Type)
	b.spec.Image.Env = envFromRecipe(recipe.Env)
	b.spec.Image.Arch = recipe.Arch
}

// FromExternalSource creates a BuildObject from an external build script or def
// file — the `-f <script>.sh` / `.def` path, which is the only way an external
// source is built. All overlays are stored in imagesDir regardless of type.
//
// The artifact name comes from the script's #TARGET: when it declares one, and
// only from the `-p` basename otherwise. The two say different things: `-p` is
// where the file goes, #TARGET: is what the payload is called, which fixes its
// /cnt/<name> prefix and, through key.Role, which dependencies count toward its
// equivalence. They are deliberately unrelated — a path artifact is addressed by
// its path, so its filename carries no naming claim.
func FromExternalSource(ctx context.Context, targetPrefix, source string, isApptainer bool, imagesDir string, update bool) (*BuildObject, error) {
	nameVersion := filepath.Base(targetPrefix)
	nameVersion = catalog.Normalize(nameVersion)

	// Determine build type from source file extension
	isDef := isApptainer || strings.HasSuffix(source, ".def")
	isShell := strings.HasSuffix(source, ".sh") || strings.HasSuffix(source, ".bash")
	externalType := "app"
	var target string
	var deps []string
	if isShell {
		parsedType, err := utils.GetTypeFromScript(source)
		if err != nil {
			return nil, fmt.Errorf("failed to parse external build type: %w", err)
		}
		externalType = parsedType

		if target, err = utils.GetTargetFromScript(source); err != nil {
			return nil, fmt.Errorf("failed to parse external build target: %w", err)
		}
		if deps, err = utils.GetDependenciesFromScript(source); err != nil {
			return nil, fmt.Errorf("failed to parse external build dependencies: %w", err)
		}
	}
	if target != "" {
		nameVersion = target
	}

	// A definition is a container root — base when it is named <>/base, os
	// otherwise; a shell build follows the declared or derived type.
	externalTyp := catalog.DeriveType(nameVersion, "", isDef, externalType)

	// The same rules a catalog recipe answers to: only data may depend on
	// anything, and an edge is a name/version rather than an overlay path.
	if err := catalog.ValidateDeps(nameVersion, externalTyp, deps); err != nil {
		return nil, err
	}
	// Checked after the type, so an app declaring #DEP: hears the more
	// fundamental refusal. Without #TARGET: the name is whatever the caller typed
	// after -p, which leaves the dependency's role — and so the artifact's
	// equivalence — decided by where the file was written. Silent history
	// classification is the trap; naming yourself is how an author opts into it.
	if target == "" && len(deps) > 0 {
		return nil, fmt.Errorf("%s declares #DEP: but no #TARGET:; an external build with dependencies must name itself", source)
	}

	targetDir := tmpRootForExternal(filepath.Dir(targetPrefix), externalTyp, isDef)
	if absDir, err := filepath.Abs(targetDir); err == nil {
		targetDir = absDir
	}
	// A definition gets no scratch image: apptainer writes its root as a sandbox
	// directory, which workspaceFor sites from isDef.
	ext := appExt3ScratchExt(externalTyp)
	if isDef {
		ext = ""
	}
	ws := workspaceFor(nameVersion, targetDir, ext, isDef)

	logging.FromContext(ctx).Debug("creating external build object",
		"nameVersion", nameVersion, "source", source,
		"targetPrefix", targetPrefix, "tmpOverlay", ws.Overlay, "cntDir", ws.CntDir)

	base := &BuildObject{
		spec:        Spec{Image: ImageSpec{Name: nameVersion, Type: externalTyp}},
		ws:          ws,
		tgt:         targetFor(targetPrefix + ".sqf"),
		buildSource: source,
		submitJob:   config.Global.SubmitJob,
		update:      update,
	}

	if err := base.clearStaleLock(ctx); err != nil {
		return nil, err
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
	base.captureLocalSourceSpec(isDef)

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
		base.spec.Image.Type = catalog.TypeApp
		base.buildType = BuildTypeConda
		base.setCondaSpec()
		return base, nil
	}

	if isContainer {
		base.asDefinitionBuild()
		return base, nil
	}

	// The entry may have declared a type the name shape does not imply, which
	// moves where the build works.
	base.retargetWorkspace()

	// It's a shell script (no extension or .sh/.bash)
	if err := base.parseScriptMetadata(ctx); err != nil {
		return nil, err
	}
	base.buildType = BuildTypeScript
	if base.spec.Source.BuildType() == "" {
		// A script given as a path, so resolution never opened a recipe.
		base.captureLocalSourceSpec(false)
	}
	return base, nil
}

// asDefinitionBuild sites the build as a definition — a sandbox directory on
// fast local scratch — and captures the Spec when resolution never opened a
// recipe. Call it after resolution: it reads the Spec's name.
func (b *BuildObject) asDefinitionBuild() {
	dir := tmpRootForDef()
	if abs, err := filepath.Abs(dir); err == nil {
		dir = abs
	}
	b.ws = workspaceFor(b.spec.Image.Name, dir, "", true)
	b.buildType = BuildTypeDef
	if b.spec.Source.BuildType() == "" {
		// A definition given as a path or a scheme:// URI, so resolution never
		// opened a recipe to take the Spec from.
		b.captureLocalSourceSpec(true)
	}
}

// resolveBuildSource resolves the module through the catalog and materializes its
// recipe. Returns (isConda, isContainer, error). A name no source provides is
// conda's, not a failure. Recipe text comes back already expanded.
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
	match, found, err := cat.Lookup(ctx, base.spec.Image.Name)
	if err != nil {
		return false, false, err
	}
	// After the lookup, which is what populates Err: an unreachable source is
	// skipped, and the next one — or conda — answers in its place.
	config.WarnUnreachableSources(ctx, cat)
	if !found {
		slog.Default().Debug("no recipe found, using conda", "name", base.spec.Image.Name)
		return true, false, nil
	}
	isContainer = strings.HasSuffix(match.Entry.Path, ".def")
	base.vars = match.Vars
	base.catalogSource = match.Source
	base.spec.Image.Type = match.Entry.Type

	// Descriptive metadata comes from the index entry, so it is available even
	// when the recipe itself is never opened (the installed-and-not-updating
	// case below returns early).
	base.spec.Image = ImageSpec{
		Name:        base.spec.Image.Name,
		Type:        match.Entry.Type,
		Description: match.Entry.Description,
		URL:         match.Entry.URL,
	}

	// Nothing to materialize for an overlay that already exists and is not
	// being updated — dependency walks reach installed nodes routinely.
	if !base.update && base.IsInstalled() {
		slog.Default().Debug("target already exists, skipping recipe fetch", "name", base.spec.Image.Name)
		return false, isContainer, nil
	}

	recipe, err := cat.Open(ctx, base.spec.Image.Name, base.vars)
	if err != nil {
		return false, false, fmt.Errorf("failed to read recipe for %s: %w", base.spec.Image.Name, err)
	}
	// Validation happens here rather than while indexing, so one bad recipe stops
	// its own build instead of taking a whole collection out of every listing.
	if err := recipe.Validate(); err != nil {
		return false, false, err
	}
	for _, warning := range recipe.Lint() {
		logging.FromContext(ctx).Warn("recipe lint", "detail", warning)
	}

	if isContainer {
		base.asDefinitionBuild()
	} else {
		base.retargetWorkspace()
	}
	path, err := writeRecipeFile(recipe, base.ws.Source)
	if err != nil {
		return false, false, err
	}
	base.buildSource = path
	base.tempSource = true
	base.inputPrompts = recipe.Inputs

	// The recipe is authoritative over the index entry, which can be stale.
	base.spec.Image.Type = recipe.Type
	base.spec.Image.Description = recipe.Description
	base.spec.Image.URL = recipe.URL
	base.spec.Image.License = recipe.License
	base.spec.Image.Redistribute = recipe.Redistributable()
	base.spec.Dependencies = recipe.Deps
	base.spec.Image.Prefix = meta.Prefix(base.spec.Image.Name, recipe.Type)
	base.spec.Image.Env = envFromRecipe(recipe.Env)
	base.spec.Image.Arch = recipe.Arch
	// recipe.Text is the template, tokens intact — what the image embeds and what
	// a rebuild starts from. The expansion went to the workspace file above and
	// is never recorded; the variant is recorded as its selected values.
	source := SourceSpec{Script: &ScriptSource{
		File:    SourceFile{Name: filepath.Base(path), Data: recipe.Text},
		Prompts: recipe.Inputs,
	}}
	if isContainer {
		source = SourceSpec{Definition: &DefinitionSource{
			File: SourceFile{Name: filepath.Base(path), Data: recipe.Text},
		}}
	}
	source.Placeholders = selectedPlaceholders(recipe)
	source.TargetTemplate = recipe.TargetTemplate
	source.RequiresInput = len(recipe.Inputs) > 0
	source.Collection = match.Source.Desc.Source
	base.spec.Source = source
	base.embedSource(SourceFile{Name: meta.RecipeFileName, Data: recipe.Text})
	slog.Default().Debug("materialized recipe", "path", path, "source", match.Source.Name)

	return false, isContainer, nil
}

// selectedPlaceholders reduces an expanded recipe's PH map to the one value each
// placeholder was resolved to. Expand leaves a single-element list per name, so
// anything else is a recipe that was never a template.
func selectedPlaceholders(recipe *catalog.Recipe) map[string]string {
	// Only an expanded template has selections. A recipe declaring #PH: without a
	// #TARGET: is not a template at all, and its PH still holds whole menus.
	if recipe.IsTemplate || recipe.TargetTemplate == "" || len(recipe.PH) == 0 {
		return nil
	}
	out := make(map[string]string, len(recipe.PH))
	for name, values := range recipe.PH {
		if len(values) == 1 {
			out[name] = values[0]
		}
	}
	if len(out) == 0 {
		return nil
	}
	return out
}

// writeRecipeFile writes the runnable recipe to a temp file the build executes.
// A template is written expanded; what the image embeds is the template itself.
func writeRecipeFile(recipe *catalog.Recipe, path string) (string, error) {
	if err := utils.MkdirAllShared(filepath.Dir(path)); err != nil {
		return "", fmt.Errorf("failed to create tmp directory: %w", err)
	}

	file, err := utils.CreateFileWritable(path)
	if err != nil {
		return "", fmt.Errorf("failed to create %s: %w", path, err)
	}
	defer file.Close()
	if _, err := file.Write(recipe.Script()); err != nil {
		return "", fmt.Errorf("failed to write %s: %w", path, err)
	}
	utils.ShareWithParentGroup(path)
	return path, nil
}
