package build

import (
	"context"
	"fmt"
	"os"
	"path/filepath"
	"strings"

	"log/slog"

	"github.com/Justype/condatainer/catalog"
	"github.com/Justype/condatainer/internal/config"
	"github.com/Justype/condatainer/internal/container"
	execpkg "github.com/Justype/condatainer/internal/exec"
	"github.com/Justype/condatainer/internal/logging"
	"github.com/Justype/condatainer/internal/scheduler"
	"github.com/Justype/condatainer/internal/utils"
)

// getAllBaseDirs returns all configured base directories that actually exist
func getAllBaseDirs() []string {
	var dirs []string

	addIfExists := func(dir string) {
		if dir == "" {
			return
		}
		if _, err := os.Stat(dir); err == nil {
			dirs = append(dirs, dir)
		}
	}

	for _, entry := range config.GetExtraImageDirs() {
		path, _ := config.ParseDirEntry(entry)
		addIfExists(path)
	}
	addIfExists(config.GetExtraRootDir())
	addIfExists(config.GetRootDir())
	addIfExists(config.GetScratchDataDir())
	addIfExists(config.GetUserDataDir())

	return dirs
}

// buildScript implements the shell script build workflow on BuildObject.
// Workflow:
//  1. Check if overlay already exists (skip if yes)
//  2. Create temporary ext3 overlay or host dirs
//  3. Check and build missing dependencies if buildDeps=true
//  4. Run shell script inside container with overlays
//  5. Create SquashFS from the output directory
//  6. Extract and save ENV variables if present
//  7. Set permissions, atomic install, and cleanup
func (b *BuildObject) buildScript(ctx context.Context, buildDeps bool) error {
	targetPath, finalPath := buildOverlayPaths(b)
	log := logging.FromContext(ctx)

	if skip, err := checkShouldBuild(b); skip || err != nil {
		return err
	}

	if err := b.createBuildLock(); err != nil {
		return err
	}
	defer b.removeBuildLock()

	log.Info("building overlay", "overlay", filepath.Base(targetPath), "mode", buildModeLabel(b))

	if err := prepareBuildWorkspace(ctx, b, config.Global.Build.UseTmpOverlay); err != nil {
		return err
	}

	log.Info("populating overlay", "overlay", filepath.Base(targetPath), "source", b.buildSource)

	if err := b.buildDependencies(ctx, buildDeps); err != nil {
		return err
	}

	if err := b.runBuildScript(ctx); err != nil {
		b.Cleanup(true)
		return err
	}

	if err := b.packOutput(ctx, finalPath); err != nil {
		return err
	}

	utils.ShareWithParentGroup(finalPath)

	if err := atomicInstall(finalPath, targetPath, b.update); err != nil {
		return err
	}

	log.Info("overlay ready", "kind", "success", "path", targetPath)
	b.Cleanup(false)
	return nil
}

// buildDependencies checks for missing dependencies and optionally builds them.
func (b *BuildObject) buildDependencies(ctx context.Context, buildDeps bool) error {
	missingDeps, err := b.GetMissingDependencies()
	if err != nil {
		return fmt.Errorf("failed to check dependencies: %w", err)
	}
	if len(missingDeps) == 0 {
		return nil
	}

	depList := strings.Join(missingDeps, ", ")

	if !buildDeps {
		logging.FromContext(ctx).Error("missing dependencies", "overlay", filepath.Base(b.targetOverlayPath), "deps", depList)
		return fmt.Errorf("missing dependencies for %s: %s. Please install them first", b.nameVersion, depList)
	}

	logging.FromContext(ctx).Info("building missing dependencies", "overlay", filepath.Base(b.targetOverlayPath), "deps", depList)

	writableImagesDir, err := config.GetWritableImagesDir()
	if err != nil {
		return fmt.Errorf("no writable images directory found: %w", err)
	}

	for _, dep := range missingDeps {
		// Strip any version constraint (e.g. "samtools/1.21>=1.16" → "samtools/1.21")
		// before passing to NewBuildObject, which expects a plain name/version.
		preferredNV := dep
		if parsed, err := catalog.ParseDep(dep); err == nil {
			preferredNV = parsed.NameVersion()
		}
		depObj, err := NewBuildObject(ctx, preferredNV, false, writableImagesDir, config.GetWritableTmpDir(), false)
		if err != nil {
			return fmt.Errorf("failed to create build object for dependency %s: %w", preferredNV, err)
		}
		if err := depObj.Build(ctx, false); err != nil {
			return fmt.Errorf("failed to build dependency %s: %w", preferredNV, err)
		}
	}

	logging.FromContext(ctx).Info("all dependencies built", "kind", "success", "overlay", filepath.Base(b.targetOverlayPath))
	return nil
}

// hostPayload reports whether the payload is written to a host directory bound
// into the container rather than into the tmp ext3 overlay. Dir mode always is;
// with a tmp overlay only data is, since a 20GB ext3 cannot hold a genome index.
func hostPayload(b *BuildObject) bool {
	return !config.Global.Build.UseTmpOverlay || b.kind == catalog.KindData
}

// buildExecOpts constructs the exec.Options for running the build script inside
// the container, plus the IO carrying any #INPUT: answers on stdin.
func (b *BuildObject) buildExecOpts() (execpkg.Options, execpkg.IO, error) {
	// Same splitter the catalog resolves with, so the two cannot disagree on
	// where the name ends and the version begins (data names carry slashes).
	dep, err := catalog.ParseDep(b.nameVersion)
	if err != nil {
		return execpkg.Options{}, execpkg.IO{}, err
	}

	effRS := buildEffectiveResourceSpec(b.scriptSpecs)
	buildRS := &scheduler.ResourceSpec{
		Nodes:        1,
		TasksPerNode: 1,
		CpusPerTask:  b.effectiveNcpus(),
		MemPerCpuMB:  effRS.MemPerCpuMB,
		MemPerNodeMB: effRS.MemPerNodeMB,
	}

	// Every build writes to /cnt/<name>/<version> — the path the artifact will be
	// mounted at, so anything the payload bakes in stays valid at run time.
	prefix := "/cnt/" + b.nameVersion

	kind := b.kind
	if kind == "" {
		kind = catalog.KindApp
	}

	// NCPUS/MEM/... come from the scheduler, already normalized across SLURM,
	// PBS and LSF, so this contract does not mint CNT_ spellings beside them.
	envSettings := append(
		scheduler.ResourceEnvVars(buildRS),
		"CNT_NAME="+dep.Name,
		"CNT_VERSION="+dep.Version,
		"CNT_KIND="+string(kind),
		"CNT_PREFIX="+prefix,
		"CNT_TMP=/ext3/tmp",
		"TMPDIR=/ext3/tmp",
		"IN_CONDATAINER=1",
	)
	bashScript := fmt.Sprintf(`
trap 'exit 130' INT TERM

mkdir -p $TMPDIR
bash -euo pipefail %s
if [ $? -ne 0 ]; then
    echo "Build script %s failed."
    exit 1
fi
# Embed the build script next to the payload for provenance (skipped if the
# script populated nothing — packOutput reports that as a build failure).
if [ -d "$CNT_PREFIX" ]; then
    cp %s "$CNT_PREFIX/%s" 2>/dev/null || true
fi
`, b.buildSource, b.buildSource, b.buildSource, utils.BuildScriptName)

	overlayArgs, err := GetOverlayArgsFromDependencies(b.dependencies)
	if err != nil {
		slog.Default().Warn("failed to get overlay args from dependencies", "err", err)
	}

	var overlays []string
	if config.Global.Build.UseTmpOverlay {
		overlays = []string{b.tmpOverlayPath}
	}
	for i := 0; i < len(overlayArgs); i += 2 {
		if overlayArgs[i] == "--overlay" && i+1 < len(overlayArgs) {
			overlays = append(overlays, overlayArgs[i+1])
		}
	}

	bindDirs := container.DeduplicateBindPaths(getAllBaseDirs())
	if hostPayload(b) {
		// Bind the leaf, not /cnt: covering /cnt would hide every dependency
		// overlay mounted beside it.
		payloadDir := filepath.Join(b.cntDirPath, b.nameVersion)
		if err := utils.MkdirAllShared(payloadDir); err != nil {
			return execpkg.Options{}, execpkg.IO{}, fmt.Errorf("failed to create payload dir %s: %w", payloadDir, err)
		}
		bindDirs = append(bindDirs, payloadDir+":"+prefix)
	}
	if !config.Global.Build.UseTmpOverlay {
		bindDirs = append(bindDirs, getBuildTmpDir(b)+":/ext3/tmp")
	}

	opts := execpkg.Options{
		BaseImage:    config.GetBaseImage(),
		ApptainerBin: config.Global.ApptainerBin,
		Overlays:     overlays,
		BindPaths:    bindDirs,
		EnvSettings:  envSettings,
		Command:      []string{"/bin/bash", "-c", bashScript},
		HidePrompt:   true,
		WritableImg:  config.Global.Build.UseTmpOverlay,
	}
	if !config.Global.Build.UseTmpOverlay {
		opts.ApptainerFlags = []string{"--writable-tmpfs"}
	}

	// #INPUT: answers go in on stdin, one line each, for the recipe to `read`.
	// Not an env var: apptainer shell-evaluates env values, so a $ or backtick
	// in a pasted URL would be mangled or executed.
	opts.PassThruStdin = true
	var ioStreams execpkg.IO
	if len(b.interactiveInputs) > 0 {
		ioStreams.Stdin = strings.NewReader(strings.Join(b.interactiveInputs, "\n") + "\n")
	}

	return opts, ioStreams, nil
}

// runBuildScript sets up a context watcher and runs the build script.
func (b *BuildObject) runBuildScript(ctx context.Context) error {
	opts, ioStreams, err := b.buildExecOpts()
	if err != nil {
		return err
	}

	logging.FromContext(ctx).Debug("running build script",
		"name", b.nameVersion, "overlays", opts.Overlays, "bindPaths", opts.BindPaths, "passThruStdin", opts.PassThruStdin)

	done := watchContext(ctx, "build script")
	defer close(done)

	if ioStreams.IsZero() {
		ioStreams = execpkg.IOFromContext(ctx)
	}
	if err := execpkg.Run(ctx, opts, ioStreams); err != nil {
		if isCancelledByUser(err) {
			return ErrBuildCancelled
		}
		return fmt.Errorf("build script %s failed: %w", b.buildSource, err)
	}

	return nil
}

// packOutput squashes the payload into the target overlay.
//
// A host payload is packed from its build directory, whose basename is cnt, so
// mksquashfs -keep-as-directory yields exactly cnt/<name>/<version>/... and
// nothing else — the build's tmp/ is a sibling and never enters the archive.
func (b *BuildObject) packOutput(ctx context.Context, finalPath string) error {
	log := logging.FromContext(ctx)
	isData := b.kind == catalog.KindData

	if hostPayload(b) {
		payloadDir := filepath.Join(b.cntDirPath, b.nameVersion)
		if entries, err := os.ReadDir(payloadDir); err != nil || len(entries) == 0 {
			b.Cleanup(true)
			return fmt.Errorf("build script did not create any files in %s", payloadDir)
		}
		log.Info("creating SquashFS", "source", b.cntDirPath, "overlay", filepath.Base(b.targetOverlayPath))
		if err := createSquashfs(ctx, b, isData, b.cntDirPath, finalPath); err != nil {
			b.Cleanup(true)
			return err
		}
		return nil
	}

	log.Info("preparing SquashFS from /cnt", "overlay", filepath.Base(b.targetOverlayPath))
	if err := createSquashfs(ctx, b, isData, "/cnt", finalPath); err != nil {
		b.Cleanup(true)
		return err
	}
	return nil
}

// saveCondaEnvFile fetches a package description from anaconda.org and writes the .env file.
// Skipped for multi-package / YAML builds where packageName is a user-chosen env name, not a real package.
func (b *BuildObject) saveCondaEnvFile(targetPath string) {
	if b.packageName == "" || b.buildSource != "" {
		return
	}
	whatis := utils.FetchCondaSummary(b.packageName, config.Global.Build.Channels)
	if err := SaveEnvFile(targetPath, map[string]EnvEntry{}, b.nameVersion, whatis); err != nil {
		slog.Default().Warn("failed to save ENV file", "err", err)
	}
}
