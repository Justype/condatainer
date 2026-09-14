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
	"github.com/Justype/condatainer/internal/logging"
	"github.com/Justype/condatainer/internal/runtime/container"
	execpkg "github.com/Justype/condatainer/internal/runtime/exec"
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

	addIfExists(config.GetExtraRootDir())
	addIfExists(config.GetRootDir())
	addIfExists(config.GetScratchDataDir())
	addIfExists(config.GetUserDataDir())

	return dirs
}

// buildScript runs a recipe as a shell script inside the base image, with its
// dependencies mounted, and packs what it wrote. buildDeps also builds any
// dependency that is missing.
func (b *BuildObject) buildScript(ctx context.Context, buildDeps bool) error {
	targetPath := b.tgt.Path
	log := logging.FromContext(ctx)

	if skip, err := checkShouldBuild(b); skip || err != nil {
		b.Cleanup(err != nil) //nolint:errcheck
		return err
	}

	if err := b.createBuildLock(); err != nil {
		b.Cleanup(true) //nolint:errcheck
		return err
	}
	defer b.removeBuildLock()
	preparedPath := b.tgt.Prepared

	// Data equivalence includes dependency equivalence, so dependencies must be
	// present before deciding whether a published artifact can substitute. This
	// still happens before creating the target's workspace or running its recipe.
	if err := b.buildDependencies(ctx, buildDeps); err != nil {
		b.Cleanup(true) //nolint:errcheck
		return err
	}
	if pulled, err := b.tryPrebuilt(ctx); err != nil || pulled {
		b.Cleanup(err != nil) //nolint:errcheck
		return err
	}

	// Dependencies are resolved by now, so the identity this build would produce
	// is fully determined — and under --store that, not the name, decides whether
	// there is anything to do.
	if b.skipIfInstalled(ctx) {
		b.Cleanup(false) //nolint:errcheck
		return nil
	}

	// A pulled artifact does not need the local build base. Resolve it only once
	// registry acquisition has declined and recipe execution will actually run.
	if err := b.resolveBase(ctx); err != nil {
		b.Cleanup(true) //nolint:errcheck
		return err
	}

	log.Info("building overlay", "kind", "note", "overlay", filepath.Base(targetPath), "mode", buildModeLabel(b))

	if err := prepareBuildWorkspace(ctx, b); err != nil {
		b.Cleanup(true) //nolint:errcheck
		return err
	}

	log.Info("populating overlay", "overlay", filepath.Base(targetPath), "source", b.buildSource)

	if err := b.runBuildScript(ctx); err != nil {
		b.Cleanup(true)
		return err
	}
	// After the recipe's container has run, so apptainer.ResolveBin has
	// already decided which binary that used and captureCommonBuildTools can
	// read it back rather than resolving a possibly different one.
	b.captureCommonBuildTools(ctx)

	// After the build, so every dependency it needed is installed and can be
	// read for the identity its records pin.
	if err := b.deriveKeys(ctx); err != nil {
		b.Cleanup(true)
		return err
	}

	metaDir, err := stageMetadata(ctx, b)
	if err != nil {
		b.Cleanup(true)
		return err
	}

	if err := b.packOutput(ctx, metaDir, preparedPath); err != nil {
		return err
	}

	utils.ShareWithParentGroup(preparedPath)

	installed, err := b.publish(preparedPath, targetPath)
	if err != nil {
		b.Cleanup(true) //nolint:errcheck
		return err
	}

	log.Info("overlay ready", "kind", "success", "path", installed)
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

	// A locked rebuild is handed exact paths, so a miss is a vanished file, not
	// something to resolve by name. Building it here would consult the catalog,
	// which is what a locked rebuild refuses to do.
	if b.locked {
		return fmt.Errorf("dependency images for %s are no longer present: %s", b.spec.Image.Name, depList)
	}

	if !buildDeps {
		logging.FromContext(ctx).Error("missing dependencies", "overlay", filepath.Base(b.tgt.Path), "deps", depList)
		return fmt.Errorf("missing dependencies for %s: %s. Please install them first", b.spec.Image.Name, depList)
	}

	logging.FromContext(ctx).Info("building missing dependencies", "kind", "note",
		"overlay", filepath.Base(b.tgt.Path), "deps", depList)

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
		depObj, err := NewBuildObject(ctx, preferredNV, false, writableImagesDir, false)
		if err != nil {
			return fmt.Errorf("failed to create build object for dependency %s: %w", preferredNV, err)
		}
		if err := depObj.Build(ctx, false); err != nil {
			return fmt.Errorf("failed to build dependency %s: %w", preferredNV, err)
		}
	}

	logging.FromContext(ctx).Info("all dependencies built", "kind", "success", "overlay", filepath.Base(b.tgt.Path))
	return nil
}

// buildExecOpts constructs the exec.Options for running the build script inside
// the container, plus the IO carrying any #INPUT: answers on stdin.
func (b *BuildObject) buildExecOpts() (execpkg.Options, execpkg.IO, error) {
	// The recipe's whole environment is a projection of Spec and Options — see
	// BuildEnv. A backend does not assemble it, so every launch path agrees.
	envSettings := buildEnv(b.spec, Options{Update: b.update, ScriptSpecs: b.scriptSpecs})

	// The payload's install prefix, as the recipe sees it via $CNT_PREFIX.
	prefix := b.spec.Image.Prefix
	if prefix == "" {
		prefix = "/cnt/" + b.spec.Image.Name
	}

	bashScript := fmt.Sprintf(`
trap 'exit 130' INT TERM

mkdir -p $TMPDIR
bash -euo pipefail %s
if [ $? -ne 0 ]; then
    echo "Build script %s failed."
    exit 1
fi
`, b.buildSource, b.buildSource)

	depOverlays, err := dependencyOverlays(b.spec.Dependencies)
	if err != nil {
		slog.Default().Warn("failed to resolve dependency overlays", "err", err)
	}

	var overlays []string
	if b.ws.UsesImage() {
		overlays = []string{b.ws.Overlay}
	}
	overlays = append(overlays, depOverlays...)

	bindDirs := container.DeduplicateBindPaths(getAllBaseDirs())
	if b.ws.HostPayload() {
		// Bind the leaf, not /cnt: covering /cnt would hide every dependency
		// overlay mounted beside it. prepareBuildWorkspace created it.
		bindDirs = append(bindDirs, filepath.Join(b.ws.CntDir, b.spec.Image.Name)+":"+prefix)
	}
	if !b.ws.UsesImage() {
		bindDirs = append(bindDirs, b.ws.TmpDir+":"+ScratchPath)
	}

	opts := execpkg.Options{
		BaseImage:   b.spec.Base,
		Overlays:    overlays,
		BindPaths:   bindDirs,
		EnvSettings: envSettings,
		Command:     []string{"/bin/bash", "-c", bashScript},
		HidePrompt:  true,
		WritableImg: b.ws.UsesImage(),
	}
	if !b.ws.UsesImage() {
		opts.ApptainerFlags = []string{"--writable-tmpfs"}
	}

	// #INPUT: answers go in on stdin, one line each, for the recipe to `read`.
	// Not an env var: apptainer shell-evaluates env values, so a $ or backtick
	// in a pasted URL would be mangled or executed.
	opts.PassThruStdin = true
	var ioStreams execpkg.IO
	if len(b.inputAnswers) > 0 {
		ioStreams.Stdin = strings.NewReader(strings.Join(b.inputAnswers, "\n") + "\n")
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
		"name", b.spec.Image.Name, "overlays", opts.Overlays, "bindPaths", opts.BindPaths, "passThruStdin", opts.PassThruStdin)

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

// packOutput squashes the payload and the staged metadata into the target image.
// A host payload is packed from its build directory, whose basename is cnt.
// metaDir is a second archive root and may be empty.
func (b *BuildObject) packOutput(ctx context.Context, metaDir, preparedPath string) error {
	log := logging.FromContext(ctx)
	isData := b.spec.Image.Type == catalog.TypeData

	if b.ws.HostPayload() {
		payloadDir := filepath.Join(b.ws.CntDir, b.spec.Image.Name)
		if entries, err := os.ReadDir(payloadDir); err != nil || len(entries) == 0 {
			b.Cleanup(true)
			return fmt.Errorf("build produced no files in %s", payloadDir)
		}
		log.Info("creating SquashFS", "source", b.ws.CntDir, "overlay", filepath.Base(b.tgt.Path))
		if err := createSquashfs(ctx, b, isData, b.ws.CntDir, metaDir, preparedPath); err != nil {
			b.Cleanup(true)
			return err
		}
		return nil
	}

	log.Info("preparing SquashFS from /cnt", "overlay", filepath.Base(b.tgt.Path))
	if err := createSquashfs(ctx, b, isData, "/cnt", metaDir, preparedPath); err != nil {
		b.Cleanup(true)
		return err
	}
	return nil
}

// describeCondaPackage fills in the description from anaconda.org, so a Conda
// image says what it is without a recipe to take a #DESC: from. Skipped for
// multi-package and YAML builds; a failed lookup leaves the description empty.
func (b *BuildObject) describeCondaPackage() {
	if b.packageName == "" || b.buildSource != "" || b.spec.Image.Description != "" {
		return
	}
	b.spec.Image.Description = utils.FetchCondaSummary(b.packageName, config.Global.Build.Channels)
}
