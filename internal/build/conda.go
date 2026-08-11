package build

import (
	"bytes"
	"context"
	"fmt"
	"path/filepath"
	"strings"

	"github.com/Justype/condatainer/internal/conda"

	"github.com/Justype/condatainer/internal/config"
	"github.com/Justype/condatainer/internal/logging"
	"github.com/Justype/condatainer/internal/runtime/exec"
	"github.com/Justype/condatainer/internal/utils"
)

// buildConda installs an environment with micromamba and packs it. The
// install -> stage -> pack sequence is the script backend's too, which is why
// both share packOutput.
func (b *BuildObject) buildConda(ctx context.Context) error {
	targetPath := b.tgt.Path
	log := logging.FromContext(ctx)

	if skip, err := checkShouldBuild(b); skip || err != nil {
		return err
	}

	// After the skip check, so an already-installed overlay never triggers a
	// base build it has no use for.
	if err := b.resolveBase(ctx); err != nil {
		return err
	}

	if err := b.createBuildLock(); err != nil {
		return err
	}
	defer b.removeBuildLock()
	preparedPath := b.tgt.Prepared

	log.Info("building overlay", "overlay", filepath.Base(targetPath), "mode", buildModeLabel(b))

	if err := prepareBuildWorkspace(ctx, b); err != nil {
		return err
	}

	if err := b.installConda(ctx); err != nil {
		b.Cleanup(true)
		return err
	}

	b.captureCondaExports(ctx)
	b.describeCondaPackage()

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

	if err := atomicInstall(preparedPath, targetPath); err != nil {
		return err
	}

	log.Info("overlay ready", "kind", "success", "path", targetPath)
	b.Cleanup(false)
	return nil
}

// buildChannelFlags builds the micromamba -c flags from the configured channels.
func buildChannelFlags() string {
	var parts []string
	for _, ch := range config.Global.Build.Channels {
		parts = append(parts, "-c "+ch)
	}
	return strings.Join(parts, " ")
}

// buildInstallCmd returns the micromamba install command string and any extra bind paths.
// Handles three modes: YAML file, comma-separated packages, and single package.
func (b *BuildObject) buildInstallCmd() (cmd string, extraBindPaths []string, err error) {
	var quietFlag string
	if utils.QuietMode {
		quietFlag = "-q"
	}
	channelFlags := buildChannelFlags()

	if utils.IsCondaFile(b.buildSource) {
		// Mode 3: env/spec file (-p prefix -f environment.yml or explicit .txt)
		absFilePath, err := filepath.Abs(b.buildSource)
		if err != nil {
			return "", nil, fmt.Errorf("failed to get absolute path for %s: %w", b.buildSource, err)
		}
		extraBindPaths = []string{filepath.Dir(absFilePath)}
		cmd = fmt.Sprintf("micromamba create -r "+ScratchPath+" %s -y %s -p /cnt/%s -f %s",
			channelFlags, quietFlag, b.spec.Image.Name, absFilePath)
	} else if b.buildSource != "" {
		// Mode 2: Multiple packages (-n name pkg1 pkg2 ...)
		packages := strings.Split(b.buildSource, ",")
		for i, pkg := range packages {
			packages[i] = strings.ReplaceAll(strings.TrimSpace(pkg), "/", "=")
		}
		cmd = fmt.Sprintf("micromamba create -r "+ScratchPath+" %s -y %s -p /cnt/%s %s",
			channelFlags, quietFlag, b.spec.Image.Name, strings.Join(packages, " "))
	} else {
		// Mode 1: Single package (name/version)
		cmd = fmt.Sprintf("micromamba create -r "+ScratchPath+" %s -y %s -p /cnt/%s %s=%s",
			channelFlags, quietFlag, b.spec.Image.Name, b.packageName, b.packageVersion)
	}

	return cmd, extraBindPaths, nil
}

// captureCondaExports records what was installed, from the environment itself
// rather than from a second solve: explicit.txt pins exact package URLs, and
// environment.yml pins names and versions. They are the Conda app's identity and
// equivalence, so `sha256sum` on either reproduces a key by hand.
//
// A failure is a warning, not a build failure. An hour of solving is worth more
// than the records, and an image without them is the same unrecorded case as
// every image built before this format.
func (b *BuildObject) captureCondaExports(ctx context.Context) {
	log := logging.FromContext(ctx)

	raw, err := b.condaExport(ctx, "--explicit", "--no-md5")
	if err == nil {
		var explicit []byte
		if explicit, err = conda.CanonicalExplicit(raw); err == nil {
			b.embedSource(SourceFile{Name: conda.ExplicitFileName, Data: explicit})
		}
	}
	if err != nil {
		log.Warn("could not record the installed package set", "name", b.spec.Image.Name, "err", err)
	}

	raw, err = b.condaExport(ctx, "--no-builds")
	if err == nil {
		var channels []string
		if b.spec.Source.Conda != nil {
			channels = b.spec.Source.Conda.Channels
		}
		var environment []byte
		if environment, err = conda.CanonicalEnvironment(raw, channels); err == nil {
			b.embedSource(SourceFile{Name: conda.EnvironmentFileName, Data: environment})
		}
	}
	if err != nil {
		log.Warn("could not record the installed environment", "name", b.spec.Image.Name, "err", err)
	}
}

// condaExport runs `micromamba env export` against the installed prefix and
// returns its stdout. It reuses the install's container setup, so the export
// reads the environment that was just built whichever workspace mode is in use.
func (b *BuildObject) condaExport(ctx context.Context, args ...string) ([]byte, error) {
	opts, err := b.condaExecOpts(fmt.Sprintf("micromamba env export -p /cnt/%s %s",
		b.spec.Image.Name, strings.Join(args, " ")), nil)
	if err != nil {
		return nil, err
	}
	opts.PassThruStdin = false

	var out bytes.Buffer
	streams := exec.IOFromContext(ctx)
	streams.Stdin = nil
	streams.Stdout = &out
	if err := exec.Run(ctx, opts, streams); err != nil {
		return nil, err
	}
	return out.Bytes(), nil
}

// condaInstallExecOpts constructs exec.Options for the micromamba run. Installs
// only — the run never sees the output path or binds the images directory.
func (b *BuildObject) condaInstallExecOpts() (exec.Options, error) {
	installCmd, extraBindPaths, err := b.buildInstallCmd()
	if err != nil {
		return exec.Options{}, err
	}

	var echoPrefix string
	if utils.QuietMode {
		echoPrefix = ": #"
	} else {
		echoPrefix = "echo"
	}

	bashScript := fmt.Sprintf(`
trap 'exit 130' INT TERM
set -e

mkdir -p $TMPDIR
%[1]s "Creating conda environment in image..."
%[2]s

if [ -z "$(ls -A /cnt 2>/dev/null)" ]; then
    %[1]s "Conda environment is empty, nothing to pack."
    exit 1
fi
`, echoPrefix, installCmd)

	return b.condaExecOpts(bashScript, extraBindPaths)
}

// condaExecOpts sites a micromamba run against this build's payload, whichever
// workspace mode is in use: bound host directories, or the scratch image the
// payload lives inside. Install and export share it so an export always reads
// the environment the install just wrote.
func (b *BuildObject) condaExecOpts(bashScript string, extraBindPaths []string) (exec.Options, error) {
	bindPaths := extraBindPaths

	if !b.ws.UsesImage() {
		bindPaths = append(bindPaths,
			b.ws.TmpDir+":"+ScratchPath,
			b.ws.CntDir+":/cnt",
		)
		return exec.Options{
			BaseImage:      b.spec.Base,
			ApptainerBin:   config.Global.ApptainerBin,
			Overlays:       []string{},
			BindPaths:      bindPaths,
			EnvSettings:    []string{"TMPDIR=" + ScratchPath},
			Command:        []string{"/bin/bash", "-c", bashScript},
			HidePrompt:     true,
			WritableImg:    false,
			ApptainerFlags: []string{"--writable-tmpfs"},
			PassThruStdin:  true,
		}, nil
	}
	return exec.Options{
		BaseImage:     b.spec.Base,
		ApptainerBin:  config.Global.ApptainerBin,
		Overlays:      []string{b.ws.Overlay},
		BindPaths:     bindPaths,
		EnvSettings:   []string{"TMPDIR=" + ScratchPath},
		Command:       []string{"/bin/bash", "-c", bashScript},
		HidePrompt:    true,
		WritableImg:   true,
		PassThruStdin: true,
	}, nil
}

// installConda populates the payload with micromamba. The caller is responsible
// for calling Cleanup(true) if an error is returned.
func (b *BuildObject) installConda(ctx context.Context) error {
	opts, err := b.condaInstallExecOpts()
	if err != nil {
		return err
	}

	logging.FromContext(ctx).Debug("installing conda environment",
		"name", b.spec.Image.Name, "overlays", opts.Overlays, "bindPaths", opts.BindPaths)

	done := watchContext(ctx, "conda install")
	defer close(done)

	if err := exec.Run(ctx, opts, exec.IOFromContext(ctx)); err != nil {
		if isCancelledByUser(err) {
			return ErrBuildCancelled
		}
		return fmt.Errorf("failed to build conda package %s: %w", b.spec.Image.Name, err)
	}

	return nil
}
