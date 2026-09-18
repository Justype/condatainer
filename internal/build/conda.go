package build

import (
	"bytes"
	"context"
	"encoding/json"
	"fmt"
	"os"
	osexec "os/exec"
	"path/filepath"
	"strings"

	"github.com/Justype/condatainer/internal/artifact/key"
	"github.com/Justype/condatainer/internal/artifact/meta"
	"github.com/Justype/condatainer/internal/conda"

	"github.com/Justype/condatainer/internal/config"
	"github.com/Justype/condatainer/internal/image/freeze"
	"github.com/Justype/condatainer/internal/libexec"
	"github.com/Justype/condatainer/internal/logging"
	"github.com/Justype/condatainer/internal/runtime/exec"
	"github.com/Justype/condatainer/internal/toolpath"
	"github.com/Justype/condatainer/internal/utils"
)

// shellQuote renders a value as a single-quoted shell word.
func shellQuote(s string) string {
	return "'" + strings.ReplaceAll(s, "'", `'\''`) + "'"
}

// buildConda installs an environment with micromamba and packs it. The
// install -> stage -> pack sequence is the script backend's too, which is why
// both share packOutput.
func (b *BuildObject) buildConda(ctx context.Context) error {
	targetPath := b.tgt.Path
	log := logging.FromContext(ctx)

	if skip, err := checkShouldBuild(b); skip || err != nil {
		b.Cleanup(err != nil) //nolint:errcheck
		return err
	}

	// After the skip check, so an already-installed overlay never triggers a
	// base build it has no use for.
	if err := b.resolveBase(ctx); err != nil {
		return err
	}

	// The solve runs inside the base, so this follows resolveBase. It is the whole
	// cost of knowing the identity in advance, and it saves creating, packing and
	// discarding an environment that is already installed.
	if b.skipIfInstalled(ctx) {
		b.Cleanup(false) //nolint:errcheck
		return nil
	}

	if err := b.createBuildLock(); err != nil {
		b.Cleanup(true) //nolint:errcheck
		return err
	}
	defer b.removeBuildLock()
	preparedPath := b.tgt.Prepared

	log.Info("building overlay", "kind", "note", "overlay", filepath.Base(targetPath), "mode", buildModeLabel(b))

	if err := prepareBuildWorkspace(ctx, b); err != nil {
		b.Cleanup(true) //nolint:errcheck
		return err
	}

	if err := b.installConda(ctx); err != nil {
		b.Cleanup(true)
		return err
	}
	// After the recipe's container has run, so apptainer.ResolveBin has
	// already decided which binary that used and captureCommonBuildTools can
	// read it back rather than resolving a possibly different one.
	b.captureCommonBuildTools(ctx)
	b.captureMicromambaVersion(ctx)

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

	installed, err := b.publish(preparedPath, targetPath)
	if err != nil {
		b.Cleanup(true) //nolint:errcheck
		return err
	}

	log.Info("overlay ready", "kind", "success", "path", installed)
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
	return b.buildCreateCmd("/cnt/" + b.spec.Image.Name)
}

// micromambaCmd names the micromamba binary a generated build script should
// invoke: the self-provisioned one's absolute path, so this never depends on
// PATH order or a same-named tool elsewhere on PATH silently shadowing it.
// Errors when libexec has not been provisioned, rather than falling back to a
// bare "micromamba" that would only fail later as an unattributed in-container
// shell error.
func micromambaCmd() (string, error) {
	path, ok := libexec.MicromambaPath()
	if !ok {
		return "", libexec.NotInstalledError("micromamba")
	}
	return path, nil
}

// buildCreateCmd returns the micromamba create command for a target prefix, and
// any extra bind paths it needs. Handles three modes: YAML file,
// comma-separated packages, and single package.
//
// The prefix is a parameter so a dry-run solve can name a throwaway one. It does
// not change what is resolved — create solves for a prefix that does not exist
// yet either way — but it does decide what has to be mounted.
//
// --no-rc keeps the solve reproducible regardless of who runs the build: without
// it, the invoking user's own ~/.condarc (channels, channel_priority, ...) is
// visible inside the container — Apptainer binds $HOME by default — and would
// silently influence a build two different users expect to produce the same
// artifact.
func (b *BuildObject) buildCreateCmd(prefix string) (cmd string, extraBindPaths []string, err error) {
	var quietFlag string
	if utils.QuietMode {
		quietFlag = "-q"
	}
	channelFlags := buildChannelFlags()
	mmCmd, err := micromambaCmd()
	if err != nil {
		return "", nil, err
	}

	if utils.IsCondaFile(b.buildSource) {
		// Mode 3: env/spec file (-p prefix -f environment.yml or explicit .txt)
		absFilePath, err := filepath.Abs(b.buildSource)
		if err != nil {
			return "", nil, fmt.Errorf("failed to get absolute path for %s: %w", b.buildSource, err)
		}
		extraBindPaths = []string{filepath.Dir(absFilePath)}
		cmd = fmt.Sprintf(mmCmd+" create -r "+ScratchPath+" --no-rc %s -y %s -p %s -f %s",
			channelFlags, quietFlag, prefix, absFilePath)
	} else if b.buildSource != "" {
		// Mode 2: Multiple packages (-n name pkg1 pkg2 ...)
		packages := strings.Split(b.buildSource, ",")
		for i, pkg := range packages {
			packages[i] = strings.ReplaceAll(strings.TrimSpace(pkg), "/", "=")
		}
		cmd = fmt.Sprintf(mmCmd+" create -r "+ScratchPath+" --no-rc %s -y %s -p %s %s",
			channelFlags, quietFlag, prefix, strings.Join(packages, " "))
	} else {
		// Mode 1: Single package (name/version)
		cmd = fmt.Sprintf(mmCmd+" create -r "+ScratchPath+" --no-rc %s -y %s -p %s %s=%s",
			channelFlags, quietFlag, prefix, b.packageName, b.packageVersion)
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
// returns its stdout, directly on host — a plain host path in dir mode, or
// the same fuse2fs mount packFromScratchImage uses to read a scratch .img in
// ext3 mode (squashfs.go).
func (b *BuildObject) condaExport(ctx context.Context, args ...string) ([]byte, error) {
	mmCmd, err := micromambaCmd()
	if err != nil {
		return nil, err
	}

	if !b.ws.UsesImage() {
		prefix := filepath.Join(b.ws.CntDir, b.spec.Image.Name)
		cmdArgs := append([]string{"env", "export", "-p", prefix}, args...)
		var out bytes.Buffer
		cmd := osexec.CommandContext(ctx, mmCmd, cmdArgs...)
		cmd.Stdout = &out
		if err := cmd.Run(); err != nil {
			return nil, err
		}
		return out.Bytes(), nil
	}

	fuse2fsBin, err := toolpath.Resolve("fuse2fs")
	if err != nil {
		return nil, err
	}
	mnt := filepath.Join(b.ws.TmpDir, "export-mnt")
	if err := os.MkdirAll(mnt, 0o755); err != nil {
		return nil, fmt.Errorf("create export mountpoint: %w", err)
	}
	defer os.RemoveAll(mnt)

	prefix := filepath.Join(mnt, freeze.UpperDir, "cnt", b.spec.Image.Name)
	outFile := filepath.Join(b.ws.TmpDir, "export-out")
	defer os.Remove(outFile)
	script := fmt.Sprintf("%s env export -p %s %s > %s",
		shellQuote(mmCmd), shellQuote(prefix), strings.Join(args, " "), shellQuote(outFile))
	if err := freeze.MountedRun(ctx, fuse2fsBin, []string{"-o", "ro", b.ws.Overlay}, mnt, script, exec.IO{}); err != nil {
		return nil, err
	}
	return os.ReadFile(outFile)
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
// payload lives inside. installConda is its only caller; solveConda builds
// its own Options, since a dry-run solve mounts no payload at all.
//
// Unlike script.go/squashfs.go, this call binds neither getAllBaseDirs() nor
// container.BindPaths(), so the toolchain's tier needs an explicit bind here.
// Checks libexec.Dir(), not Ensure — a build must never trigger a first-time
// download — and errors rather than silently omitting the bind when nothing
// is provisioned, since micromambaCmd's absolute path would not resolve
// in-container without it.
func (b *BuildObject) condaExecOpts(bashScript string, extraBindPaths []string) (exec.Options, error) {
	bindPaths := extraBindPaths

	dir, ok := libexec.Dir()
	if !ok {
		return exec.Options{}, libexec.ErrNotProvisioned
	}
	bindPaths = append(bindPaths, dir)

	if !b.ws.UsesImage() {
		bindPaths = append(bindPaths,
			b.ws.TmpDir+":"+ScratchPath,
			b.ws.CntDir+":/cnt",
		)
		return exec.Options{
			BaseImage:      b.spec.Base,
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

// predictCondaIdentity derives the identity this Conda build will produce, by
// solving without installing.
//
// The solve is the only part that decides the answer: identity hashes the
// canonical explicit.txt, whose whole content is the resolved package URLs, and
// `micromamba create --dry-run --json` reports exactly those. Creating the
// environment, packing it and reading it back would give the same key for
// considerably more work.
//
// Both files are written through conda.ExplicitFrom and conda.EnvironmentFrom,
// which end in the same canonicalizers captureCondaExports uses. That shared
// ending is what keeps one artifact from having two identities.
func (b *BuildObject) predictCondaIdentity(ctx context.Context) (meta.KeyRef, error) {
	packages, err := b.solveConda(ctx)
	if err != nil {
		return meta.KeyRef{}, fmt.Errorf("%w: %v", ErrNoPrediction, err)
	}
	explicit, err := conda.ExplicitFrom(packages)
	if err != nil {
		return meta.KeyRef{}, fmt.Errorf("%w: %v", ErrNoPrediction, err)
	}
	var channels []string
	if b.spec.Source.Conda != nil {
		channels = b.spec.Source.Conda.Channels
	}
	environment, err := conda.EnvironmentFrom(packages, channels)
	if err != nil {
		return meta.KeyRef{}, fmt.Errorf("%w: %v", ErrNoPrediction, err)
	}

	// A manifest copy carrying the sources this build would embed. b is left
	// untouched: the real capture happens against the installed environment.
	manifest := b.Manifest()
	manifest.Keys = meta.Keys{}
	manifest.Source.Files = append(manifest.Source.Files,
		conda.ExplicitFileName, conda.EnvironmentFileName)
	sources := b.keySources()
	sources[conda.ExplicitFileName] = explicit
	sources[conda.EnvironmentFileName] = environment

	derived, err := key.Generate(manifest, sources)
	if err != nil {
		return meta.KeyRef{}, fmt.Errorf("%w: %v", ErrNoPrediction, err)
	}
	return derived.Identity.Ref, nil
}

// solveConda runs the build's own create command as a dry run and returns what
// it would install.
//
// It mounts none of the build workspace. A dry run writes nothing, so the
// payload directory it would install into need not exist — which matters,
// because this runs before the workspace is prepared, and preparing one for a
// build that is about to be skipped is the cost this whole path avoids. Only the
// producer's scratch root is bound, as micromamba's root prefix, so the repodata
// it downloads is reused by the build that follows.
//
// Nothing it prints reaches the user. A solve is a question asked on the way to
// a decision, and Apptainer's mount chatter or Micromamba's progress would read
// as a build that had started. Output is kept and reported only if it fails.
func (b *BuildObject) solveConda(ctx context.Context) ([]conda.Package, error) {
	// mkdir only — never prepareBuildWorkspace, which would create the scratch
	// image this deliberately does without.
	if err := ensureWorkspaceRoot(b); err != nil {
		return nil, err
	}
	cmd, extraBindPaths, err := b.buildCreateCmd(ScratchPath + "/solve")
	if err != nil {
		return nil, err
	}

	opts := exec.Options{
		BaseImage:      b.spec.Base,
		Overlays:       []string{},
		BindPaths:      append(extraBindPaths, b.ws.Root+":"+ScratchPath),
		EnvSettings:    []string{"TMPDIR=" + ScratchPath},
		Command:        []string{"/bin/bash", "-c", cmd + " --dry-run --json"},
		HidePrompt:     true,
		WritableImg:    false,
		ApptainerFlags: []string{"--writable-tmpfs"},
		PassThruStdin:  false,
	}

	// Every stream is the caller's own buffer, so exec.Run stays silent: it writes
	// to a terminal only through writers it is handed.
	var out, errOut bytes.Buffer
	if err := exec.Run(ctx, opts, exec.IO{Stdout: &out, Stderr: &errOut}); err != nil {
		if detail := strings.TrimSpace(errOut.String()); detail != "" {
			return nil, fmt.Errorf("%w: %s", err, detail)
		}
		return nil, err
	}

	// Micromamba prefixes the JSON with progress on some versions, so start at
	// the document rather than assuming the stream is clean.
	data := out.Bytes()
	if i := bytes.IndexByte(data, '{'); i > 0 {
		data = data[i:]
	}
	var report conda.DryRun
	if err := json.Unmarshal(data, &report); err != nil {
		return nil, fmt.Errorf("cannot read the solve: %w", err)
	}
	packages := report.Resolved()
	if len(packages) == 0 {
		return nil, fmt.Errorf("the solve resolved no packages")
	}
	return packages, nil
}
