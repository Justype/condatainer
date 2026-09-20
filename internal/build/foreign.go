package build

import (
	"context"
	"encoding/json"
	"errors"
	"fmt"
	"io/fs"
	"os"
	"path/filepath"
	"runtime"
	"strconv"

	"github.com/Justype/condatainer/catalog"
	"github.com/Justype/condatainer/internal/artifact/meta"
	"github.com/Justype/condatainer/internal/config"
	"github.com/Justype/condatainer/internal/image/freeze"
	"github.com/Justype/condatainer/internal/image/sif"
	"github.com/Justype/condatainer/internal/image/tool"
	"github.com/Justype/condatainer/internal/logging"
	execpkg "github.com/Justype/condatainer/internal/runtime/exec"
	"github.com/Justype/condatainer/internal/toolpath"
	"github.com/Justype/condatainer/internal/utils"
)

// foreignRoot is a container root condatainer did not build: a .sif or an
// Apptainer sandbox directory, handed to -f already built. buildForeign packs
// it into a real .sqf instead of running apptainer.Build.
type foreignRoot struct {
	path    string // the .sif file or sandbox directory passed to -f
	sandbox bool   // true when path is already a plain directory
	boot    bootstrap
	digest  string // labels.json's base-image digest, or meta.Unrecorded
}

// FromForeignRoot creates a BuildObject that packs an already-built .sif or
// Apptainer sandbox directory into a real .sqf, carrying identity and
// equivalence derived from the root's own embedded build record rather than
// left absent. See this package's README, Importing a foreign root, for why
// this is a separate constructor from FromExternalSource rather than a third
// case of it.
func FromForeignRoot(ctx context.Context, targetPrefix, source, imagesDir string, update bool) (*BuildObject, error) {
	sandbox := utils.IsSandboxDir(source)
	if !sandbox && !utils.IsSif(source) {
		return nil, fmt.Errorf("%s is neither a .sif nor an Apptainer sandbox directory", source)
	}

	boot, digest, err := readForeignBootstrap(source, sandbox)
	if err != nil {
		return nil, fmt.Errorf("%s: %w", source, err)
	}
	if err := allowedForeignBootstrap(boot.Agent); err != nil {
		return nil, fmt.Errorf("%s: %w", source, err)
	}

	if !sandbox {
		if arch, err := sif.Arch(source); err == nil && arch != runtime.GOARCH {
			return nil, fmt.Errorf("%s was built for %s, this host is %s: it would pack but never run here",
				source, arch, runtime.GOARCH)
		}
	}

	nameVersion := catalog.Normalize(filepath.Base(targetPrefix))
	externalTyp := catalog.DeriveType(nameVersion, "", true, "")

	targetDir := tmpRootForExternal(filepath.Dir(targetPrefix), externalTyp, true)
	if absDir, err := filepath.Abs(targetDir); err == nil {
		targetDir = absDir
	}
	// isDef=false here, deliberately unlike targetDir's selection above: this
	// build never asks Apptainer to write a sandbox, so ws.Sandbox must stay
	// empty — packSources reads straight from the foreign root instead.
	ws := workspaceFor(nameVersion, targetDir, false)

	logging.FromContext(ctx).Debug("creating foreign-root build object",
		"nameVersion", nameVersion, "source", source, "sandbox", sandbox, "targetPrefix", targetPrefix)

	base := &BuildObject{
		spec:        Spec{Image: ImageSpec{Name: nameVersion, Type: externalTyp}},
		ws:          ws,
		tgt:         targetFor(targetPrefix + ".sqf"),
		buildSource: source,
		submitJob:   config.Global.SubmitJob,
		update:      update,
		buildType:   BuildTypeDef,
		foreignRoot: &foreignRoot{path: source, sandbox: sandbox, boot: boot, digest: digest},
	}
	if err := base.clearStaleLock(ctx); err != nil {
		return nil, err
	}
	return base, nil
}

// readForeignBootstrap extracts the Bootstrap/From directives and, when
// present, the upstream digest Apptainer recorded, from a foreign root's own
// /.singularity.d metadata — never by resolving the reference again, which
// would follow a tag that has moved since this root was built.
func readForeignBootstrap(path string, sandbox bool) (bootstrap, string, error) {
	read := sif.ReadFile
	if sandbox {
		read = readSandboxFile
	}

	def, err := read(path, ".singularity.d/Singularity")
	if err != nil {
		return bootstrap{}, "", fmt.Errorf("no build record at .singularity.d/Singularity: %w", err)
	}
	boot := parseBootstrap(def)

	digest := meta.Unrecorded
	if raw, err := read(path, ".singularity.d/labels.json"); err == nil {
		digest = digestFromLabels(raw)
	}
	return boot, digest, nil
}

// readSandboxFile reads one file out of an unpacked container root, the same
// shape sif.ReadFile has, so readForeignBootstrap can use either through one
// variable.
func readSandboxFile(root, inner string) ([]byte, error) {
	data, err := os.ReadFile(filepath.Join(root, inner))
	if errors.Is(err, fs.ErrNotExist) {
		return nil, fmt.Errorf("%w: %s", tool.ErrFileNotFound, inner)
	}
	return data, err
}

// digestFromLabels reads the base-image digest Apptainer stamped into
// labels.json at build time. Absent on many real .sif files — an older
// Apptainer/Singularity version never wrote it — so a missing or unparsable
// label is not an error, just no digest.
func digestFromLabels(raw []byte) string {
	var labels map[string]string
	if err := json.Unmarshal(raw, &labels); err != nil {
		return meta.Unrecorded
	}
	if d := labels["org.opencontainers.image.base.digest"]; d != "" {
		return d
	}
	return meta.Unrecorded
}

// allowedForeignBootstrap refuses every bootstrap agent that cannot be
// rebuilt from a fresh machine — shub is dead, localimage and scratch name
// nothing reproducible elsewhere, and yum/zypper/debootstrap have no single
// pinnable reference the identity model represents.
func allowedForeignBootstrap(agent string) error {
	switch agent {
	case "docker", "oras", "library":
		return nil
	case "":
		return errors.New("no Bootstrap: agent recorded; cannot determine what this root was built from")
	default:
		return fmt.Errorf("built from bootstrap agent %q, which cannot be rebuilt from a fresh machine; "+
			"only docker://, oras://, and library:// sources can be imported", agent)
	}
}

// buildForeign packs a foreign root (see FromForeignRoot) into a real .sqf.
// It shares every generic stage with buildDef — the lock, the prebuilt check,
// key derivation, metadata staging, publish — and replaces only the two
// things a live build does that an already-built root cannot: there is no
// apptainer.Build, and the upstream digest comes from the root's own record
// instead of resolveUpstream.
func (b *BuildObject) buildForeign(ctx context.Context) error {
	fr := b.foreignRoot
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
	if pulled, err := b.tryPrebuilt(ctx); err != nil || pulled {
		b.Cleanup(err != nil) //nolint:errcheck
		return err
	}
	if b.skipIfInstalled(ctx) {
		b.Cleanup(false) //nolint:errcheck
		return nil
	}

	log.Info("importing image", "kind", "note", "image", filepath.Base(targetPath), "source", fr.path)

	done := watchContext(ctx, "foreign import")
	defer close(done)

	if err := ensureWorkspaceRoot(b); err != nil {
		b.Cleanup(true) //nolint:errcheck
		return err
	}

	// The same directive lines a live build from this source would produce, so
	// the recipe hash matches an equivalent fresh `create --from`.
	uri := fr.boot.Agent + "://" + fr.boot.From
	defPath, err := synthesizeDefFromURI(uri, b.ws.Root)
	if err != nil {
		b.Cleanup(true) //nolint:errcheck
		return err
	}
	defData, err := os.ReadFile(defPath)
	if err != nil {
		b.Cleanup(true) //nolint:errcheck
		return err
	}
	b.captureSynthesizedRecipe(defPath, defData)
	b.spec.Source.Definition.From = &meta.From{Bootstrap: fr.boot.Agent, Ref: fr.boot.From, Digest: fr.digest}

	if err := b.deriveKeys(ctx); err != nil {
		b.Cleanup(true) //nolint:errcheck
		return err
	}

	metaDir, err := stageMetadata(ctx, b)
	if err != nil {
		b.Cleanup(true) //nolint:errcheck
		return err
	}

	if fr.sandbox {
		if err := checkBash(fr.path); err != nil {
			b.Cleanup(true) //nolint:errcheck
			return err
		}
		if err := createSquashfs(ctx, b, false, fr.path, metaDir, preparedPath); err != nil {
			os.Remove(preparedPath) //nolint:errcheck
			b.Cleanup(true)
			return err
		}
	} else {
		if err := sif.RequireBash(fr.path); err != nil {
			b.Cleanup(true) //nolint:errcheck
			return err
		}
		if err := packFromSIF(ctx, b, fr.path, metaDir, preparedPath); err != nil {
			os.Remove(preparedPath) //nolint:errcheck
			b.Cleanup(true)
			return err
		}
	}
	utils.ShareWithParentGroup(preparedPath)

	installed, err := b.publish(preparedPath, targetPath)
	if err != nil {
		b.Cleanup(true) //nolint:errcheck
		return err
	}

	log.Info("image ready", "kind", "success", "path", installed)
	b.Cleanup(false)
	return nil
}

// packFromSIF mounts a .sif's primary SquashFS partition read-only via
// squashfuse, at the partition's own byte offset, inside freeze.MountedRun's
// unprivileged namespace, and packs it plus metaDir directly.
func packFromSIF(ctx context.Context, b *BuildObject, sifPath, metaDir, targetPath string) error {
	if absTarget, err := filepath.Abs(targetPath); err == nil {
		targetPath = absTarget
	}
	done := watchContext(ctx, "SquashFS creation")
	defer close(done)

	mksquashfsBin, err := toolpath.Resolve("mksquashfs")
	if err != nil {
		return err
	}
	squashfuseBin, err := freeze.FindSquashfuse()
	if err != nil {
		return err
	}
	part, err := sif.PrimarySystemPartition(sifPath)
	if err != nil {
		return err
	}

	mnt := filepath.Join(b.ws.TmpDir, "import-mnt")
	if err := os.MkdirAll(mnt, 0o755); err != nil {
		return fmt.Errorf("create import mountpoint: %w", err)
	}
	defer os.RemoveAll(mnt)

	ncpus := b.effectiveNcpus()
	compressArgs := config.Global.Build.CompressArgs
	blockSize := config.Global.Build.BlockSize
	io := execpkg.IOFromContext(ctx)

	logging.FromContext(ctx).Info("packing SquashFS", "source", sifPath, "target", targetPath)
	script := squashfsScript(mksquashfsBin, []string{mnt, metaDir}, targetPath, ncpus, blockSize, compressArgs, true)
	offsetArg := "offset=" + strconv.FormatInt(part.Offset, 10)

	runErr := freeze.MountedRun(ctx, squashfuseBin, []string{"-o", offsetArg, sifPath}, mnt, script, io)
	if runErr != nil {
		if isCancelledByUser(runErr) {
			return ErrBuildCancelled
		}
		return fmt.Errorf("failed to create SquashFS: %w", runErr)
	}
	return nil
}
