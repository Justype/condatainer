package build

import (
	"context"
	"fmt"
	"path/filepath"
	"strings"

	"github.com/Justype/condatainer/catalog"
	"github.com/Justype/condatainer/internal/artifact/meta"
	"github.com/Justype/condatainer/internal/config"
	"github.com/Justype/condatainer/internal/logging"
	"github.com/Justype/condatainer/internal/runtime/container"
	execpkg "github.com/Justype/condatainer/internal/runtime/exec"
)

// metaMountPath is where the staged metadata directory is bound while packing.
// mksquashfs names an archive root after its source's basename and cannot rename
// one, so the directory must already be called .cnt when mksquashfs sees it.
const metaMountPath = "/" + meta.DirName

// sandboxMountPath is where a definition's sandbox is bound while packing. It is
// an empty directory in every base, so binding over it hides nothing.
const sandboxMountPath = "/mnt"

// createSquashfs runs mksquashfs inside a container to pack sourceDir into
// targetPath. metaDir is the staged .cnt directory, packed as a second archive
// root; empty packs no metadata. A sandbox build passes no metaDir: its .cnt is
// already inside the tree. The caller must Cleanup(true) on an error.
func createSquashfs(ctx context.Context, b *BuildObject, isData bool, sourceDir, metaDir, targetPath string) error {
	if absTarget, err := filepath.Abs(targetPath); err == nil {
		targetPath = absTarget
	}

	done := watchContext(ctx, "SquashFS creation")
	defer close(done)

	bashScript, packOverlays, packBindDirs := buildSquashfsOpts(b, isData, sourceDir, metaDir, targetPath)

	// A base packs itself; everything else is packed by a base.
	baseImage := b.spec.Base
	if baseImage == "" && b.spec.Image.Type == catalog.TypeBase {
		baseImage = b.ws.Sandbox
	}

	opts := execpkg.Options{
		BaseImage:    baseImage,
		ApptainerBin: config.Global.ApptainerBin,
		Overlays:     packOverlays,
		BindPaths:    packBindDirs,
		Command:      []string{"/bin/bash", "-c", bashScript},
		HidePrompt:   true,
		WritableImg:  b.ws.UsesImage(),
	}
	if !b.ws.UsesImage() {
		opts.ApptainerFlags = []string{"--writable-tmpfs"}
	}

	// The sandbox is reported by its host path: sourceDir is where it is bound
	// inside the packing container, which names nothing the user has.
	source := sourceDir
	if b.ws.UsesSandbox() {
		source = b.ws.Sandbox
	}

	log := logging.FromContext(ctx)
	log.Debug("creating SquashFS", "name", b.spec.Image.Name, "overlays", opts.Overlays, "bindPaths", opts.BindPaths)
	log.Info("packing SquashFS", "source", source, "target", targetPath)

	if err := execpkg.Run(ctx, opts, execpkg.IOFromContext(ctx)); err != nil {
		if isCancelledByUser(err) {
			return ErrBuildCancelled
		}
		return fmt.Errorf("failed to create SquashFS: %w", err)
	}

	return nil
}

// buildSquashfsOpts constructs the bash script, overlay list, and bind dirs for
// mksquashfs. The three cases are the workspace's three modes and nothing else:
// a sandbox packs itself, an ext3 scratch image is mounted, dir mode binds a
// host path.
func buildSquashfsOpts(b *BuildObject, isData bool, sourceDir, metaDir, targetPath string) (bashScript string, overlays, bindDirs []string) {
	ncpus := b.effectiveNcpus()
	compressArgs := config.Global.Build.CompressArgs
	blockSize := config.Global.Build.BlockSize
	if isData {
		blockSize = config.Global.Build.DataBlockSize
	}

	// One archive root per source, which is what puts /cnt/... and /.cnt/... at
	// the same level without copying the metadata into the payload.
	sources := []string{sourceDir}
	bindDirs = container.DeduplicateBindPaths(getAllBaseDirs())

	if b.ws.UsesSandbox() {
		// The sandbox is the container root, so its own contents are at / along
		// with apptainer's runtime binds. Bind it a second time and pack that,
		// or mksquashfs would descend into the host.
		overlays = []string{}
		bindDirs = append(bindDirs, b.ws.Sandbox+":"+sandboxMountPath, filepath.Dir(targetPath))
		bashScript = squashfsScript([]string{sandboxMountPath}, targetPath, ncpus, blockSize, compressArgs, false)
		return bashScript, overlays, bindDirs
	}

	if b.ws.UsesImage() {
		// sourceDir is a path inside the image, so only the staged metadata has
		// to come in from the host.
		overlays = []string{b.ws.Overlay}
		if metaDir != "" {
			sources = append(sources, metaMountPath)
			bindDirs = append(bindDirs, metaDir+":"+metaMountPath)
		}
	} else {
		// Dir mode: sourceDir is a host path; bind it and the output dir.
		overlays = []string{}
		bindDirs = append(bindDirs, sourceDir, filepath.Dir(targetPath))
		if metaDir != "" {
			// A host path already named .cnt, so it is its own source.
			sources = append(sources, metaDir)
			bindDirs = append(bindDirs, metaDir)
		}
	}

	bashScript = squashfsScript(sources, targetPath, ncpus, blockSize, compressArgs, true)
	return
}

// squashfsScript renders the mksquashfs invocation.
//
// keepAsDirectory only affects a lone source, where it keeps that directory
// instead of unwrapping it into the archive root; with two sources mksquashfs
// already does that. A sandbox turns it off: its contents belong at the archive
// root, dotfiles and all, not nested under the directory's name.
//
// -no-xattrs because nothing here reads them back. A shared filesystem hands
// mksquashfs attributes it cannot store (`system.nfs4_acl`) and unsquashfs
// attributes it cannot restore unprivileged (`security.selinux`), so carrying
// them buys warnings on both sides of the round trip and nothing else. Whatever
// an image needs at run time comes from its mode bits and its metadata.
//
// -quiet suppresses the final filesystem statistics but leaves the progress bar
// enabled. Do not pair it with -no-progress: progress is the useful build output.
func squashfsScript(sources []string, targetPath string, ncpus int, blockSize, compressArgs string, keepAsDirectory bool) string {
	keep := ""
	if keepAsDirectory {
		keep = "-keep-as-directory "
	}
	return fmt.Sprintf(`
trap 'exit 130' INT TERM
echo "Packing overlay to SquashFS..."
mksquashfs %s %s -processors %d -b %s %s-all-root -no-xattrs -quiet %s
`, strings.Join(sources, " "), targetPath, ncpus, blockSize, keep, compressArgs)
}
