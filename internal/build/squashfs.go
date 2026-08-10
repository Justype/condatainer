package build

import (
	"context"
	"fmt"
	"path/filepath"
	"strings"

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

// createSquashfs runs mksquashfs inside a container to pack sourceDir into
// targetPath. metaDir is the staged .cnt directory, packed as a second archive
// root; empty packs no metadata. The caller must Cleanup(true) on an error.
func createSquashfs(ctx context.Context, b *BuildObject, isData bool, sourceDir, metaDir, targetPath string) error {
	if absTarget, err := filepath.Abs(targetPath); err == nil {
		targetPath = absTarget
	}

	done := watchContext(ctx, "SquashFS creation")
	defer close(done)

	bashScript, packOverlays, packBindDirs := buildSquashfsOpts(b, isData, sourceDir, metaDir, targetPath)

	opts := execpkg.Options{
		BaseImage:    b.spec.Base,
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

	log := logging.FromContext(ctx)
	log.Debug("creating SquashFS", "name", b.spec.Image.Name, "overlays", opts.Overlays, "bindPaths", opts.BindPaths)
	log.Info("packing SquashFS", "source", sourceDir, "target", targetPath)

	if err := execpkg.Run(ctx, opts, execpkg.IOFromContext(ctx)); err != nil {
		if isCancelledByUser(err) {
			return ErrBuildCancelled
		}
		return fmt.Errorf("failed to create SquashFS: %w", err)
	}

	return nil
}

// buildSquashfsOpts constructs the bash script, overlay list, and bind dirs for
// mksquashfs. The two cases are the workspace's two modes and nothing else.
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

	bashScript = squashfsScript(sources, targetPath, ncpus, blockSize, compressArgs)
	return
}

// squashfsScript renders the mksquashfs invocation.
//
// -keep-as-directory only affects a lone source, where it keeps that directory
// instead of unwrapping it into the archive root; with two sources mksquashfs
// already does that.
//
// -no-xattrs because nothing here reads them back. A shared filesystem hands
// mksquashfs attributes it cannot store (`system.nfs4_acl`) and unsquashfs
// attributes it cannot restore unprivileged (`security.selinux`), so carrying
// them buys warnings on both sides of the round trip and nothing else. Whatever
// an image needs at run time comes from its mode bits and its metadata.
func squashfsScript(sources []string, targetPath string, ncpus int, blockSize, compressArgs string) string {
	return fmt.Sprintf(`
trap 'exit 130' INT TERM
echo "Packing overlay to SquashFS..."
mksquashfs %s %s -processors %d -b %s -keep-as-directory -all-root -no-xattrs %s
`, strings.Join(sources, " "), targetPath, ncpus, blockSize, compressArgs)
}
