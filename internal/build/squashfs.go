package build

import (
	"context"
	"fmt"
	"os"
	"os/exec"
	"path/filepath"
	"strings"

	"github.com/Justype/condatainer/internal/config"
	"github.com/Justype/condatainer/internal/image/freeze"
	"github.com/Justype/condatainer/internal/logging"
	execpkg "github.com/Justype/condatainer/internal/runtime/exec"
	"github.com/Justype/condatainer/internal/toolpath"
)

// createSquashfs packs sourceDir (and metaDir, if any) into targetPath with
// mksquashfs, run directly on the host — no container, no Apptainer. Every
// source is already a plain host path (sandbox, dir mode) or is made one by
// mounting the scratch .img with fuse2fs (ext3 mode), the same
// apptainer-free mechanism internal/image/freeze uses to read a .img.
// isData picks the data block size over the general one. The caller must
// Cleanup(true) on an error.
func createSquashfs(ctx context.Context, b *BuildObject, isData bool, sourceDir, metaDir, targetPath string) error {
	if absTarget, err := filepath.Abs(targetPath); err == nil {
		targetPath = absTarget
	}

	done := watchContext(ctx, "SquashFS creation")
	defer close(done)

	mksquashfsBin, err := toolpath.Resolve("mksquashfs")
	if err != nil {
		return err
	}

	ncpus := b.effectiveNcpus()
	compressArgs := config.Global.Build.CompressArgs
	blockSize := config.Global.Build.BlockSize
	if isData {
		blockSize = config.Global.Build.DataBlockSize
	}

	log := logging.FromContext(ctx)
	io := execpkg.IOFromContext(ctx)

	var runErr error
	if b.ws.UsesImage() {
		log.Info("packing SquashFS", "source", "/cnt", "target", targetPath)
		runErr = packFromScratchImage(ctx, b, metaDir, targetPath, mksquashfsBin, ncpus, blockSize, compressArgs, io)
	} else {
		source := sourceDir
		if b.ws.UsesSandbox() {
			source = b.ws.Sandbox
		}
		sources, keepAsDirectory := packSources(b, sourceDir, metaDir, "")
		log.Debug("creating SquashFS", "name", b.spec.Image.Name, "sources", sources)
		log.Info("packing SquashFS", "source", source, "target", targetPath)
		script := squashfsScript(mksquashfsBin, sources, targetPath, ncpus, blockSize, compressArgs, keepAsDirectory)
		runErr = runHostScript(ctx, script, io)
	}

	if runErr != nil {
		if isCancelledByUser(runErr) {
			return ErrBuildCancelled
		}
		return fmt.Errorf("failed to create SquashFS: %w", runErr)
	}

	return nil
}

// packSources decides mksquashfs's source list, and whether a lone source
// keeps its own directory name in the archive, for the workspace's mode.
// ext3Base is the fuse2fs mount's upper/ directory and is only consulted for
// ext3 mode; the other two modes' paths are already host paths.
func packSources(b *BuildObject, sourceDir, metaDir, ext3Base string) (sources []string, keepAsDirectory bool) {
	if b.ws.UsesSandbox() {
		// The sandbox is the container root apptainer wrote: its own contents
		// belong at the archive root, dotfiles and all, not nested under the
		// directory's name.
		return []string{b.ws.Sandbox}, false
	}
	if b.ws.UsesImage() {
		sources = []string{filepath.Join(ext3Base, "cnt")}
	} else {
		sources = []string{sourceDir}
	}
	if metaDir != "" {
		// A host path already named .cnt (staged by stageMetadata), so it is
		// its own source — mksquashfs names an archive root after a source's
		// basename and cannot rename one.
		sources = append(sources, metaDir)
	}
	return sources, true
}

// packFromScratchImage mounts the build's scratch .img read-only with
// fuse2fs, inside freeze.MountedRun's unprivileged namespace, and packs from
// there. Read-only because packing only reads; the mount is torn down (and
// its mountpoint removed) the moment mksquashfs finishes.
func packFromScratchImage(ctx context.Context, b *BuildObject, metaDir, targetPath, mksquashfsBin string, ncpus int, blockSize, compressArgs string, io execpkg.IO) error {
	fuse2fsBin, err := toolpath.Resolve("fuse2fs")
	if err != nil {
		return err
	}
	mnt := filepath.Join(b.ws.TmpDir, "pack-mnt")
	if err := os.MkdirAll(mnt, 0o755); err != nil {
		return fmt.Errorf("create pack mountpoint: %w", err)
	}
	defer os.RemoveAll(mnt)

	sources, keepAsDirectory := packSources(b, "", metaDir, filepath.Join(mnt, freeze.UpperDir))
	script := squashfsScript(mksquashfsBin, sources, targetPath, ncpus, blockSize, compressArgs, keepAsDirectory)
	return freeze.MountedRun(ctx, fuse2fsBin, []string{"-o", "ro", b.ws.Overlay}, mnt, script, io)
}

// runHostScript runs script with /bin/bash directly on the host, wiring the
// caller's own IO — no container involved.
func runHostScript(ctx context.Context, script string, io execpkg.IO) error {
	cmd := exec.CommandContext(ctx, "/bin/bash", "-c", script)
	cmd.Stdin, cmd.Stdout, cmd.Stderr = io.Stdin, io.Stdout, io.Stderr
	return cmd.Run()
}

// squashfsScript renders the mksquashfs invocation. mksquashfsBin is already
// resolved (toolpath.Resolve, in createSquashfs) — packing must not trigger a
// first-time libexec download, and by the time this runs resolution has
// already succeeded or createSquashfs has already returned its error.
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
func squashfsScript(mksquashfsBin string, sources []string, targetPath string, ncpus int, blockSize, compressArgs string, keepAsDirectory bool) string {
	keep := ""
	if keepAsDirectory {
		keep = "-keep-as-directory "
	}
	return fmt.Sprintf(`
trap 'exit 130' INT TERM
echo "Packing overlay to SquashFS..."
%s %s %s -processors %d -b %s %s-all-root -no-xattrs -quiet %s
`, mksquashfsBin, strings.Join(sources, " "), targetPath, ncpus, blockSize, keep, compressArgs)
}
