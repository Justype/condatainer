package build

import (
	"context"
	"fmt"
	"path/filepath"

	"github.com/Justype/condatainer/internal/config"
	"github.com/Justype/condatainer/internal/container"
	execpkg "github.com/Justype/condatainer/internal/exec"
	"github.com/Justype/condatainer/internal/utils"
)

// createSquashfs runs mksquashfs inside a container to pack sourceDir into targetPath.
// It handles three modes depending on useTmpOverlay and sourceDir:
//   - Dir mode (!useTmpOverlay): sourceDir is a host path; no overlay needed.
//   - Ext3 app mode (useTmpOverlay, sourceDir=="/cnt"): pack from /cnt inside the overlay.
//   - Ext3 ref mode (useTmpOverlay, other sourceDir): pack from a host path with no overlay.
//
// The caller is responsible for calling b.Cleanup(true) if an error is returned.
func createSquashfs(ctx context.Context, b *BuildObject, isRef bool, sourceDir, targetPath string) error {
	if absTarget, err := filepath.Abs(targetPath); err == nil {
		targetPath = absTarget
	}

	done := watchContext(ctx, "SquashFS creation")
	defer close(done)

	bashScript, packOverlays, packBindDirs := buildSquashfsOpts(b, isRef, sourceDir, targetPath)

	opts := execpkg.Options{
		BaseImage:    config.GetBaseImage(),
		ApptainerBin: config.Global.ApptainerBin,
		Overlays:     packOverlays,
		BindPaths:    packBindDirs,
		Command:      []string{"/bin/bash", "-c", bashScript},
		HidePrompt:   true,
		WritableImg:  config.Global.Build.UseTmpOverlay,
	}
	if !config.Global.Build.UseTmpOverlay {
		opts.ApptainerFlags = []string{"--writable-tmpfs"}
	}

	utils.PrintDebug("[BUILD] Creating SquashFS for %s with options: overlays=%v, bindPaths=%v",
		b.nameVersion, opts.Overlays, opts.BindPaths)
	utils.PrintMessage("Packing %s into %s", utils.StylePath(sourceDir), utils.StylePath(targetPath))

	if err := execpkg.Run(ctx, opts); err != nil {
		if isCancelledByUser(err) {
			return ErrBuildCancelled
		}
		return fmt.Errorf("failed to create SquashFS: %w", err)
	}

	return nil
}

// buildSquashfsOpts constructs the bash script, overlay list, and bind dirs for mksquashfs.
func buildSquashfsOpts(b *BuildObject, isRef bool, sourceDir, targetPath string) (bashScript string, overlays, bindDirs []string) {
	ncpus := b.effectiveNcpus()
	compressArgs := config.Global.Build.CompressArgs

	if !config.Global.Build.UseTmpOverlay {
		// Dir mode: sourceDir is a host path; bind it and the output dir.
		blockSize := config.Global.Build.BlockSize
		if isRef {
			blockSize = config.Global.Build.DataBlockSize
		}
		bashScript = fmt.Sprintf(`
trap 'exit 130' INT TERM
echo "Packing overlay to SquashFS..."
mksquashfs %s %s -processors %d -b %s -keep-as-directory -all-root %s
`, sourceDir, targetPath, ncpus, blockSize, compressArgs)
		overlays = []string{}
		bindDirs = append(container.DeduplicateBindPaths(getAllBaseDirs()), sourceDir, filepath.Dir(targetPath))
	} else if sourceDir == "/cnt" {
		// Ext3 mode, app overlays: pack from /cnt inside the overlay.
		bashScript = fmt.Sprintf(`
trap 'exit 130' INT TERM
echo "Packing overlay to SquashFS..."
mksquashfs /cnt %s -processors %d -b %s -keep-as-directory -all-root %s
`, targetPath, ncpus, config.Global.Build.BlockSize, compressArgs)
		overlays = []string{b.tmpOverlayPath}
		bindDirs = container.DeduplicateBindPaths(getAllBaseDirs())
	} else {
		// Ext3 mode, ref overlays: pack from host path (no overlay needed).
		bashScript = fmt.Sprintf(`
trap 'exit 130' INT TERM
echo "Packing overlay to SquashFS..."
mksquashfs %s %s -processors %d -b %s -keep-as-directory -all-root %s
`, sourceDir, targetPath, ncpus, config.Global.Build.DataBlockSize, compressArgs)
		overlays = []string{}
		bindDirs = container.DeduplicateBindPaths(getAllBaseDirs())
	}
	return
}
