package exec

import (
	"context"
	"fmt"
	"os"
	"strings"

	"github.com/Justype/condatainer/internal/image/ext3"
	"github.com/Justype/condatainer/internal/logging"
	"github.com/Justype/condatainer/internal/runtime/container"
	"github.com/Justype/condatainer/internal/utils"
)

// DescribeInitialCondaPackages returns user-facing text for the initial package list.
func DescribeInitialCondaPackages(pkgs []string) string {
	if len(pkgs) == 0 {
		return "empty conda environment"
	}
	return strings.Join(pkgs, " ")
}

// CreateCondaOverlay creates a new user-owned ext3 overlay at opts.Path. When
// pkgs is non-empty it initializes the conda environment before moving and
// allocating the final file; otherwise initialization remains lazy.
//
// postInstallCmd is run inside the overlay after conda init (empty = skip).
// fakeroot should be false for normal user-owned overlays.
func CreateCondaOverlay(ctx context.Context, opts *ext3.CreateOptions, pkgs []string, postInstallCmd string, fakeroot bool, io IO) error {
	if opts.UID == 0 && opts.GID == 0 {
		opts.UID = os.Getuid()
		opts.GID = os.Getgid()
	}
	tmpPath, err := ext3.CreateInTmp(ctx, opts)
	if err != nil {
		return err
	}
	cleanup := func() {
		os.Remove(tmpPath)
		utils.RemoveDirIfEmpty(tmpPath)
	}

	// tmpPath is a scratch path, disconnected from opts.Path — container.Setup's
	// own autoload looks beside the .img it's given, so it can never find a
	// snapshot that lives beside the *final* destination instead. Looked up
	// here and passed through explicitly, so install sees it and writes only
	// the incremental diff rather than reinstalling what the snapshot already has.
	snapshot := container.LookupSnapshot(opts.Path).Path

	if err := InitCondaEnv(ctx, tmpPath, snapshot, pkgs, fakeroot, io); err != nil {
		cleanup()
		return err
	}

	if postInstallCmd != "" {
		if err := RunPostInstall(ctx, tmpPath, snapshot, postInstallCmd, fakeroot, io); err != nil {
			cleanup()
			return fmt.Errorf("post-install failed: %w", err)
		}
	}

	logging.FromContext(ctx).Info(fmt.Sprintf("moving overlay to %s", opts.Path))
	copied, err := ext3.MoveOverlayCopied(ctx, tmpPath, opts.Path, opts.Sparse)
	if err != nil {
		cleanup()
		return err
	}
	if !opts.Sparse && !copied {
		ext3.AllocateOverlay(ctx, opts.Path, opts.SizeMB)
	}
	return nil
}

// InitCondaEnv creates a new conda environment inside imgPath using `condatainer env install`,
// then cleans the micromamba package cache to reduce overlay size.
// Use for initial environment creation on a fresh image.
// For adding packages to an existing environment, use InstallPackages.
//
// snapshot, if non-empty, is mounted read-only beneath imgPath — pass the
// result of container.LookupSnapshot when imgPath is a scratch path that
// won't autoload one on its own (see CreateCondaOverlay); "" elsewhere.
func InitCondaEnv(ctx context.Context, imgPath, snapshot string, pkgs []string, fakeroot bool, io IO) error {
	if len(pkgs) == 0 {
		return nil
	}
	overlays := condaScratchOverlays(imgPath, snapshot)
	cmd := append([]string{container.BoundExecPath, "env", "install", "-y"}, pkgs...)
	if err := Run(ctx, Options{
		Overlays:    overlays,
		WritableImg: true,
		Fakeroot:    fakeroot,
		Command:     cmd,
		HidePrompt:  true,
	}, io); err != nil {
		return err
	}
	// Non-fatal cache clean to reduce overlay size.
	_ = Run(ctx, Options{
		Overlays:    overlays,
		WritableImg: true,
		Fakeroot:    fakeroot,
		Command:     []string{container.BoundExecPath, "env", "clean", "-a", "-y", "-q"},
		HidePrompt:  true,
	}, IO{})
	return nil
}

// condaScratchOverlays is the overlay list for a conda operation on imgPath,
// with snapshot appended when given.
func condaScratchOverlays(imgPath, snapshot string) []string {
	if snapshot == "" {
		return []string{imgPath}
	}
	return []string{imgPath, snapshot}
}

// InstallPackages installs additional conda packages into an existing base
// environment inside imgPath using `condatainer env install` (respects the overlay's .condarc channels).
func InstallPackages(ctx context.Context, imgPath string, pkgs []string, fakeroot bool, io IO) error {
	return Run(ctx, Options{
		Overlays:    []string{imgPath},
		WritableImg: true,
		Fakeroot:    fakeroot,
		Command:     append([]string{container.BoundExecPath, "env", "install", "-y"}, pkgs...),
		HidePrompt:  true,
	}, io)
}

// RemovePackages removes conda packages from the base environment inside imgPath.
func RemovePackages(ctx context.Context, imgPath string, pkgs []string, fakeroot bool, io IO) error {
	return Run(ctx, Options{
		Overlays:    []string{imgPath},
		WritableImg: true,
		Fakeroot:    fakeroot,
		Command:     append([]string{container.BoundExecPath, "env", "remove", "-y"}, pkgs...),
		HidePrompt:  true,
	}, io)
}

// RunPostInstall runs an arbitrary command inside imgPath after package
// installation. snapshot is as in InitCondaEnv.
func RunPostInstall(ctx context.Context, imgPath, snapshot, postCmd string, fakeroot bool, io IO) error {
	return Run(ctx, Options{
		Overlays:    condaScratchOverlays(imgPath, snapshot),
		WritableImg: true,
		Fakeroot:    fakeroot,
		Command:     []string{"bash", "-c", postCmd},
		HidePrompt:  true,
	}, io)
}
