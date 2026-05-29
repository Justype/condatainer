package exec

import (
	"context"
	"fmt"
	"os"
	"strings"

	"github.com/Justype/condatainer/internal/logging"
	"github.com/Justype/condatainer/internal/overlay"
	"github.com/Justype/condatainer/internal/utils"
)

const minimalCondaPackage = "zlib"

// ResolveInitialCondaPackages returns the package list used to initialize a
// fresh conda overlay. An empty list gets a tiny package so mm-create creates
// the environment instead of receiving no specs.
func ResolveInitialCondaPackages(pkgs []string) []string {
	if len(pkgs) == 0 {
		return []string{minimalCondaPackage}
	}
	return pkgs
}

// DescribeInitialCondaPackages returns user-facing text for the resolved initial package list.
func DescribeInitialCondaPackages(pkgs []string) string {
	if len(pkgs) == 0 {
		return "minimal conda environment with small package (" + minimalCondaPackage + ")"
	}
	return strings.Join(pkgs, " ")
}

// CreateCondaOverlay creates a new user-owned ext3 overlay at opts.Path,
// initializes a conda environment with pkgs using mm-create, then moves and
// allocates the final file. This mirrors the `condatainer overlay create -- <pkgs>` flow.
//
// postInstallCmd is run inside the overlay after conda init (empty = skip).
// fakeroot should be false for normal user-owned overlays.
func CreateCondaOverlay(ctx context.Context, opts *overlay.CreateOptions, pkgs []string, postInstallCmd string, fakeroot bool, io IO) error {
	if opts.UID == 0 && opts.GID == 0 {
		opts.UID = os.Getuid()
		opts.GID = os.Getgid()
	}
	pkgs = ResolveInitialCondaPackages(pkgs)

	tmpPath, err := overlay.CreateInTmp(ctx, opts)
	if err != nil {
		return err
	}
	cleanup := func() {
		os.Remove(tmpPath)
		utils.RemoveDirIfEmpty(tmpPath)
	}

	if err := InitCondaEnv(ctx, tmpPath, pkgs, fakeroot, io); err != nil {
		cleanup()
		return err
	}

	if postInstallCmd != "" {
		if err := RunPostInstall(ctx, tmpPath, postInstallCmd, fakeroot, io); err != nil {
			cleanup()
			return fmt.Errorf("post-install failed: %w", err)
		}
	}

	logging.FromContext(ctx).Info(fmt.Sprintf("moving overlay to %s", opts.Path))
	copied, err := overlay.MoveOverlayCopied(ctx, tmpPath, opts.Path, opts.Sparse)
	if err != nil {
		cleanup()
		return err
	}
	if !opts.Sparse && !copied {
		overlay.AllocateOverlay(ctx, opts.Path, opts.SizeMB)
	}
	return nil
}

// InitCondaEnv creates a new conda environment inside imgPath using mm-create,
// then cleans the micromamba package cache to reduce overlay size.
// Use for initial environment creation on a fresh overlay.
// For adding packages to an existing environment, use InstallPackages.
func InitCondaEnv(ctx context.Context, imgPath string, pkgs []string, fakeroot bool, io IO) error {
	pkgs = ResolveInitialCondaPackages(pkgs)
	cmd := append([]string{"mm-create", "-y"}, pkgs...)
	if err := Run(ctx, Options{
		Overlays:    []string{imgPath},
		WritableImg: true,
		Fakeroot:    fakeroot,
		Command:     cmd,
		HidePrompt:  true,
	}, io); err != nil {
		return err
	}
	// Non-fatal cache clean to reduce overlay size.
	_ = Run(ctx, Options{
		Overlays:    []string{imgPath},
		WritableImg: true,
		Fakeroot:    fakeroot,
		Command:     []string{"mm-clean", "-a", "-y", "-q"},
		HidePrompt:  true,
	}, IO{})
	return nil
}

// InstallPackages installs additional conda packages into an existing base
// environment inside imgPath using mm-install (respects overlay's .condarc channels).
func InstallPackages(ctx context.Context, imgPath string, pkgs []string, fakeroot bool, io IO) error {
	return Run(ctx, Options{
		Overlays:    []string{imgPath},
		WritableImg: true,
		Fakeroot:    fakeroot,
		Command:     append([]string{"mm-install", "-y"}, pkgs...),
		HidePrompt:  true,
	}, io)
}

// RemovePackages removes conda packages from the base environment inside imgPath.
func RemovePackages(ctx context.Context, imgPath string, pkgs []string, fakeroot bool, io IO) error {
	return Run(ctx, Options{
		Overlays:    []string{imgPath},
		WritableImg: true,
		Fakeroot:    fakeroot,
		Command:     append([]string{"mm-remove", "-y"}, pkgs...),
		HidePrompt:  true,
	}, io)
}

// RunPostInstall runs an arbitrary command inside imgPath after package installation.
func RunPostInstall(ctx context.Context, imgPath, postCmd string, fakeroot bool, io IO) error {
	return Run(ctx, Options{
		Overlays:    []string{imgPath},
		WritableImg: true,
		Fakeroot:    fakeroot,
		Command:     []string{"bash", "-c", postCmd},
		HidePrompt:  true,
	}, io)
}
