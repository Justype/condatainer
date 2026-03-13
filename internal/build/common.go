package build

import (
	"context"
	"errors"
	"fmt"
	"os"
	"path/filepath"

	"github.com/Justype/condatainer/internal/overlay"
	"github.com/Justype/condatainer/internal/utils"
)

// checkShouldBuild returns (skip=true, nil) if the overlay already exists and update=false.
// In update mode, if the overlay is locked by a running container, returns an error.
func checkShouldBuild(b *BaseBuildObject) (skip bool, err error) {
	styledOverlay := utils.StyleName(filepath.Base(b.targetOverlayPath))
	if !b.update {
		if _, err := os.Stat(b.targetOverlayPath); err == nil {
			utils.PrintMessage("Overlay %s already exists at %s. Skipping creation.",
				styledOverlay, utils.StylePath(b.targetOverlayPath))
			return true, nil
		}
	}
	if b.update && utils.FileExists(b.targetOverlayPath) {
		if lock, err := overlay.AcquireLock(b.targetOverlayPath, true); err != nil {
			return false, fmt.Errorf("cannot update %s: %w", b.nameVersion, err)
		} else {
			lock.Close()
		}
	}
	return false, nil
}

// watchContext starts a goroutine that logs a warning when ctx is cancelled.
// The caller must close(done) when the protected region exits to stop the goroutine.
func watchContext(ctx context.Context, label string) (done chan struct{}) {
	done = make(chan struct{})
	go func() {
		select {
		case <-ctx.Done():
			utils.PrintWarning("Build cancelled. Interrupting %s...", label)
			// Cleanup is the caller's responsibility after exec returns.
		case <-done:
			return
		}
	}()
	return done
}

// buildFinalPath returns targetPath+".new" if update=true, else targetPath.
// In update mode, builds write to a .new file for atomic replacement on success.
func buildFinalPath(targetPath string, update bool) string {
	if update {
		return targetPath + ".new"
	}
	return targetPath
}

// atomicInstall atomically installs finalPath as targetPath.
// In update mode: removes old target, renames finalPath → targetPath.
// In non-update mode: no-op (finalPath == targetPath already).
func atomicInstall(finalPath, targetPath string, update bool) error {
	if !update {
		return nil
	}
	os.Remove(targetPath) //nolint:errcheck
	if err := os.Rename(finalPath, targetPath); err != nil {
		os.Remove(finalPath) //nolint:errcheck
		return fmt.Errorf("failed to replace overlay %s: %w", targetPath, err)
	}
	return nil
}

// prepareBuildWorkspace creates the build workspace (tmp overlay or host dirs).
// Handles stale-artifact detection with one retry after cleanup.
func prepareBuildWorkspace(ctx context.Context, b *BaseBuildObject, useTmpOverlay bool) error {
	if !useTmpOverlay {
		if err := b.CreateBuildDirs(ctx, false); err != nil {
			if !errors.Is(err, ErrTmpOverlayExists) {
				return fmt.Errorf("failed to create build dirs: %w", err)
			}
			utils.PrintWarning("Stale build directory found for %s. Cleaning up...", b.nameVersion)
			b.Cleanup(true)
			if err := b.CreateBuildDirs(ctx, false); err != nil {
				return fmt.Errorf("failed to create build dirs: %w", err)
			}
		}
	} else {
		if err := b.CreateTmpOverlay(ctx, false); err != nil {
			if !errors.Is(err, ErrTmpOverlayExists) {
				return fmt.Errorf("failed to create temporary overlay: %w", err)
			}
			utils.PrintWarning("Stale temporary overlay found for %s. Cleaning up...", b.nameVersion)
			b.Cleanup(true)
			if err := b.CreateTmpOverlay(ctx, false); err != nil {
				return fmt.Errorf("failed to create temporary overlay: %w", err)
			}
		}
	}
	return nil
}

// buildModeLabel returns the build mode string for display ("local" or "sbatch").
func buildModeLabel(b *BaseBuildObject) string {
	if b.RequiresScheduler() {
		return "sbatch"
	}
	return "local"
}

// buildOverlayPaths returns (targetPath, finalPath) where finalPath is the atomic-write
// destination (.new suffix in update mode).
func buildOverlayPaths(b *BaseBuildObject) (targetPath, finalPath string) {
	targetPath = b.targetOverlayPath
	if abs, err := filepath.Abs(targetPath); err == nil {
		targetPath = abs
	}
	finalPath = buildFinalPath(targetPath, b.update)
	return
}
