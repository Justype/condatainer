package build

import (
	"context"
	"encoding/json"
	"errors"
	"fmt"
	"os"
	"os/exec"
	"path/filepath"
	"runtime"
	"strings"
	"syscall"
	"time"

	"github.com/Justype/condatainer/internal/config"
	"github.com/Justype/condatainer/internal/container"
	"github.com/Justype/condatainer/internal/overlay"
	"github.com/Justype/condatainer/internal/scheduler"
	"github.com/Justype/condatainer/internal/utils"
)

// cachedInstalledOverlays is the set of "name/version" strings for all overlays found
// across all image search paths. Nil means the cache is cold or has been invalidated.
var cachedInstalledOverlays map[string]bool

// getInstalledOverlays returns the cached installed-overlay set, scanning all image
// search paths on a cold cache. Follows the same pattern as cachedLocalScripts in fetch.go.
func getInstalledOverlays() map[string]bool {
	if cachedInstalledOverlays != nil {
		return cachedInstalledOverlays
	}
	installed := make(map[string]bool)
	for _, imagesDir := range config.GetImageSearchPaths() {
		if !utils.DirExists(imagesDir) {
			continue
		}
		entries, err := os.ReadDir(imagesDir)
		if err != nil {
			utils.PrintWarning("Failed to read directory %s: %v", imagesDir, err)
			continue
		}
		for _, entry := range entries {
			if entry.IsDir() || !utils.IsOverlay(entry.Name()) {
				continue
			}
			nameVersion := strings.TrimSuffix(entry.Name(), filepath.Ext(entry.Name()))
			installed[strings.ReplaceAll(nameVersion, "--", "/")] = true
		}
	}
	cachedInstalledOverlays = installed
	return installed
}

// checkShouldBuild returns (skip=true, nil) if the overlay already exists and update=false.
// In update mode, if the overlay is locked by a running container, returns an error.
func checkShouldBuild(b *BuildObject) (skip bool, err error) {
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
// In non-update mode: finalPath == targetPath already; only invalidates caches.
func atomicInstall(finalPath, targetPath string, update bool) error {
	if update {
		os.Remove(targetPath) //nolint:errcheck
		if err := os.Rename(finalPath, targetPath); err != nil {
			os.Remove(finalPath) //nolint:errcheck
			return fmt.Errorf("failed to replace overlay %s: %w", targetPath, err)
		}
	}
	cachedInstalledOverlays = nil                // invalidate so next dep-check sees the new overlay
	container.InvalidateInstalledOverlaysCache() // invalidate container resolve cache too
	return nil
}

// prepareBuildWorkspace creates the build workspace (tmp overlay or host dirs).
// Stale artifacts are detected by Create*; on detection a warning is printed and
// the workspace is re-created with force=true, which removes only build dirs/overlays
// and leaves the remote build source intact.
func prepareBuildWorkspace(ctx context.Context, b *BuildObject, useTmpOverlay bool) error {
	if !useTmpOverlay {
		if err := b.CreateBuildDirs(ctx, false); err != nil {
			if !errors.Is(err, ErrTmpOverlayExists) {
				return fmt.Errorf("failed to create build dirs: %w", err)
			}
			utils.PrintWarning("Stale build directory found for %s. Cleaning up...", b.nameVersion)
			if err := b.CreateBuildDirs(ctx, true); err != nil {
				return fmt.Errorf("failed to create build dirs: %w", err)
			}
		}
	} else {
		if err := b.CreateTmpOverlay(ctx, false); err != nil {
			if !errors.Is(err, ErrTmpOverlayExists) {
				return fmt.Errorf("failed to create temporary overlay: %w", err)
			}
			utils.PrintWarning("Stale temporary overlay found for %s. Cleaning up...", b.nameVersion)
			if err := b.CreateTmpOverlay(ctx, true); err != nil {
				return fmt.Errorf("failed to create temporary overlay: %w", err)
			}
		}
	}
	return nil
}

// buildModeLabel returns the build mode string for display ("local" or "sbatch").
func buildModeLabel(b *BuildObject) string {
	if b.RequiresScheduler() {
		return "sbatch"
	}
	return "local"
}

// buildOverlayPaths returns (targetPath, finalPath) where finalPath is the atomic-write
// destination (.new suffix in update mode).
func buildOverlayPaths(b *BuildObject) (targetPath, finalPath string) {
	targetPath = b.targetOverlayPath
	if abs, err := filepath.Abs(targetPath); err == nil {
		targetPath = abs
	}
	finalPath = buildFinalPath(targetPath, b.update)
	return
}

// tryDownloadPrebuilt attempts to download a prebuilt asset from the given prebuiltLink base URL.
// ext is the file extension without a dot (e.g. "sif", "sqf").
// downloadFn is called with (ctx, url, destPath) and should return an error on failure.
// Returns false immediately if prebuiltLink is empty (no prebuilt source for this script).
func tryDownloadPrebuilt(ctx context.Context, nameVersion, destPath, ext, prebuiltLink string, downloadFn func(context.Context, string, string) error) bool {
	if prebuiltLink == "" {
		return false
	}
	archMap := map[string]string{
		"amd64": "x86_64",
		"arm64": "aarch64",
	}
	archName, ok := archMap[runtime.GOARCH]
	if !ok {
		return false
	}
	normalized := utils.NormalizeNameVersion(nameVersion)
	parts := strings.SplitN(normalized, "/", 2)
	if len(parts) != 2 {
		return false
	}
	url := fmt.Sprintf("%s/%s/%s_%s.%s", prebuiltLink, parts[0], parts[1], archName, ext)
	if !utils.URLExists(ctx, url) {
		return false
	}
	utils.PrintMessage("Found pre-built %s. Downloading...", utils.StyleName(normalized))
	if err := downloadFn(ctx, url, destPath); err != nil {
		utils.PrintWarning("Download failed. Falling back to local build.")
		return false
	}
	utils.PrintSuccess("Pre-built %s downloaded.", utils.StyleName(normalized))
	return true
}

// isCancelledByUser checks if the error is due to user cancellation (Ctrl+C)
// Exit code 130 = 128 + SIGINT(2), checks for "signal: killed/interrupt" or context errors
func isCancelledByUser(err error) bool {
	if errors.Is(err, context.Canceled) || errors.Is(err, context.DeadlineExceeded) {
		return true
	}
	errMsg := err.Error()
	if strings.Contains(errMsg, "signal: killed") || strings.Contains(errMsg, "signal: interrupt") {
		return true
	}
	var exitErr *exec.ExitError
	if errors.As(err, &exitErr) {
		// 130 (SIGINT) or -1 (signal killed)
		return exitErr.ExitCode() == 130 || exitErr.ExitCode() == -1
	}
	return false
}

// shortHostname returns the unqualified hostname (strips domain suffix).
// os.Hostname() may return "cn001" or "cn001.cluster.edu" depending on system
// configuration; always store/compare the short form to avoid false mismatches.
func shortHostname() string {
	h, _ := os.Hostname()
	if idx := strings.Index(h, "."); idx > 0 {
		return h[:idx]
	}
	return h
}

// buildDefaults holds resource defaults for build operations.
// Set from config at CLI startup via SetBuildDefaults.
var buildDefaults = scheduler.ResourceSpec{
	Nodes:        1,
	TasksPerNode: 1,
	CpusPerTask:  4, // conservative default for local builds
	MemPerNodeMB: 8192,
	Time:         2 * time.Hour,
}

// SetBuildDefaults sets the resource defaults used for build job submissions.
func SetBuildDefaults(d scheduler.ResourceSpec) { buildDefaults = d }

// GetBuildDefaults returns the current build resource defaults.
func GetBuildDefaults() scheduler.ResourceSpec { return buildDefaults }

// buildEffectiveResourceSpec resolves resources for build using the priority chain:
//
//	buildDefaults → scriptSpecs.Spec (when HasDirectives=true) → scheduler job resources
func buildEffectiveResourceSpec(specs *scheduler.ScriptSpecs) *scheduler.ResourceSpec {
	var jobRes *scheduler.ResourceSpec
	if sched := scheduler.ActiveScheduler(); sched != nil {
		jobRes = sched.GetJobResources()
	}
	return scheduler.ResolveResourceSpecFrom(buildDefaults, jobRes, specs)
}

// acquireBuildLockFile is a package-level helper that creates a lock file
// atomically (O_CREATE|O_EXCL) and writes JSON metadata.
// Used by both BuildObject and graph.go's submitJob.
func acquireBuildLockFile(path string, info BuildLockInfo) error {
	data, err := json.Marshal(info)
	if err != nil {
		return fmt.Errorf("failed to marshal build lock: %w", err)
	}
	f, err := os.OpenFile(path, os.O_CREATE|os.O_EXCL|os.O_WRONLY, utils.PermFile)
	if err != nil {
		return err // caller checks os.IsExist
	}
	defer f.Close()
	_, err = f.Write(data)
	return err
}

// overwriteBuildLockFile overwrites an existing lock file with new JSON metadata.
// The caller must already hold the lock (i.e. have created it via acquireBuildLockFile).
func overwriteBuildLockFile(path string, info BuildLockInfo) error {
	data, err := json.Marshal(info)
	if err != nil {
		return fmt.Errorf("failed to marshal build lock: %w", err)
	}
	return os.WriteFile(path, data, utils.PermFile)
}

// readBuildLockFile reads and parses a lock file at the given path.
func readBuildLockFile(path string) (BuildLockInfo, error) {
	data, err := os.ReadFile(path)
	if err != nil {
		return BuildLockInfo{}, err
	}
	if len(data) == 0 {
		// Old empty-lock format — treat as stale.
		return BuildLockInfo{}, nil
	}
	var info BuildLockInfo
	if err := json.Unmarshal(data, &info); err != nil {
		return BuildLockInfo{}, fmt.Errorf("corrupt lock file: %w", err)
	}
	return info, nil
}

// isBuildLockStale returns whether the lock is stale, the job's current status, and any
// uncertainty error. Returns (true, Unknown, nil) when definitely stale, (false, status, nil)
// when definitely alive, or (false, Unknown, err) when the state cannot be verified.
func isBuildLockStale(info BuildLockInfo) (stale bool, status scheduler.JobStatus, err error) {
	// Empty type means old empty-lock format → treat as stale for backward compat.
	if info.Type == "" {
		return true, scheduler.JobStatusUnknown, nil
	}

	if info.Type != "local" {
		if info.JobID == "" {
			// Lock was written before submit returned — treat as stale.
			return true, scheduler.JobStatusUnknown, nil
		}
		sched := scheduler.ActiveScheduler()
		if sched == nil {
			// Can't check without a scheduler — be conservative.
			return false, scheduler.JobStatusUnknown, fmt.Errorf("scheduler unavailable, cannot verify job %s", info.JobID)
		}
		st, err := sched.GetJobStatus(info.JobID)
		if err != nil {
			return false, scheduler.JobStatusUnknown, fmt.Errorf("cannot check job %s: %w", info.JobID, err)
		}
		if st == scheduler.JobStatusUnknown {
			// Can't determine state — be conservative (treat as alive).
			return false, st, fmt.Errorf("cannot determine status of job %s", info.JobID)
		}
		return !st.IsAlive(), st, nil
	}

	// Local lock: compare node + PID.
	if info.Node != shortHostname() {
		return false, scheduler.JobStatusUnknown, fmt.Errorf("lock held by node %q (current: %q); cannot verify remotely", info.Node, shortHostname())
	}
	// Same node: check if the PID is still alive via signal 0.
	// EPERM means process exists (different owner); ESRCH means no such process.
	proc, err := os.FindProcess(info.PID)
	if err != nil {
		return true, scheduler.JobStatusUnknown, nil // process not found → stale
	}
	if err := proc.Signal(syscall.Signal(0)); err != nil {
		if err == syscall.EPERM {
			return false, scheduler.JobStatusRunning, nil // alive, different owner
		}
		return true, scheduler.JobStatusUnknown, nil // ESRCH → process gone → stale
	}
	return false, scheduler.JobStatusRunning, nil // process alive
}
