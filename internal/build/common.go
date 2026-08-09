package build

import (
	"context"
	"encoding/json"
	"errors"
	"fmt"
	"os"
	"os/exec"
	"path/filepath"
	"strconv"
	"strings"
	"syscall"

	"log/slog"

	"github.com/Justype/condatainer/internal/config"
	"github.com/Justype/condatainer/internal/image"
	"github.com/Justype/condatainer/internal/logging"
	"github.com/Justype/condatainer/internal/runtime/container"
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
	// No aliases: a #DEP: names an image exactly, so a bare name must not be
	// satisfied by the distro overlay that happens to share it.
	scan, err := image.ScanOverlays(image.ScanOptions{})
	if err != nil {
		slog.Default().Warn("failed to read image directory", "err", err)
	}
	cachedInstalledOverlays = image.Names(scan)
	return cachedInstalledOverlays
}

// checkShouldBuild returns (skip=true, nil) if the overlay already exists and update=false.
// In update mode, if the overlay is locked by a running container, returns an error.
func checkShouldBuild(b *BuildObject) (skip bool, err error) {
	// IsInstalled, not a stat of the target: for a base that means every image
	// search path, so one supplied by a shared install is not rebuilt into the
	// user's own directory. Identical to a stat for every other type.
	if !b.update && b.IsInstalled() {
		slog.Default().Info("overlay already exists, skipping",
			"overlay", filepath.Base(b.tgt.Path), "path", b.tgt.Path)
		return true, nil
	}
	if b.update && utils.FileExists(b.tgt.Path) {
		if lock, err := image.AcquireLock(b.tgt.Path, true); err != nil {
			return false, fmt.Errorf("cannot update %s: %w", b.spec.Image.Name, err)
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
			logging.FromContext(ctx).Warn("build cancelled, interrupting", "step", label)
			// Cleanup is the caller's responsibility after exec returns.
		case <-done:
			return
		}
	}()
	return done
}

// preparedSuffix marks a build's in-progress output.
const preparedSuffix = ".part"

// preparedPathFor derives where a build writes its output: beside the target,
// tagged with the lock owner so stale-lock cleanup can recompute it.
// See the README's Prepared output.
func preparedPathFor(targetPath string, info BuildLockInfo) string {
	owner := info.Runner
	if owner == "" {
		owner = "local"
	}
	tag := owner
	if info.JobID != "" {
		tag += "-" + info.JobID
	} else {
		tag += "-" + info.Node + "-" + strconv.Itoa(info.PID)
	}
	return targetPath + "." + sanitizeTag(tag) + preparedSuffix
}

// sanitizeTag keeps a lock owner usable as a filename component.
func sanitizeTag(tag string) string {
	return strings.Map(func(r rune) rune {
		switch {
		case r >= 'a' && r <= 'z', r >= 'A' && r <= 'Z', r >= '0' && r <= '9', r == '-', r == '_':
			return r
		default:
			return '-'
		}
	}, tag)
}

// atomicInstall renames preparedPath over targetPath and invalidates the
// installed-overlay caches. The installed image is never removed first — see
// the README's Prepared output.
func atomicInstall(preparedPath, targetPath string) error {
	if preparedPath != targetPath {
		if err := os.Rename(preparedPath, targetPath); err != nil {
			os.Remove(preparedPath) //nolint:errcheck
			return fmt.Errorf("failed to install overlay %s: %w", targetPath, err)
		}
	}
	cachedInstalledOverlays = nil                // invalidate so next dep-check sees the new overlay
	container.InvalidateInstalledOverlaysCache() // invalidate container resolve cache too
	return nil
}

// prepareBuildWorkspace creates the build workspace: an ext3 scratch image, or
// host directories. Script and Conda builds only. A stale workspace is warned
// about and re-created, leaving a fetched build source intact.
func prepareBuildWorkspace(ctx context.Context, b *BuildObject) error {
	if !b.ws.UsesImage() {
		if err := b.CreateBuildDirs(ctx, false); err != nil {
			if !errors.Is(err, ErrTmpOverlayExists) {
				return fmt.Errorf("failed to create build dirs: %w", err)
			}
			logging.FromContext(ctx).Warn("stale build directory found, cleaning up", "name", b.spec.Image.Name)
			if err := b.CreateBuildDirs(ctx, true); err != nil {
				return fmt.Errorf("failed to create build dirs: %w", err)
			}
		}
	} else {
		if err := b.CreateTmpOverlay(ctx, false); err != nil {
			if !errors.Is(err, ErrTmpOverlayExists) {
				return fmt.Errorf("failed to create temporary overlay: %w", err)
			}
			logging.FromContext(ctx).Warn("stale temporary overlay found, cleaning up", "name", b.spec.Image.Name)
			if err := b.CreateTmpOverlay(ctx, true); err != nil {
				return fmt.Errorf("failed to create temporary overlay: %w", err)
			}
		}
	}

	// A host payload is bound over the install prefix, so its leaf has to exist
	// before the container starts. In ext3 mode CreateTmpOverlay makes no host
	// directories at all, which is why this is here and not in either branch.
	if b.ws.HostPayload() {
		payloadDir := filepath.Join(b.ws.CntDir, b.spec.Image.Name)
		if err := utils.MkdirAllShared(payloadDir); err != nil {
			return fmt.Errorf("failed to create payload dir %s: %w", payloadDir, err)
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

// buildDefaults holds resource defaults for build operations. Set from config at
// CLI startup via SetBuildDefaults; these values are what a caller that skips the
// CLI sees, so they track the config defaults rather than restating them.
var buildDefaults = scheduler.ResourceSpec{
	Nodes:        1,
	TasksPerNode: 1,
	CpusPerTask:  config.DefaultNcpus,
	MemPerNodeMB: config.DefaultMemMB,
	Time:         config.DefaultBuildDuration,
}

// SetBuildDefaults sets the resource defaults used for build job submissions.
func SetBuildDefaults(d scheduler.ResourceSpec) { buildDefaults = d }

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
	if _, err = f.Write(data); err != nil {
		return err
	}
	// Lives next to the image in the images/ data dir, which may be a shared install;
	// share with the parent group so other members can clear a stale lock.
	utils.ShareWithParentGroup(path)
	return nil
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
	if info.Runner == "" {
		return true, scheduler.JobStatusUnknown, nil
	}

	if info.Runner != "local" {
		if info.JobID == "" {
			// Lock was written before submit returned — treat as stale.
			return true, scheduler.JobStatusUnknown, nil
		}
		sched := scheduler.ActiveScheduler()
		if sched == nil {
			// Can't check without a scheduler — be conservative.
			return false, scheduler.JobStatusUnknown, fmt.Errorf("scheduler unavailable, cannot verify job %s", info.JobID)
		}
		st, err := sched.GetJobStatus(context.Background(), info.JobID)
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
