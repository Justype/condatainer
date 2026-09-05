package build

import (
	"context"
	"errors"
	"fmt"
	"os"
	"os/exec"
	"path/filepath"
	"strings"

	"log/slog"

	"github.com/Justype/condatainer/internal/artifact/compare"
	"github.com/Justype/condatainer/internal/config"
	"github.com/Justype/condatainer/internal/image"
	"github.com/Justype/condatainer/internal/image/producer"
	"github.com/Justype/condatainer/internal/logging"
	"github.com/Justype/condatainer/internal/runtime/container"
	"github.com/Justype/condatainer/internal/scheduler"
	"github.com/Justype/condatainer/internal/store"
	"github.com/Justype/condatainer/internal/utils"
)

// cachedInstalledOverlays is the set of "name/version" strings for all overlays found
// across all image search paths. Nil means the cache is cold or has been invalidated.
var cachedInstalledOverlays map[string]bool

// Kept as the build package's spelling for tests and documentation; producer
// owns the value and the path construction.
const preparedSuffix = producer.PreparedSuffix

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

func invalidateInstalledOverlays() {
	cachedInstalledOverlays = nil
	container.InvalidateInstalledOverlaysCache()
}

// checkShouldBuild returns (skip=true, nil) if the overlay already exists and update=false.
// In update mode, if the overlay is locked by a running container, returns an error.
//
// storeOverflow never skips: the point of it is to produce a second identity of
// a name that is already installed, and whether this build is that second one
// cannot be known until it has been built and its keys read.
func checkShouldBuild(b *BuildObject) (skip bool, err error) {
	// IsInstalled, not a stat of the target: for a base that means every image
	// search path, so one supplied by a shared install is not rebuilt into the
	// user's own directory. Identical to a stat for every other type.
	if !b.update && !b.storeOverflow && b.IsInstalled() {
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

// preparedPathFor derives where a build writes its output: beside the target,
// tagged with the lock owner so stale-lock cleanup can recompute it.
// See the README's Prepared output.
func preparedPathFor(targetPath string, info BuildLockInfo) string {
	return producer.PreparedPath(targetPath, info)
}

// ensureWorkspaceRoot creates the selected scratch root with its configured
// permissions, then the nested producer-private directory beneath it. The scratch
// root is also what a definition build hands apptainer as APPTAINER_TMPDIR, so
// creating it here is what makes that directory exist.
func ensureWorkspaceRoot(b *BuildObject) error {
	if err := utils.EnsureTmpSubdir(b.ws.BaseRoot); err != nil {
		return fmt.Errorf("failed to create tmp dir %s: %w", b.ws.BaseRoot, err)
	}
	if err := utils.MkdirAllShared(b.ws.Root); err != nil {
		return fmt.Errorf("failed to create producer workspace %s: %w", b.ws.Root, err)
	}
	return nil
}

// publish installs a finished build and reports where it landed.
//
// Ordinarily that is the flat target and one rename. Under storeOverflow the
// artifact is filed by its identity instead, which is what lets a second build
// of a name exist beside the first rather than being skipped.
func (b *BuildObject) publish(preparedPath, targetPath string) (string, error) {
	if !b.storeOverflow {
		if err := atomicInstall(preparedPath, targetPath); err != nil {
			return "", err
		}
		return targetPath, nil
	}
	return b.publishToStore(preparedPath, targetPath)
}

// publishToStore files the finished build under its own identity.
//
// StoreOnly is not a choice here. The build holds the producer lock on the flat
// target for its whole duration (Target.Lock is Path + ".lock"), which is the
// very lock Begin would take to claim the bare name — so letting it consider
// flat placement would deadlock the build against itself. The store is where a
// second identity belongs anyway; `store use` is how one becomes the default.
//
// Keys are regenerated from the packed output rather than read from what the
// build believes it produced, so Commit refuses anything that does not match
// the address it reserved.
func (b *BuildObject) publishToStore(preparedPath, targetPath string) (string, error) {
	artifact, err := compare.Read(preparedPath)
	if err != nil {
		os.Remove(preparedPath) //nolint:errcheck
		return "", fmt.Errorf("cannot file %s by identity: %w", filepath.Base(targetPath), err)
	}

	tx, err := store.Begin(artifact.Name, artifact.IdentityRef(), store.BeginOptions{
		ImagesDir: filepath.Dir(targetPath),
		Equiv:     artifact.EquivRef(),
		StoreOnly: true,
	})
	if err != nil {
		os.Remove(preparedPath) //nolint:errcheck
		return "", err
	}
	// An identical build is already installed. It cost the build but not the
	// disk, and the caller is told where the artifact actually is.
	if tx.Adopted != nil {
		os.Remove(preparedPath) //nolint:errcheck
		return tx.Adopted.Path, nil
	}
	defer tx.Abort()

	// Both paths are in one images root, so this is a rename and not a copy.
	if err := os.Rename(preparedPath, tx.Prepared); err != nil {
		os.Remove(preparedPath) //nolint:errcheck
		return "", fmt.Errorf("failed to stage %s in the store: %w", filepath.Base(targetPath), err)
	}
	candidate, err := tx.Commit()
	if err != nil {
		return "", err
	}
	invalidateInstalledOverlays()
	return candidate.Path, nil
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
	invalidateInstalledOverlays()
	return nil
}

// prepareBuildWorkspace creates the build workspace: an ext3 scratch image, or
// host directories. Script and Conda builds only. A stale workspace is warned
// about and re-created, leaving a fetched build source intact.
func prepareBuildWorkspace(ctx context.Context, b *BuildObject) error {
	if err := ensureWorkspaceRoot(b); err != nil {
		return err
	}
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
	return producer.ShortHostname()
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

// EffectiveResourceSpec resolves the resources a build runs with, in priority
// order:
//
//	buildDefaults → scriptSpecs.Spec (when HasDirectives=true) → scheduler job resources
func EffectiveResourceSpec(specs *scheduler.ScriptSpecs) *scheduler.ResourceSpec {
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
	return producer.Acquire(path, info)
}

// overwriteBuildLockFile overwrites an existing lock file with new JSON metadata.
// The caller must already hold the lock (i.e. have created it via acquireBuildLockFile).
func overwriteBuildLockFile(path string, info BuildLockInfo) error {
	return producer.Overwrite(path, info)
}

// readBuildLockFile reads and parses a lock file at the given path.
func readBuildLockFile(path string) (BuildLockInfo, error) {
	return producer.Read(path)
}

// isBuildLockStale returns whether the lock is stale, the job's current status, and any
// uncertainty error. Returns (true, Unknown, nil) when definitely stale, (false, status, nil)
// when definitely alive, or (false, Unknown, err) when the state cannot be verified.
func isBuildLockStale(info BuildLockInfo) (stale bool, status scheduler.JobStatus, err error) {
	return producer.IsStale(info)
}
