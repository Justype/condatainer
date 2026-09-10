package libexec

import (
	"bytes"
	"context"
	"fmt"
	"os"
	"os/exec"
	"path/filepath"
	"regexp"
	"runtime"
	"strconv"
	"strings"

	"github.com/Justype/condatainer/internal/config"
	"github.com/Justype/condatainer/internal/logging"
	"github.com/Justype/condatainer/internal/utils"
)

// commandContext builds a context-bound command for a host-side tool
// invocation (micromamba, or a staged apptainer/mksquashfs/squashfuse being
// verified) — never apptainer-mediated, these all run directly on the host.
func commandContext(ctx context.Context, name string, args ...string) *exec.Cmd {
	return exec.CommandContext(ctx, name, args...)
}

// packages is the full self-provisioned toolchain, installed together in one
// solve. micromamba is included so the standalone bootstrap binary that
// creates the environment is immediately superseded by a tracked, updatable
// package inside it, the same way the official installer manages itself.
var packages = []string{"micromamba", "squashfs-tools", "squashfuse", "apptainer"}

// zstdFloor is the minimum apptainer version this toolchain must provide. See
// the README, Verification, for why this is independent of
// internal/runtime/apptainer.CheckZstdSupport despite the same threshold.
const zstdFloorMajor, zstdFloorMinor = 1, 4

// writableTarget picks where a fresh libexec should be created: the first
// writable tier, as a path that does not exist yet (see the README, Bootstrap
// sequence, for why create requires that).
func writableTarget() (string, error) {
	dir, err := config.GetWritableLibexecDir()
	if err != nil {
		return "", fmt.Errorf("no writable location for the toolchain: %w", err)
	}
	os.Remove(dir) //nolint:errcheck // best-effort; create fails loudly below if this didn't clear it
	return dir, nil
}

// Update refreshes the toolchain: create a fresh generation at a staging
// path, verify it, move the outgoing generation aside and activate the new
// one, then remove what was moved aside — all under the live generation's
// exclusive lock (see the README, Locking).
func Update(ctx context.Context) error {
	live, hadLive := Dir()

	var target string
	var liveLock *utils.FlockHandle
	if hadLive {
		lockPath := filepath.Join(live, lockFileName)
		if !utils.FileExists(lockPath) {
			// Self-heal a generation that predates this mechanism (see README).
			if f, err := utils.CreateFileWritable(lockPath); err == nil {
				f.Close()
			}
		}
		lock, err := utils.AcquireFlock(lockPath, true)
		if err != nil {
			return fmt.Errorf("condatainer is currently running (toolchain is locked); stop all running condatainer sessions before updating: %w", err)
		}
		liveLock = lock
		target = live
	} else {
		var err error
		target, err = writableTarget()
		if err != nil {
			return err
		}
	}
	defer func() {
		if liveLock != nil {
			liveLock.Close()
		}
	}()

	staging := target + ".new"
	os.RemoveAll(staging) // leftover from a previously-failed attempt

	// Reuse the live toolchain's own micromamba when there is one — it is a
	// separate binary from the staging path being created, so this is the
	// ordinary "binary at path A creates an environment at path B" case, not
	// the self-referential one. Only first-time provisioning has nothing to
	// borrow from.
	var mmBin string
	if hadLive {
		if p := filepath.Join(live, "bin", "micromamba"); utils.FileExists(p) {
			mmBin = p
		}
	}
	var cleanup func()
	if mmBin == "" {
		var err error
		mmBin, cleanup, err = downloadBootstrap(ctx)
		if err != nil {
			return err
		}
	}
	if cleanup != nil {
		defer cleanup()
	}

	if err := provisionAt(ctx, mmBin, staging); err != nil {
		os.RemoveAll(staging)
		return err
	}

	if err := verifyToolchain(ctx, staging); err != nil {
		os.RemoveAll(staging)
		return fmt.Errorf("newly provisioned toolchain failed verification, discarded: %w", err)
	}

	// old is where the outgoing generation is moved aside: a rename succeeds
	// unconditionally even while a file inside is open or executing, unlike
	// removal, so it never blocks activation on anything being in use.
	old := ""
	if hadLive {
		old = target + ".old"
		os.RemoveAll(old) // leftover from a previous Update that couldn't remove it
		if err := os.Rename(target, old); err != nil {
			return fmt.Errorf("failed to move the outgoing toolchain aside: %w", err)
		}
		// liveLock's own file just moved with the rest of old; release it
		// now rather than at the end, or removal below would be unlinking a
		// file this same process still has open.
		liveLock.Close()
		liveLock = nil
	}
	if err := os.Rename(staging, target); err != nil {
		return fmt.Errorf("failed to activate the new toolchain: %w", err)
	}
	if old != "" {
		if err := os.RemoveAll(old); err != nil {
			logging.FromContext(ctx).Warn("could not remove the outgoing toolchain generation; left on disk",
				"path", old, "err", err)
		}
	}
	return nil
}

// provisionAt runs the full bootstrap sequence into a not-yet-existing
// prefix: create, clean, strip conda-meta (see the README, Bootstrap
// sequence), then write this generation's own lock sentinel.
func provisionAt(ctx context.Context, mmBin, prefix string) error {
	if err := createAt(ctx, mmBin, prefix); err != nil {
		return fmt.Errorf("failed to provision the toolchain: %w", err)
	}
	if err := cleanAt(ctx, mmBin, prefix); err != nil {
		return fmt.Errorf("failed to clean the toolchain cache: %w", err)
	}
	if err := os.RemoveAll(filepath.Join(prefix, "conda-meta")); err != nil {
		return fmt.Errorf("failed to remove conda-meta: %w", err)
	}
	// micromamba creates every file under prefix itself, bypassing this
	// codebase's own MkdirAllShared/CreateFileWritable entirely, so nothing
	// below the top level picks up "2775" group-write on its own — confirmed
	// against a real provisioned tree, where only prefix itself came out
	// group-writable and every file and directory under it stayed at
	// micromamba's own umask-derived mode. Sharing the whole tree here is
	// what makes Update's later os.RemoveAll (a different group member than
	// whoever provisioned this generation) actually able to remove it.
	if err := utils.ShareTreeWithParentGroup(prefix); err != nil {
		return fmt.Errorf("failed to share the toolchain with the group: %w", err)
	}
	lockFile, err := utils.CreateFileWritable(filepath.Join(prefix, lockFileName))
	if err != nil {
		return fmt.Errorf("failed to create the lock sentinel: %w", err)
	}
	lockFile.Close()
	return nil
}

// createAt runs `micromamba create` into prefix from a binary that lives
// outside it (see the README, Bootstrap sequence, for the exact flags).
func createAt(ctx context.Context, mmBin, prefix string) error {
	args := append([]string{"-r", prefix, "--no-rc", "create", "-y", "-p", prefix, "-c", "conda-forge"}, packages...)
	return runMicromamba(ctx, mmBin, args...)
}

// cleanAt reclaims disk after provisioning (see the README, Bootstrap
// sequence, for why -f/--force-pkgs-dirs is required).
func cleanAt(ctx context.Context, mmBin, prefix string) error {
	return runMicromamba(ctx, mmBin, "-r", prefix, "--no-rc", "clean", "-a", "-f", "-y")
}

// runMicromamba runs one micromamba invocation with the caller's ambient
// conda/mamba env vars scrubbed, so an inherited CONDA_PREFIX or
// MAMBA_ROOT_PREFIX from wherever condatainer itself is running cannot leak
// into the environment being provisioned.
func runMicromamba(ctx context.Context, mmBin string, args ...string) error {
	managed := map[string]bool{
		"CONDA_PREFIX": true, "CONDARC": true, "MAMBA_NO_RC": true,
		"MAMBA_ROOT_PREFIX": true, "MAMBA_TARGET_PREFIX": true,
	}
	env := make([]string, 0, len(os.Environ()))
	for _, entry := range os.Environ() {
		name, _, _ := strings.Cut(entry, "=")
		if !managed[name] {
			env = append(env, entry)
		}
	}

	cmd := commandContext(ctx, mmBin, args...)
	cmd.Env = env
	var out bytes.Buffer
	cmd.Stdout = &out
	cmd.Stderr = &out
	if err := cmd.Run(); err != nil {
		return fmt.Errorf("%s %s: %w: %s", mmBin, strings.Join(args, " "), err, strings.TrimSpace(out.String()))
	}
	return nil
}

// verifyToolchain sanity-checks a staged (not yet live) toolchain before it
// is allowed to replace the running one: each tool must run (checked by
// output content, not exit status — see the README, Verification), and the
// staged apptainer must clear the zstd floor.
func verifyToolchain(ctx context.Context, prefix string) error {
	for _, name := range []string{"apptainer", "mksquashfs", "squashfuse"} {
		path := filepath.Join(prefix, "bin", name)
		out, _ := commandContext(ctx, path, "--version").CombinedOutput()
		if !strings.Contains(strings.ToLower(string(out)), name) {
			return fmt.Errorf("%s does not run: unexpected output: %s", name, strings.TrimSpace(string(out)))
		}
	}

	out, err := commandContext(ctx, filepath.Join(prefix, "bin", "apptainer"), "--version").Output()
	if err != nil {
		return fmt.Errorf("could not read apptainer's version: %w", err)
	}
	version := parseVersion(string(out))
	if version == "" {
		return fmt.Errorf("could not parse apptainer's version from %q", strings.TrimSpace(string(out)))
	}
	if !meetsZstdFloor(version) {
		return fmt.Errorf("provisioned apptainer %s is below the required %d.%d (zstd support)",
			version, zstdFloorMajor, zstdFloorMinor)
	}
	return nil
}

var versionPattern = regexp.MustCompile(`(\d+)\.(\d+)(\.\d+)?`)

func parseVersion(output string) string {
	return versionPattern.FindString(strings.TrimSpace(output))
}

// ToolVersion is one provisioned binary's name and the version its own
// --version reported.
type ToolVersion struct {
	Name    string
	Version string
}

// versionedBins is the one binary per package this package installs whose
// version is worth reporting — not unsquashfs/squashfuse_ll, since each
// always shares its sibling's version.
var versionedBins = []string{"apptainer", "mksquashfs", "squashfuse", "micromamba"}

// Versions runs --version against the live toolchain's own binaries. A
// binary that fails to run or prints nothing parseVersion recognizes reports
// "unknown" rather than failing the rest.
func Versions(ctx context.Context) ([]ToolVersion, error) {
	dir, ok := Dir()
	if !ok {
		return nil, ErrNotProvisioned
	}
	versions := make([]ToolVersion, 0, len(versionedBins))
	for _, name := range versionedBins {
		out, _ := commandContext(ctx, filepath.Join(dir, "bin", name), "--version").CombinedOutput()
		version := parseVersion(string(out))
		if version == "" {
			version = "unknown"
		}
		versions = append(versions, ToolVersion{Name: name, Version: version})
	}
	return versions, nil
}

func meetsZstdFloor(version string) bool {
	parts := strings.Split(version, ".")
	if len(parts) < 2 {
		return false
	}
	major, _ := strconv.Atoi(parts[0])
	minor, _ := strconv.Atoi(parts[1])
	return major > zstdFloorMajor || (major == zstdFloorMajor && minor >= zstdFloorMinor)
}

// downloadBootstrap fetches a standalone micromamba binary to a transient
// location outside any libexec prefix. The caller must call the returned
// cleanup once done with it.
func downloadBootstrap(ctx context.Context) (path string, cleanup func(), err error) {
	asset, err := micromambaAssetName()
	if err != nil {
		return "", nil, err
	}

	root := utils.GetTmpDir()
	if err := utils.EnsureTmpSubdir(root); err != nil {
		return "", nil, fmt.Errorf("failed to create tmp dir %s: %w", root, err)
	}
	dir := filepath.Join(root, "libexec-bootstrap")
	if err := utils.MkdirAllShared(dir); err != nil {
		return "", nil, fmt.Errorf("failed to create bootstrap dir %s: %w", dir, err)
	}

	binPath := filepath.Join(dir, "micromamba")
	url := "https://github.com/mamba-org/micromamba-releases/releases/latest/download/" + asset
	if err := utils.DownloadExecutable(ctx, url, binPath); err != nil {
		os.RemoveAll(dir)
		return "", nil, fmt.Errorf("failed to download micromamba: %w", err)
	}
	return binPath, func() { os.RemoveAll(dir) }, nil
}

// micromambaAssetName maps the running OS/arch to the standalone release
// asset name mamba-org/micromamba-releases publishes. Apptainer itself is
// Linux-only, so this tool has no reason to support anything else.
func micromambaAssetName() (string, error) {
	if runtime.GOOS != "linux" {
		return "", fmt.Errorf("the self-provisioned toolchain is only available on linux (GOOS=%s)", runtime.GOOS)
	}
	switch runtime.GOARCH {
	case "amd64":
		return "micromamba-linux-64", nil
	case "arm64":
		return "micromamba-linux-aarch64", nil
	default:
		return "", fmt.Errorf("no micromamba release for linux/%s", runtime.GOARCH)
	}
}
