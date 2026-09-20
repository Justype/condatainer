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
	"sync"

	"github.com/Justype/condatainer/internal/config"
	"github.com/Justype/condatainer/internal/utils"
)

// commandContext builds a context-bound command for a host-side tool
// invocation (micromamba, or a staged apptainer/mksquashfs/squashfuse being
// verified) — never apptainer-mediated, these all run directly on the host.
var ensureMu sync.Mutex

func commandContext(ctx context.Context, name string, args ...string) *exec.Cmd {
	return exec.CommandContext(ctx, name, args...)
}

// pkg is one conda package this toolchain can hold.
type pkg struct {
	name   string   // conda package name
	bins   []string // binaries it installs; bins[0] marks it installed
	verify []string // binaries whose version output must name themselves
	report string   // binary whose version is reported and floor-checked ("" = none)

	floorMajor, floorMinor int    // minimum version of report; 0 = no floor
	floorNeed              string // what the floor buys, for a failure message
}

// packages lists what a prefix may hold, in install order. micromamba is the
// base every prefix carries: it installs and updates the rest, and its binary
// marks a tier as provisioned.
var packages = []pkg{
	{name: "micromamba", bins: []string{"micromamba"}, report: "micromamba"},
	{
		name: "squashfs-tools", bins: []string{"mksquashfs", "unsquashfs"},
		verify: []string{"mksquashfs", "unsquashfs"}, report: "mksquashfs",
		floorMajor: 4, floorMinor: 4, floorNeed: "zstd and -offset support",
	},
	{name: "squashfuse", bins: []string{"squashfuse", "squashfuse_ll"}, verify: []string{"squashfuse"}, report: "squashfuse"},
	{
		name: "apptainer", bins: []string{"apptainer"}, verify: []string{"apptainer"}, report: "apptainer",
		floorMajor: 1, floorMinor: 4, floorNeed: "zstd support",
	},
}

// packageFor finds the package that is called name or installs a binary called name.
func packageFor(name string) (pkg, bool) {
	for _, p := range packages {
		if p.name == name {
			return p, true
		}
		for _, bin := range p.bins {
			if bin == name {
				return p, true
			}
		}
	}
	return pkg{}, false
}

// installedAt reports whether p's binary is present in prefix's bin/.
func (p pkg) installedAt(prefix string) bool {
	return utils.FileExists(filepath.Join(prefix, "bin", p.bins[0]))
}

// installedPackages returns the packages present in prefix, in install order.
func installedPackages(prefix string) []pkg {
	var have []pkg
	for _, p := range packages {
		if p.installedAt(prefix) {
			have = append(have, p)
		}
	}
	return have
}

// versionFlag is the flag a provisioned binary prints its own version for.
// squashfs-tools (mksquashfs, unsquashfs) predates GNU-style long options and
// only recognizes the single-dash form; everything else this package
// provisions takes "--version".
func versionFlag(name string) string {
	switch name {
	case "mksquashfs", "unsquashfs":
		return "-version"
	default:
		return "--version"
	}
}

// libexecLockName is the sibling lock that serializes Update calls on a tier,
// including the first one, when there is no prefix yet to hold a lock.
const libexecLockName = ".libexec.lock"

// resolveTools maps requested package names to packages.
func resolveTools(names []string) ([]pkg, error) {
	var out []pkg
	for _, name := range names {
		var found bool
		for _, p := range packages {
			if p.name == name {
				out = append(out, p)
				found = true
				break
			}
		}
		if !found {
			valid := make([]string, len(packages))
			for i, p := range packages {
				valid[i] = p.name
			}
			return nil, fmt.Errorf("unknown toolchain package %q (available: %s)", name, strings.Join(valid, ", "))
		}
	}
	return out, nil
}

// Update creates or updates the toolchain prefix in place, under an exclusive
// lock. With no tools it updates what is installed, or creates a prefix holding
// micromamba on a tier with none; naming tools also installs any that are
// missing. A prefix without conda-meta/ is recreated.
func Update(ctx context.Context, tools ...string) error {
	return run(ctx, tools, false, false)
}

// EnsureMicromamba creates a toolchain holding micromamba alone when no tier
// has one, and does nothing otherwise. Callers are serialized in this process
// and, through the tier lock, across processes: one that loses the race finds
// the toolchain provisioned and leaves it alone.
func EnsureMicromamba(ctx context.Context) error {
	ensureMu.Lock()
	defer ensureMu.Unlock()
	if Installed("micromamba") {
		return nil
	}
	return run(ctx, nil, false, true)
}

// Sync completes the toolchain: it creates the prefix if there is none, and
// otherwise installs any package in install that is missing and updates every
// installed package.
func Sync(ctx context.Context, install ...string) error {
	return run(ctx, install, true, false)
}

// run is Update and Sync: updateAll updates everything installed rather than
// only the named packages. ifMissing returns once the tier lock is held if
// micromamba is by then installed.
func run(ctx context.Context, tools []string, updateAll, ifMissing bool) error {
	requested, err := resolveTools(tools)
	if err != nil {
		return err
	}

	live, hadLive := Dir()
	target := live
	if !hadLive {
		dir, err := config.GetWritableLibexecDir()
		if err != nil {
			return fmt.Errorf("no writable location for the toolchain: %w", err)
		}
		target = dir
	}
	parent := filepath.Dir(target)
	if err := utils.MkdirAllShared(parent); err != nil {
		return fmt.Errorf("failed to create %s: %w", parent, err)
	}

	tierLock, err := acquireTierLock(filepath.Join(parent, libexecLockName))
	if err != nil {
		return err
	}
	defer tierLock.Close()

	if ifMissing && Installed("micromamba") {
		return nil
	}

	var useLock *utils.FileLock
	if hadLive {
		lockPath := filepath.Join(target, lockFileName)
		if !utils.FileExists(lockPath) {
			if f, err := utils.CreateFileWritable(lockPath); err == nil {
				f.Close()
			}
		}
		useLock, err = utils.AcquireFileLock(lockPath, true)
		if err != nil {
			return fmt.Errorf("condatainer is currently running (toolchain is locked); stop all running condatainer sessions before updating: %w", err)
		}
		defer func() {
			if useLock != nil {
				useLock.Close()
			}
		}()
	}

	inPlace := hadLive && utils.DirExists(filepath.Join(target, "conda-meta"))
	if inPlace {
		return updateInPlace(ctx, target, requested, updateAll)
	}

	// A prefix with no conda-meta/, or a half-built one, cannot be updated by
	// micromamba; recreate it holding the same tools plus any requested.
	want := requested
	if hadLive {
		want = mergePackages(installedPackages(target), requested)
		// The in-prefix lock file goes with the directory, so release it first.
		useLock.Close()
		useLock = nil
	}
	if err := utils.RemoveAllWritable(target); err != nil {
		return fmt.Errorf("failed to remove the old toolchain: %w", err)
	}
	return createPrefix(ctx, target, want)
}

// mergePackages returns the packages in a and b, once each, in install order.
func mergePackages(a, b []pkg) []pkg {
	var out []pkg
	for _, p := range packages {
		for _, list := range [][]pkg{a, b} {
			if containsPkg(list, p) {
				out = append(out, p)
				break
			}
		}
	}
	return out
}

func containsPkg(list []pkg, p pkg) bool {
	for _, q := range list {
		if q.name == p.name {
			return true
		}
	}
	return false
}

// acquireTierLock takes the tier's exclusive Update lock, creating its file.
func acquireTierLock(path string) (*utils.FileLock, error) {
	if !utils.FileExists(path) {
		if f, err := utils.CreateFileWritable(path); err == nil {
			f.Close()
		}
	}
	lock, err := utils.AcquireFileLock(path, true)
	if err != nil {
		return nil, fmt.Errorf("another condatainer update of the toolchain is running: %w", err)
	}
	return lock, nil
}

// createPrefix builds a new prefix at target (which must not exist) holding
// micromamba plus want, from a downloaded bootstrap binary, and removes the
// prefix again if any step fails so nothing half-built looks provisioned.
func createPrefix(ctx context.Context, target string, want []pkg) error {
	mmBin, cleanup, err := downloadBootstrap(ctx)
	if err != nil {
		return err
	}
	defer cleanup()

	names := packageNames(mergePackages(want, packages[:1]))
	args := append([]string{"create", "-y", "-p", target, "-c", "conda-forge"}, names...)
	if err := runMicromamba(ctx, mmBin, target, args...); err != nil {
		utils.RemoveAllWritable(target) //nolint:errcheck
		return fmt.Errorf("failed to provision the toolchain: %w", err)
	}
	if err := finishPrefix(ctx, mmBin, target); err != nil {
		utils.RemoveAllWritable(target) //nolint:errcheck
		return err
	}
	return nil
}

// updateInPlace installs any requested package the prefix lacks, then updates
// everything installed when updateAll or none were named, otherwise only the
// requested packages, using the prefix's own micromamba.
func updateInPlace(ctx context.Context, target string, requested []pkg, updateAll bool) error {
	mmBin := filepath.Join(target, "bin", "micromamba")
	if !utils.FileExists(mmBin) {
		return fmt.Errorf("%s is missing; remove %s and run `condatainer update --libexec`", mmBin, target)
	}

	var missing []pkg
	for _, p := range requested {
		if !p.installedAt(target) {
			missing = append(missing, p)
		}
	}
	if len(missing) > 0 {
		args := append([]string{"install", "-y", "-p", target, "-c", "conda-forge"}, packageNames(missing)...)
		if err := runMicromamba(ctx, mmBin, target, args...); err != nil {
			return fmt.Errorf("failed to install into the toolchain: %w", err)
		}
	}

	args := []string{"update", "-y", "-p", target, "-c", "conda-forge"}
	if updateAll || len(requested) == 0 {
		args = append(args, "-a")
	} else {
		args = append(args, packageNames(requested)...)
	}
	if err := runMicromamba(ctx, mmBin, target, args...); err != nil {
		return fmt.Errorf("failed to update the toolchain; run `condatainer update --libexec` again: %w", err)
	}
	return finishPrefix(ctx, mmBin, target)
}

func packageNames(list []pkg) []string {
	names := make([]string, len(list))
	for i, p := range list {
		names[i] = p.name
	}
	return names
}

// finishPrefix runs after a successful create or update: link fusermount3,
// clean the package cache, share the tree with the group, ensure the lock
// sentinel, and verify what is installed.
func finishPrefix(ctx context.Context, mmBin, prefix string) error {
	if err := ensureFusermount3InBin(prefix); err != nil {
		return fmt.Errorf("failed to link fusermount3: %w", err)
	}
	if err := runMicromamba(ctx, mmBin, prefix, "clean", "-a", "-f", "-y"); err != nil {
		return fmt.Errorf("failed to clean the toolchain cache: %w", err)
	}
	os.Remove(filepath.Join(prefix, ".cache")) //nolint:errcheck // an empty leftover of the contained cache
	// micromamba creates every file under prefix itself, bypassing this
	// codebase's own MkdirAllShared/CreateFileWritable, so nothing below the
	// top level picks up group-write on its own.
	if err := utils.ShareTreeWithParentGroup(prefix); err != nil {
		return fmt.Errorf("failed to share the toolchain with the group: %w", err)
	}
	lockPath := filepath.Join(prefix, lockFileName)
	if !utils.FileExists(lockPath) {
		lockFile, err := utils.CreateFileWritable(lockPath)
		if err != nil {
			return fmt.Errorf("failed to create the lock sentinel: %w", err)
		}
		lockFile.Close()
	}
	if err := verifyToolchain(ctx, prefix); err != nil {
		return fmt.Errorf("the toolchain failed verification; run `condatainer update --libexec` again: %w", err)
	}
	return nil
}

// ensureFusermount3InBin symlinks bin/fusermount3 to ../sbin/fusermount3 when
// the fuse3 package installed the setuid helper there instead of bin/ — some
// builds do. libfuse3 bakes an absolute <prefix>/bin/fusermount3 path into
// itself at install time; if the real binary only exists in sbin/, that baked
// path never resolves regardless of $PATH.
func ensureFusermount3InBin(prefix string) error {
	sbinPath := filepath.Join(prefix, "sbin", "fusermount3")
	if !utils.FileExists(sbinPath) {
		return nil // this build put it in bin/ already, or doesn't ship it
	}
	binPath := filepath.Join(prefix, "bin", "fusermount3")
	if utils.FileExists(binPath) {
		return nil
	}
	return os.Symlink(filepath.Join("..", "sbin", "fusermount3"), binPath)
}

// runMicromamba runs one micromamba subcommand rooted at prefix, with the
// package cache, repodata cache and HOME kept inside or apart from it so
// nothing lands in the invoking user's home, and the caller's ambient
// conda/mamba variables scrubbed.
func runMicromamba(ctx context.Context, mmBin, prefix string, sub ...string) error {
	home, err := os.MkdirTemp("", "cnt-libexec-home-")
	if err != nil {
		return fmt.Errorf("failed to create a scratch home: %w", err)
	}
	defer os.RemoveAll(home)

	managed := map[string]bool{
		"CONDA_PREFIX": true, "CONDARC": true, "MAMBA_NO_RC": true,
		"MAMBA_ROOT_PREFIX": true, "MAMBA_TARGET_PREFIX": true,
		"CONDA_PKGS_DIRS": true, "MAMBA_PKGS_DIRS": true,
		"XDG_CACHE_HOME": true, "HOME": true,
	}
	env := make([]string, 0, len(os.Environ())+3)
	for _, entry := range os.Environ() {
		name, _, _ := strings.Cut(entry, "=")
		if !managed[name] {
			env = append(env, entry)
		}
	}
	env = append(env,
		"CONDA_PKGS_DIRS="+filepath.Join(prefix, "pkgs"),
		"XDG_CACHE_HOME="+filepath.Join(prefix, ".cache"),
		"HOME="+home,
	)

	args := append([]string{"-r", prefix, "--no-rc"}, sub...)
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

// verifyToolchain sanity-checks the packages installed in prefix: each listed
// binary must run (checked by output content, not exit status — see the
// README, Verification) and each floor must hold.
func verifyToolchain(ctx context.Context, prefix string) error {
	for _, p := range installedPackages(prefix) {
		for _, name := range p.verify {
			path := filepath.Join(prefix, "bin", name)
			out, _ := commandContext(ctx, path, versionFlag(name)).CombinedOutput()
			if !strings.Contains(strings.ToLower(string(out)), name) {
				return fmt.Errorf("%s does not run: unexpected output: %s", name, strings.TrimSpace(string(out)))
			}
		}
		if p.floorMajor > 0 {
			if err := verifyFloor(ctx, prefix, p.report, p.floorMajor, p.floorMinor, p.floorNeed); err != nil {
				return err
			}
		}
	}
	return nil
}

// verifyFloor reads name's own version from prefix and checks it against
// floorMajor.floorMinor, naming what that floor buys in a failure.
func verifyFloor(ctx context.Context, prefix, name string, floorMajor, floorMinor int, need string) error {
	return checkFloorAt(ctx, filepath.Join(prefix, "bin", name), name, floorMajor, floorMinor, need)
}

// CheckFloor checks the binary at path against the minimum version of the
// package that provides name, for a binary found outside this toolchain. A name
// with no floor passes.
func CheckFloor(ctx context.Context, path, name string) error {
	p, ok := packageFor(name)
	if !ok || p.floorMajor == 0 {
		return nil
	}
	return checkFloorAt(ctx, path, name, p.floorMajor, p.floorMinor, p.floorNeed)
}

func checkFloorAt(ctx context.Context, path, name string, floorMajor, floorMinor int, need string) error {
	// Exit status is ignored: unsquashfs prints its version and exits 1.
	out, runErr := commandContext(ctx, path, versionFlag(name)).CombinedOutput()
	version := parseVersion(string(out))
	if version == "" {
		if runErr != nil {
			return fmt.Errorf("could not read %s's version: %w", name, runErr)
		}
		return fmt.Errorf("could not parse %s's version from %q", name, strings.TrimSpace(string(out)))
	}
	if !meetsFloor(version, floorMajor, floorMinor) {
		return fmt.Errorf("%s %s is below the required %d.%d (%s)", name, version, floorMajor, floorMinor, need)
	}
	return nil
}

var versionPattern = regexp.MustCompile(`(\d+)\.(\d+)(\.\d+)?`)

func parseVersion(output string) string {
	return versionPattern.FindString(strings.TrimSpace(output))
}

// ToolVersion is one provisioned binary's name and the version its own
// version flag (see versionFlag) reported.
type ToolVersion struct {
	Name    string
	Version string
}

// Versions runs each installed package's reported binary's own version flag
// (see versionFlag) against the live toolchain. A binary that fails to run or
// prints nothing parseVersion recognizes reports "unknown" rather than failing
// the rest.
func Versions(ctx context.Context) ([]ToolVersion, error) {
	dir, ok := Dir()
	if !ok {
		return nil, ErrNotProvisioned
	}
	var versions []ToolVersion
	for _, p := range installedPackages(dir) {
		if p.report == "" {
			continue
		}
		out, _ := commandContext(ctx, filepath.Join(dir, "bin", p.report), versionFlag(p.report)).CombinedOutput()
		version := parseVersion(string(out))
		if version == "" {
			version = "unknown"
		}
		versions = append(versions, ToolVersion{Name: p.report, Version: version})
	}
	return versions, nil
}

// meetsFloor reports whether version is at or above floorMajor.floorMinor,
// shared by every version-floor check this package makes (apptainer's zstd
// floor, squashfs-tools' -offset floor).
func meetsFloor(version string, floorMajor, floorMinor int) bool {
	parts := strings.Split(version, ".")
	if len(parts) < 2 {
		return false
	}
	major, _ := strconv.Atoi(parts[0])
	minor, _ := strconv.Atoi(parts[1])
	return major > floorMajor || (major == floorMajor && minor >= floorMinor)
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
