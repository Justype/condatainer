package freeze

import (
	"context"
	"fmt"
	"os"
	"os/exec"
	"syscall"

	"github.com/Justype/condatainer/internal/image/tool"
	execpkg "github.com/Justype/condatainer/internal/runtime/exec"
)

// executablePath resolves condatainer's own binary path, for MountedRun to
// re-exec into the sentinel. Overridden in tests: os.Executable() inside `go
// test` returns the test binary, which has no _mount_sentinel command.
var executablePath = os.Executable

// SetExecutablePathForTest points MountedRun's self re-exec at path instead
// of the real condatainer binary, and returns a func that restores the
// previous value. For tests outside this package (e.g. internal/build's,
// whose packFromScratchImage calls MountedRun too) that build
// testdata/mountharness as a stand-in condatainer binary; tests inside this
// package can set the unexported executablePath directly instead.
func SetExecutablePathForTest(path string) (restore func()) {
	old := executablePath
	executablePath = func() (string, error) { return path, nil }
	return func() { executablePath = old }
}

// MountedRun mounts fuseBin at mnt inside a fresh, unprivileged mount+user
// namespace, waits for the mount to appear, runs work (a bash script fragment
// that sees mnt as an ordinary directory), then ends the mount by killing the
// FUSE process. fuseArgs holds every flag the tool needs; mnt is appended.
// Exported so internal/build's own SquashFS packing can read a scratch .img
// the same apptainer-free way, not just this package's own freeze/unfreeze.
//
// `unshare --mount --user --map-root-user` is what lets an ordinary user call
// mount() at all: it maps namespace-uid 0 to the real caller, which is enough
// capability for a FUSE mount and nothing else — real root has no mapping
// inside it. That is also why the mount is never ended with umount/fusermount:
// those are normally setuid-root, and a setuid binary's owning uid (real root)
// isn't mapped in this namespace, so some kernels refuse to exec them at all
// rather than quietly ignoring the bit. Killing the foregrounded FUSE process
// tears down its own session instead, and the whole namespace — mount
// included — disappears the moment nothing is left running in it, so a killed
// or crashed run leaks nothing.
//
// MountedRun doesn't call unshare directly: it re-execs itself into the
// hidden `_mount_sentinel` command (RunSentinel) instead, so that something
// outside the namespace stays reachable if condatainer itself dies mid-mount.
// See RunSentinel's doc comment and this package's README for why.
func MountedRun(ctx context.Context, fuseBin string, fuseArgs []string, mnt string, work string, io execpkg.IO) error {
	if err := tool.CheckDependencies([]string{"unshare"}); err != nil {
		return err
	}

	self, err := executablePath()
	if err != nil {
		return fmt.Errorf("locating condatainer binary: %w", err)
	}

	sentinelArgs := append([]string{"_mount_sentinel", fuseBin, mnt}, fuseArgs...)
	cmd := exec.CommandContext(ctx, self, sentinelArgs...)
	cmd.Env = append(os.Environ(), EnvSentinelWork+"="+work)
	cmd.Stdin = io.Stdin
	cmd.Stdout = io.Stdout
	cmd.Stderr = io.Stderr

	// Setpgid: true (no Pgid) makes the sentinel its own group leader, and
	// RunSentinel joins bash to that same group — so a single group kill,
	// from either direction below, reaches the sentinel, bash, and the
	// backgrounded FUSE process together.
	//
	// Cancel covers a still-running condatainer choosing to cancel: it kills
	// the group directly, same as before. Pdeathsig covers the other
	// failure mode — condatainer dying with no chance to run any Go code at
	// all (SIGKILL, OOM-kill, crash): the kernel delivers SIGTERM to the
	// sentinel directly (Pdeathsig survives here because the sentinel never
	// enters the user namespace itself), and RunSentinel's own trap does the
	// group kill from the inside instead.
	cmd.SysProcAttr = &syscall.SysProcAttr{Setpgid: true, Pdeathsig: syscall.SIGTERM}
	cmd.Cancel = func() error {
		return syscall.Kill(-cmd.Process.Pid, syscall.SIGKILL)
	}

	return cmd.Run()
}
