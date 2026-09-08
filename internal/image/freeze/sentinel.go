package freeze

import (
	"fmt"
	"os"
	"os/exec"
	"os/signal"
	"strings"
	"syscall"
)

// EnvSentinelWork is the env var MountedRun uses to hand the work script to
// the re-exec'd sentinel process; cmd's hidden `_mount_sentinel` command
// reads the same name. fuseBin, mnt, and fuseArgs travel as plain CLI args
// instead (os/exec rejects any env var value containing a NUL byte, which
// ruled out joining fuseArgs into one; args need no joining or escaping at
// all, so everything that isn't a multi-line script goes there).
const EnvSentinelWork = "CNT_MOUNT_WORK"

// RunSentinel is the body of the hidden `_mount_sentinel` command that
// MountedRun re-execs itself into, rather than calling `unshare` directly.
//
// It exists as a separate process for one reason: entering the user
// namespace (the next step here, via `unshare --user`) clears Pdeathsig on
// whatever process does it — the same kernel rule that clears it across a
// setuid exec, since becoming root-in-namespace is exactly that kind of
// privilege-elevating credential change. So nothing inside the namespace can
// ever be told "your supervisor just died" by the kernel. RunSentinel itself
// never escalates privilege, so its own Pdeathsig (SIGTERM, set by
// MountedRun) keeps working normally, and on receiving it, it kills its own
// process group by hand — taking the namespaced unshare/bash/FUSE tree down
// with it before it can be orphaned. See this package's README for the
// failure mode this closes (condatainer itself dying abruptly, e.g. an
// out-of-memory kill, mid-mount) and the one it doesn't (the sentinel itself
// being killed directly leaves the same gap, just in a far smaller and
// shorter-lived process than condatainer).
func RunSentinel(fuseBin, mnt, work string, fuseArgs []string) error {
	pgid, err := syscall.Getpgid(0)
	if err != nil {
		return fmt.Errorf("reading own process group: %w", err)
	}

	var quoted []string
	for _, a := range append(append([]string{}, fuseArgs...), mnt) {
		quoted = append(quoted, shellQuote(a))
	}
	// See mount.go's former inline copy of this script for why work runs in a
	// subshell and the cleanup at the end exists.
	script := fmt.Sprintf(`set -o pipefail
%s -f %s &
fpid=$!
for i in $(seq 1 50); do grep -qF ' %s ' /proc/mounts && break; sleep 0.1; done
grep -qF ' %s ' /proc/mounts || { echo "mount never appeared" >&2; kill $fpid 2>/dev/null; wait $fpid 2>/dev/null; exit 1; }
( %s )
rc=$?
kill $fpid 2>/dev/null
wait $fpid 2>/dev/null
exit $rc
`, shellQuote(fuseBin), strings.Join(quoted, " "), mnt, mnt, work)

	cmd := exec.Command("unshare", "--mount", "--user", "--map-root-user", "--", "/bin/bash", "-c", script)
	cmd.Stdin = os.Stdin
	cmd.Stdout = os.Stdout
	cmd.Stderr = os.Stderr
	// Join the sentinel's own group instead of starting a fresh one: both the
	// SIGTERM path below and MountedRun's own cmd.Cancel (a still-alive
	// condatainer choosing to cancel) kill by group, and need bash/FUSE in it.
	cmd.SysProcAttr = &syscall.SysProcAttr{Setpgid: true, Pgid: pgid}

	if err := cmd.Start(); err != nil {
		return fmt.Errorf("starting unshare: %w", err)
	}

	sigc := make(chan os.Signal, 1)
	signal.Notify(sigc, syscall.SIGTERM)
	defer signal.Stop(sigc)

	done := make(chan error, 1)
	go func() { done <- cmd.Wait() }()

	select {
	case <-sigc:
		// condatainer is gone; nothing further to report to. Best-effort:
		// take the whole group down, including this process.
		syscall.Kill(-pgid, syscall.SIGKILL) //nolint:errcheck
		<-done
		return fmt.Errorf("mount sentinel: condatainer exited before the mount finished")
	case err := <-done:
		return err
	}
}
