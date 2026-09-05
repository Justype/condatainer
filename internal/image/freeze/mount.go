package freeze

import (
	"context"
	"fmt"
	"os/exec"
	"strings"

	"github.com/Justype/condatainer/internal/image/tool"
	execpkg "github.com/Justype/condatainer/internal/runtime/exec"
)

// mountedRun mounts fuseBin at mnt inside a fresh, unprivileged mount+user
// namespace, waits for the mount to appear, runs work (a bash script fragment
// that sees mnt as an ordinary directory), then ends the mount by killing the
// FUSE process. fuseArgs holds every flag the tool needs; mnt is appended.
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
func mountedRun(ctx context.Context, fuseBin string, fuseArgs []string, mnt string, work string, io execpkg.IO) error {
	if err := tool.CheckDependencies([]string{"unshare"}); err != nil {
		return err
	}

	var quoted []string
	for _, a := range append(append([]string{}, fuseArgs...), mnt) {
		quoted = append(quoted, shellQuote(a))
	}
	// work runs in a subshell: it commonly sets its own `set -e` (every pack and
	// mke2fs script does), and without the subshell that setting would outlive
	// the substitution point and apply to the kill/wait cleanup below — where
	// wait's exit status (the just-killed FUSE daemon's own, often non-zero for
	// SIGTERM) would then abort the script silently, before it ever reaches the
	// exit $rc that was supposed to report the real result.
	script := fmt.Sprintf(`set -o pipefail
%s -f %s &
fpid=$!
for i in $(seq 1 50); do grep -qF ' %s ' /proc/mounts && break; sleep 0.1; done
grep -qF ' %s ' /proc/mounts || { echo "mount never appeared" >&2; exit 1; }
( %s )
rc=$?
kill $fpid 2>/dev/null
wait $fpid 2>/dev/null
exit $rc
`, shellQuote(fuseBin), strings.Join(quoted, " "), mnt, mnt, work)

	cmd := exec.CommandContext(ctx, "unshare", "--mount", "--user", "--map-root-user", "--", "/bin/bash", "-c", script)
	cmd.Stdin = io.Stdin
	cmd.Stdout = io.Stdout
	cmd.Stderr = io.Stderr
	return cmd.Run()
}
