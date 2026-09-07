// Package toolpath finds a runnable path for an external tool by name — the
// one place that decision is made. internal/libexec's self-provisioned
// toolchain is checked first and required over whatever the host has, when
// it exists at all; nothing else in this package or its caller decides that
// preference. A name libexec never provisions (debugfs, e2fsck, ...) simply
// isn't found there and falls through to PATH, then the FHS fallback
// directories, unaffected.
package toolpath

import (
	"context"
	"errors"
	"fmt"
	"os"
	"os/exec"
	"path/filepath"
	"strings"

	"github.com/Justype/condatainer/internal/libexec"
	"github.com/Justype/condatainer/internal/logging"
)

// ErrToolMissing reports that name could not be found anywhere Resolve looks.
var ErrToolMissing = errors.New("required tool not available")

// fhsFallbackDirs is Apptainer's own BinaryPath default search order for
// this class of external tool (apptainer's pkg/util/apptainerconf,
// unexported there): PATH, then the standard FHS directories. Matching it
// exactly means a tool findable by Apptainer's own mount is findable by
// CondaTainer's direct freeze/unfreeze path too, rather than a shorter,
// invented list.
var fhsFallbackDirs = []string{
	"/usr/local/sbin", "/usr/local/bin", "/usr/sbin", "/usr/bin", "/sbin", "/bin",
}

// Resolve finds name: libexec's self-provisioned copy first — required, not
// merely tried first among equals, so it wins even over a same-named binary
// already on PATH — then the caller's PATH, then the FHS directories
// directly, regardless of PATH. Returns ErrToolMissing, wrapping
// libexec.NotProvisionedMessage(name), if none of them has it — naming the
// `condatainer update --libexec` fix only when libexec.Provides(name), so a
// tool it never installs (debugfs, e2fsck, ...) gets a plain message instead
// of a fix that would not help.
//
// This is for a host-side invocation only, one with no bind mount to
// construct. A container-bound call (packing inside a build, a conda
// install) needs the containing directory too, so it uses
// libexec.ApptainerPath/MksquashfsPath/MicromambaPath directly instead.
func Resolve(name string) (string, error) {
	if p, ok := libexec.Path(name); ok {
		if info, err := os.Stat(p); err == nil && !info.IsDir() && info.Mode()&0111 != 0 {
			return p, nil
		}
	}
	if p, err := exec.LookPath(name); err == nil {
		return p, nil
	}
	for _, dir := range fhsFallbackDirs {
		p := filepath.Join(dir, name)
		if info, err := os.Stat(p); err == nil && !info.IsDir() && info.Mode()&0111 != 0 {
			return p, nil
		}
	}
	return "", fmt.Errorf("%w: %s", ErrToolMissing, NotFoundMessage(name))
}

// NotFoundMessage explains why name could not be resolved — libexec's own
// wording, re-exported here so a caller that already resolves through this
// package (unsquashfsBin's callers, CheckDependencies, findSquashfuse, ...)
// does not also need to import internal/libexec just to build one. Naming
// the `condatainer update --libexec` fix is libexec's call, not this
// package's — Provides(name) is what decides it, and that decision lives
// where the provisioning spec does.
func NotFoundMessage(name string) string {
	return libexec.NotProvisionedMessage(name)
}

// Command resolves name (see Resolve) and returns a context-bound *exec.Cmd
// ready to configure and run, with the invocation logged at debug level.
// The caller wires stdout/stderr and wraps a failure in whatever error
// shape its own domain needs — this package has no opinion on that.
func Command(ctx context.Context, name string, args ...string) (*exec.Cmd, error) {
	resolved, err := Resolve(name)
	if err != nil {
		return nil, err
	}
	logging.FromContext(ctx).Debug("running "+name, "args", strings.Join(args, " "))
	return exec.CommandContext(ctx, resolved, args...), nil
}
