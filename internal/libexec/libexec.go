// Package libexec resolves and provisions CondaTainer's self-provisioned
// toolchain — mksquashfs, squashfuse, and an ordinary (non-fakeroot)
// apptainer, installed via micromamba into one of the four data-directory
// tiers. See the package README for the design.
package libexec

import (
	"errors"
	"fmt"
	"os"
	"path/filepath"

	"github.com/Justype/condatainer/internal/config"
	"github.com/Justype/condatainer/internal/utils"
)

// binMarker is the executable that must exist for a tier's libexec/ to count
// as provisioned, rather than an empty or half-written directory.
const binMarker = "apptainer"

// ErrNotProvisioned is what a container-bound caller returns when nothing is
// provisioned at any tier. Condatainer never auto-provisions on its own — an
// ordinary exec/run or a build must never trigger a first-time download — so
// every message directs the user to the one command that does.
var ErrNotProvisioned = errors.New("the self-provisioned toolchain is not installed; run `condatainer update --libexec` first")

// providedBins are the exact binary names provision.go's package spec
// installs. Path/BinDir build a path for any name regardless of whether
// this package actually installs it (see Path's own doc comment — that is
// what lets toolpath.Resolve ask about names like "debugfs" it never
// provisions), so this is the one place a caller can ask "is name actually
// one of mine" before pointing someone at `condatainer update --libexec`.
var providedBins = map[string]bool{
	"apptainer":     true, // apptainer
	"mksquashfs":    true, // squashfs-tools
	"unsquashfs":    true, // squashfs-tools
	"squashfuse":    true, // squashfuse
	"squashfuse_ll": true, // squashfuse
	"micromamba":    true, // micromamba
}

// Provides reports whether name is one of the binaries this package's own
// spec installs, independent of whether a toolchain is actually provisioned
// anywhere right now.
func Provides(name string) bool {
	return providedBins[name]
}

// NotProvisionedMessage explains why name could not be found — for embedding
// in a generated shell script's own failure message, where a Go error value
// cannot reach. It names the toolchain-install fix only when Provides(name):
// pointing someone at `condatainer update --libexec` for a tool this package
// never installs (e2fsprogs, say) would not help. Single-quoted, not
// backtick-quoted like ErrNotProvisioned's own text: this string is meant to
// sit inside a double-quoted shell echo, where a backtick would attempt
// command substitution instead of printing literally.
func NotProvisionedMessage(name string) string {
	if Provides(name) {
		return fmt.Sprintf("%s not found; run 'condatainer update --libexec' to install it", name)
	}
	return name + " not found"
}

// lockFileName is the per-generation lock sentinel a reader holds LOCK_SH on
// and Update takes LOCK_EX on. See the README, Locking.
const lockFileName = ".lock"

// Dir returns the nearest tier's libexec directory that is actually
// provisioned (has bin/apptainer), resolved to its real (symlink-free) path,
// and true. Returns "", false if none is.
func Dir() (string, bool) {
	for _, dir := range config.GetLibexecSearchPaths() {
		if _, err := os.Stat(filepath.Join(dir, "bin", binMarker)); err == nil {
			if real, err := filepath.EvalSymlinks(dir); err == nil {
				dir = real
			}
			return dir, true
		}
	}
	return "", false
}

// BinDir returns the nearest provisioned tier's bin/ directory.
func BinDir() (string, bool) {
	dir, ok := Dir()
	if !ok {
		return "", false
	}
	return filepath.Join(dir, "bin"), true
}

// ApptainerPath returns the path to the self-provisioned apptainer binary,
// the one ordinary exec/run and conda/script builds use — never the
// fakeroot-capable system/module one an os/base .def build needs.
func ApptainerPath() (string, bool) {
	return Path("apptainer")
}

// MksquashfsPath returns the path to the self-provisioned mksquashfs binary.
func MksquashfsPath() (string, bool) {
	return Path("mksquashfs")
}

// UnsquashfsPath returns the path to the self-provisioned unsquashfs binary —
// the squashfs-tools package installs it alongside mksquashfs.
func UnsquashfsPath() (string, bool) {
	return Path("unsquashfs")
}

// SquashfusePath returns the path to the self-provisioned squashfuse binary.
func SquashfusePath() (string, bool) {
	return Path("squashfuse")
}

// MicromambaPath returns the path to the self-provisioned micromamba binary.
func MicromambaPath() (string, bool) {
	return Path("micromamba")
}

// Path returns the path a binary called name would have in the provisioned
// bin/, and true if a toolchain is provisioned at all — regardless of
// whether name is actually one of the tools it installs. The named
// accessors above are for a container-bound caller that already knows which
// tool it wants; internal/toolpath.Resolve calls this one directly for an
// arbitrary name, since a tool this package never provisions simply isn't
// found on the host-side stat check that follows.
func Path(name string) (string, bool) {
	bin, ok := BinDir()
	if !ok {
		return "", false
	}
	return filepath.Join(bin, name), true
}

// LockPath returns the live generation's lock sentinel, and true if a
// toolchain is provisioned at all.
func LockPath() (string, bool) {
	dir, ok := Dir()
	if !ok {
		return "", false
	}
	return filepath.Join(dir, lockFileName), true
}

// AcquireUse takes a shared lock on the live toolchain for the duration of
// one use — one apptainer subprocess call. Returns (nil, nil) if nothing is
// provisioned yet, in which case there is nothing to protect and callers must
// treat a nil lock as "no lock held," not an error. A non-nil error means
// Update holds the exclusive lock right now.
func AcquireUse() (*utils.FlockHandle, error) {
	path, ok := LockPath()
	if !ok {
		return nil, nil
	}
	if !utils.FileExists(path) {
		// Self-heal a generation that predates this mechanism (see README).
		if f, err := utils.CreateFileWritable(path); err == nil {
			f.Close()
		}
	}
	lock, err := utils.AcquireFlock(path, false)
	if err != nil {
		return nil, fmt.Errorf("the toolchain is being updated right now; try again in a moment: %w", err)
	}
	return lock, nil
}
