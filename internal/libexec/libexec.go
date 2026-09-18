// Package libexec resolves and provisions CondaTainer's self-provisioned
// toolchain — micromamba plus, on request, mksquashfs, squashfuse and an
// ordinary (non-fakeroot) apptainer, installed via micromamba into one of the
// four data-directory tiers. See the package README for the design.
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
// as provisioned, rather than an empty or half-written directory. micromamba
// is the one package every prefix holds.
const binMarker = "micromamba"

// ErrNotProvisioned is what a container-bound caller returns when nothing is
// provisioned at any tier. Condatainer never auto-provisions on its own — an
// ordinary exec/run or a build must never trigger a first-time download — so
// every message directs the user to the one command that does.
var ErrNotProvisioned = errors.New("the self-provisioned toolchain is not installed; run `condatainer update --libexec` first")

// Provides reports whether name is a binary or package this package can
// install, independent of whether it is installed anywhere right now. Path
// and BinDir build a path for any name (that is what lets toolpath.Resolve
// ask about names like "debugfs" this package never provisions), so this is
// the one place a caller can ask "is name actually one of mine" before
// pointing someone at `condatainer update --libexec`.
func Provides(name string) bool {
	_, ok := packageFor(name)
	return ok
}

// NotProvisionedMessage explains why name could not be found — for embedding
// in a generated shell script's own failure message, where a Go error value
// cannot reach. It names the install command only when Provides(name):
// pointing someone at `condatainer update --libexec` for a tool this package
// never installs (e2fsprogs, say) would not help. Single-quoted, not
// backtick-quoted like ErrNotProvisioned's own text: this string is meant to
// sit inside a double-quoted shell echo, where a backtick would attempt
// command substitution instead of printing literally.
func NotProvisionedMessage(name string) string {
	if p, ok := packageFor(name); ok {
		return fmt.Sprintf("%s not found; run 'condatainer update --libexec %s' to install it", name, p.name)
	}
	return name + " not found"
}

// NotInstalledError is what a container-bound caller returns when name is not
// in the toolchain: ErrNotProvisioned when no tier is provisioned at all,
// otherwise a message naming the package that installs it.
func NotInstalledError(name string) error {
	if _, ok := Dir(); !ok {
		return ErrNotProvisioned
	}
	p, ok := packageFor(name)
	if !ok {
		return fmt.Errorf("%s is not part of the self-provisioned toolchain", name)
	}
	return fmt.Errorf("%s is not installed in the self-provisioned toolchain; run `condatainer update --libexec %s`", name, p.name)
}

// lockFileName is the per-generation lock sentinel a reader holds LOCK_SH on
// and Update takes LOCK_EX on. See the README, Locking.
const lockFileName = ".lock"

// Dir returns the nearest tier's libexec directory that is actually
// provisioned (has bin/micromamba), resolved to its real (symlink-free) path,
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

// Path returns the path of the binary called name in the provisioned bin/,
// and true only if it is installed there. internal/toolpath.Resolve calls this
// for an arbitrary name, and a name this package never provisions is simply
// not found.
func Path(name string) (string, bool) {
	bin, ok := BinDir()
	if !ok {
		return "", false
	}
	path := filepath.Join(bin, name)
	if !utils.FileExists(path) {
		return "", false
	}
	return path, true
}

// Installed reports whether the nearest provisioned tier has the binary name.
func Installed(name string) bool {
	_, ok := Path(name)
	return ok
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
