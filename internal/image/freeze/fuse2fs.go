package freeze

import (
	"errors"
	"fmt"

	"github.com/Justype/condatainer/internal/toolpath"
)

// ErrNoFuse2fs reports that no fuse2fs could be found to read the overlay with.
var ErrNoFuse2fs = errors.New("no fuse2fs available to read the overlay")

// ErrNoSquashfuse reports that no squashfuse could be found to read an artifact.
var ErrNoSquashfuse = errors.New("no squashfuse available to read the artifact")

// FindSquashfuse locates the squashfuse that mounts a frozen artifact, for
// MountedRun to run directly — no container involved to bind it into.
// toolpath.Resolve prefers internal/libexec's own provisioned copy, so this
// finds one even on a host with no squashfuse of its own. squashfuse_ll is
// tried first (Apptainer's own preference) but the final message names
// "squashfuse": libexec provisions both from the same package, so either
// name works.
//
// Exported so internal/build's sif-import path can mount a .sif's own
// partition the same apptainer-free way this package mounts a .sqf.
func FindSquashfuse() (string, error) {
	for _, name := range []string{"squashfuse_ll", "squashfuse"} {
		if p, err := toolpath.Resolve(name); err == nil {
			return p, nil
		}
	}
	return "", fmt.Errorf("%w: %s", ErrNoSquashfuse, toolpath.NotFoundMessage("squashfuse"))
}

// findFuse2fs locates the fuse2fs that will read the image, the same way
// FindSquashfuse locates its counterpart — except libexec never provisions
// it (internal/image/README.md), so NotFoundMessage never suggests
// `condatainer update --libexec` for it, correctly: that command would not
// help.
func findFuse2fs() (string, error) {
	if p, err := toolpath.Resolve("fuse2fs"); err == nil {
		return p, nil
	}
	return "", fmt.Errorf("%w: %s", ErrNoFuse2fs, toolpath.NotFoundMessage("fuse2fs"))
}
