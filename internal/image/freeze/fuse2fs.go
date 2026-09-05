package freeze

import (
	"errors"
	"fmt"
	"os/exec"
)

// ErrNoFuse2fs reports that no fuse2fs could be found to read the overlay with.
var ErrNoFuse2fs = errors.New("no fuse2fs available to read the overlay")

// ErrNoSquashfuse reports that no squashfuse could be found to read an artifact.
var ErrNoSquashfuse = errors.New("no squashfuse available to read the artifact")

// findSquashfuse locates the squashfuse that mounts a frozen artifact, for
// mountedRun to run directly — no container involved to bind it into.
func findSquashfuse() (string, error) {
	for _, name := range []string{"squashfuse_ll", "squashfuse"} {
		if p, err := exec.LookPath(name); err == nil {
			return p, nil
		}
	}
	return "", fmt.Errorf("%w: looked on PATH", ErrNoSquashfuse)
}

// findFuse2fs locates the fuse2fs that will read the image, the same way
// findSquashfuse locates its counterpart.
func findFuse2fs() (string, error) {
	if p, err := exec.LookPath("fuse2fs"); err == nil {
		return p, nil
	}
	return "", fmt.Errorf("%w: looked on PATH", ErrNoFuse2fs)
}
