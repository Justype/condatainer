package ext3

import (
	"context"
	"fmt"
	"os/exec"
	"path/filepath"
	"strings"

	"github.com/Justype/condatainer/internal/image"
	"github.com/Justype/condatainer/internal/image/tool"
	"github.com/Justype/condatainer/internal/logging"
	"github.com/Justype/condatainer/internal/toolpath"
	"github.com/Justype/condatainer/internal/utils"
)

// CheckIntegrity runs a filesystem check (e2fsck) on the overlay image.
// force: If true, adds '-f' to force the check even if the filesystem appears clean.
// Returns an error immediately if the image is currently in use (mounted writable).
func CheckIntegrity(ctx context.Context, path string, force bool) error {
	if err := image.CheckAvailable(path, true); err != nil {
		return fmt.Errorf("%s is currently in use — stop any running jobs using it first", utils.StylePath(path))
	}

	e2fsckPath, err := toolpath.Resolve("e2fsck")
	if err != nil {
		return err
	}

	args := []string{"-p"}
	if force {
		args = append(args, "-f")
	}
	args = append(args, path)

	log := logging.FromContext(ctx)
	log.Info(fmt.Sprintf("checking integrity of %s", utils.StylePath(filepath.Base(path))))
	log.Debug("e2fsck " + strings.Join(args, " "))

	cmd := exec.CommandContext(ctx, e2fsckPath, args...)
	out, err := cmd.CombinedOutput()

	// e2fsck exit codes: 0 = clean, 1 = errors corrected, 2+ = critical failure.
	if err != nil {
		exitErr, ok := err.(*exec.ExitError)
		if !ok || exitErr.ExitCode() > 1 {
			return &tool.Error{
				Op:      "check integrity",
				Path:    path,
				Tool:    "e2fsck",
				Output:  string(out),
				BaseErr: err,
			}
		}
	}

	if len(out) > 0 {
		log.Debug("e2fsck output: " + string(out))
	}

	log.Info(fmt.Sprintf("filesystem check completed for %s", utils.StylePath(filepath.Base(path))), "kind", "success")
	return nil
}
