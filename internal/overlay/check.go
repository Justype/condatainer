package overlay

import (
	"os/exec"

	"github.com/Justype/condatainer/internal/utils"
)

// CheckIntegrity runs a filesystem check (e2fsck) on the overlay image.
// force: If true, adds '-f' to force the check even if the filesystem appears clean.
func CheckIntegrity(path string, force bool) error {
	// 0. Check Dependencies
	if err := checkDependencies([]string{"e2fsck"}); err != nil {
		return err
	}

	// 1. Prepare Arguments
	// -p: Automatically repair (preen) without asking questions
	args := []string{"-p"}

	if force {
		args = append(args, "-f")
	}

	args = append(args, path)

	utils.PrintMessage("Checking integrity of %s (Force: %v)...", utils.StylePath(path), force)

	// 2. Run e2fsck
	cmd := exec.Command("e2fsck", args...)
	out, err := cmd.CombinedOutput()

	// 3. Handle Exit Codes
	// e2fsck uses non-standard exit codes:
	// 0 = No errors
	// 1 = File system errors corrected (Success for our purposes)
	// 2+ = Critical errors or reboot required
	if err != nil {
		exitErr, ok := err.(*exec.ExitError)
		// If casting failed (not an exit error) OR exit code > 1, treat as failure
		if !ok || exitErr.ExitCode() > 1 {
			return &Error{
				Op:      "check integrity",
				Path:    path,
				Tool:    "e2fsck",
				Output:  string(out),
				BaseErr: err,
			}
		}
	}

	// If we got here, it's either 0 (clean) or 1 (fixed).
	if len(out) > 0 {
		utils.PrintDebug("e2fsck output: %s", string(out))
	}

	utils.PrintSuccess("Filesystem is clean.")
	return nil
}
