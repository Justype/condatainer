package overlay

import (
	"context"
	"fmt"
	"os/exec"
	"strings"

	"github.com/Justype/condatainer/internal/utils"
)

// CheckIntegrity runs a filesystem check (e2fsck) on the overlay image.
// force: If true, adds '-f' to force the check even if the filesystem appears clean.
// Returns an error immediately if the image is currently in use (mounted writable).
func CheckIntegrity(ctx context.Context, path string, force bool) error {
	// 0. Refuse to run e2fsck on a live-mounted image.
	if err := CheckAvailable(path, true); err != nil {
		return fmt.Errorf("%s is currently in use — stop any running jobs using it first", utils.StylePath(path))
	}

	// 1. Check Dependencies
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

	utils.PrintMessage("Checking integrity of %s (force=%v)", utils.StylePath(path), force)
	utils.PrintDebug("[CHECK] e2fsck %s", strings.Join(args, " "))

	// 2. Run e2fsck
	cmd := exec.CommandContext(ctx, "e2fsck", args...)
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

	utils.PrintSuccess("Filesystem check completed for %s", utils.StylePath(path))
	return nil
}

// imgPathExists checks whether `entry` exists inside an ext3 image.
// It uses debugfs to stat the 'upper/<entry>' path in the image.
func imgPathExists(imgPath, entry string) bool {
	entry = strings.TrimPrefix(entry, "/")
	dbg, err := exec.LookPath("debugfs")
	if err != nil {
		return false
	}
	statArg := "stat upper/" + entry
	cmd := exec.Command(dbg, "-R", statArg, imgPath)
	out, err := cmd.CombinedOutput()
	if err != nil {
		outStr := string(out)
		if strings.Contains(outStr, "File not found") || strings.Contains(outStr, "No such file") {
			return false
		}
		return false
	}
	return true
}

// PathExists dispatches to the appropriate backend check depending on the overlay type.
func PathExists(overlayPath, entry string) bool {
	if utils.IsSqf(overlayPath) {
		return sqfPathExists(overlayPath, entry)
	}
	if utils.IsImg(overlayPath) {
		return imgPathExists(overlayPath, entry)
	}
	return false
}

// ListCondaPackages reads the conda-meta directory from a writable ext3 overlay
// image and returns a map of package name → installed version for every conda
// package found. The conda environment is expected at /ext3/env inside the
// container, which maps to upper/ext3/env/conda-meta inside the .img file.
//
// Returns (nil, nil) for non-.img files or when no conda-meta directory exists.
func ListCondaPackages(imgPath string) (map[string]string, error) {
	if !utils.IsImg(imgPath) {
		return nil, nil
	}
	dbg, err := exec.LookPath("debugfs")
	if err != nil {
		return nil, fmt.Errorf("debugfs not found: %w", err)
	}
	// ls -p outputs one entry per line: /<inode>/<mode>/<uid>/<gid>/<name>/<size>/
	cmd := exec.Command(dbg, "-R", "ls -p upper/ext3/env/conda-meta", imgPath)
	out, err := cmd.Output()
	if err != nil {
		// conda-meta directory missing means no conda environment installed yet.
		return nil, nil
	}
	pkgs := make(map[string]string)
	for _, line := range strings.Split(string(out), "\n") {
		line = strings.TrimSpace(line)
		// format: /<inode>/<mode>/<uid>/<gid>/<name>/<size>/
		parts := strings.Split(line, "/")
		if len(parts) < 6 {
			continue
		}
		filename := parts[5]
		if !strings.HasSuffix(filename, ".json") {
			continue
		}
		name, ver := parseCondaMetaFilename(filename)
		if name != "" {
			pkgs[name] = ver
		}
	}
	return pkgs, nil
}

// parseCondaMetaFilename parses a conda-meta filename such as
// "jupyterlab-4.2.1-pyhd8ed1ab_0.json" into (name, version).
// Conda versions always start with a digit, so the name ends at the first
// '-' that is immediately followed by a digit.
// Returns ("", "") if the filename does not match the expected pattern.
func parseCondaMetaFilename(filename string) (name, version string) {
	s := strings.TrimSuffix(filename, ".json")
	for i := 0; i < len(s)-1; i++ {
		if s[i] == '-' && s[i+1] >= '0' && s[i+1] <= '9' {
			name = s[:i]
			rest := s[i+1:]
			// version is the segment up to the next '-' (build string follows)
			if j := strings.Index(rest, "-"); j >= 0 {
				version = rest[:j]
			} else {
				version = rest
			}
			return
		}
	}
	return "", ""
}
