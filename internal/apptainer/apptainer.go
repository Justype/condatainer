package apptainer

import (
	"bytes"
	"errors"
	"fmt"
	"io"
	"os"
	"os/exec"
	"regexp"
	"strconv"
	"strings"

	"github.com/Justype/condatainer/internal/utils"
)

// apptainerCmd holds the resolved, absolute path to the binary.
var apptainerCmd string

// SetBin configures and validates the Apptainer binary path.
func SetBin(path string) error {
	target := path
	if target == "" {
		target = "apptainer"
	}

	fullPath, err := exec.LookPath(target)
	if err != nil {
		baseErr := &ApptainerNotFoundError{Path: target}

		// Use StyleHint for the label and StyleAction for the command
		hint := fmt.Sprintf(
			"%s On HPC systems, you likely need to load the module:\n\t%s",
			utils.StyleHint("Hint:"),
			utils.StyleAction("module load apptainer"),
		)
		return fmt.Errorf("%w\n%s", baseErr, hint)
	}

	apptainerCmd = fullPath
	utils.PrintDebug("Apptainer binary resolved to: %s", utils.StylePath(apptainerCmd))
	return nil
}

// Which returns the currently configured Apptainer executable path.
func Which() string {
	return apptainerCmd
}

// GetVersion returns the version of the currently loaded Apptainer binary.
func GetVersion() (string, error) {
	// We use exec.Command directly here because we need the output string
	// and runApptainer (below) is designed for void/error returns.
	cmd := exec.Command(apptainerCmd, "--version")
	output, err := cmd.Output()

	if err != nil {
		// Use the new ApptainerError struct with Op/Cmd
		return "", &ApptainerError{
			Op:      "check version",
			Cmd:     fmt.Sprintf("%s --version", apptainerCmd),
			Output:  "", // Output is often empty on binary execution failure
			BaseErr: err,
		}
	}

	outStr := strings.TrimSpace(string(output))
	re := regexp.MustCompile(`(\d+\.\d+(\.\d+)?)`)
	match := re.FindString(outStr)

	if match == "" {
		return "", fmt.Errorf("could not parse version from output: %s", outStr)
	}

	utils.PrintDebug("Detected Apptainer Version: %s", utils.StyleNumber(match))
	return match, nil
}

// CheckZstdSupport checks if the *current* version supports ZSTD.
func CheckZstdSupport(currentVersion string) bool {
	parts := strings.Split(currentVersion, ".")
	if len(parts) < 2 {
		return false
	}

	major, _ := strconv.Atoi(parts[0])
	minor, _ := strconv.Atoi(parts[1])

	// Apptainer supported ZSTD starting around version 1.4
	if major > 1 || (major == 1 && minor >= 4) {
		return true
	}

	utils.PrintDebug("ZSTD compression %s (requires Apptainer >= 1.4, found %s)",
		utils.StyleInfo("disabled"),
		utils.StyleNumber(currentVersion),
	)
	return false
}

// ---------------------------------------------------------
// Internal Helper
// ---------------------------------------------------------

// runApptainer executes an apptainer command and wraps errors with specific context.
// Use this for any command where you don't strictly need to parse stdout.
//
// op: "exec", "pull", "build", "instance start", etc.
// imagePath: optional, helpful for logging which container failed
// capture: whether to capture stdout/stderr for storing on the error.
func runApptainer(op string, imagePath string, capture bool, args ...string) error {
	cmd := exec.Command(apptainerCmd, args...)
	cmd.Stdin = os.Stdin

	var stdoutBuf, stderrBuf bytes.Buffer
	var stdoutWriter, stderrWriter io.Writer
	if capture {
		stdoutWriter = &stdoutBuf
		stderrWriter = &stderrBuf
	} else {
		stdoutWriter = io.MultiWriter(os.Stdout, &stdoutBuf)
		stderrWriter = io.MultiWriter(os.Stderr, &stderrBuf)
	}

	cmd.Stdout = stdoutWriter
	cmd.Stderr = stderrWriter

	err := cmd.Run()
	if err != nil {
		fullCmd := fmt.Sprintf("%s %s", apptainerCmd, strings.Join(args, " "))
		return &ApptainerError{
			Op:         op,
			Cmd:        fullCmd,
			Path:       imagePath,
			Output:     captureOutput(&stdoutBuf, &stderrBuf),
			HideOutput: !capture,
			BaseErr:    err,
		}
	}
	return nil
}

func captureOutput(stdoutBuf, stderrBuf *bytes.Buffer) string {
	combined := strings.TrimSpace(stdoutBuf.String() + stderrBuf.String())
	return combined
}

// IsBuildCancelled returns true when the Apptainer error looks like the
// user declined to continue after being asked about overwriting an existing
// build target.
func IsBuildCancelled(err error) bool {
	var apErr *ApptainerError
	if !errors.As(err, &apErr) {
		return false
	}

	output := strings.ToLower(apErr.Output)
	if output == "" {
		return false
	}

	return strings.Contains(output, "build target") &&
		strings.Contains(output, "do you want to continue")
}
