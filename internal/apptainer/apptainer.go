package apptainer

import (
	"bytes"
	"context"
	"errors"
	"fmt"
	"io"
	"os"
	"os/exec"
	"regexp"
	"strconv"
	"strings"
	"syscall"
	"time"

	"github.com/Justype/condatainer/internal/utils"
)

// apptainerCmd holds the resolved, absolute path to the binary.
var apptainerCmd string

// cachedVersion holds the cached apptainer version to avoid repeated calls
var cachedVersion string

// SetBin configures and validates the Apptainer binary path.
// If already set to the same path, this is a no-op.
func SetBin(path string) error {
	target := path
	if target == "" {
		target = "apptainer"
	}

	fullPath, err := exec.LookPath(target)
	if err != nil {
		// Return structured error without hints - let caller decide what to suggest
		return &ApptainerNotFoundError{Path: target}
	}

	// Skip if already set to the same path
	if apptainerCmd == fullPath {
		return nil
	}

	apptainerCmd = fullPath
	cachedVersion = "" // Clear cached version when binary changes
	utils.PrintDebug("Apptainer binary resolved to: %s", utils.StylePath(apptainerCmd))
	return nil
}

// Which returns the currently configured Apptainer executable path.
func Which() string {
	return apptainerCmd
}

// GetVersion returns the version of the currently loaded Apptainer binary.
// Results are cached to avoid repeated calls.
func GetVersion() (string, error) {
	// Return cached version if available
	if cachedVersion != "" {
		return cachedVersion, nil
	}

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

	cachedVersion = match
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
func runApptainer(ctx context.Context, op string, imagePath string, capture bool, args ...string) error {
	return runApptainerWithOutput(ctx, op, imagePath, capture, false, nil, args...)
}

// runApptainerWithOutput executes an apptainer command with control over output handling.
//
// op: "exec", "pull", "build", "instance start", etc.
// imagePath: optional, helpful for logging which container failed
// capture: whether to capture stdout/stderr for storing on the error.
// hideOutput: whether to redirect stdout/stderr to /dev/null
func runApptainerWithOutput(ctx context.Context, op string, imagePath string, capture bool, hideOutput bool, stdin io.Reader, args ...string) error {
	cmd := exec.CommandContext(ctx, apptainerCmd, args...)

	// Set stdin - use provided stdin or default to os.Stdin
	if stdin != nil {
		cmd.Stdin = stdin
	} else {
		cmd.Stdin = os.Stdin
	}

	// For build operations, unset SINGULARITY_BIND and APPTAINER_BIND to prevent
	// mount conflicts during container build (e.g., when %post tries to access bound paths)
	if op == "build" {
		env := os.Environ()
		filteredEnv := make([]string, 0, len(env))
		for _, e := range env {
			if !strings.HasPrefix(e, "SINGULARITY_BIND=") && !strings.HasPrefix(e, "APPTAINER_BIND=") {
				filteredEnv = append(filteredEnv, e)
			}
		}
		cmd.Env = filteredEnv
	}

	var stdoutBuf, stderrBuf bytes.Buffer
	var stdoutWriter, stderrWriter io.Writer

	if hideOutput {
		// Redirect to /dev/null - suppress all output
		devNull, err := os.OpenFile(os.DevNull, os.O_WRONLY, 0)
		if err != nil {
			return fmt.Errorf("failed to open /dev/null: %w", err)
		}
		defer devNull.Close()
		stdoutWriter = devNull
		stderrWriter = devNull
	} else if capture {
		// Only buffer output when explicitly requested for capture
		// This prevents memory exhaustion on long-running builds in SLURM
		stdoutWriter = io.MultiWriter(os.Stdout, &stdoutBuf)
		stderrWriter = io.MultiWriter(os.Stderr, &stderrBuf)
	} else {
		// Normal streaming output
		stdoutWriter = os.Stdout
		stderrWriter = os.Stderr
	}

	cmd.Stdout = stdoutWriter
	cmd.Stderr = stderrWriter

	cmd.WaitDelay = 5 * time.Second
	cmd.Cancel = func() error {
		// Use SIGTERM first for graceful shutdown
		return cmd.Process.Signal(syscall.SIGTERM)
	}

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
