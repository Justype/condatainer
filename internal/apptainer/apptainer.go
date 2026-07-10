package apptainer

import (
	"bytes"
	"context"
	"errors"
	"fmt"
	"io"
	"os"
	"os/exec"
	"path/filepath"
	"regexp"
	"strconv"
	"strings"
	"syscall"
	"time"

	"github.com/Justype/condatainer/internal/config"
)

// apptainerCmd holds the resolved, absolute path to the binary.
var apptainerCmd string

// cachedVersion holds the cached apptainer version to avoid repeated calls
var cachedVersion string

// SetBin configures and validates the Apptainer binary path.
// If already set to the same path, this is a no-op.
// When path is empty, tries "apptainer" then "singularity" in PATH order.
func SetBin(path string) error {
	if path == "" {
		// Try apptainer first, then singularity as fallback
		for _, name := range []string{"apptainer", "singularity"} {
			if fullPath, err := exec.LookPath(name); err == nil {
				if apptainerCmd == fullPath {
					return nil
				}
				apptainerCmd = fullPath
				cachedVersion = ""
				return nil
			}
		}
		return &ApptainerNotFoundError{Path: "apptainer/singularity"}
	}

	fullPath, err := exec.LookPath(path)
	if err != nil {
		// Return structured error without hints - let caller decide what to suggest
		return &ApptainerNotFoundError{Path: path}
	}

	// Skip if already set to the same path
	if apptainerCmd == fullPath {
		return nil
	}

	apptainerCmd = fullPath
	cachedVersion = "" // Clear cached version when binary changes
	return nil
}

// EnsureApptainer checks if apptainer binary is available and configured.
// Returns an error if apptainer cannot be found.
func EnsureApptainer() error {
	return SetBin(config.Global.ApptainerBin)
}

// IsSingularity returns true if the configured binary is Singularity (not Apptainer).
// Singularity defaults to gzip compression for SquashFS images.
func IsSingularity() bool {
	return strings.Contains(strings.ToLower(filepath.Base(apptainerCmd)), "singularity")
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

	return false
}

// ---------------------------------------------------------
// Internal Helper
// ---------------------------------------------------------

// runApptainerWithOutput executes an apptainer command with control over output handling.
//
// op: "exec", "pull", "build", etc.
// imagePath: optional, helpful for logging which container failed
// capture: whether to capture stdout/stderr for storing on the error.
// stdin: custom stdin reader (optional; nil means no stdin)
// stdout/stderr: redirect streams to these writers (optional; nil means discard).
//
//	Both are always teed with an internal buffer so ApptainerError.Output is
//	populated on failure regardless of redirection.
func runApptainerWithOutput(ctx context.Context, op string, imagePath string, capture bool, stdin io.Reader, stdout, stderr io.Writer, args ...string) error {
	cmd := exec.CommandContext(ctx, apptainerCmd, args...)

	cmd.Stdin = stdin

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

	if capture {
		stdoutWriter = &stdoutBuf
		stderrWriter = &stderrBuf
	} else {
		// Tee custom writers with error buffers so ApptainerError.Output
		// is populated on failure when doing so does not hide a real TTY from
		// the child process. Interactive CLI callers pass *os.File streams, and
		// those must be attached directly so shells and tools can detect width,
		// colors, job control, etc.
		if stdout != nil {
			if f, ok := stdout.(*os.File); ok {
				stdoutWriter = f
			} else {
				stdoutWriter = io.MultiWriter(stdout, &stdoutBuf)
			}
		} else {
			stdoutWriter = io.Discard
		}
		if stderr != nil {
			if f, ok := stderr.(*os.File); ok {
				stderrWriter = f
			} else {
				stderrWriter = io.MultiWriter(stderr, &stderrBuf)
			}
		} else {
			stderrWriter = io.Discard
		}
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
// build target, or when the build was interrupted by the user (Ctrl+C).
func IsBuildCancelled(err error) bool {
	// Check for context cancellation (Ctrl+C)
	if errors.Is(err, context.Canceled) {
		return true
	}

	// Check for signal interrupts in error message
	errMsg := err.Error()
	if strings.Contains(errMsg, "signal: interrupt") || strings.Contains(errMsg, "signal: killed") {
		return true
	}

	// Check for SIGINT exit code (130 = 128 + 2)
	var exitErr *exec.ExitError
	if errors.As(err, &exitErr) {
		if exitErr.ExitCode() == 130 || exitErr.ExitCode() == -1 {
			return true
		}
	}

	// Check for interactive cancellation (user declining to overwrite)
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
