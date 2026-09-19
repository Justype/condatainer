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
	"sync"
	"syscall"
	"time"

	"github.com/Justype/condatainer/internal/config"
	"github.com/Justype/condatainer/internal/libexec"
	"github.com/Justype/condatainer/internal/toolpath"
)

// Bin is one resolved apptainer (or singularity) binary. A launch carries the
// Bin its caller resolved, so nothing about which binary runs is package state.
type Bin struct {
	Path    string
	Libexec bool // the self-provisioned copy in libexec/
}

// ErrNeedsHost reports an action that needs the host's apptainer: a container
// has no setuid starter to escalate with.
var ErrNeedsHost = errors.New("needs a host apptainer, which is not available inside a container")

// insideContainer is config.IsInsideContainer, replaceable by tests.
var insideContainer = config.IsInsideContainer

// last is the Bin the latest resolver returned, kept only so Current can report
// it: a build runs its container through exec.Run, which does not hand the Bin back.
var (
	lastMu sync.Mutex
	last   Bin
)

func remember(b Bin) Bin {
	lastMu.Lock()
	last = b
	lastMu.Unlock()
	return b
}

// systemBin returns the system/module apptainer: build.system_apptainer, which
// startup filled from PATH when unset.
func systemBin() (Bin, error) {
	path := config.Global.Build.SystemApptainer
	if path == "" {
		return Bin{}, &ApptainerNotFoundError{Path: "apptainer/singularity"}
	}
	full, err := exec.LookPath(path)
	if err != nil {
		return Bin{}, &ApptainerNotFoundError{Path: path}
	}
	return Bin{Path: full}, nil
}

// Normal returns the apptainer for an ordinary launch: the one installed in
// libexec, else the system/module one, which must be apptainer >= 1.4.
func Normal() (Bin, error) {
	if path, ok := libexec.ApptainerPath(); ok {
		return remember(Bin{Path: path, Libexec: true}), nil
	}
	bin, err := systemBin()
	if err == nil {
		err = bin.requireZstd("exec")
	}
	if err != nil {
		if insideContainer() {
			return Bin{}, fmt.Errorf("no usable apptainer in this container: %w; nested running needs the apptainer from libexec/ or an apptainer overlay (see nested_run)", err)
		}
		return Bin{}, fmt.Errorf("no usable apptainer: %w; install one with `condatainer update --libexec apptainer`, or load an apptainer module (>= 1.4)", err)
	}
	return remember(bin), nil
}

// Fakeroot returns the apptainer for a fakeroot exec: the system/module one,
// apptainer >= 1.4, since only its setuid starter can escalate. Never libexec's.
func Fakeroot() (Bin, error) {
	if insideContainer() {
		return Bin{}, fmt.Errorf("fakeroot exec %w", ErrNeedsHost)
	}
	bin, err := systemBin()
	if err == nil {
		err = bin.requireZstd("fakeroot exec")
	}
	if err != nil {
		return Bin{}, err
	}
	return remember(bin), nil
}

// ForBuild returns the apptainer for a .def build: the system/module one, with
// no version check, since the output is a sandbox and no zstd overlay is
// mounted during the build.
func ForBuild() (Bin, error) {
	if insideContainer() {
		return Bin{}, fmt.Errorf("a definition build %w", ErrNeedsHost)
	}
	bin, err := systemBin()
	if err != nil {
		return Bin{}, err
	}
	return remember(bin), nil
}

// requireZstd refuses singularity and an apptainer that cannot mount
// condatainer's zstd-compressed overlays. what names the caller in the message.
func (b Bin) requireZstd(what string) error {
	if b.IsSingularity() {
		return fmt.Errorf("%s requires apptainer (found singularity), which cannot mount condatainer's zstd-compressed overlays", what)
	}
	version, err := b.Version()
	if err != nil {
		return fmt.Errorf("could not determine the system apptainer's version: %w", err)
	}
	if !CheckZstdSupport(version) {
		return fmt.Errorf("system apptainer %s does not support zstd (>= 1.4 required); a %s must mount condatainer's zstd-compressed overlays with this binary — upgrade apptainer or the loaded module", version, what)
	}
	return nil
}

// ErrNotResolved reports that no apptainer binary has been resolved yet.
var ErrNotResolved = errors.New("no apptainer binary has been resolved yet")

// Current reports the implementation and version of the binary the latest
// resolver returned, performing no resolution of its own. A caller that already
// resolved one to launch a container reads the result back through this rather
// than deciding it again.
func Current() (implementation, version string, err error) {
	lastMu.Lock()
	bin := last
	lastMu.Unlock()
	if bin.Path == "" {
		return "", "", ErrNotResolved
	}
	version, err = bin.Version()
	return bin.Implementation(), version, err
}

// IsSingularity returns true if the binary is Singularity (not Apptainer).
// Singularity defaults to gzip compression for SquashFS images.
func (b Bin) IsSingularity() bool {
	return strings.Contains(strings.ToLower(filepath.Base(b.Path)), "singularity")
}

// Implementation returns the compatible implementation name. The manifest needs
// this beside the version because Apptainer and Singularity can use overlapping
// version numbers.
func (b Bin) Implementation() string {
	if b.IsSingularity() {
		return "singularity"
	}
	return "apptainer"
}

// versions caches each binary's version for the process, by path.
var versions sync.Map

var versionPattern = regexp.MustCompile(`(\d+\.\d+(\.\d+)?)`)

// Version returns the binary's version. It runs the binary only when this
// process and the per-user cache (unchanged size and modification time) have
// no answer for its path.
func (b Bin) Version() (string, error) {
	if cached, ok := versions.Load(b.Path); ok {
		return cached.(string), nil
	}
	version, err := toolpath.Remember("version", b.Path, b.runVersion)
	if err != nil {
		return "", err
	}
	versions.Store(b.Path, version)
	return version, nil
}

// runVersion runs the binary's --version and parses the number out.
func (b Bin) runVersion() (string, error) {
	// exec.Command directly: we need the output string, and runApptainer
	// (below) is designed for void/error returns.
	cmd := exec.Command(b.Path, "--version")
	output, err := cmd.Output()
	if err != nil {
		return "", &ApptainerError{
			Op:      "check version",
			Cmd:     fmt.Sprintf("%s --version", b.Path),
			Output:  "", // Output is often empty on binary execution failure
			BaseErr: err,
		}
	}

	outStr := strings.TrimSpace(string(output))
	match := versionPattern.FindString(outStr)
	if match == "" {
		return "", fmt.Errorf("could not parse version from output: %s", outStr)
	}
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
//
// procEnv: extra KEY=VALUE settings for apptainer's own environment, on top of
// the parent's. This is how APPTAINERENV_* vars reach the container.
func runApptainerWithOutput(ctx context.Context, bin Bin, op string, imagePath string, capture bool, stdin io.Reader, stdout, stderr io.Writer, procEnv []string, args ...string) error {
	var cmd *exec.Cmd
	if dir, ok := libexec.Dir(); bin.Libexec && ok {
		cmd = exec.CommandContext(ctx, "bash", append([]string{"-c", libexecActivationScript(dir), bin.Path}, args...)...)
	} else {
		cmd = exec.CommandContext(ctx, bin.Path, args...)
	}

	cmd.Stdin = stdin

	env := os.Environ()

	// For build operations, unset SINGULARITY_BIND and APPTAINER_BIND to prevent
	// mount conflicts during container build (e.g., when %post tries to access bound paths)
	if op == "build" {
		filteredEnv := make([]string, 0, len(env))
		for _, e := range env {
			if !strings.HasPrefix(e, "SINGULARITY_BIND=") && !strings.HasPrefix(e, "APPTAINER_BIND=") {
				filteredEnv = append(filteredEnv, e)
			}
		}
		env = filteredEnv
	}

	// Apptainer needs unsquashfs/mksquashfs in PATH (e.g. to extract a .sqf
	// into a sandbox). When the binary is the self-provisioned libexec copy,
	// add its bin/ to PATH so it finds them there.
	if bin.Libexec {
		env = append(env, "PATH="+filepath.Dir(bin.Path)+string(os.PathListSeparator)+os.Getenv("PATH"))
	}

	cmd.Env = append(env, procEnv...)

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
		fullCmd := fmt.Sprintf("%s %s", bin.Path, strings.Join(args, " "))
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

func libexecActivationScript(prefix string) string {
	q := "'" + strings.ReplaceAll(prefix, "'", `'\''`) + "'"
	return fmt.Sprintf(`if [ -d %[1]s/etc/conda/activate.d ]; then
  for __cnt_f in %[1]s/etc/conda/activate.d/*.sh; do
    [ -e "$__cnt_f" ] || continue
    CONDA_PREFIX=%[1]s . "$__cnt_f"
  done
fi
exec "$0" "$@"
`, q)
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
