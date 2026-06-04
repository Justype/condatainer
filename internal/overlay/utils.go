package overlay

import (
	"bytes"
	"context"
	"fmt"
	"io"
	"os/exec"
	"strings"

	"github.com/Justype/condatainer/internal/logging"
	"github.com/Justype/condatainer/internal/utils"
)

// checkDependencies verifies that all tools in the provided list are available in the system PATH.
// It returns a consolidated error listing all missing tools, or nil if all are present.
func checkDependencies(tools []string) error {
	var missing []string

	for _, tool := range tools {
		if _, err := exec.LookPath(tool); err != nil {
			missing = append(missing, tool)
		}
	}

	if len(missing) > 0 {
		return fmt.Errorf("missing required system tools: %s",
			utils.StyleError(strings.Join(missing, ", ")))
	}

	return nil
}

// runCommand executes a shell command and wraps failures in overlay.Error.
// When ctx carries a writer (via logging.WithWriter), stdout/stderr are streamed
// there in real time (web terminal). Otherwise output is buffered and only
// surfaced on error (CLI behaviour unchanged).
func runCommand(ctx context.Context, op, path, tool string, args ...string) error {
	cmd := exec.CommandContext(ctx, tool, args...)
	logging.FromContext(ctx).Debug("running "+tool, "op", op, "args", strings.Join(args, " "))

	var errBuf bytes.Buffer
	if w := logging.WriterFromCtx(ctx); w != nil {
		out := io.MultiWriter(w, &errBuf)
		cmd.Stdout = out
		cmd.Stderr = out
	} else {
		cmd.Stdout = &errBuf
		cmd.Stderr = &errBuf
	}

	if err := cmd.Run(); err != nil {
		return &Error{
			Op:      op,
			Path:    path,
			Tool:    tool,
			Output:  errBuf.String(),
			BaseErr: err,
		}
	}
	return nil
}
