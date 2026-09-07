package tool

import (
	"bytes"
	"context"
	"fmt"
	"io"
	"strings"

	"github.com/Justype/condatainer/internal/logging"
	"github.com/Justype/condatainer/internal/toolpath"
	"github.com/Justype/condatainer/internal/utils"
)

// CheckDependencies verifies that all tools in the provided list can be
// resolved (toolpath.Resolve: libexec, then PATH, then the FHS fallback
// directories). It returns a consolidated error listing all missing tools,
// or nil if all are present. Each entry is toolpath.NotFoundMessage(tool)
// rather than the bare name, so a caller that lists a libexec-provisioned
// tool (unfreeze.go's "unsquashfs", say) gets pointed at
// `condatainer update --libexec`, while an e2fsprogs/core-OS tool in the same
// list does not.
func CheckDependencies(tools []string) error {
	var missing []string

	for _, tool := range tools {
		if _, err := toolpath.Resolve(tool); err != nil {
			missing = append(missing, toolpath.NotFoundMessage(tool))
		}
	}

	if len(missing) > 0 {
		return fmt.Errorf("missing required system tools: %s",
			utils.StyleError(strings.Join(missing, "; ")))
	}

	return nil
}

// RunCommand executes a shell command and wraps failures in Error. tool is
// resolved via toolpath.Command (libexec, then PATH, then the FHS fallback
// directories) rather than handed to exec as a bare name.
// When ctx carries a writer (via logging.WithWriter), stdout/stderr are streamed
// there in real time (web terminal). Otherwise output is buffered and only
// surfaced on error (CLI behaviour unchanged).
func RunCommand(ctx context.Context, op, path, tool string, args ...string) error {
	cmd, err := toolpath.Command(ctx, tool, args...)
	if err != nil {
		return &Error{Op: op, Path: path, Tool: tool, BaseErr: err}
	}

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
