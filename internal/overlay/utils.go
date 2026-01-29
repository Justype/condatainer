package overlay

import (
	"fmt"
	"os/exec"
	"strings"

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
func runCommand(op, path, tool string, args ...string) error {
	cmd := exec.Command(tool, args...)
	utils.PrintDebug("[%s] Running %s %s", op, tool, strings.Join(args, " "))
	output, err := cmd.CombinedOutput()

	if err != nil {
		return &Error{
			Op:      op,
			Path:    path,
			Tool:    tool,
			Output:  string(output),
			BaseErr: err,
		}
	}
	return nil
}
