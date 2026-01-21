package apptainer

import (
	"fmt"
	"strings"

	"github.com/Justype/condatainer/internal/utils"
)

// ApptainerNotFoundError is returned when the apptainer binary cannot be found
type ApptainerNotFoundError struct {
	Path string // The path or command name that was searched for
}

func (e *ApptainerNotFoundError) Error() string {
	return fmt.Sprintf("apptainer binary not found: %s", e.Path)
}

// ApptainerError represents errors from apptainer command execution
type ApptainerError struct {
	Op         string // Operation being performed (e.g., "build", "exec", "pull")
	Cmd        string // Full command that was executed
	Path       string // Container/image path (optional)
	Output     string // Command output (stdout/stderr)
	BaseErr    error  // Underlying error
	HideOutput bool   // Suppress output in the formatted error
}

func (e *ApptainerError) Error() string {
	hint := e.analyzeError()

	var msg strings.Builder
	if e.Path != "" {
		msg.WriteString(fmt.Sprintf("Apptainer failed to %s %s.\n",
			e.Op, utils.StylePath(e.Path)))
	} else {
		msg.WriteString(fmt.Sprintf("Apptainer %s failed.\n", e.Op))
	}

	if e.Output != "" && !e.HideOutput {
		cleanOut := strings.TrimSpace(e.Output)
		if cleanOut != "" {
			msg.WriteString(fmt.Sprintf("\t%s: %s\n",
				utils.StyleHint("Output"),
				utils.StyleError(cleanOut)))
		}
	}

	if hint != "" {
		msg.WriteString(fmt.Sprintf("\t%s %s\n",
			utils.StyleHint("Hint:"), hint))
	}

	msg.WriteString(fmt.Sprintf("\t%s: %v",
		utils.StyleHint("Error"), e.BaseErr))
	return msg.String()
}

func (e *ApptainerError) analyzeError() string {
	out := e.Output
	if strings.Contains(out, "FATAL:") && strings.Contains(out, "container creation failed") {
		return "Container creation failed. Check if the base image path is correct."
	}
	if strings.Contains(out, "not a squashfs image") {
		return "The SIF image appears corrupted or invalid."
	}
	if strings.Contains(out, "No space left on device") {
		return "Host storage is full."
	}
	if strings.Contains(out, "Permission denied") {
		return "You do not have permission to access the container image or bind path."
	}
	if strings.Contains(out, "not found") && strings.Contains(out, "command") {
		return "The command was not found inside the container."
	}
	return ""
}
