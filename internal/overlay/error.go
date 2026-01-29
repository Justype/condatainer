package overlay

import (
	"fmt"
	"strings"

	"github.com/Justype/condatainer/internal/utils"
)

// Error represents a failure in filesystem tools (dd, mkfs, debugfs, etc.).
// Usage: if err, ok := err.(*overlay.Error); ok { ... }
type Error struct {
	Op      string // high level intent: "create", "resize"
	Tool    string // low level tool: "mke2fs", "dd", "debugfs"
	Path    string // the file being manipulated
	Output  string // Captured Stderr/Stdout
	BaseErr error  // The underlying execution error
}

func (e *Error) Error() string {
	hint := e.analyze()
	var msg strings.Builder

	msg.WriteString(fmt.Sprintf("Overlay filesystem operation '%s' failed.\n", utils.StyleAction(e.Op)))
	msg.WriteString(fmt.Sprintf("\tTarget:  %s\n", utils.StylePath(e.Path)))
	msg.WriteString(fmt.Sprintf("\tTool:    %s\n", utils.StyleCommand(e.Tool)))

	if e.Output != "" {
		cleanOut := strings.TrimSpace(e.Output)
		if len(cleanOut) > 0 {
			// Indent output slightly for better readability
			msg.WriteString(fmt.Sprintf("\tOutput:  %s\n", utils.StyleError(cleanOut)))
		}
	}

	if hint != "" {
		msg.WriteString(fmt.Sprintf("\t%s    %s\n", utils.StyleHint("Hint:"), hint))
	}

	msg.WriteString(fmt.Sprintf("\tError:   %v", e.BaseErr))

	return msg.String()
}

// Unwrap allows errors.Is/As to see the underlying BaseErr
func (e *Error) Unwrap() error {
	return e.BaseErr
}

func (e *Error) analyze() string {
	out := e.Output

	// --- General System Errors ---
	if strings.Contains(out, "No space left on device") {
		return "Host storage is full. Cannot allocate image file."
	}
	if strings.Contains(out, "Permission denied") {
		return "Check file permissions. You may not have write access to the destination."
	}
	if strings.Contains(out, "Read-only file system") {
		return "Destination filesystem is Read-Only."
	}
	// VERY Common: User tries to resize an image while container is running
	if strings.Contains(out, "Device or resource busy") || strings.Contains(out, "Text file busy") {
		return "The overlay image is currently in use. Stop any running containers using it first."
	}

	// --- Tool Specific: Resize2fs ---
	if strings.Contains(out, "New size smaller than minimum") {
		return "Cannot shrink image below current usage. Try a larger size."
	}

	// --- Tool Specific: E2fsck / Tune2fs ---
	if strings.Contains(out, "Bad magic number") || strings.Contains(out, "Not a valid filesystem") {
		return "Target file is not a valid ext3 filesystem. Is this an overlay image?"
	}
	if strings.Contains(out, "needs human intervention") {
		return fmt.Sprintf("Filesystem corrupted. Run 'e2fsck -y %s' manually.", e.Path)
	}
	if strings.Contains(out, "is mounted") {
		return "Cannot perform this operation while the image is mounted."
	}

	// --- Tool Specific: Debugfs ---
	if strings.Contains(out, "File not found") && strings.Contains(e.Tool, "debugfs") {
		return "Internal overlay structure is damaged or missing."
	}

	return ""
}
