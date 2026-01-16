package apptainer

import (
	"fmt"
	"strings"

	"github.com/Justype/condatainer/internal/utils"
)

// ---------------------------------------------------------
// ApptainerError: For Container Runtime Failures
// (Exec, Build, Pull, SIF extraction)
// ---------------------------------------------------------

type ApptainerError struct {
	Op      string 
	Path    string 
	Cmd     string 
	Output  string 
	BaseErr error  
}

func (e *ApptainerError) Error() string {
	hint := e.analyze()

	var msg strings.Builder
	msg.WriteString(fmt.Sprintf("Apptainer failed to %s at %s.\n", e.Op, utils.StylePath(e.Path)))
	
	if e.Output != "" {
		// Use RedText for container errors as they are often critical execution failures
		msg.WriteString(fmt.Sprintf("\tOutput: %s\n", utils.RedText(strings.TrimSpace(e.Output))))
	}
	
	if hint != "" {
		msg.WriteString(fmt.Sprintf("\t%s %s\n", utils.YellowText("Hint:"), hint))
	}
	
	msg.WriteString(fmt.Sprintf("\tError: %v", e.BaseErr))
	return msg.String()
}

func (e *ApptainerError) analyze() string {
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
	return ""
}


// ---------------------------------------------------------
// OverlayError: For Filesystem/Image Tool Failures
// (Resize2fs, Debugfs, E2fsck, DD, Overlay Create)
// ---------------------------------------------------------

type OverlayError struct {
	Op      string 
	Path    string 
	Tool    string // e.g. "resize2fs", "debugfs", "mkfs"
	Output  string 
	BaseErr error  
}

func (e *OverlayError) Error() string {
	hint := e.analyze()

	var msg strings.Builder
	msg.WriteString(fmt.Sprintf("Overlay operation '%s' failed on %s.\n", e.Op, utils.StylePath(e.Path)))
	
	if e.Tool != "" {
		msg.WriteString(fmt.Sprintf("\tTool: %s\n", e.Tool))
	}
	
	if e.Output != "" {
		// Clean up output (remove excessive newlines common in fs tools)
		cleanOut := strings.TrimSpace(e.Output)
		if len(cleanOut) > 0 {
			msg.WriteString(fmt.Sprintf("\tOutput: %s\n", utils.RedText(cleanOut)))
		}
	}
	
	if hint != "" {
		msg.WriteString(fmt.Sprintf("\t%s %s\n", utils.YellowText("Hint:"), hint))
	}
	
	msg.WriteString(fmt.Sprintf("\tError: %v", e.BaseErr))
	return msg.String()
}

func (e *OverlayError) analyze() string {
	out := e.Output
	
	// Resize2fs errors
	if strings.Contains(out, "New size smaller than minimum") {
		return "You cannot shrink the image below its current usage. Try a larger size."
	}
	if strings.Contains(out, "Filesystem has unsupported features") {
		return "The overlay format is too new for this kernel version."
	}
	
	// E2fsck errors
	if strings.Contains(out, "needs human intervention") || strings.Contains(out, "Run e2fsck MANUALLY") {
		return "The filesystem is corrupted. Run 'e2fsck -y " + e.Path + "' manually to fix it."
	}

	// Debugfs errors
	if strings.Contains(out, "File not found") && strings.Contains(e.Tool, "debugfs") {
		return "The target directory does not exist inside the overlay image."
	}
	if strings.Contains(out, "Bad magic number") {
		return "This file is not a valid ext3 image. It might be corrupted."
	}

	// Common System errors
	if strings.Contains(out, "No space left on device") {
		return "Storage full. Cannot expand/create image."
	}
	if strings.Contains(out, "Read-only file system") {
		return "The destination filesystem is Read-Only."
	}
	
	return ""
}
