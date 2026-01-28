package utils

import (
	"fmt"
	"os"

	"github.com/fatih/color"
)

// DebugMode controls whether PrintDebug output is visible.
var DebugMode = false

// QuietMode controls whether verbose messages are suppressed (errors/warnings still shown)
var QuietMode = false

// projectPrefix is the standard tag for all logs.
const projectPrefix = "[CNT]"

// ---------------------------------------------------------
// 1. Private Color Definitions
//    (We hide these so we don't use raw colors in logic)
// ---------------------------------------------------------

var (
	red         = color.New(color.FgRed).SprintFunc()
	green       = color.New(color.FgGreen).SprintFunc()
	yellow      = color.New(color.FgYellow).SprintFunc()
	blueBold    = color.New(color.FgBlue, color.Bold).SprintFunc()
	magenta     = color.New(color.FgMagenta).SprintFunc()
	magentaBold = color.New(color.FgMagenta, color.Bold).SprintFunc()
	cyan        = color.New(color.FgCyan).SprintFunc()
	cyanBold    = color.New(color.FgCyan, color.Bold).SprintFunc()
	gray        = color.New(color.FgWhite).SprintFunc() // FgWhite = Gray in ANSI
	bold        = color.New(color.Bold).SprintFunc()
)

// ---------------------------------------------------------
// 2. Semantic Styles (The "Style..." API)
//    Use these for formatting specific types of data.
// ---------------------------------------------------------

// StyleError formats critical failure messages (Red).
func StyleError(msg string) string { return red(msg) }

// StyleSuccess formats success messages (Green).
func StyleSuccess(msg string) string { return green(msg) }

// StyleWarning formats non-critical warnings (Yellow).
func StyleWarning(msg string) string { return yellow(msg) }

// StyleHint formats helpful tips or suggestions (Cyan).
func StyleHint(msg string) string { return cyan(msg) }

// StyleNote formats neutral notes or annotations (Magenta).
func StyleNote(msg string) string { return magenta(msg) }

// StyleInfo formats status labels or properties (Magenta)
func StyleInfo(msg string) string { return magenta(msg) }

// StyleDebug formats low-level technical info (Gray).
func StyleDebug(msg string) string { return gray(msg) }

// StyleCommand formats shell commands or flags (Gray/Faint).
func StyleCommand(cmd string) string { return gray(cmd) }

// StyleAction formats verbs or active operations (Yellow).
func StyleAction(act string) string { return yellow(act) }

// StyleTitle
func StyleTitle(title string) string { return bold(cyan(title)) }

// StyleNumber formats counts, sizes, or IDs (Magenta).
func StyleNumber(num interface{}) string {
	return magenta(fmt.Sprintf("%v", num))
}

// StylePath formats file paths with context-aware coloring.
func StylePath(path string) string {
	// 1. Writable Overlay Images -> Bold Magenta
	if IsImg(path) {
		return magentaBold(path)
	}
	// 2. Read-Only Containers (SIF/SquashFS) -> Bold Cyan
	if IsSif(path) || IsSqf(path) {
		return cyanBold(path)
	}
	// 3. Standard System Paths -> Bold Blue
	return blueBold(path)
}

// StyleName formats names, identifiers, or keys (Yellow).
func StyleName(name string) string { return yellow(name) }

// StyleHighlight formats search matches or highlighted text (Yellow Bold).
func StyleHighlight(text string) string { return bold(yellow(text)) }

// FormatBytes formats bytes into human-readable format (e.g., "1.50 GB").
func FormatBytes(bytes int64) string {
	const (
		KB = 1024
		MB = KB * 1024
		GB = MB * 1024
		TB = GB * 1024
	)

	switch {
	case bytes >= TB:
		return fmt.Sprintf("%.2f TB", float64(bytes)/TB)
	case bytes >= GB:
		return fmt.Sprintf("%.2f GB", float64(bytes)/GB)
	case bytes >= MB:
		return fmt.Sprintf("%.2f MB", float64(bytes)/MB)
	case bytes >= KB:
		return fmt.Sprintf("%.2f KB", float64(bytes)/KB)
	default:
		return fmt.Sprintf("%d B", bytes)
	}
}

// ---------------------------------------------------------
// 3. Log Printers
//    High-level functions that print entire lines with tags.
// ---------------------------------------------------------

// PrintMessage prints a standard info message.
// Output: [CNT] Message...
func PrintMessage(format string, a ...interface{}) {
	if QuietMode {
		return
	}
	msg := fmt.Sprintf(format, a...)
	fmt.Fprintf(os.Stdout, "%s %s\n", projectPrefix, msg)
}

// PrintSuccess prints a success message with a Green tag.
// Output: [CNT][SUCCESS] Operation complete.
func PrintSuccess(format string, a ...interface{}) {
	if QuietMode {
		return
	}
	msg := fmt.Sprintf(format, a...)
	tag := StyleSuccess("[PASS]")
	fmt.Fprintf(os.Stdout, "%s%s %s\n", projectPrefix, tag, msg)
}

// PrintError prints an error message with a Red tag to Stderr.
// Output: [CNT][ERR] Something failed.
func PrintError(format string, a ...interface{}) {
	msg := fmt.Sprintf(format, a...)
	tag := StyleError("[ERR] ")
	fmt.Fprintf(os.Stderr, "%s%s %s\n", projectPrefix, tag, msg)
}

// PrintWarning prints a warning with a Yellow tag to Stderr.
// Output: [CNT][WARN] Disk is almost full.
func PrintWarning(format string, a ...interface{}) {
	msg := fmt.Sprintf(format, a...)
	tag := StyleWarning("[WARN]")
	fmt.Fprintf(os.Stderr, "%s%s %s\n", projectPrefix, tag, msg)
}

// PrintHint prints a helpful hint with a Cyan tag.
// Output: [CNT][HINT] Try running with --force.
func PrintHint(format string, a ...interface{}) {
	if QuietMode {
		return
	}
	msg := fmt.Sprintf(format, a...)
	tag := StyleHint("[HINT]")
	// Added extra spacing for alignment
	fmt.Fprintf(os.Stdout, "%s%s %s\n", projectPrefix, tag, msg)
}

// PrintNote prints a note with a Magenta tag.
// Output: [CNT][NOTE] This might take a while.
func PrintNote(format string, a ...interface{}) {
	if QuietMode {
		return
	}
	msg := fmt.Sprintf(format, a...)
	tag := StyleNote("[NOTE]")
	fmt.Fprintf(os.Stdout, "%s%s %s\n", projectPrefix, tag, msg)
}

// PrintDebug prints a debug message with a Gray tag (only if DebugMode is true).
// Output: [CNT] [DBG] Executing: rm -rf /tmp/foo
func PrintDebug(format string, a ...interface{}) {
	if DebugMode {
		msg := fmt.Sprintf(format, a...)
		tag := StyleDebug("[DBG] ")
		fmt.Fprintf(os.Stderr, "%s%s %s\n", projectPrefix, tag, msg)
	}
}

// ---------------------------------------------------------
// 4. Terminal Detection
// ---------------------------------------------------------

// IsInteractiveShell checks if stdout is connected to a TTY (interactive terminal).
// Returns true if the program is running in an interactive shell, false otherwise.
func IsInteractiveShell() bool {
	fileInfo, err := os.Stdout.Stat()
	if err != nil {
		return false
	}
	return (fileInfo.Mode() & os.ModeCharDevice) != 0
}
