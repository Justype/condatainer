package utils

import (
	"fmt"
	"os"

	"github.com/fatih/color"
)

// Global Debug flag
var DebugMode bool = false

// --- 1. Base Color Functions ---
// Use these when you need raw colored strings
var (
	BlueText    = color.New(color.FgHiBlue).SprintFunc()
	RedText     = color.New(color.FgHiRed).SprintFunc()
	YellowText  = color.New(color.FgHiYellow).SprintFunc()
	GreenText   = color.New(color.FgHiGreen).SprintFunc()
	CyanText    = color.New(color.FgHiCyan).SprintFunc()
	MagentaText = color.New(color.FgHiMagenta).SprintFunc()
)

// --- 2. Semantic Helpers (The "Business Logic" of Colors) ---
// Use these in your logic code to ensure consistency

// StylePath formats file paths, URLs, and shell commands (Blue)
func StylePath(text string) string {
	return BlueText(text)
}

// StyleName formats app names, versions, and key variables (Yellow)
func StyleName(text string) string {
	return YellowText(text)
}

// StyleNumber formats counts and integers (Magenta)
func StyleNumber(n interface{}) string {
	// Handles int, float, etc.
	return MagentaText(fmt.Sprintf("%v", n))
}

// --- 3. Printing Functions ---

// PrintMessage prints standard info: [CNT] message
func PrintMessage(format string, a ...interface{}) {
	msg := fmt.Sprintf(format, a...)
	fmt.Printf("[CNT] %s\n", msg)
}

// PrintDebug prints yellow debug info: [CNT][DEBUG] message
func PrintDebug(format string, a ...interface{}) {
	if DebugMode {
		msg := fmt.Sprintf(format, a...)
		fmt.Printf("[CNT][%s] %s\n", YellowText("DEBUG"), msg)
	}
}

// PrintNote prints blue note info: [CNT][NOTE] message
func PrintNote(format string, a ...interface{}) {
	msg := fmt.Sprintf(format, a...)
	fmt.Printf("[CNT][%s] %s\n", BlueText("NOTE"), msg)
}

// PrintWarning prints yellow warning info to Stderr: [CNT][WARNING] message
func PrintWarning(format string, a ...interface{}) {
	msg := fmt.Sprintf(format, a...)
	fmt.Fprintf(os.Stderr, "[CNT][%s] %s\n", YellowText("WARNING"), msg)
}

// PrintError prints red error info to Stderr: [CNT][ERROR] message
func PrintError(format string, a ...interface{}) {
	msg := fmt.Sprintf(format, a...)
	fmt.Fprintf(os.Stderr, "[CNT][%s] %s\n", RedText("ERROR"), msg)
}

// PrintSuccess prints green success info: [CNT][SUCCESS] message
func PrintSuccess(format string, a ...interface{}) {
	msg := fmt.Sprintf(format, a...)
	fmt.Printf("[CNT][%s] %s\n", GreenText("SUCCESS"), msg)
}
