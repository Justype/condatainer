package cmd

import (
	"os"
	"strings"

	"github.com/spf13/cobra"
)

// CommonFlags holds the common flags used by exec, instance start, and instance exec
type CommonFlags struct {
	Overlays    []string
	WritableImg bool
	EnvSettings []string
	BaseImage   string
	BindPaths   []string
	Fakeroot    bool
}

// RegisterCommonFlags registers common flags on a cobra command
func RegisterCommonFlags(cmd *cobra.Command, flags *CommonFlags) {
	cmd.Flags().StringSliceVarP(&flags.Overlays, "overlay", "o", nil, "overlay file to mount (can be used multiple times)")
	cmd.Flags().BoolVarP(&flags.WritableImg, "writable", "w", false, "mount .img overlays as writable (default: read-only)")
	cmd.Flags().BoolVar(&flags.WritableImg, "writable-img", false, "Alias for --writable")
	cmd.Flags().StringSliceVar(&flags.EnvSettings, "env", nil, "set environment variable 'KEY=VALUE' (can be used multiple times)")
	cmd.Flags().StringVarP(&flags.BaseImage, "base-image", "b", "", "base image to use instead of default")
	cmd.Flags().StringSliceVar(&flags.BindPaths, "bind", nil, "bind path 'HOST:CONTAINER' (can be used multiple times)")
	cmd.Flags().BoolVarP(&flags.Fakeroot, "fakeroot", "f", false, "run container with fakeroot privileges")

	// Register completions
	cmd.RegisterFlagCompletionFunc("overlay", overlayFlagCompletion(true, true))
	cmd.RegisterFlagCompletionFunc("base-image", baseImageFlagCompletion())

	// Stop flag parsing after the first positional argument
	cmd.Flags().SetInterspersed(false)

	// Allow unknown flags to pass through (for apptainer)
	cmd.FParseErrWhitelist.UnknownFlags = true
}

// KnownFlags returns a map of all known condatainer flags
func KnownFlags() map[string]bool {
	return map[string]bool{
		"--overlay": true, "-o": true,
		"--writable": true, "--writable-img": true, "-w": true,
		"--env":        true,
		"--bind":       true,
		"--base-image": true, "-b": true,
		"--fakeroot": true, "-f": true,
		"--debug": true,
		"--local": true,
		"--quiet": true, "-q": true,
		"--yes": true, "-y": true,
	}
}

// ParseCommandArgs parses arguments from os.Args after a given subcommand
// Returns: commands (positional args), apptainerFlags (unknown flags for apptainer)
func ParseCommandArgs(subcommand string) ([]string, []string) {
	var commands, apptainerFlags []string

	// Find where subcommand appears in os.Args
	cmdIdx := -1
	for i, arg := range os.Args {
		if arg == subcommand {
			cmdIdx = i
			break
		}
	}

	if cmdIdx == -1 {
		// Fallback: no args found
		return commands, apptainerFlags
	}

	knownFlags := KnownFlags()

	// Parse from after subcommand in os.Args
	commandStarted := false
	for i := cmdIdx + 1; i < len(os.Args); i++ {
		arg := os.Args[i]

		// Once we hit the command, everything after is part of the command
		if commandStarted {
			commands = append(commands, arg)
			continue
		}

		// Skip known flags (already handled by cobra)
		if knownFlags[arg] || isKnownFlagWithEquals(knownFlags, arg) {
			// Value flags (space-separated)
			if knownFlags[arg] && needsValue(arg) && i+1 < len(os.Args) {
				i++ // Skip value
			}
			continue
		}

		// Unknown flag → must use --flag=value format for apptainer pass-through
		if strings.HasPrefix(arg, "-") {
			apptainerFlags = append(apptainerFlags, arg)
			continue
		}

		// First non-flag argument starts the command
		commandStarted = true
		commands = append(commands, arg)
	}

	return commands, apptainerFlags
}

// ParseInstanceArgs parses arguments for instance start/run commands
// Returns: instanceName, apptainerFlags
func ParseInstanceArgs(subcommand string) (string, []string) {
	var instanceName string
	var apptainerFlags []string

	// Find where "instance <subcommand>" appears in os.Args
	cmdIdx := -1
	for i := 0; i < len(os.Args)-1; i++ {
		if os.Args[i] == "instance" && os.Args[i+1] == subcommand {
			cmdIdx = i + 1 // Point to the subcommand
			break
		}
	}

	if cmdIdx == -1 {
		// Fallback: no args found
		return instanceName, apptainerFlags
	}

	knownFlags := KnownFlags()

	// Parse from after subcommand in os.Args
	for i := cmdIdx + 1; i < len(os.Args); i++ {
		arg := os.Args[i]

		// Skip known flags (already handled by cobra)
		if knownFlags[arg] || isKnownFlagWithEquals(knownFlags, arg) {
			// Value flags (space-separated)
			if knownFlags[arg] && needsValue(arg) && i+1 < len(os.Args) {
				i++ // Skip value
			}
			continue
		}

		// Unknown flag → pass through to apptainer
		if strings.HasPrefix(arg, "-") {
			apptainerFlags = append(apptainerFlags, arg)
			continue
		}

		// First (and only) non-flag argument is the instance name
		if instanceName == "" {
			instanceName = arg
		}
		// Ignore any additional positional arguments - overlays must use -o flag
	}

	return instanceName, apptainerFlags
}

// ParseInstanceExecArgs parses arguments for instance exec command
// Returns: instanceName, command, envSettings, apptainerFlags
func ParseInstanceExecArgs() (string, []string, []string, []string) {
	var instanceName string
	var command []string
	var envSettings []string
	var apptainerFlags []string

	// Parse os.Args directly - find where "instance exec" appears
	execIdx := -1
	for i := 0; i < len(os.Args)-1; i++ {
		if os.Args[i] == "instance" && os.Args[i+1] == "exec" {
			execIdx = i + 1 // Point to "exec"
			break
		}
	}

	if execIdx == -1 || execIdx+1 >= len(os.Args) {
		// Fallback: no args found
		return instanceName, command, envSettings, apptainerFlags
	}

	// Known condatainer flags
	knownFlags := map[string]bool{
		"--debug": true,
		"--local": true,
		"--quiet": true, "-q": true,
		"--yes": true, "-y": true,
		"--env": true,
	}

	// Parse from after "exec" in os.Args
	foundInstanceName := false
	commandStarted := false
	for i := execIdx + 1; i < len(os.Args); i++ {
		arg := os.Args[i]

		// Once command started, everything else is part of the command
		if commandStarted {
			command = append(command, arg)
			continue
		}

		// Handle --env flag
		if arg == "--env" && i+1 < len(os.Args) {
			i++
			envSettings = append(envSettings, os.Args[i])
			continue
		}

		// Handle --env=VALUE format
		if strings.HasPrefix(arg, "--env=") {
			envValue := strings.TrimPrefix(arg, "--env=")
			envSettings = append(envSettings, envValue)
			continue
		}

		// Skip other known global flags before instance name
		if knownFlags[arg] {
			continue
		}

		// Handle unknown flags (apptainer pass-through)
		if strings.HasPrefix(arg, "-") && !foundInstanceName {
			// Unknown flag before instance name - pass to apptainer
			apptainerFlags = append(apptainerFlags, arg)
			continue
		}

		// Once we have the instance name, next non-flag starts the command
		if foundInstanceName {
			// This is the first command argument
			commandStarted = true
			command = append(command, arg)
			continue
		}

		// First non-flag argument is the instance name
		if !strings.HasPrefix(arg, "-") {
			instanceName = arg
			foundInstanceName = true
		}
	}

	return instanceName, command, envSettings, apptainerFlags
}
