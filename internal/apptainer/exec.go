package apptainer

import (
	"strings"

	"github.com/Justype/condatainer/internal/utils"
)

// ExecOptions contains options for executing commands in a container
type ExecOptions struct {
	Bind          []string // Bind mounts (format: "/host/path:/container/path")
	Overlay       []string // Overlay images to use
	Writable      bool     // Mount image as writable (requires sandbox or overlay)
	Fakeroot      bool     // Run with fakeroot
	NoHome        bool     // Do not mount home directory
	Home          string   // Custom home directory inside container
	ContainAll    bool     // Contain all user environments
	CleanEnv      bool     // Clean environment before running
	Env           []string // Environment variables to set (format: "KEY=VALUE")
	Pwd           string   // Working directory inside container
	CaptureOutput bool     // Capture output instead of streaming
	Additional    []string // Additional flags to pass to apptainer exec
}

// Exec executes a command inside a container
func Exec(imagePath string, command []string, opts *ExecOptions) error {
	if opts == nil {
		opts = &ExecOptions{}
	}

	args := []string{"exec"}

	// Add bind mounts
	for _, bind := range opts.Bind {
		args = append(args, "--bind", bind)
	}

	// Add overlays
	for _, overlay := range opts.Overlay {
		args = append(args, "--overlay", overlay)
	}

	// Add flags
	if opts.Writable {
		args = append(args, "--writable")
	}
	if opts.Fakeroot {
		args = append(args, "--fakeroot")
	}
	if opts.NoHome {
		args = append(args, "--no-home")
	}
	if opts.Home != "" {
		args = append(args, "--home", opts.Home)
	}
	if opts.ContainAll {
		args = append(args, "--containall")
	}
	if opts.CleanEnv {
		args = append(args, "--cleanenv")
	}

	// Add environment variables
	for _, env := range opts.Env {
		args = append(args, "--env", env)
	}

	// Add working directory
	if opts.Pwd != "" {
		args = append(args, "--pwd", opts.Pwd)
	}

	// Add detected GPU flags first, then any additional user flags.
	args = append(args, DetectGPUFlags()...)
	args = append(args, opts.Additional...)

	// Add image path
	args = append(args, imagePath)

	// Add command
	args = append(args, command...)

	utils.PrintDebug("Executing in container %s: %s",
		utils.StylePath(imagePath),
		utils.StyleAction(strings.Join(command, " ")))

	return runApptainer("exec", imagePath, opts.CaptureOutput, args...)
}
