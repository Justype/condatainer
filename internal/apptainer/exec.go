package apptainer

import (
	"strings"

	"github.com/Justype/condatainer/internal/utils"
)

// ExecOptions contains options for executing commands in a container
type ExecOptions struct {
	Bind       []string // Bind mounts (format: "/host/path:/container/path")
	Overlay    []string // Overlay images to use
	Fakeroot   bool     // Run with fakeroot
	Env        []string // Environment variables to set (format: "KEY=VALUE")
	HideOutput bool     // Hide output by redirecting to /dev/null
	Additional []string // Additional flags to pass to apptainer exec
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

	// Add fakeroot flag
	if opts.Fakeroot {
		args = append(args, "--fakeroot")
	}

	// Add environment variables
	for _, env := range opts.Env {
		args = append(args, "--env", env)
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

	return runApptainerWithOutput("exec", imagePath, false, opts.HideOutput, args...)
}
