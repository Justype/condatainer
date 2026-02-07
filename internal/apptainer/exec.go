package apptainer

import (
	"context"
	"io"
	"strings"

	"github.com/Justype/condatainer/internal/utils"
)

// ExecOptions contains options for executing commands in a container
type ExecOptions struct {
	Bind       []string  // Bind mounts (format: "/host/path:/container/path")
	Overlay    []string  // Overlay images to use
	Fakeroot   bool      // Run with fakeroot
	Env        []string  // Environment variables to set (format: "KEY=VALUE")
	HideOutput bool      // Hide output by redirecting to /dev/null
	Additional []string  // Additional flags to pass to apptainer exec
	Stdin      io.Reader // Custom stdin reader (optional, defaults to os.Stdin)
}

// Exec executes a command inside a container
// Note: Bind paths and GPU flags should be processed by the caller (e.g., via container.Setup())
func Exec(ctx context.Context, imagePath string, command []string, opts *ExecOptions) error {
	if opts == nil {
		opts = &ExecOptions{}
	}

	args := []string{"exec"}

	for _, bind := range opts.Bind {
		args = append(args, "--bind", bind)
	}
	for _, overlay := range opts.Overlay {
		args = append(args, "--overlay", overlay)
	}
	if opts.Fakeroot {
		args = append(args, "--fakeroot")
	}
	for _, env := range opts.Env {
		args = append(args, "--env", env)
	}

	// Add additional flags (caller should have already included GPU flags if needed)
	args = append(args, opts.Additional...)

	args = append(args, imagePath)
	args = append(args, command...)

	utils.PrintDebug("Executing in container %s: %s",
		utils.StylePath(imagePath),
		utils.StyleAction(strings.Join(command, " ")))

	return runApptainerWithOutput(ctx, "exec", imagePath, false, opts.HideOutput, opts.Stdin, args...)
}
