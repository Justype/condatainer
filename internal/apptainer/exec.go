package apptainer

import (
	"path/filepath"
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

	// Deduplicate and resolve bind mounts
	seen := make(map[string]bool)
	for _, bind := range opts.Bind {
		// Parse bind format: "/host/path", "/host/path:/container/path", or "/host/path:/container/path:ro"
		parts := strings.SplitN(bind, ":", 3)
		hostPath := parts[0]

		// Resolve host path to absolute and follow symlinks
		absHostPath, err := filepath.Abs(hostPath)
		if err != nil {
			absHostPath = hostPath
		}
		if realPath, err := filepath.EvalSymlinks(absHostPath); err == nil {
			absHostPath = realPath
		}

		// Reconstruct the bind string with resolved host path
		var resolvedBind string
		switch len(parts) {
		case 3:
			resolvedBind = absHostPath + ":" + parts[1] + ":" + parts[2]
		case 2:
			resolvedBind = absHostPath + ":" + parts[1]
		default:
			resolvedBind = absHostPath
		}

		// Skip if already seen
		if seen[resolvedBind] {
			continue
		}
		seen[resolvedBind] = true
		args = append(args, "--bind", resolvedBind)
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
