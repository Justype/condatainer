package apptainer

import (
	"io"
	"path/filepath"
	"strings"

	"github.com/Justype/condatainer/internal/utils"
)

// DeduplicateBindPaths resolves, deduplicates, and filters child paths from bind directories.
// It handles formats: "/path", "/path:/container", "/path:/container:ro"
func DeduplicateBindPaths(paths []string) []string {
	// First pass: resolve all paths and deduplicate by host path
	seen := make(map[string]bool)
	resolved := make([]string, 0, len(paths))
	for _, bind := range paths {
		if bind == "" {
			continue
		}
		// Parse bind format
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

		// Skip if host path already seen
		if seen[absHostPath] {
			continue
		}
		seen[absHostPath] = true

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
		resolved = append(resolved, resolvedBind)
	}

	// Second pass: filter out child paths covered by parent paths
	filtered := make([]string, 0, len(resolved))
	for _, bind := range resolved {
		parts := strings.SplitN(bind, ":", 3)
		hostPath := parts[0]

		isChild := false
		for _, otherBind := range resolved {
			otherParts := strings.SplitN(otherBind, ":", 3)
			otherHostPath := otherParts[0]
			if hostPath == otherHostPath {
				continue
			}
			rel, err := filepath.Rel(otherHostPath, hostPath)
			if err == nil && rel != "." && !strings.HasPrefix(rel, "..") {
				isChild = true
				break
			}
		}
		if !isChild {
			filtered = append(filtered, bind)
		}
	}

	return filtered
}

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
func Exec(imagePath string, command []string, opts *ExecOptions) error {
	if opts == nil {
		opts = &ExecOptions{}
	}

	args := []string{"exec"}

	// Deduplicate and resolve bind mounts
	for _, bind := range DeduplicateBindPaths(opts.Bind) {
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

	return runApptainerWithOutput("exec", imagePath, false, opts.HideOutput, opts.Stdin, args...)
}
