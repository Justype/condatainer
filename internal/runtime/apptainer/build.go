package apptainer

import (
	"context"
	"os"

	"github.com/Justype/condatainer/internal/logging"
)

// BuildOptions contains options for building a container image
type BuildOptions struct {
	Force      bool     // Force overwrite of existing image
	NoCleanup  bool     // Do not clean up bundle after failed build
	Additional []string // Additional flags to pass to apptainer build
}

// Build builds a container image from a definition file
// Always uses --fakeroot (required for non-root users on HPC systems)
func Build(ctx context.Context, imagePath, defFile string, opts *BuildOptions) error {
	if opts == nil {
		opts = &BuildOptions{}
	}

	args := []string{"build"}

	args = append(args, "--fakeroot")

	// Add optional flags
	if opts.Force {
		args = append(args, "--force")
	}
	if opts.NoCleanup {
		args = append(args, "--no-cleanup")
	}

	// Add user-provided extras.
	// args = append(args, DetectGPUFlags()...) // Do not detect the gpu when building
	args = append(args, opts.Additional...)

	// Add image path and definition file
	args = append(args, imagePath, defFile)

	logging.FromContext(ctx).Debug("building container", "image", imagePath, "definition", defFile)

	return runApptainerWithOutput(ctx, "build", imagePath, false, os.Stdin, os.Stdout, os.Stderr, nil, args...)
}
