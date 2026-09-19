package apptainer

import (
	"context"
	"os"

	"github.com/Justype/condatainer/internal/logging"
)

// BuildOptions contains options for building a container image
type BuildOptions struct {
	Bin        Bin      // The apptainer to run, from ForBuild
	Force      bool     // Force overwrite of existing image
	NoCleanup  bool     // Do not clean up bundle after failed build
	Sandbox    bool     // Write a sandbox directory instead of a SIF
	TmpDir     string   // APPTAINER_TMPDIR for this build; "" leaves the environment alone
	Additional []string // Additional flags to pass to apptainer build
}

// Build builds a container image from a definition file
// Always uses --fakeroot (required for non-root users on HPC systems)
func Build(ctx context.Context, imagePath, defFile string, opts *BuildOptions) error {
	if opts == nil {
		opts = &BuildOptions{}
	}
	if opts.Bin.Path == "" {
		return ErrNotResolved
	}

	args := []string{"build"}

	args = append(args, "--fakeroot")

	if opts.Sandbox {
		args = append(args, "--sandbox")
	}

	// Add optional flags
	if opts.Force {
		args = append(args, "--force")
	}
	if opts.NoCleanup {
		args = append(args, "--no-cleanup")
	}

	// Add user-provided extras. GPU is never detected for a build, requested or not.
	args = append(args, opts.Additional...)

	// Add image path and definition file
	args = append(args, imagePath, defFile)

	logging.FromContext(ctx).Debug("building container", "image", imagePath, "definition", defFile)

	// Left unset, APPTAINER_TMPDIR inherits TMPDIR, which on a scheduler is
	// routinely network scratch. A sandbox assembles beside its destination
	// rather than here, so this is the smaller scratch, and it must exist.
	// procEnv is appended to the inherited environment, so a later entry wins.
	var procEnv []string
	if opts.TmpDir != "" {
		procEnv = append(procEnv, "APPTAINER_TMPDIR="+opts.TmpDir)
	}

	return runApptainerWithOutput(ctx, opts.Bin, "build", imagePath, false, os.Stdin, os.Stdout, os.Stderr, procEnv, 0, args...)
}
