package apptainer

import (
	"context"
	"io"
	"strings"

	"github.com/Justype/condatainer/internal/logging"
)

// ExecOptions contains options for executing commands in a container
type ExecOptions struct {
	Bin        Bin       // The apptainer to run, from Normal or Fakeroot
	Bind       []string  // Bind mounts (format: "/host/path:/container/path")
	Overlay    []string  // Overlay images to use
	Fakeroot   bool      // Run with fakeroot
	Env        []string  // Environment variables to set (format: "KEY=VALUE")
	Additional []string  // Additional flags to pass to apptainer exec
	Stdin      io.Reader // Custom stdin reader (optional, defaults to os.Stdin)
	Stdout     io.Writer // Redirect stdout (optional; nil = discard)
	Stderr     io.Writer // Redirect stderr (optional; nil = discard)
}

// Exec executes a command inside a container
// Note: Bind paths and GPU flags should be processed by the caller (e.g., via container.Setup())
func Exec(ctx context.Context, imagePath string, command []string, opts *ExecOptions) error {
	if opts == nil {
		opts = &ExecOptions{}
	}
	if opts.Bin.Path == "" {
		return ErrNotResolved
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

	// Add additional flags (caller should have already included GPU flags if needed)
	args = append(args, opts.Additional...)

	args = append(args, imagePath)
	args = append(args, command...)

	logging.FromContext(ctx).Debug("executing in container", "image", imagePath, "command", strings.Join(command, " "))

	return runApptainerWithOutput(ctx, opts.Bin, "exec", imagePath, false, opts.Stdin, opts.Stdout, opts.Stderr, envPrefixed(opts.Bin, opts.Env), args...)
}

// envPrefixed rewrites KEY=VALUE settings as APPTAINERENV_KEY=VALUE for the
// apptainer process's own environment.
//
// Not --env: that flag is parsed with CSV rules, so a comma or a quote in a
// value — routine in the signed download URL an #INPUT: asks for — fails the
// launch outright. The prefixed form passes values through verbatim, and keeps
// them out of the command line where any user's ps would see them.
func envPrefixed(bin Bin, settings []string) []string {
	if len(settings) == 0 {
		return nil
	}
	prefix := "APPTAINERENV_"
	if bin.IsSingularity() {
		prefix = "SINGULARITYENV_"
	}
	out := make([]string, 0, len(settings))
	for _, setting := range settings {
		key, value, ok := strings.Cut(setting, "=")
		if !ok || key == "" {
			continue
		}
		out = append(out, prefix+key+"="+value)
	}
	return out
}
