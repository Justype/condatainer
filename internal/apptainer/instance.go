package apptainer

import (
	"context"
	"fmt"
	"io"
)

// InstanceStartOptions contains options for starting an instance
type InstanceStartOptions struct {
	Bind       []string // Bind mounts (format: "/host/path:/container/path")
	Overlay    []string // Overlay images to use
	Fakeroot   bool     // Run with fakeroot
	Additional []string // Additional flags to pass to apptainer instance start
}

// InstanceStart starts a named Apptainer instance
func InstanceStart(ctx context.Context, baseImage string, instanceName string, opts *InstanceStartOptions) error {
	if opts == nil {
		opts = &InstanceStartOptions{}
	}

	args := []string{"instance", "start"}

	// Add fakeroot flag
	if opts.Fakeroot {
		args = append(args, "--fakeroot")
	}

	// Add bind mounts
	for _, bind := range opts.Bind {
		args = append(args, "--bind", bind)
	}

	// Add overlays
	for _, overlay := range opts.Overlay {
		args = append(args, "--overlay", overlay)
	}

	// Add additional flags (including GPU flags from caller)
	args = append(args, opts.Additional...)

	// Add base image and instance name
	args = append(args, baseImage, instanceName)

	return runApptainer(ctx, "instance start", baseImage, false, args...)
}

// InstanceStop stops a running instance
// If name is empty, all arguments (including instance name) should be in additionalFlags
func InstanceStop(ctx context.Context, name string, additionalFlags []string) error {
	args := []string{"instance", "stop"}

	// Add any additional flags first
	args = append(args, additionalFlags...)

	// Only add instance name if it's not empty (when name is empty, it's already in additionalFlags)
	if name != "" {
		args = append(args, name)
	}

	return runApptainer(ctx, "instance stop", "", false, args...)
}

// InstanceList lists all running instances
func InstanceList(ctx context.Context) error {
	args := []string{"instance", "list"}
	return runApptainer(ctx, "instance list", "", false, args...)
}

// InstanceStats shows statistics for a running instance
func InstanceStats(ctx context.Context, name string) error {
	args := []string{"instance", "stats", name}
	return runApptainer(ctx, "instance stats", "", false, args...)
}

// InstanceExecOptions contains options for executing commands in an instance
type InstanceExecOptions struct {
	Env        []string  // Environment variables to set (format: "KEY=VALUE")
	Additional []string  // Additional flags to pass to apptainer exec
	Stdin      io.Reader // Custom stdin reader (optional, defaults to os.Stdin)
}

// InstanceExec executes a command in a running instance
func InstanceExec(ctx context.Context, instanceName string, command []string, opts *InstanceExecOptions) error {
	if opts == nil {
		opts = &InstanceExecOptions{}
	}

	args := []string{"exec"}

	// Add environment variables
	for _, env := range opts.Env {
		args = append(args, "--env", env)
	}

	// Add additional flags
	args = append(args, opts.Additional...)

	// Add instance URI
	instanceURI := fmt.Sprintf("instance://%s", instanceName)
	args = append(args, instanceURI)

	// Add command
	args = append(args, command...)

	return runApptainerWithOutput(ctx, "instance exec", instanceURI, false, false, opts.Stdin, args...)
}
