package apptainer

import (
	"context"
	"fmt"
	"os/exec"

	"github.com/Justype/condatainer/internal/utils"
)

// InstanceStartOptions contains options for starting a container instance.
// These options control how an apptainer instance is started and what
// environment is provided to it.
type InstanceStartOptions struct {
	Bind       []string // Bind mounts (format: "/host/path:/container/path")
	Overlay    []string // Overlay images to use
	Writable   bool     // Mount image as writable
	Fakeroot   bool     // Run with fakeroot
	NoHome     bool     // Do not mount home directory
	Home       string   // Custom home directory inside container
	ContainAll bool     // Contain all user environments
	Env        []string // Environment variables to set (format: "KEY=VALUE")
	Additional []string // Additional flags to pass to apptainer instance start
}

// InstanceStart starts a named container instance
func InstanceStart(ctx context.Context, imagePath, instanceName string, opts *InstanceStartOptions) error {
	if opts == nil {
		opts = &InstanceStartOptions{}
	}

	args := []string{"instance", "start"}

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

	// Add environment variables
	for _, env := range opts.Env {
		args = append(args, "--env", env)
	}

	// Add additional flags
	args = append(args, opts.Additional...)

	// Add image path and instance name
	args = append(args, imagePath, instanceName)

	utils.PrintDebug("Starting instance %s from %s",
		utils.StyleInfo(instanceName),
		utils.StylePath(imagePath))

	return runApptainer(ctx, "instance start", imagePath, false, args...)
}

// InstanceStop stops a running container instance
func InstanceStop(ctx context.Context, instanceName string, force bool) error {
	args := []string{"instance", "stop"}

	if force {
		args = append(args, "--force")
	}

	args = append(args, instanceName)

	utils.PrintDebug("Stopping instance %s", utils.StyleInfo(instanceName))

	return runApptainer(ctx, "instance stop", instanceName, false, args...)
}

// InstanceStopAll stops all running container instances
func InstanceStopAll(ctx context.Context, force bool) error {
	args := []string{"instance", "stop"}

	if force {
		args = append(args, "--force")
	}

	args = append(args, "--all")

	utils.PrintDebug("Stopping all instances")

	return runApptainer(ctx, "instance stop", "all", false, args...)
}

// InstanceList lists running container instances
// Returns the output from apptainer instance list
func InstanceList(ctx context.Context) (string, error) {
	cmd := exec.CommandContext(ctx, apptainerCmd, "instance", "list")
	output, err := cmd.Output()

	if err != nil {
		return "", &ApptainerError{
			Op:      "instance list",
			Cmd:     fmt.Sprintf("%s instance list", apptainerCmd),
			Output:  string(output),
			BaseErr: err,
		}
	}

	return string(output), nil
}
