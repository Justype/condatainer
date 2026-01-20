package apptainer

import (
	"fmt"
	"os/exec"

	"github.com/Justype/condatainer/internal/utils"
)

// InstanceStartOptions contains options for starting a container instance
type InstanceStartOptions struct {
	Bind       []string // Bind mounts (format: "/host/path:/container/path")
	Overlay    []string // Overlay images to use
	Writable   bool     // Mount image as writable
	Fakeroot   bool     // Run with fakeroot
	NoHome     bool     // Do not mount home directory
	ContainAll bool     // Contain all user environments
	Env        []string // Environment variables to set (format: "KEY=VALUE")
	Additional []string // Additional flags to pass to apptainer instance start
}

// InstanceStart starts a named container instance
func InstanceStart(imagePath, instanceName string, opts *InstanceStartOptions) error {
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
	if opts.ContainAll {
		args = append(args, "--contain")
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

	return runApptainer("instance start", imagePath, false, args...)
}

// InstanceStop stops a running container instance
func InstanceStop(instanceName string, force bool) error {
	args := []string{"instance", "stop"}

	if force {
		args = append(args, "--force")
	}

	args = append(args, instanceName)

	utils.PrintDebug("Stopping instance %s", utils.StyleInfo(instanceName))

	return runApptainer("instance stop", instanceName, false, args...)
}

// InstanceStopAll stops all running container instances
func InstanceStopAll(force bool) error {
	args := []string{"instance", "stop"}

	if force {
		args = append(args, "--force")
	}

	args = append(args, "--all")

	utils.PrintDebug("Stopping all instances")

	return runApptainer("instance stop", "all", false, args...)
}

// InstanceList lists running container instances
// Returns the output from apptainer instance list
func InstanceList() (string, error) {
	cmd := exec.Command(apptainerCmd, "instance", "list")
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
