package instance

import (
	"context"
	"fmt"

	"log/slog"

	"github.com/Justype/condatainer/internal/apptainer"
	"github.com/Justype/condatainer/internal/config"
	"github.com/Justype/condatainer/internal/container"
	"github.com/Justype/condatainer/internal/logging"
	"github.com/Justype/condatainer/internal/utils"
)

// Options holds configuration for starting an instance
type Options struct {
	Name           string   // Instance name
	Overlays       []string // Overlay paths
	WritableImg    bool     // Whether .img overlays should be writable
	EnvSettings    []string // Environment variables
	BindPaths      []string // Bind mount paths
	ApptainerFlags []string // Additional apptainer flags
	Fakeroot       bool     // Run with fakeroot
	BaseImage      string   // Base image path
	ApptainerBin   string   // Path to apptainer binary
}

// Start starts a named instance with the given options
func Start(ctx context.Context, options Options) error {
	if options.ApptainerBin == "" {
		options.ApptainerBin = config.Global.ApptainerBin
	}

	// Ensure apptainer binary is set
	if err := apptainer.SetBin(options.ApptainerBin); err != nil {
		return err
	}

	baseImage := options.BaseImage
	if baseImage == "" {
		baseImage = config.GetBaseImage()
	}

	if !utils.FileExists(baseImage) {
		return fmt.Errorf("base image not found: %s", baseImage)
	}

	// Use shared container setup logic
	setupResult, err := container.Setup(container.SetupConfig{
		Overlays:       options.Overlays,
		WritableImg:    options.WritableImg,
		EnvSettings:    options.EnvSettings,
		BindPaths:      options.BindPaths,
		Fakeroot:       options.Fakeroot,
		ApptainerFlags: options.ApptainerFlags,
	})
	if err != nil {
		return err
	}

	// Auto-enable fakeroot if needed for writable .img
	fakeroot, _ := container.AutoEnableFakeroot(setupResult.LastImg, options.WritableImg, setupResult.Fakeroot)

	// Build apptainer instance start options
	opts := &apptainer.InstanceStartOptions{
		Bind:       setupResult.BindPaths,
		Overlay:    setupResult.OverlayArgs,
		Fakeroot:   fakeroot,
		Additional: setupResult.ApptainerFlags,
	}

	// Debug output
	slog.Default().Debug("starting instance",
		"name", options.Name, "overlays", setupResult.Overlays,
		"bindPaths", setupResult.BindPaths, "envList", setupResult.EnvList)

	// Start the instance using apptainer package
	if err := apptainer.InstanceStart(ctx, baseImage, options.Name, opts); err != nil {
		return err
	}

	// Save instance state for use with instance exec
	state := &State{
		Name:      options.Name,
		Env:       setupResult.EnvList,
		Overlays:  setupResult.Overlays,
		BindPaths: setupResult.BindPaths,
		BaseImage: baseImage,
		Notes:     setupResult.EnvNotes,
	}

	if err := saveState(state); err != nil {
		logging.FromContext(ctx).Warn("failed to save instance state", "err", err)
		logging.FromContext(ctx).Info("environment variables may not work correctly with 'instance exec'", "kind", "note")
	}

	return nil
}
