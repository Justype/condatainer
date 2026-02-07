package instance

import (
	"context"
	"fmt"

	"github.com/Justype/condatainer/internal/apptainer"
	"github.com/Justype/condatainer/internal/config"
	"github.com/Justype/condatainer/internal/container"
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
	fakeroot := container.AutoEnableFakeroot(setupResult.LastImg, options.WritableImg, setupResult.Fakeroot)

	// Build apptainer instance start options
	opts := &apptainer.InstanceStartOptions{
		Bind:       setupResult.BindPaths,
		Overlay:    setupResult.OverlayArgs,
		Fakeroot:   fakeroot,
		Additional: setupResult.ApptainerFlags,
	}

	// Debug output
	if utils.DebugMode {
		utils.PrintDebug("[INSTANCE]Instance name: %s", options.Name)
		utils.PrintDebug("[INSTANCE]Overlays: %v", setupResult.Overlays)
		utils.PrintDebug("[INSTANCE]Bind paths: %v", setupResult.BindPaths)
		utils.PrintDebug("[INSTANCE]Env list: %v", setupResult.EnvList)
	}

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
		utils.PrintWarning("Failed to save instance state: %v", err)
		utils.PrintNote("Environment variables may not work correctly with 'instance exec'")
	}

	return nil
}
