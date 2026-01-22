package apptainer

import (
	"fmt"
	"os"
	"path/filepath"
	"runtime"

	"github.com/Justype/condatainer/internal/config"
	"github.com/Justype/condatainer/internal/utils"
)

// PrebuiltBaseImagePlatforms lists architectures with prebuilt base images available
var PrebuiltBaseImagePlatforms = map[string]bool{
	"x86_64": true,
	// "aarch64": false, // Not available yet
}

// EnsureApptainer checks if apptainer binary is available and configured.
// Returns an error if apptainer cannot be found.
func EnsureApptainer() error {
	return SetBin(config.Global.ApptainerBin)
}

// EnsureBaseImage ensures the base image exists, downloading or building it if necessary.
// This is the main entry point for base image initialization.
//
// The function follows this logic:
//  1. Create images directory if it doesn't exist
//  2. Check if apptainer is available
//  3. If base image already exists, return immediately
//  4. Try to download prebuilt base image (unless debug mode)
//  5. If download fails, try to build locally from definition file
func EnsureBaseImage() error {
	// Step 1: Ensure images directory exists
	if err := os.MkdirAll(config.Global.ImagesDir, utils.PermDir); err != nil {
		return fmt.Errorf("failed to create images directory: %w", err)
	}

	// Step 2: Ensure apptainer is available
	if err := EnsureApptainer(); err != nil {
		return err
	}

	// Step 3: Check if base image already exists
	if utils.FileExists(config.Global.BaseImage) {
		utils.PrintDebug("Base image already exists at %s", utils.StylePath(config.Global.BaseImage))
		return nil
	}

	// Step 4: Try to download prebuilt base image
	// In debug mode, always build locally
	if !config.Global.Debug {
		utils.PrintMessage("Base image not found. Downloading base image...")
		if err := tryDownloadPrebuiltBaseImage(); err == nil {
			return nil
		} else {
			utils.PrintDebug("Failed to download prebuilt base image: %v", err)
		}
	} else {
		utils.PrintDebug("Base image not found. Proceeding to build locally due to debug mode.")
	}

	// Step 5: Try to build base image locally
	utils.PrintMessage("Trying to build base image locally...")

	if err := ensureBaseDef(); err != nil {
		return err
	}

	baseDefPath := filepath.Join(config.Global.BuildScriptsDir, "base_image.def")
	if err := Build(config.Global.BaseImage, baseDefPath, nil); err != nil {
		return fmt.Errorf("building base image failed: %w", err)
	}

	utils.PrintSuccess("Base image built.")

	// Set permissions
	if err := os.Chmod(config.Global.BaseImage, utils.PermDir); err != nil {
		utils.PrintWarning("Failed to set permissions on base image: %v", err)
	}

	utils.PrintNote("You can run %s to free up space.", utils.StyleAction("apptainer cache clean"))

	return nil
}

// tryDownloadPrebuiltBaseImage attempts to download a prebuilt base image for the current architecture.
// Returns nil on success, error otherwise.
func tryDownloadPrebuiltBaseImage() error {
	// Map Go arch names to the naming convention used in releases
	arch := runtime.GOARCH
	switch arch {
	case "amd64":
		arch = "x86_64"
	case "arm64":
		arch = "aarch64"
	}

	if !PrebuiltBaseImagePlatforms[arch] {
		return fmt.Errorf("pre-built base image not available for architecture: %s", arch)
	}

	url := fmt.Sprintf("https://github.com/%s/releases/download/v%s/base_image_%s.sif",
		config.GitHubRepo, config.VERSION, arch)

	utils.PrintMessage("Attempting to download pre-built base image ...")

	if err := utils.DownloadExecutable(url, config.Global.BaseImage); err != nil {
		return err
	}

	utils.PrintSuccess("Pre-built base image downloaded.")
	return nil
}

// ensureBaseDef ensures the base image definition file exists, downloading it if necessary.
func ensureBaseDef() error {
	baseDefPath := filepath.Join(config.Global.BuildScriptsDir, "base_image.def")

	if utils.FileExists(baseDefPath) {
		return nil
	}

	utils.PrintMessage("Base image definition file not found. Downloading...")

	// Ensure build-scripts directory exists
	if err := os.MkdirAll(config.Global.BuildScriptsDir, utils.PermDir); err != nil {
		return fmt.Errorf("failed to create build-scripts directory: %w", err)
	}

	url := fmt.Sprintf("https://raw.githubusercontent.com/%s/main/build-scripts/base_image.def", config.GitHubRepo)

	if err := utils.DownloadFile(url, baseDefPath); err != nil {
		return fmt.Errorf("failed to download base image definition file: %w", err)
	}

	utils.PrintSuccess("Base image definition file downloaded.")
	return nil
}
