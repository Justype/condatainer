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
//  1. Ensure apptainer is available
//  2. Search all image paths for existing base image
//  3. If not found, try to download prebuilt base image (unless debug mode)
//  4. If download fails, try to build locally from definition file
func EnsureBaseImage() error {
	// Step 1: Ensure apptainer is available
	if err := EnsureApptainer(); err != nil {
		return err
	}

	// Step 2: Search all image paths for existing base image
	if existingPath := config.FindBaseImage(); existingPath != "" {
		utils.PrintDebug("Base image found at %s", utils.StylePath(existingPath))
		return nil
	}

	// Get the path where we'll write the new base image
	baseImagePath, err := config.GetBaseImageWritePath()
	if err != nil {
		return fmt.Errorf("no writable directory for base image: %w", err)
	}

	// Ensure the parent directory exists
	if err := os.MkdirAll(filepath.Dir(baseImagePath), utils.PermDir); err != nil {
		return fmt.Errorf("failed to create images directory: %w", err)
	}

	// Step 3: Try to download prebuilt base image
	// In debug mode, always build locally
	if !config.Global.Debug {
		utils.PrintMessage("Base image not found. Downloading base image...")
		if err := tryDownloadPrebuiltBaseImage(baseImagePath); err == nil {
			return nil
		} else {
			utils.PrintDebug("Failed to download prebuilt base image: %v", err)
		}
	} else {
		utils.PrintDebug("Base image not found. Proceeding to build locally due to debug mode.")
	}

	// Step 4: Try to build base image locally
	utils.PrintMessage("Trying to build base image locally...")

	baseDefPath, err := ensureBaseDef()
	if err != nil {
		return err
	}

	if err := Build(baseImagePath, baseDefPath, nil); err != nil {
		return fmt.Errorf("building base image failed: %w", err)
	}

	utils.PrintSuccess("Base image built at %s", utils.StylePath(baseImagePath))

	// Set permissions
	if err := os.Chmod(baseImagePath, utils.PermDir); err != nil {
		utils.PrintWarning("Failed to set permissions on base image: %v", err)
	}

	utils.PrintNote("You can run %s to free up space.", utils.StyleAction("apptainer cache clean"))

	return nil
}

// tryDownloadPrebuiltBaseImage attempts to download a prebuilt base image for the current architecture.
// Downloads to the specified destPath. Returns nil on success, error otherwise.
func tryDownloadPrebuiltBaseImage(destPath string) error {
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

	url := fmt.Sprintf("%s/base_image_%s.sif", config.PrebuiltBaseURL, arch)

	utils.PrintMessage("Attempting to download pre-built base image ...")

	if err := utils.DownloadExecutable(url, destPath); err != nil {
		return err
	}

	utils.PrintSuccess("Pre-built base image downloaded to %s", utils.StylePath(destPath))
	return nil
}

// ensureBaseDef ensures the base image definition file exists, downloading it if necessary.
// Searches all build-scripts paths first, then downloads to the first writable directory.
// Returns the path to the definition file.
func ensureBaseDef() (string, error) {
	// Search all build-scripts paths for existing def file
	if existingPath := config.FindBaseImageDef(); existingPath != "" {
		utils.PrintDebug("Base image def found at %s", utils.StylePath(existingPath))
		return existingPath, nil
	}

	// Get the path where we'll write the new def file
	baseDefPath, err := config.GetBaseImageDefWritePath()
	if err != nil {
		return "", fmt.Errorf("no writable directory for base image def: %w", err)
	}

	utils.PrintMessage("Base image definition file not found. Downloading...")

	// Ensure build-scripts directory exists
	if err := os.MkdirAll(filepath.Dir(baseDefPath), utils.PermDir); err != nil {
		return "", fmt.Errorf("failed to create build-scripts directory: %w", err)
	}

	url := fmt.Sprintf("https://raw.githubusercontent.com/%s/main/build-scripts/base_image.def", config.GitHubRepo)

	if err := utils.DownloadFile(url, baseDefPath); err != nil {
		return "", fmt.Errorf("failed to download base image definition file: %w", err)
	}

	utils.PrintSuccess("Base image definition file downloaded to %s", utils.StylePath(baseDefPath))
	return baseDefPath, nil
}
