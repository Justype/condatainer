package apptainer

import (
	"fmt"
	"io"
	"net/http"
	"os"
	"path/filepath"
	"runtime"
	"time"

	"github.com/Justype/condatainer/internal/config"
	"github.com/Justype/condatainer/internal/utils"
)

// PrebuiltBaseImagePlatforms lists architectures with prebuilt base images available
var PrebuiltBaseImagePlatforms = map[string]bool{
	"x86_64": true,
	// "aarch64": false, // Not available yet
}

// PrebuiltBaseImageVersion is the version tag for prebuilt images
const PrebuiltBaseImageVersion = "v1.0.5"

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

	// Step 4: Try to download prebuilt base image (unless in debug mode)
	if !config.Global.Debug {
		utils.PrintMessage("Base image not found. Downloading base image...")
		if err := tryDownloadPrebuiltBaseImage(); err == nil {
			// Download succeeded
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
	arch := runtime.GOARCH
	// Map Go arch names to the naming convention used in releases
	if arch == "amd64" {
		arch = "x86_64"
	} else if arch == "arm64" {
		arch = "aarch64"
	}

	if !PrebuiltBaseImagePlatforms[arch] {
		return fmt.Errorf("pre-built base image not available for architecture: %s", arch)
	}

	url := fmt.Sprintf("https://github.com/%s/releases/download/%s/base_image_%s.sif",
		config.GitHubRepo, PrebuiltBaseImageVersion, arch)

	utils.PrintMessage("Attempting to download pre-built base image for %s...", utils.StyleInfo(arch))

	if err := downloadExecutable(url, config.Global.BaseImage); err != nil {
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

	if err := downloadFile(url, baseDefPath); err != nil {
		return fmt.Errorf("failed to download base image definition file: %w", err)
	}

	utils.PrintSuccess("Base image definition file downloaded.")
	return nil
}

// downloadFile downloads a file from url to destPath.
func downloadFile(url, destPath string) error {
	client := &http.Client{
		Timeout: 5 * time.Minute,
	}

	resp, err := client.Get(url)
	if err != nil {
		return fmt.Errorf("failed to fetch %s: %w", url, err)
	}
	defer resp.Body.Close()

	if resp.StatusCode != http.StatusOK {
		return fmt.Errorf("failed to download: HTTP %d", resp.StatusCode)
	}

	// Create temporary file first
	tmpPath := destPath + ".tmp"
	file, err := os.Create(tmpPath)
	if err != nil {
		return fmt.Errorf("failed to create file: %w", err)
	}

	_, err = io.Copy(file, resp.Body)
	file.Close()
	if err != nil {
		os.Remove(tmpPath)
		return fmt.Errorf("failed to write file: %w", err)
	}

	// Rename temp file to final destination
	if err := os.Rename(tmpPath, destPath); err != nil {
		os.Remove(tmpPath)
		return fmt.Errorf("failed to rename file: %w", err)
	}

	return nil
}

// downloadExecutable downloads a file and sets it as executable (0755).
func downloadExecutable(url, destPath string) error {
	if err := downloadFile(url, destPath); err != nil {
		return err
	}

	// Set executable permissions
	if err := os.Chmod(destPath, 0755); err != nil {
		return fmt.Errorf("failed to set executable permissions: %w", err)
	}

	return nil
}
