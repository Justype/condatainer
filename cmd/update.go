package cmd

import (
	"context"
	"encoding/json"
	"fmt"
	"io"
	"net/http"
	"os"
	"path/filepath"
	"runtime"
	"sort"
	"strconv"
	"strings"

	"github.com/Justype/condatainer/internal/apptainer"
	"github.com/Justype/condatainer/internal/build"
	"github.com/Justype/condatainer/internal/config"
	"github.com/Justype/condatainer/internal/utils"
	"github.com/spf13/cobra"
)

var (
	updateForce bool
	updateDev   bool
)

var updateCmd = &cobra.Command{
	Use:     "self-update",
	Aliases: []string{"update"},
	Short:   "Update condatainer to the latest version from GitHub",
	Long: `Download and replace the current condatainer binary with the latest version from GitHub.

This command downloads the latest binary from the GitHub repository
and replaces the current executable. A backup of the current version is not created.`,
	Example: `  condatainer self-update       # Update to latest stable version
  condatainer self-update --yes # Update without confirmation
  condatainer self-update -f    # Force update even if already on latest version
  condatainer self-update --dev # Include pre-release versions`,
	SilenceUsage: true, // Runtime errors should not show usage
	RunE:         runUpdate,
}

func init() {
	rootCmd.AddCommand(updateCmd)
	updateCmd.Flags().BoolVarP(&updateForce, "force", "f", false, "Force update even if already on latest version")
	updateCmd.Flags().BoolVar(&updateDev, "dev", false, "Include pre-release versions (also enabled if config branch is 'dev')")
}

func runUpdate(cmd *cobra.Command, args []string) error {
	// Capture old version before updating
	oldVersion := config.VERSION

	// Get current executable path
	exePath, err := os.Executable()
	if err != nil {
		return fmt.Errorf("failed to get executable path: %w", err)
	}

	// Resolve symlinks
	exePath, err = filepath.EvalSymlinks(exePath)
	if err != nil {
		return fmt.Errorf("failed to resolve symlink: %w", err)
	}

	// Detect OS and architecture
	osName := runtime.GOOS
	arch := runtime.GOARCH

	// Map Go arch names to common names
	archMap := map[string]string{
		"amd64": "x86_64",
		"arm64": "aarch64",
		"386":   "i386",
	}
	if mappedArch, ok := archMap[arch]; ok {
		arch = mappedArch
	}

	// Check if dev mode is enabled (via flag or config branch)
	devMode := updateDev || config.Global.Branch == "dev"

	if devMode {
		utils.PrintNote("Dev mode enabled, including pre-release versions")
	}

	utils.PrintMessage("Fetching latest release information...")

	// Get release from GitHub API
	var releaseURL string
	if devMode {
		// Fetch all releases to find the latest (including pre-releases)
		releaseURL = fmt.Sprintf("https://api.github.com/repos/%s/releases", config.GITHUB_REPO)
	} else {
		// Fetch only the latest stable release
		releaseURL = fmt.Sprintf("https://api.github.com/repos/%s/releases/latest", config.GITHUB_REPO)
	}

	resp, err := http.Get(releaseURL)
	if err != nil {
		return fmt.Errorf("failed to fetch release information: %w", err)
	}
	defer resp.Body.Close()

	if resp.StatusCode != http.StatusOK {
		return fmt.Errorf("failed to fetch release: HTTP %d", resp.StatusCode)
	}

	type releaseAsset struct {
		Name               string `json:"name"`
		BrowserDownloadURL string `json:"browser_download_url"`
	}

	type releaseInfo struct {
		TagName    string         `json:"tag_name"`
		Prerelease bool           `json:"prerelease"`
		Assets     []releaseAsset `json:"assets"`
	}

	var release releaseInfo

	if devMode {
		// Parse array of releases and find the latest
		var releases []releaseInfo
		if err := json.NewDecoder(resp.Body).Decode(&releases); err != nil {
			return fmt.Errorf("failed to parse release information: %w", err)
		}

		if len(releases) == 0 {
			return fmt.Errorf("no releases found")
		}

		// The releases are already sorted by creation date (newest first)
		// Take the first one (latest release, stable or pre-release)
		release = releases[0]
	} else {
		// Parse single latest stable release
		if err := json.NewDecoder(resp.Body).Decode(&release); err != nil {
			return fmt.Errorf("failed to parse release information: %w", err)
		}
	}

	// Check if already on latest version
	currentVersion := "v" + config.VERSION
	latestVersion := strings.TrimSpace(release.TagName)

	// Compare versions semantically (not just string comparison)
	cmp := compareVersions(currentVersion, latestVersion)

	if cmp >= 0 && !updateForce {
		// Current version >= latest version (equal or newer)
		if cmp == 0 {
			utils.PrintSuccess("Already on the latest version %s!", utils.StyleNumber(latestVersion))
		} else {
			utils.PrintSuccess("Already on a newer version %s (latest: %s)",
				utils.StyleNumber(currentVersion), utils.StyleNumber(latestVersion))
		}
		return nil
	}

	utils.PrintMessage("Current version: %s", utils.StyleNumber(currentVersion))
	if release.Prerelease {
		utils.PrintMessage("Latest pre-release: %s", utils.StyleNumber(latestVersion))
	} else {
		utils.PrintMessage("Latest version: %s", utils.StyleNumber(latestVersion))
	}

	// Ask for confirmation unless global --yes flag or --force flag is provided
	if !utils.ShouldAnswerYes() && !updateForce {
		fmt.Print("Are you sure you want to download and replace the current binary from GitHub releases? [y/N]: ")
		var confirm string
		fmt.Scanln(&confirm)
		confirm = strings.ToLower(strings.TrimSpace(confirm))

		if confirm != "y" && confirm != "yes" {
			utils.PrintNote("Update cancelled by user.")
			return nil
		}
	}

	// Find matching binary for current OS/arch
	// Expected format: condatainer_{os}_{arch} (e.g., condatainer_linux_x86_64)
	binaryName := fmt.Sprintf("condatainer_%s_%s", osName, arch)
	var downloadURL string

	for _, asset := range release.Assets {
		if asset.Name == binaryName {
			downloadURL = asset.BrowserDownloadURL
			break
		}
	}

	if downloadURL == "" {
		return fmt.Errorf("no binary found for %s/%s in release %s", osName, arch, release.TagName)
	}

	utils.PrintMessage("Downloading condatainer %s for %s/%s...", utils.StyleNumber(release.TagName), osName, arch)

	// Download to temporary file
	tempPath := exePath + ".tmp"
	if err := downloadFile(downloadURL, tempPath); err != nil {
		return fmt.Errorf("failed to download latest version: %w", err)
	}

	// Make executable
	if err := os.Chmod(tempPath, utils.PermExec); err != nil {
		os.Remove(tempPath)
		return fmt.Errorf("failed to set executable permissions: %w", err)
	}

	// Replace current executable
	// On Unix systems, we can replace the file while it's running
	if err := os.Rename(tempPath, exePath); err != nil {
		os.Remove(tempPath)
		return fmt.Errorf("failed to replace executable: %w", err)
	}

	utils.PrintSuccess("condatainer updated to %s!", utils.StyleNumber(release.TagName))

	// Update base image only if minor or major version changed
	if shouldUpdateBaseImage(oldVersion, release.TagName) {
		fmt.Println()
		utils.PrintMessage("Downloading base image for version %s...", release.TagName)

		// Download new base image (update=true, downloadOnly=true)
		// This downloads to .new and replaces old if successful
		if err := apptainer.EnsureBaseImage(context.Background(), true, true); err != nil {
			utils.PrintWarning("Failed to update base image: %v", err)
			utils.PrintWarning("You may need to manually rebuild the base image.")
			// Don't fail - continue with version check
		} else {
			utils.PrintSuccess("Base image updated successfully")
		}
	} else {
		utils.PrintDebug("Base image update skipped (patch version change only)")
	}

	// Check for major version upgrade
	oldMajor, err1 := parseMajorVersion(oldVersion)
	newMajor, err2 := parseMajorVersion(release.TagName)

	if err1 == nil && err2 == nil && oldMajor != newMajor {
		fmt.Println()
		utils.PrintWarning("Major version upgrade detected: %d.x → %d.x", oldMajor, newMajor)
		fmt.Println()

		// Get def-built containers
		containers := getDefBuiltContainers()

		if len(containers) > 0 {
			fmt.Println("The following def-built containers may need to be rebuilt:")
			for _, name := range containers {
				fmt.Printf("  - %s\n", utils.StyleName(name))
			}
		}
	}

	return nil
}

// parseVersion extracts major, minor, patch from "v1.2.3" or "1.2.3"
// Handles semantic versioning: strips pre-release tags and build metadata
// Examples: "v1.2.3-alpha+build.123" → (1, 2, 3)
// Returns (major, minor, patch, error)
func parseVersion(version string) (int, int, int, error) {
	// Strip "v" prefix
	version = strings.TrimPrefix(version, "v")

	// Strip build metadata (everything after "+")
	// e.g., "1.2.3+build.123" → "1.2.3"
	if idx := strings.Index(version, "+"); idx != -1 {
		version = version[:idx]
	}

	// Strip pre-release tag (everything after "-")
	// e.g., "1.2.3-alpha" → "1.2.3"
	if idx := strings.Index(version, "-"); idx != -1 {
		version = version[:idx]
	}

	// Split by "."
	parts := strings.Split(version, ".")
	if len(parts) < 2 {
		return 0, 0, 0, fmt.Errorf("invalid version format")
	}

	major, err := strconv.Atoi(parts[0])
	if err != nil {
		return 0, 0, 0, err
	}

	minor, err := strconv.Atoi(parts[1])
	if err != nil {
		return 0, 0, 0, err
	}

	patch := 0
	if len(parts) >= 3 {
		patch, err = strconv.Atoi(parts[2])
		if err != nil {
			return 0, 0, 0, err
		}
	}

	return major, minor, patch, nil
}

// parseMajorVersion extracts major version from "v1.2.3" or "1.2.3"
func parseMajorVersion(version string) (int, error) {
	major, _, _, err := parseVersion(version)
	return major, err
}

// compareVersions compares two semantic versions
// Returns: -1 if v1 < v2, 0 if v1 == v2, 1 if v1 > v2
func compareVersions(v1, v2 string) int {
	major1, minor1, patch1, err1 := parseVersion(v1)
	major2, minor2, patch2, err2 := parseVersion(v2)

	// If parsing fails, treat as different versions
	if err1 != nil || err2 != nil {
		return -1 // Assume v1 < v2 to trigger update
	}

	// Compare major
	if major1 < major2 {
		return -1
	}
	if major1 > major2 {
		return 1
	}

	// Major equal, compare minor
	if minor1 < minor2 {
		return -1
	}
	if minor1 > minor2 {
		return 1
	}

	// Major and minor equal, compare patch
	if patch1 < patch2 {
		return -1
	}
	if patch1 > patch2 {
		return 1
	}

	// All equal
	return 0
}

// shouldUpdateBaseImage checks if base image should be updated based on version change
// Returns true if major or minor version changed (not just patch)
func shouldUpdateBaseImage(oldVersion, newVersion string) bool {
	oldMajor, oldMinor, _, err1 := parseVersion(oldVersion)
	newMajor, newMinor, _, err2 := parseVersion(newVersion)

	if err1 != nil || err2 != nil {
		// If we can't parse, assume update is needed
		return true
	}

	// Update if major or minor changed
	return oldMajor != newMajor || oldMinor != newMinor
}

// getDefBuiltContainers returns list of def-built containers (no "/" in name)
func getDefBuiltContainers() []string {
	defList := build.GetDefBuiltList()
	containers := []string{}
	for name := range defList {
		if !strings.Contains(name, "/") {
			containers = append(containers, name)
		}
	}
	sort.Strings(containers)
	return containers
}

// downloadFile downloads a file from a URL to a local path
func downloadFile(url, destPath string) error {
	resp, err := http.Get(url)
	if err != nil {
		return err
	}
	defer resp.Body.Close()

	if resp.StatusCode != http.StatusOK {
		return fmt.Errorf("HTTP %d: %s", resp.StatusCode, resp.Status)
	}

	// Create parent directory
	if err := os.MkdirAll(filepath.Dir(destPath), utils.PermDir); err != nil {
		return err
	}

	// Create destination file
	out, err := os.Create(destPath)
	if err != nil {
		return err
	}
	defer out.Close()

	// Copy data
	_, err = io.Copy(out, resp.Body)
	return err
}
