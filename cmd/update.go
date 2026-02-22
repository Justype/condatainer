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
	"strconv"
	"strings"

	"golang.org/x/mod/semver"

	"github.com/Justype/condatainer/internal/apptainer"
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

	// Compare versions semantically; if parsing fails we assume the
	// latest version is newer (so we update).
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
	if getMajorMinor(currentVersion) != getMajorMinor(latestVersion) {
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

	return nil
}

// getMajorNumber returns the major version as an integer. The input may
// include or omit a leading 'v'. If parsing fails the error is returned.
func getMajorNumber(version string) (int, error) {
	if !strings.HasPrefix(version, "v") {
		version = "v" + version
	}
	c := semver.Canonical(version)
	if c == "" {
		return 0, fmt.Errorf("invalid version %q", version)
	}
	maj := strings.TrimPrefix(semver.Major(c), "v")
	return strconv.Atoi(maj)
}

// getMajorMinor returns the "vMAJOR.MINOR" string for a version, ignoring
// patch and any suffixes. Returns empty string on failure.
func getMajorMinor(version string) string {
	if !strings.HasPrefix(version, "v") {
		version = "v" + version
	}
	c := semver.Canonical(version)
	if c == "" {
		return ""
	}
	return semver.MajorMinor(c)
}

func isMinorChange(oldVersion, newVersion string) bool {
	return getMajorMinor(oldVersion) != getMajorMinor(newVersion)
}

// isMajorChange reports whether the major version component has changed.
// It returns false if either version cannot be parsed.
func isMajorChange(oldVersion, newVersion string) bool {
	m1, err1 := getMajorNumber(oldVersion)
	m2, err2 := getMajorNumber(newVersion)
	return err1 == nil && err2 == nil && m1 != m2
}

// compareVersions compares two semantic versions. It returns:
//
//	-1 if v1 < v2, 0 if v1 == v2, 1 if v1 > v2.
//
// Pre-release data is taken into account according to semver rules
// (e.g. "1.2.3-alpha" < "1.2.3"). Build metadata is used only as a
// secondary lexicographic tie‑breaker.
func compareVersions(v1, v2 string) int {
	// semver package requires a leading 'v'; add it if missing so that
	// canonicalization succeeds for numeric-only tags.
	if !strings.HasPrefix(v1, "v") {
		v1 = "v" + v1
	}
	if !strings.HasPrefix(v2, "v") {
		v2 = "v" + v2
	}
	c1 := semver.Canonical(v1)
	c2 := semver.Canonical(v2)
	if c1 == "" || c2 == "" {
		// If we can't parse a version, assume the first is older so an update
		// will be attempted.
		return -1
	}
	res := semver.Compare(c1, c2)
	if res != 0 {
		return res
	}
	b1 := semver.Build(v1)
	b2 := semver.Build(v2)
	if b1 != b2 {
		if b1 < b2 {
			return -1
		}
		return 1
	}
	return 0
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
