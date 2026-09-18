package cmd

import (
	"encoding/json"
	"fmt"
	"io"
	"net/http"
	"os"
	"path/filepath"
	"runtime"
	"strings"

	"golang.org/x/mod/semver"

	"github.com/Justype/condatainer/internal/config"
	"github.com/Justype/condatainer/internal/helper"
	"github.com/Justype/condatainer/internal/libexec"
	"github.com/Justype/condatainer/internal/utils"
	"github.com/spf13/cobra"
)

// ── update command ──────────────────────────────────────────────────────────

var (
	updateBuild       bool
	updateHelpScripts bool
	updateLibexec     bool
)

var updateCmd = &cobra.Command{
	Use:   "update",
	Short: "Update script metadata caches or the toolchain",
	Long: `Update build script metadata, helper script metadata, or the self-provisioned toolchain.

With no flags, refreshes both the build and helper script metadata caches.

With --libexec, updates the tools installed in the toolchain, or creates it with
micromamba when there is none. Name packages after it to install them too:
apptainer, squashfs-tools, squashfuse.`,
	Example: `  condatainer update                 # Refresh build + helper metadata (default)
  condatainer update --build         # Build script metadata only
  condatainer update --helper        # Helper script metadata only
  condatainer update --libexec       # Update the installed toolchain
  condatainer update --libexec apptainer squashfs-tools  # Install and update these`,
	Args: func(cmd *cobra.Command, args []string) error {
		if len(args) > 0 && !updateLibexec {
			return fmt.Errorf("packages can only be named with --libexec")
		}
		return nil
	},
	SilenceUsage: true,
	RunE:         runUpdate,
}

func init() {
	rootCmd.AddCommand(updateCmd)
	updateCmd.Flags().BoolVar(&updateBuild, "build", false, "Refresh build script metadata cache")
	updateCmd.Flags().BoolVar(&updateHelpScripts, "helper", false, "Refresh helper script metadata cache")
	updateCmd.Flags().BoolVar(&updateLibexec, "libexec", false, "Update the self-provisioned toolchain; name packages to install them")
}

func runUpdate(cmd *cobra.Command, args []string) error {
	// Default: both --build and --help-scripts when no content flags given
	if !updateBuild && !updateHelpScripts && !updateLibexec {
		updateBuild = true
		updateHelpScripts = true
	}

	if updateBuild {
		if err := config.RefreshCatalogCache(); err != nil {
			return fmt.Errorf("failed to clear the recipe cache: %w", err)
		}
		cat, err := config.OpenCatalog(cmd.Context())
		if err != nil {
			return fmt.Errorf("failed to open recipe sources: %w", err)
		}
		for _, src := range cat {
			utils.PrintMessage("Fetching recipes from %s (%s) ...", src.Name, src.Base)
		}
		entries := cat.Entries(cmd.Context())
		config.WarnUnreachableSources(cmd.Context(), cat)
		utils.PrintSuccess("%d recipes available.", len(entries))
	}

	if updateHelpScripts {
		if _, err := helper.RefreshRemoteMetadata(cmd.Context(), true, os.Stdout); err != nil {
			utils.PrintWarning("Failed to update helper script metadata: %v", err)
		} else {
			utils.PrintSuccess("Helper script metadata updated.")
		}
	}

	if updateLibexec {
		// libexec.Update holds its own lock and refuses internally if the
		// toolchain is in use; no separate check needed here.
		utils.PrintMessage("Updating the self-provisioned toolchain...")
		if err := libexec.Update(cmd.Context(), args...); err != nil {
			return fmt.Errorf("failed to update the toolchain: %w", err)
		}
		utils.PrintSuccess("Toolchain updated successfully.")
		if versions, err := libexec.Versions(cmd.Context()); err == nil {
			for _, v := range versions {
				utils.PrintMessage("  %-10s %s", v.Name, v.Version)
			}
		}
	}

	return nil
}

// ── self-update command ──────────────────────────────────────────────────────

var (
	selfUpdateForce bool
	selfUpdateDev   bool
)

var selfUpdateCmd = &cobra.Command{
	Use:   "self-update",
	Short: "Update condatainer to the latest version",
	Long: `Download the latest condatainer from GitHub and replace the current binary.

Note: No backup of the current version is kept.`,
	Example: `  condatainer self-update        # Update to latest stable version
  condatainer self-update --yes  # Update without confirmation
  condatainer self-update -f     # Force update even if already on latest version
  condatainer self-update --dev  # Include pre-release versions`,
	Args:         cobra.NoArgs,
	SilenceUsage: true,
	RunE:         runSelfUpdate,
}

func init() {
	rootCmd.AddCommand(selfUpdateCmd)
	selfUpdateCmd.Flags().BoolVarP(&selfUpdateForce, "force", "f", false, "Force update even if already on latest version")
	selfUpdateCmd.Flags().BoolVar(&selfUpdateDev, "dev", false, "Include pre-release versions")
}

func runSelfUpdate(cmd *cobra.Command, args []string) error {
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

	if selfUpdateDev {
		utils.PrintNote("Dev mode enabled, including pre-release versions")
	}

	utils.PrintMessage("Fetching latest release information...")

	// Get release from GitHub API
	var releaseURL string
	if selfUpdateDev {
		// Fetch all releases to find the latest (including pre-releases)
		releaseURL = fmt.Sprintf("https://api.github.com/repos/%s/releases", config.GitHubRepo)
	} else {
		// Fetch only the latest stable release
		releaseURL = fmt.Sprintf("https://api.github.com/repos/%s/releases/latest", config.GitHubRepo)
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

	if selfUpdateDev {
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

	if cmp >= 0 && !selfUpdateForce {
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
	if !utils.ShouldAnswerYes() && !selfUpdateForce {
		fmt.Print("Are you sure to update to the latest version? [y/N]: ")
		confirm, err := utils.ReadLineContext(cmd.Context())
		if err != nil || (confirm != "y" && confirm != "yes") {
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
	if err := utils.MakeExecutable(tempPath); err != nil {
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

	return nil
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
	if err := utils.MkdirAllShared(filepath.Dir(destPath)); err != nil {
		return err
	}

	// Create destination file
	out, err := utils.CreateFileWritable(destPath)
	if err != nil {
		return err
	}
	defer out.Close()

	// Copy data. Permissions were already set by CreateFileWritable (umask-subject,
	// shared with the parent group); io.Copy only writes content.
	_, err = io.Copy(out, resp.Body)
	return err
}
