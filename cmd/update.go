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

	"github.com/Justype/condatainer/internal/config"
	"github.com/Justype/condatainer/internal/utils"
	"github.com/spf13/cobra"
)

var (
	updateYes bool
)

var updateCmd = &cobra.Command{
	Use:     "self-update",
	Aliases: []string{"update"},
	Short:   "Update condatainer to the latest version from GitHub",
	Long: `Download and replace the current condatainer binary with the latest version from GitHub.

This command downloads the latest Go binary (condatainer_go) from the GitHub repository
and replaces the current executable. A backup of the current version is not created.`,
	Example: `  condatainer self-update       # Update with confirmation prompt
  condatainer self-update -y    # Update without confirmation`,
	RunE: runUpdate,
}

func init() {
	rootCmd.AddCommand(updateCmd)
	updateCmd.Flags().BoolVarP(&updateYes, "yes", "y", false, "Skip confirmation prompt")
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

	// Ask for confirmation unless --yes flag is provided
	if !updateYes {
		fmt.Print("Are you sure you want to download and replace the current binary from GitHub releases? [y/N]: ")
		var confirm string
		fmt.Scanln(&confirm)
		confirm = strings.ToLower(strings.TrimSpace(confirm))

		if confirm != "y" && confirm != "yes" {
			utils.PrintNote("Update cancelled by user.")
			return nil
		}
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

	utils.PrintMessage("Fetching latest release information...")

	// Get latest release from GitHub API
	releaseURL := fmt.Sprintf("https://api.github.com/repos/%s/releases/latest", config.GITHUB_REPO)
	resp, err := http.Get(releaseURL)
	if err != nil {
		return fmt.Errorf("failed to fetch release information: %w", err)
	}
	defer resp.Body.Close()

	if resp.StatusCode != http.StatusOK {
		return fmt.Errorf("failed to fetch release: HTTP %d", resp.StatusCode)
	}

	var release struct {
		TagName string `json:"tag_name"`
		Assets  []struct {
			Name               string `json:"name"`
			BrowserDownloadURL string `json:"browser_download_url"`
		} `json:"assets"`
	}

	if err := json.NewDecoder(resp.Body).Decode(&release); err != nil {
		return fmt.Errorf("failed to parse release information: %w", err)
	}

	// Find matching binary for current OS/arch
	// Expected format: condatainer_go_{os}_{arch} (e.g., condatainer_go_linux_x86_64)
	binaryName := fmt.Sprintf("condatainer_go_%s_%s", osName, arch)
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
	if err := os.Chmod(tempPath, 0755); err != nil {
		os.Remove(tempPath)
		return fmt.Errorf("failed to set executable permissions: %w", err)
	}

	// Replace current executable
	// On Unix systems, we can replace the file while it's running
	if err := os.Rename(tempPath, exePath); err != nil {
		os.Remove(tempPath)
		return fmt.Errorf("failed to replace executable: %w", err)
	}

	utils.PrintSuccess("condatainer updated successfully to %s!", utils.StyleNumber(release.TagName))
	utils.PrintNote("The new version will be used on the next run.")

	return nil
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
	if err := os.MkdirAll(filepath.Dir(destPath), 0755); err != nil {
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
