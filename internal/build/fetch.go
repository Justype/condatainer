package build

import (
	"compress/gzip"
	"encoding/json"
	"fmt"
	"io"
	"net/http"
	"os"
	"path/filepath"
	"strings"
	"time"

	"github.com/Justype/condatainer/internal/config"
	"github.com/Justype/condatainer/internal/utils"
)

// RemoteMetadataURL is the URL for the remote build scripts metadata
var RemoteMetadataURL = fmt.Sprintf("https://raw.githubusercontent.com/%s/main/metadata/build-scripts.json.gz", config.GitHubRepo)

// ScriptInfo holds information about a build script
type ScriptInfo struct {
	Name        string // name/version format (e.g., "samtools/1.21")
	Path        string // full path to build script (local) or relative path (remote)
	IsContainer bool   // true if .def file
	IsRemote    bool   // true if from remote repository
}

// RemoteScriptEntry represents a single script entry in the metadata
type RemoteScriptEntry struct {
	RelativePath string   `json:"relative_path"`
	Deps         []string `json:"deps"`
	Sbatch       bool     `json:"sbatch"`
	Whatis       string   `json:"whatis"`
}

// GetLocalBuildScripts scans all build-scripts directories and returns a map of name -> ScriptInfo.
// Searches user → scratch → legacy → system directories.
// First match wins (user scripts shadow system ones).
func GetLocalBuildScripts() (map[string]ScriptInfo, error) {
	scripts := make(map[string]ScriptInfo)

	// Scan all build script directories in priority order
	for _, buildScriptsDir := range config.GetBuildScriptSearchPaths() {
		if !utils.DirExists(buildScriptsDir) {
			continue
		}

		err := filepath.Walk(buildScriptsDir, func(path string, info os.FileInfo, err error) error {
			if err != nil {
				return nil // Skip errors, continue with other files
			}

			// Skip directories
			if info.IsDir() {
				return nil
			}

			// Skip .py and .sh files (helper scripts)
			if strings.HasSuffix(info.Name(), ".py") || strings.HasSuffix(info.Name(), ".sh") {
				return nil
			}

			// Skip template files
			if strings.Contains(path, "template") {
				return nil
			}

			// Generate key (relative path)
			relPath, err := filepath.Rel(buildScriptsDir, path)
			if err != nil {
				return nil
			}

			// Handle .def files
			isContainer := strings.HasSuffix(relPath, ".def")
			if isContainer {
				// Skip base image/overlay scripts
				if strings.HasPrefix(relPath, "base_image") || strings.HasPrefix(relPath, "base-overlay") {
					return nil
				}
				relPath = strings.TrimSuffix(relPath, ".def")
			}

			// First match wins - don't overwrite if already found in higher-priority dir
			if _, exists := scripts[relPath]; !exists {
				scripts[relPath] = ScriptInfo{
					Name:        relPath,
					Path:        path,
					IsContainer: isContainer,
					IsRemote:    false,
				}
			}
			return nil
		})

		if err != nil {
			// Log but continue with other directories
			utils.PrintDebug("Error scanning %s: %v", buildScriptsDir, err)
		}
	}

	return scripts, nil
}

// GetRemoteBuildScripts fetches the build scripts metadata from GitHub
// Returns a map of name -> ScriptInfo
func GetRemoteBuildScripts() (map[string]ScriptInfo, error) {
	scripts := make(map[string]ScriptInfo)

	// Create HTTP client with timeout
	client := &http.Client{
		Timeout: 30 * time.Second,
	}

	// Fetch the gzipped metadata
	resp, err := client.Get(RemoteMetadataURL)
	if err != nil {
		return nil, fmt.Errorf("failed to fetch remote metadata: %w", err)
	}
	defer resp.Body.Close()

	if resp.StatusCode != http.StatusOK {
		return nil, fmt.Errorf("failed to fetch remote metadata: HTTP %d", resp.StatusCode)
	}

	// Decompress gzip
	gzReader, err := gzip.NewReader(resp.Body)
	if err != nil {
		return nil, fmt.Errorf("failed to decompress metadata: %w", err)
	}
	defer gzReader.Close()

	// Read all data
	data, err := io.ReadAll(gzReader)
	if err != nil {
		return nil, fmt.Errorf("failed to read metadata: %w", err)
	}

	// Parse JSON - it's a flat map of name -> entry
	var metadata map[string]RemoteScriptEntry
	if err := json.Unmarshal(data, &metadata); err != nil {
		return nil, fmt.Errorf("failed to parse metadata: %w", err)
	}

	// Convert to ScriptInfo map
	for name, entry := range metadata {
		// Determine if it's a container (.def file)
		isContainer := strings.HasSuffix(entry.RelativePath, ".def")
		scripts[name] = ScriptInfo{
			Name:        name,
			Path:        entry.RelativePath,
			IsContainer: isContainer,
			IsRemote:    true,
		}
	}

	return scripts, nil
}

// GetAllBuildScripts returns both local and remote build scripts
// Local scripts take precedence over remote scripts with the same name
func GetAllBuildScripts(includeRemote bool) (map[string]ScriptInfo, error) {
	// Get local scripts first
	scripts, err := GetLocalBuildScripts()
	if err != nil {
		return nil, err
	}

	if !includeRemote {
		return scripts, nil
	}

	// Get remote scripts
	remoteScripts, err := GetRemoteBuildScripts()
	if err != nil {
		// Log warning but don't fail - remote is optional
		utils.PrintDebug("Failed to fetch remote scripts: %v", err)
		return scripts, nil
	}

	// Merge remote scripts (local takes precedence)
	for name, info := range remoteScripts {
		if _, exists := scripts[name]; !exists {
			scripts[name] = info
		}
	}

	return scripts, nil
}

// FindBuildScript looks for a build script by name/version
// First checks local, then remote if not found locally
// Returns the ScriptInfo and a boolean indicating if it was found
func FindBuildScript(nameVersion string) (ScriptInfo, bool) {
	normalized := utils.NormalizeNameVersion(nameVersion)

	// Check local first
	localScripts, err := GetLocalBuildScripts()
	if err == nil {
		if info, found := localScripts[normalized]; found {
			return info, true
		}
	}

	// Check remote
	remoteScripts, err := GetRemoteBuildScripts()
	if err == nil {
		if info, found := remoteScripts[normalized]; found {
			return info, true
		}
	}

	return ScriptInfo{}, false
}

// DownloadRemoteScript downloads a remote build script to the local tmp directory
// Returns the local path to the downloaded script
func DownloadRemoteScript(info ScriptInfo, tmpDir string) (string, error) {
	if !info.IsRemote {
		return info.Path, nil
	}

	// Build the raw GitHub URL
	rawURL := fmt.Sprintf("https://raw.githubusercontent.com/%s/main/build-scripts/%s", config.GitHubRepo, info.Path)
	if info.IsContainer {
		rawURL += ".def"
	}

	// Create HTTP client with timeout
	client := &http.Client{
		Timeout: 60 * time.Second,
	}

	// Fetch the script
	resp, err := client.Get(rawURL)
	if err != nil {
		return "", fmt.Errorf("failed to download script: %w", err)
	}
	defer resp.Body.Close()

	if resp.StatusCode != http.StatusOK {
		return "", fmt.Errorf("failed to download script: HTTP %d", resp.StatusCode)
	}

	// Create local path
	localPath := filepath.Join(tmpDir, "remote_scripts", info.Name)
	if info.IsContainer {
		localPath += ".def"
	}

	// Ensure directory exists
	if err := os.MkdirAll(filepath.Dir(localPath), 0755); err != nil {
		return "", fmt.Errorf("failed to create directory: %w", err)
	}

	// Write to file
	file, err := os.Create(localPath)
	if err != nil {
		return "", fmt.Errorf("failed to create file: %w", err)
	}
	defer file.Close()

	if _, err := io.Copy(file, resp.Body); err != nil {
		return "", fmt.Errorf("failed to write script: %w", err)
	}

	return localPath, nil
}
