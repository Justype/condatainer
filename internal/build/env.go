package build

import (
	"bufio"
	"fmt"
	"os"
	"path/filepath"
	"strings"

	"github.com/Justype/condatainer/internal/container"
	"github.com/Justype/condatainer/internal/utils"
)

// env.go provides functionality for reading and managing environment variables
// from build scripts and dependencies.
//
// Key features:
// - Parse #ENV: and #ENVNOTE: directives from build scripts
// - Generate PATH environment from dependencies
// - Collect overlay mount arguments from dependencies
// - Load .env files from dependency overlays
// - Save environment variables to .env files for overlays
//
// This mirrors the Python version's behavior for dependency and ENV handling.

// GetEnvDictFromBuildScript parses #ENV:KEY=VALUE and #ENVNOTE:KEY description lines from build script.
// Returns a map of environment variables with their values and notes.
//
// Example lines in script:
//
//	#ENV:STAR_INDEX_DIR=$app_root/star
//	#ENVNOTE:STAR_INDEX_DIR STAR index dir
func GetEnvDictFromBuildScript(scriptPath string) (map[string]EnvEntry, error) {
	envDict := make(map[string]EnvEntry)

	if !utils.FileExists(scriptPath) {
		return envDict, fmt.Errorf("build script not found at %s", scriptPath)
	}

	file, err := os.Open(scriptPath)
	if err != nil {
		return envDict, fmt.Errorf("failed to open build script: %w", err)
	}
	defer file.Close()

	// Read all lines into memory for lookahead
	var lines []string
	scanner := bufio.NewScanner(file)
	for scanner.Scan() {
		lines = append(lines, scanner.Text())
	}
	if err := scanner.Err(); err != nil {
		return envDict, fmt.Errorf("failed to read build script: %w", err)
	}

	// Find all #ENV: lines
	for i, line := range lines {
		line = strings.TrimSpace(line)
		if !strings.HasPrefix(line, "#ENV:") {
			continue
		}

		// Parse KEY=VALUE (strip inline comments)
		content := utils.StripInlineComment(line[len("#ENV:"):])
		if !strings.Contains(content, "=") {
			continue
		}

		parts := strings.SplitN(content, "=", 2)
		key := strings.TrimSpace(parts[0])
		value := strings.TrimSpace(parts[1])

		entry := EnvEntry{
			Value: value,
			Note:  "",
		}

		// Check if the next line is an #ENVNOTE:
		if i+1 < len(lines) {
			nextLine := strings.TrimSpace(lines[i+1])
			if strings.HasPrefix(nextLine, "#ENVNOTE:") {
				note := strings.TrimSpace(nextLine[len("#ENVNOTE:"):])
				entry.Note = note
			}
		}

		envDict[key] = entry
	}

	return envDict, nil
}

// EnvEntry holds an environment variable value and its note
type EnvEntry struct {
	Value string
	Note  string
}

// SaveEnvFile saves environment variables to a .env file next to the overlay.
// The $app_root placeholder is replaced with the actual overlay path.
func SaveEnvFile(overlayPath string, envDict map[string]EnvEntry, relativePath string) error {
	if len(envDict) == 0 {
		return nil
	}

	envFilePath := overlayPath + ".env"
	file, err := os.Create(envFilePath)
	if err != nil {
		return fmt.Errorf("failed to create env file: %w", err)
	}
	defer file.Close()

	for key, entry := range envDict {
		// Replace $app_root placeholder with actual path
		value := strings.ReplaceAll(entry.Value, "$app_root", fmt.Sprintf("/cnt/%s", relativePath))

		// Write KEY=VALUE
		if _, err := fmt.Fprintf(file, "%s=%s\n", key, value); err != nil {
			return fmt.Errorf("failed to write env entry: %w", err)
		}

		// Write #ENVNOTE:KEY=Note if present
		if entry.Note != "" {
			if _, err := fmt.Fprintf(file, "#ENVNOTE:%s=%s\n", key, entry.Note); err != nil {
				return fmt.Errorf("failed to write env note: %w", err)
			}
		}
	}

	if err := os.Chmod(envFilePath, 0o664); err != nil {
		utils.PrintDebug("Failed to set permissions on %s: %v", envFilePath, err)
	}

	utils.PrintMessage("ENV file created at %s", utils.StylePath(envFilePath))
	return nil
}

// GetPathEnvFromDependencies generates the PATH environment variable from dependencies.
// It prepends /cnt/<dep>/bin paths for each dependency overlay.
func GetPathEnvFromDependencies(dependencies []string) (string, error) {
	if len(dependencies) == 0 {
		return "", nil
	}

	// Resolve dependency paths
	depPaths, err := container.ResolveOverlayPaths(dependencies)
	if err != nil {
		return "", fmt.Errorf("failed to resolve dependency paths: %w", err)
	}

	// Build PATH components
	paths := []string{}
	for _, depPath := range depPaths {
		// Strip :ro or :rw suffix for path checking
		cleanPath := strings.TrimSuffix(strings.TrimSuffix(depPath, ":ro"), ":rw")
		name := strings.TrimSuffix(filepath.Base(cleanPath), filepath.Ext(cleanPath))
		normalized := utils.NormalizeNameVersion(name)
		if normalized == "" {
			continue
		}
		// Skip ref overlays (more than one slash)
		if strings.Count(normalized, "/") > 1 {
			continue
		}

		var binPath string
		if utils.IsImg(cleanPath) {
			binPath = "/ext3/env/bin"
		} else if utils.IsSqf(cleanPath) {
			binPath = fmt.Sprintf("/cnt/%s/bin", normalized)
		} else {
			utils.PrintWarning("Unknown overlay file extension for %s. Skipping PATH addition.", utils.StylePath(depPath))
			continue
		}
		paths = append(paths, binPath)
	}

	if len(paths) == 0 {
		return "", nil
	}

	// Return as PATH-style string
	return strings.Join(paths, ":"), nil
}

// GetOverlayArgsFromDependencies generates overlay mount arguments for Apptainer from dependencies.
// Each dependency is mounted as read-only.
func GetOverlayArgsFromDependencies(dependencies []string) ([]string, error) {
	if len(dependencies) == 0 {
		return nil, nil
	}

	// Resolve dependency paths
	depPaths, err := container.ResolveOverlayPaths(dependencies)
	if err != nil {
		return nil, fmt.Errorf("failed to resolve dependency paths: %w", err)
	}

	// Build overlay arguments
	args := []string{}
	for _, depPath := range depPaths {
		args = append(args, "--overlay", depPath+":ro")
	}

	return args, nil
}

// GetOverlayEnvConfigsFromDependencies reads .env files from dependencies and returns environment settings.
// Returns a list of ["--env", "KEY=VALUE", ...] arguments.
func GetOverlayEnvConfigsFromDependencies(dependencies []string) ([]string, error) {
	if len(dependencies) == 0 {
		return nil, nil
	}

	// Resolve dependency paths
	depPaths, err := container.ResolveOverlayPaths(dependencies)
	if err != nil {
		return nil, fmt.Errorf("failed to resolve dependency paths: %w", err)
	}

	configs := make(map[string]string)

	// Read .env files from each dependency
	for _, depPath := range depPaths {
		// Strip :ro or :rw suffix for file path
		cleanPath := strings.TrimSuffix(strings.TrimSuffix(depPath, ":ro"), ":rw")
		envPath := cleanPath + ".env"
		if !utils.FileExists(envPath) {
			continue
		}

		file, err := os.Open(envPath)
		if err != nil {
			utils.PrintWarning("Unable to read overlay env %s: %v", utils.StylePath(envPath), err)
			continue
		}

		scanner := bufio.NewScanner(file)
		for scanner.Scan() {
			line := strings.TrimSpace(scanner.Text())
			if line == "" || strings.HasPrefix(line, "#") {
				continue
			}

			parts := strings.SplitN(line, "=", 2)
			if len(parts) != 2 {
				continue
			}

			key := strings.TrimSpace(parts[0])
			value := strings.TrimSpace(parts[1])
			if key == "" {
				continue
			}

			if _, exists := configs[key]; exists {
				utils.PrintMessage("Environment variable %s is defined in multiple overlays. Using value from %s.",
					utils.StyleName(key), utils.StylePath(filepath.Base(depPath)))
			}
			configs[key] = value
		}

		if err := scanner.Err(); err != nil {
			utils.PrintWarning("Failed to scan %s: %v", utils.StylePath(envPath), err)
		}
		file.Close()
	}

	// Convert to KEY=VALUE format (not --env prefixed, exec.Options.EnvSettings expects raw KEY=VALUE)
	args := []string{}
	for key, value := range configs {
		args = append(args, fmt.Sprintf("%s=%s", key, value))
	}

	return args, nil
}
