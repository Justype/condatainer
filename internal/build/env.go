package build

import (
	"fmt"
	"strings"

	"log/slog"

	"github.com/Justype/condatainer/internal/runtime/container"
	"github.com/Justype/condatainer/internal/utils"
)

// env.go writes the sidecar .env file and collects overlay mount arguments from
// dependencies. Recipe #ENV: declarations are read from the embedded build
// script at load time — see internal/runtime/container/env.go.

// EnvEntry holds an environment variable value and its note
type EnvEntry struct {
	Value string
	Note  string
}

// SaveEnvFile saves environment variables to a .env file next to the image.
// The {prefix} placeholder is replaced with the overlay's mount root.
// description is written as a #DESCRIPTION: comment at the top if non-empty.
func SaveEnvFile(overlayPath string, envDict map[string]EnvEntry, relativePath string, description string) error {
	if len(envDict) == 0 && description == "" {
		return nil
	}

	envFilePath := overlayPath + ".env"
	file, err := utils.CreateFileWritable(envFilePath)
	if err != nil {
		return fmt.Errorf("failed to create env file: %w", err)
	}
	defer file.Close()

	if description != "" {
		if _, err := fmt.Fprintf(file, "#DESCRIPTION:%s\n", description); err != nil {
			return fmt.Errorf("failed to write description: %w", err)
		}
	}

	for key, entry := range envDict {
		// Replace {prefix} placeholder with actual path
		value := strings.ReplaceAll(entry.Value, "{prefix}", fmt.Sprintf("/cnt/%s", relativePath))

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

	utils.ShareWithParentGroup(envFilePath)

	slog.Default().Info("ENV file created", "path", envFilePath)
	return nil
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
