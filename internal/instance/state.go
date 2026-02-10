package instance

import (
	"encoding/json"
	"fmt"
	"os"
	"path/filepath"

	"github.com/Justype/condatainer/internal/config"
	"github.com/Justype/condatainer/internal/utils"
)

// State holds the runtime state of an instance
// apptainer exec instance://NAME cannot use the ENV used with instance start, so we have to save them here to be loaded by exec
type State struct {
	Name      string            `json:"name"`
	Env       []string          `json:"env"`        // Environment variables
	Overlays  []string          `json:"overlays"`   // Overlay paths
	BindPaths []string          `json:"bind_paths"` // Bind paths
	BaseImage string            `json:"base_image"` // Base image path
	Notes     map[string]string `json:"notes"`      // Environment variable notes for display
}

// getStateDir returns the directory where instance state files are stored
func getStateDir() (string, error) {
	userStateDir := config.GetUserStateDir()
	if userStateDir == "" {
		return "", fmt.Errorf("unable to determine user state directory")
	}
	stateDir := filepath.Join(userStateDir, "instance")
	return stateDir, nil
}

// getStatePath returns the path to the state file for a named instance
func getStatePath(name string) (string, error) {
	stateDir, err := getStateDir()
	if err != nil {
		return "", err
	}
	return filepath.Join(stateDir, name+".json"), nil
}

// saveState saves the instance state to disk
func saveState(state *State) error {
	stateDir, err := getStateDir()
	if err != nil {
		return err
	}

	// Create state directory if it doesn't exist
	if err := os.MkdirAll(stateDir, utils.PermDir); err != nil {
		return fmt.Errorf("failed to create state directory: %w", err)
	}

	statePath, err := getStatePath(state.Name)
	if err != nil {
		return err
	}

	data, err := json.MarshalIndent(state, "", "  ")
	if err != nil {
		return fmt.Errorf("failed to marshal state: %w", err)
	}

	if err := os.WriteFile(statePath, data, utils.PermFile); err != nil {
		return fmt.Errorf("failed to write state file: %w", err)
	}

	utils.PrintDebug("[INSTANCE]Saved state to %s", statePath)
	return nil
}

// loadState loads the instance state from disk
func loadState(name string) (*State, error) {
	statePath, err := getStatePath(name)
	if err != nil {
		return nil, err
	}

	data, err := os.ReadFile(statePath)
	if err != nil {
		if os.IsNotExist(err) {
			return nil, fmt.Errorf("no state found for instance %s", name)
		}
		return nil, fmt.Errorf("failed to read state file: %w", err)
	}

	var state State
	if err := json.Unmarshal(data, &state); err != nil {
		return nil, fmt.Errorf("failed to unmarshal state: %w", err)
	}

	utils.PrintDebug("[INSTANCE]Loaded state from %s", statePath)
	return &state, nil
}

// deleteState removes the instance state file
func deleteState(name string) error {
	statePath, err := getStatePath(name)
	if err != nil {
		return err
	}

	if err := os.Remove(statePath); err != nil && !os.IsNotExist(err) {
		return fmt.Errorf("failed to remove state file: %w", err)
	}

	utils.PrintDebug("[INSTANCE]Deleted state file %s", statePath)
	return nil
}

// deleteAllStates removes all instance state files
func deleteAllStates() error {
	stateDir, err := getStateDir()
	if err != nil {
		return err
	}

	// Read all state files
	entries, err := os.ReadDir(stateDir)
	if err != nil {
		if os.IsNotExist(err) {
			return nil // No state directory means nothing to clean
		}
		return fmt.Errorf("failed to read state directory: %w", err)
	}

	// Delete each .json file
	for _, entry := range entries {
		if entry.IsDir() {
			continue
		}
		if filepath.Ext(entry.Name()) == ".json" {
			statePath := filepath.Join(stateDir, entry.Name())
			if err := os.Remove(statePath); err != nil {
				utils.PrintDebug("[INSTANCE]Failed to delete state file %s: %v", statePath, err)
			} else {
				utils.PrintDebug("[INSTANCE]Deleted state file %s", statePath)
			}
		}
	}

	return nil
}

// deleteStatesMatching removes state files matching a pattern
func deleteStatesMatching(pattern string) error {
	stateDir, err := getStateDir()
	if err != nil {
		return err
	}

	// Read all state files
	entries, err := os.ReadDir(stateDir)
	if err != nil {
		if os.IsNotExist(err) {
			return nil
		}
		return fmt.Errorf("failed to read state directory: %w", err)
	}

	// Delete matching .json files
	for _, entry := range entries {
		if entry.IsDir() {
			continue
		}
		if filepath.Ext(entry.Name()) != ".json" {
			continue
		}

		// Extract instance name from filename (remove .json extension)
		instanceName := entry.Name()[:len(entry.Name())-5]

		// Check if instance name matches the pattern
		matched, err := filepath.Match(pattern, instanceName)
		if err != nil {
			utils.PrintDebug("[INSTANCE]Invalid pattern %s: %v", pattern, err)
			continue
		}

		if matched {
			statePath := filepath.Join(stateDir, entry.Name())
			if err := os.Remove(statePath); err != nil {
				utils.PrintDebug("[INSTANCE]Failed to delete state file %s: %v", statePath, err)
			} else {
				utils.PrintDebug("[INSTANCE]Deleted state file %s", statePath)
			}
		}
	}

	return nil
}
