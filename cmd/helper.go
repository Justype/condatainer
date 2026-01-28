package cmd

import (
	"encoding/json"
	"fmt"
	"io"
	"net/http"
	"os"
	"os/exec"
	"path/filepath"

	"github.com/Justype/condatainer/internal/config"
	"github.com/Justype/condatainer/internal/scheduler"
	"github.com/Justype/condatainer/internal/utils"
	"github.com/spf13/cobra"
)

var (
	helperPath   bool
	helperUpdate bool
)

var helperCmd = &cobra.Command{
	Use:   "helper [script-name] [args...]",
	Short: "Manage and run helper scripts",
	Long: `Manage helper scripts for running services like code-server, rstudio-server, etc.

Helper scripts come in two categories:
  - headless: Direct execution scripts
  - slurm: SLURM job submission scripts (auto-selected when scheduler available)

Common helper scripts:
  - code-server: VS Code server
  - rstudio-server: RStudio server
  - vscode-tunnel: VS Code tunnel`,
	Example: `  condatainer helper --update           # Update all helper scripts
  condatainer helper --update code-server  # Update specific script
  condatainer helper --path              # Show helper scripts directory
  condatainer helper code-server 8080    # Run code-server on port 8080`,
	RunE: runHelper,
}

func init() {
	rootCmd.AddCommand(helperCmd)
	helperCmd.Flags().BoolVar(&helperPath, "path", false, "Show path to helper scripts directory")
	helperCmd.Flags().BoolVarP(&helperUpdate, "update", "u", false, "Update helper scripts from remote")
}

func runHelper(cmd *cobra.Command, args []string) error {
	helperScriptsDir := filepath.Join(config.Global.BaseDir, "helper-scripts")

	// --- Path Mode ---
	if helperPath {
		if len(args) > 0 {
			scriptName := args[0]
			localPath := filepath.Join(helperScriptsDir, scriptName)
			if utils.FileExists(localPath) {
				absPath, _ := filepath.Abs(localPath)
				fmt.Println(absPath)
				return nil
			}
			return fmt.Errorf("helper script '%s' not found locally", scriptName)
		}

		// Show directory path
		if err := os.MkdirAll(helperScriptsDir, 0775); err != nil {
			return fmt.Errorf("failed to create helper scripts directory: %w", err)
		}
		absPath, _ := filepath.Abs(helperScriptsDir)
		fmt.Println(absPath)
		return nil
	}

	// --- Update Mode ---
	if helperUpdate {
		return updateHelperScripts(args, helperScriptsDir)
	}

	// --- Run Mode ---
	if len(args) == 0 {
		return fmt.Errorf("no helper script name provided - use --update to update scripts")
	}

	scriptName := args[0]
	scriptArgs := args[1:]

	localPath := filepath.Join(helperScriptsDir, scriptName)
	if !utils.FileExists(localPath) {
		return fmt.Errorf("helper script '%s' not found locally at %s\nRun '%s' to fetch available helper scripts",
			scriptName, localPath, utils.StyleAction("condatainer helper --update"))
	}

	// Ensure executable
	if err := os.Chmod(localPath, 0755); err != nil {
		utils.PrintDebug("Failed to chmod helper script: %v", err)
	}

	// Run the helper script
	cmdArgs := append([]string{localPath}, scriptArgs...)
	utils.PrintDebug("Running helper command: %v", cmdArgs)

	helperCmd := exec.Command(cmdArgs[0], cmdArgs[1:]...)
	helperCmd.Stdin = os.Stdin
	helperCmd.Stdout = os.Stdout
	helperCmd.Stderr = os.Stderr

	return helperCmd.Run()
}

// updateHelperScripts downloads and updates helper scripts from remote metadata
func updateHelperScripts(args []string, helperScriptsDir string) error {
	// Fetch remote metadata
	metadata, err := fetchRemoteHelperMetadata()
	if err != nil {
		return fmt.Errorf("failed to fetch remote helper metadata: %w", err)
	}

	// Choose category based on scheduler availability
	category := "headless"
	if config.Global.SubmitJob {
		if sched, err := scheduler.DetectSchedulerWithBinary(config.Global.SchedulerBin); err == nil {
			info := sched.GetInfo()
			if info.Available && !info.InJob {
				category = "slurm"
			}
		}
	}

	utils.PrintDebug("Helper category selected: %s", category)

	entries, ok := metadata[category]
	if !ok {
		return fmt.Errorf("no helper scripts found for category '%s'", category)
	}

	// Create helper scripts directory
	if err := os.MkdirAll(helperScriptsDir, 0775); err != nil {
		return fmt.Errorf("failed to create helper scripts directory: %w", err)
	}

	// Filter to specific script if provided
	if len(args) > 0 {
		scriptName := args[0]
		if entry, ok := entries[scriptName]; ok {
			entries = map[string]HelperScriptEntry{scriptName: entry}
		} else {
			return fmt.Errorf("helper script '%s' not found in remote metadata for category '%s'", scriptName, category)
		}
	} else {
		utils.PrintMessage("Updating all helper scripts for %s...", category)
	}

	// Download each helper script
	for name, entry := range entries {
		if entry.Path == "" {
			continue
		}

		url := fmt.Sprintf("https://raw.githubusercontent.com/%s/main/%s",
			config.GITHUB_REPO, entry.Path)
		destName := filepath.Base(entry.Path)
		dest := filepath.Join(helperScriptsDir, destName)

		utils.PrintMessage("Updating %s/%s", category, name)
		if err := downloadExecutable(url, dest); err != nil {
			utils.PrintWarning("Failed to update %s/%s: %v", category, name, err)
			continue
		}
	}

	utils.PrintSuccess("Helper update finished.")
	return nil
}

// HelperScriptEntry represents a helper script metadata entry
type HelperScriptEntry struct {
	Path string `json:"path"`
}

// fetchRemoteHelperMetadata fetches the helper scripts metadata from GitHub
func fetchRemoteHelperMetadata() (map[string]map[string]HelperScriptEntry, error) {
	url := fmt.Sprintf("https://raw.githubusercontent.com/%s/main/metadata/helper-scripts.json",
		config.GITHUB_REPO)

	resp, err := http.Get(url)
	if err != nil {
		return nil, err
	}
	defer resp.Body.Close()

	if resp.StatusCode != http.StatusOK {
		return nil, fmt.Errorf("failed to fetch metadata: HTTP %d", resp.StatusCode)
	}

	data, err := io.ReadAll(resp.Body)
	if err != nil {
		return nil, err
	}

	var metadata map[string]map[string]HelperScriptEntry
	if err := json.Unmarshal(data, &metadata); err != nil {
		return nil, err
	}

	return metadata, nil
}

// downloadExecutable downloads a file and makes it executable
func downloadExecutable(url, destPath string) error {
	resp, err := http.Get(url)
	if err != nil {
		return err
	}
	defer resp.Body.Close()

	if resp.StatusCode != http.StatusOK {
		return fmt.Errorf("HTTP %d", resp.StatusCode)
	}

	// Create parent directory
	if err := os.MkdirAll(filepath.Dir(destPath), 0775); err != nil {
		return err
	}

	// Write file
	out, err := os.Create(destPath)
	if err != nil {
		return err
	}
	defer out.Close()

	if _, err := io.Copy(out, resp.Body); err != nil {
		return err
	}

	// Make executable
	return os.Chmod(destPath, 0755)
}
