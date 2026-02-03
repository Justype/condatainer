package cmd

import (
	"encoding/json"
	"fmt"
	"io"
	"net/http"
	"os"
	"os/exec"
	"path/filepath"
	"strings"

	"github.com/Justype/condatainer/internal/apptainer"
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
	Use:   "helper [flags] [script-name] [script-args...]",
	Short: "Manage and run helper scripts",
	Long: `Manage helper scripts for running services like code-server, rstudio-server, etc.

Available helper scripts:
  - code-server: code server
  - rstudio-server: RStudio server
  - vscode-tunnel: VS Code tunnel

Note: Helper is not available inside a container or a scheduler job.`,
	Example: `  condatainer helper --update                # Update all helper scripts
  condatainer helper --path                  # Show helper scripts directory
  condatainer helper code-server -p 11908    # Run code-server and forward flags`,
	RunE:              runHelper,
	ValidArgsFunction: completeHelperScripts,
}

func init() {
	rootCmd.AddCommand(helperCmd)
	helperCmd.Flags().BoolVar(&helperPath, "path", false, "Show path to helper scripts directory")
	helperCmd.Flags().BoolVarP(&helperUpdate, "update", "u", false, "Update helper scripts from remote")

	// Stop flag parsing after the first positional argument so script flags (like -p or -v)
	// are passed through to the helper script rather than being interpreted by cobra.
	helperCmd.Flags().SetInterspersed(false)
}

// completeHelperScripts provides shell completion for helper script names
func completeHelperScripts(cmd *cobra.Command, args []string, toComplete string) ([]string, cobra.ShellCompDirective) {
	// Only complete the first argument (script name), then return to default file completion
	if len(args) > 0 {
		return nil, cobra.ShellCompDirectiveDefault
	}

	// Collect scripts from all search paths (user → system)
	seen := make(map[string]bool)
	var scripts []string

	for _, dir := range config.GetHelperScriptSearchPaths() {
		entries, err := os.ReadDir(dir)
		if err != nil {
			continue
		}

		for _, entry := range entries {
			name := entry.Name()
			// Skip hidden files and directories
			if !entry.IsDir() && !strings.HasPrefix(name, ".") {
				if !seen[name] {
					seen[name] = true
					scripts = append(scripts, name)
				}
			}
		}
	}

	return scripts, cobra.ShellCompDirectiveNoFileComp
}

func runHelper(cmd *cobra.Command, args []string) error {
	// Prevent helper commands when inside a container or scheduler job (except --path)
	if config.IsInsideContainer() && !helperPath {
		cmd.SilenceUsage = true
		utils.PrintError("helper commands are not available inside a container")
		os.Exit(1)
	}
	if scheduler.IsInsideJob() && !helperPath {
		cmd.SilenceUsage = true
		utils.PrintError("helper commands are not available inside a scheduler job")
		os.Exit(1)
	}

	// --- Path Mode ---
	if helperPath {
		if len(args) > 0 {
			// Find specific script in all search paths
			scriptName := args[0]
			scriptPath, err := config.FindHelperScript(scriptName)
			if err != nil {
				cmd.SilenceUsage = true
				utils.PrintError("helper script '%s' not found (searched: %v)", scriptName, config.GetHelperScriptSearchPaths())
				os.Exit(1)
			}
			fmt.Println(scriptPath)
			return nil
		}

		// Show all helper script search paths
		fmt.Println("Helper script search paths (priority order):")
		for i, dir := range config.GetHelperScriptSearchPaths() {
			exists := ""
			if !config.DirExists(dir) {
				exists = " (not found)"
			}
			fmt.Printf("  %d. %s%s\n", i+1, dir, exists)
		}

		// Show writable directory
		if writableDir, err := config.GetWritableHelperScriptsDir(); err == nil {
			fmt.Printf("\nWritable directory: %s\n", writableDir)
		}
		return nil
	}

	// --- Update Mode ---
	if helperUpdate {
		// Get writable directory for updates
		helperScriptsDir, err := config.GetWritableHelperScriptsDir()
		if err != nil {
			cmd.SilenceUsage = true
			utils.PrintError("failed to find writable helper scripts directory: %v", err)
			os.Exit(1)
		}
		if err := updateHelperScripts(args, helperScriptsDir); err != nil {
			cmd.SilenceUsage = true
			utils.PrintError("%v", err)
			os.Exit(1)
		}
		return nil
	}

	// --- Run Mode ---
	if len(args) == 0 {
		// Parameter error - show usage
		return fmt.Errorf("no helper script name provided - use --update to update scripts")
	}

	scriptName := args[0]
	scriptArgs := args[1:]

	// Find script in all search paths
	scriptPath, err := config.FindHelperScript(scriptName)
	if err != nil {
		cmd.SilenceUsage = true
		utils.PrintError("helper script '%s' not found\nSearched: %v\nRun '%s' to fetch available helper scripts",
			scriptName, config.GetHelperScriptSearchPaths(), utils.StyleAction("condatainer helper --update"))
		os.Exit(1)
	}

	// Ensure executable
	if err := os.Chmod(scriptPath, 0755); err != nil {
		utils.PrintDebug("Failed to chmod helper script: %v", err)
	}

	// Check apptainer is available before running helper script
	if err := apptainer.EnsureApptainer(); err != nil {
		cmd.SilenceUsage = true
		utils.PrintError("apptainer is required to run helper scripts: %v", err)
		os.Exit(1)
	}

	// Run the helper script
	cmdArgs := append([]string{scriptPath}, scriptArgs...)
	utils.PrintDebug("Running helper command: %v", cmdArgs)

	helperCmd := exec.Command(cmdArgs[0], cmdArgs[1:]...)
	helperCmd.Stdin = os.Stdin
	helperCmd.Stdout = os.Stdout
	helperCmd.Stderr = os.Stderr

	if err := helperCmd.Run(); err != nil {
		// Propagate the exit code from the helper script
		if exitErr, ok := err.(*exec.ExitError); ok {
			os.Exit(exitErr.ExitCode())
		}
		return err
	}
	return nil
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
