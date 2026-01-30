package cmd

import (
	"fmt"
	"os"

	"github.com/Justype/condatainer/internal/apptainer"
	"github.com/Justype/condatainer/internal/build"
	"github.com/Justype/condatainer/internal/config"
	"github.com/Justype/condatainer/internal/utils"
	"github.com/spf13/cobra"
)

var (
	checkAutoInstall bool
)

var scriptCheckCmd = &cobra.Command{
	Use:   "check [script]",
	Short: "Check if the dependencies of a script are installed",
	Long: `Check dependencies declared in a build script using #DEP: tags.

The script can be a local file path or a package name (e.g., samtools/1.22).
If the script is not found locally, it will be downloaded from remote metadata.

Dependencies are declared in scripts using comments like:
  #DEP: package/version
  #DEP: another-package/1.0`,
	Example: `  condatainer check script.sh           # Check local script
  condatainer check samtools/1.22       # Check build script by name
  condatainer check script.sh -a        # Check and auto-install missing deps`,
	Args:         cobra.ExactArgs(1),
	SilenceUsage: true, // Runtime errors should not show usage
	RunE:         runCheck,
}

func init() {
	rootCmd.AddCommand(scriptCheckCmd)
	scriptCheckCmd.Flags().BoolVarP(&checkAutoInstall, "auto-install", "a", false, "Automatically install missing dependencies")
	scriptCheckCmd.Flags().BoolP("install", "i", false, "Alias for --auto-install")
}

func runCheck(cmd *cobra.Command, args []string) error {
	scriptPathOrName := args[0]

	// Check if -i flag was used (alias for -a)
	if installFlag, _ := cmd.Flags().GetBool("install"); installFlag {
		checkAutoInstall = true
	}

	// Resolve script path
	scriptPath, isRemote, err := resolveScriptPath(scriptPathOrName)
	if err != nil {
		return err
	}

	// Clean up remote script at the end
	if isRemote {
		defer os.Remove(scriptPath)
	}

	// Get dependencies from script
	deps, err := utils.GetDependenciesFromScript(scriptPath)
	if err != nil {
		return fmt.Errorf("failed to parse dependencies: %w", err)
	}

	if len(deps) == 0 {
		utils.PrintMessage("No dependencies found in script.")
		return nil
	}

	// Check which overlays are installed
	installedOverlays, err := getInstalledOverlaysMap()
	if err != nil {
		return err
	}

	missingDeps := []string{}
	for _, dep := range deps {
		normalized := utils.NormalizeNameVersion(dep)
		if _, ok := installedOverlays[normalized]; ok {
			utils.PrintMessage("Dependency: %s", utils.StyleName(dep))
		} else {
			utils.PrintMessage("Dependency: %s %s", utils.StyleName(dep), utils.StyleError("(missing)"))
			missingDeps = append(missingDeps, dep)
		}
	}

	if len(missingDeps) == 0 {
		utils.PrintSuccess("All dependencies are installed.")
		return nil
	}

	// Handle missing dependencies
	if !checkAutoInstall {
		utils.PrintMessage("Run the command again with %s or %s to install missing dependencies.",
			utils.StyleAction("-a"), utils.StyleAction("--auto-install"))
		return nil
	}

	// Auto-install missing dependencies
	utils.PrintMessage("Attempting to auto-install missing dependencies...")

	// Ensure base image exists
	if err := apptainer.EnsureBaseImage(); err != nil {
		return err
	}

	// Check for custom overlays (files) that can't be auto-installed
	for _, dep := range missingDeps {
		if utils.IsOverlay(dep) {
			return fmt.Errorf("custom overlay %s is missing - cannot proceed with auto-installation", utils.StyleName(dep))
		}
	}

	// Create build objects
	buildObjects := make([]build.BuildObject, 0, len(missingDeps))
	for _, pkg := range missingDeps {
		bo, err := build.NewBuildObject(pkg, false, "", "")
		if err != nil {
			return fmt.Errorf("failed to create build object for %s: %w", pkg, err)
		}
		buildObjects = append(buildObjects, bo)
	}

	// Build graph and execute
	graph, err := build.NewBuildGraph(buildObjects, "", "", config.Global.SubmitJob)
	if err != nil {
		return fmt.Errorf("failed to create build graph: %w", err)
	}

	if err := graph.Run(); err != nil {
		return fmt.Errorf("some overlays failed to install: %w", err)
	}

	utils.PrintSuccess("All selected overlays installed/submitted.")
	return nil
}

// resolveScriptPath resolves a script path or name to an actual file path
// Returns (path, isRemote, error)
func resolveScriptPath(scriptPathOrName string) (string, bool, error) {
	// If it's a file, use it directly
	if utils.FileExists(scriptPathOrName) {
		return scriptPathOrName, false, nil
	}

	// Try to find as build script
	normalized := utils.NormalizeNameVersion(scriptPathOrName)
	utils.PrintDebug("[CHECK] Checking for build script %s locally and remotely...", normalized)

	// Check local build scripts
	localScripts, err := build.GetLocalBuildScripts()
	if err == nil {
		if info, ok := localScripts[normalized]; ok {
			utils.PrintMessage("Found local build script %s", utils.StylePath(info.Path))
			return info.Path, false, nil
		}
	}

	// Check remote build scripts
	remoteScripts, err := build.GetRemoteBuildScripts()
	if err != nil {
		return "", false, fmt.Errorf("build script for %s not found", normalized)
	}

	if info, ok := remoteScripts[normalized]; ok {
		utils.PrintMessage("Downloading build script for %s from remote metadata...", normalized)

		// Download to temp location
		tempPath := fmt.Sprintf("/tmp/%s.sh", utils.NormalizeNameVersion(normalized))
		if err := utils.DownloadFile(info.Path, tempPath); err != nil {
			return "", false, fmt.Errorf("failed to download build script: %w", err)
		}

		// Make executable
		if err := os.Chmod(tempPath, 0755); err != nil {
			os.Remove(tempPath)
			return "", false, fmt.Errorf("failed to make script executable: %w", err)
		}

		utils.PrintMessage("Downloaded build script to %s", utils.StylePath(tempPath))
		return tempPath, true, nil
	}

	return "", false, fmt.Errorf("build script for %s not found", normalized)
}
