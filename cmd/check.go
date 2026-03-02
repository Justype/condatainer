package cmd

import (
	"context"
	"errors"
	"fmt"
	"os"
	"path/filepath"
	"strings"

	"github.com/Justype/condatainer/internal/build"
	"github.com/Justype/condatainer/internal/config"
	"github.com/Justype/condatainer/internal/utils"
	"github.com/spf13/cobra"
)

var (
	checkAutoInstall     bool
	checkParseModuleLoad bool
)

var scriptCheckCmd = &cobra.Command{
	Use:   "check <script|dir> [script|dir ...]",
	Short: "Check if the dependencies of a script are installed",
	Long: `Check dependencies declared in build scripts using #DEP: tags.

Each argument can be:
  - A local .sh file path
  - A directory (all .sh files under it are checked recursively)
  - A package name (e.g., samtools/1.22) resolved from local or remote build scripts

Dependencies from all scripts are deduplicated and checked together.

Note: If creation jobs are submitted to a scheduler, exits with code 3.`,
	Example: `  condatainer check script.sh             # Check single local script
  condatainer check samtools/1.22         # Check build script by name
  condatainer check script1.sh script2.sh # Check multiple scripts
  condatainer check ./my-scripts/         # Check all .sh files in a directory
  condatainer check -a ./src/*sh          # Check and auto-install missing dependencies from multiple scripts`,
	Args:         cobra.MinimumNArgs(1),
	SilenceUsage: true, // Runtime errors should not show usage
	RunE:         runCheck,
}

func init() {
	rootCmd.AddCommand(scriptCheckCmd)
	scriptCheckCmd.Flags().BoolVarP(&checkAutoInstall, "auto-install", "a", false, "Automatically install missing dependencies")
	scriptCheckCmd.Flags().BoolP("install", "i", false, "Alias for --auto-install")
	scriptCheckCmd.Flags().BoolVar(&checkParseModuleLoad, "module", false, "Also parse 'module load' / 'ml' lines as dependencies")
}

func runCheck(cmd *cobra.Command, args []string) error {
	// Check if -i flag was used (alias for -a)
	if installFlag, _ := cmd.Flags().GetBool("install"); installFlag {
		checkAutoInstall = true
	}

	// Resolve all args to concrete script paths
	scriptPaths, remotePaths, err := resolveAllScriptPaths(args)
	if err != nil {
		return err
	}

	// Clean up any downloaded remote scripts at the end
	for _, rp := range remotePaths {
		defer os.Remove(rp)
	}

	if len(scriptPaths) == 0 {
		utils.PrintWarning("No scripts found to check.")
		return nil
	}

	// Collect dependencies from all scripts, deduplicating across scripts.
	// All deps are aggregated into a single list before any checking or installing.
	multiScript := len(scriptPaths) > 1
	parseModuleLoad := config.Global.ParseModuleLoad || checkParseModuleLoad

	seen := make(map[string]bool)
	var deps []string

	for _, scriptPath := range scriptPaths {
		if multiScript {
			utils.PrintMessage("Checking script: %s", utils.StylePath(scriptPath))
		}
		scriptDeps, parseErr := utils.GetDependenciesFromScript(scriptPath, parseModuleLoad)
		if parseErr != nil {
			return fmt.Errorf("failed to parse dependencies from %s: %w", scriptPath, parseErr)
		}
		for _, dep := range scriptDeps {
			key := dep
			if !utils.IsOverlay(dep) {
				key = utils.NormalizeNameVersion(dep)
			}
			if !seen[key] {
				seen[key] = true
				deps = append(deps, dep)
			}
		}
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
		// Check if it's an external overlay file by extension
		isExternalOverlay := utils.IsOverlay(dep)

		if isExternalOverlay {
			// External overlay file - check if file exists
			if utils.FileExists(dep) {
				utils.PrintMessage("Dependency: %s", utils.StyleName(dep))
			} else {
				utils.PrintMessage("Dependency: %s %s", utils.StyleName(dep), utils.StyleError("(missing local file)"))
				missingDeps = append(missingDeps, dep)
			}
		} else {
			// Package name - check if installed
			normalized := utils.NormalizeNameVersion(dep)
			if _, ok := installedOverlays[normalized]; ok {
				utils.PrintMessage("Dependency: %s", utils.StyleName(dep))
			} else {
				utils.PrintMessage("Dependency: %s %s", utils.StyleName(dep), utils.StyleError("(missing)"))
				missingDeps = append(missingDeps, dep)
			}
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
	if err := ensureBaseImage(cmd.Context()); err != nil {
		return err
	}

	// Check for custom overlays (files) that can't be auto-installed
	externalFiles := []string{}
	for _, dep := range missingDeps {
		if utils.IsOverlay(dep) {
			externalFiles = append(externalFiles, dep)
		}
	}

	if len(externalFiles) > 0 {
		fileList := strings.Join(externalFiles, ", ")
		utils.PrintWarning("External overlay(s) %s not found - cannot be auto-installed", utils.StyleName(fileList))
		// Remove external files from the install list
		nonExternalDeps := []string{}
		for _, dep := range missingDeps {
			if !utils.IsOverlay(dep) {
				nonExternalDeps = append(nonExternalDeps, dep)
			}
		}
		missingDeps = nonExternalDeps
	}

	if len(missingDeps) == 0 {
		if len(externalFiles) > 0 {
			return nil
		}
	}

	// Create build objects
	imagesDir, err := config.GetWritableImagesDir()
	if err != nil {
		return fmt.Errorf("no writable images directory found: %w", err)
	}
	buildObjects := make([]build.BuildObject, 0, len(missingDeps))
	for _, pkg := range missingDeps {
		bo, err := build.NewBuildObject(cmd.Context(), pkg, false, imagesDir, config.GetWritableTmpDir(), false)
		if err != nil {
			return fmt.Errorf("failed to create build object for %s: %w", pkg, err)
		}
		buildObjects = append(buildObjects, bo)
	}

	// Build graph and execute
	graph, err := build.NewBuildGraph(cmd.Context(), buildObjects, imagesDir, config.GetWritableTmpDir(), config.Global.SubmitJob, false)
	if err != nil {
		return fmt.Errorf("failed to create build graph: %w", err)
	}

	if err := graph.Run(cmd.Context()); err != nil {
		if errors.Is(err, build.ErrBuildCancelled) ||
			errors.Is(cmd.Context().Err(), context.Canceled) ||
			strings.Contains(err.Error(), "signal: killed") ||
			strings.Contains(err.Error(), "context canceled") {
			utils.PrintWarning("Installation cancelled.")
			return nil
		}
		return fmt.Errorf("some overlays failed to install: %w", err)
	}

	// If jobs were submitted to the scheduler, exit with the same distinct code
	// so downstream tooling can detect overlays will be created asynchronously.
	ExitIfJobsSubmitted(graph)

	utils.PrintSuccess("All selected overlays installed.")
	return nil
}

// findScriptsInDir recursively finds all .sh files under dir.
// If no scripts are found, a warning is printed but an empty slice is returned without error.
func findScriptsInDir(dir string) ([]string, error) {
	var found []string
	err := filepath.WalkDir(dir, func(path string, d os.DirEntry, err error) error {
		if err != nil {
			utils.PrintWarning("Skipping inaccessible path %s: %v", path, err)
			return nil // continue walking
		}
		if !d.IsDir() && strings.HasSuffix(d.Name(), ".sh") {
			found = append(found, path)
		}
		return nil
	})
	if err != nil {
		return nil, fmt.Errorf("failed to walk directory %s: %w", dir, err)
	}
	if len(found) == 0 {
		utils.PrintWarning("No .sh files found in directory %s", utils.StylePath(dir))
	}
	return found, nil
}

// resolveAllScriptPaths resolves all positional arguments to concrete file paths.
// Directories are expanded to all .sh files found recursively.
// Other args are resolved via resolveScriptPath (file, local build script, or remote download).
// Returns scriptPaths (deduplicated), remotePaths (to defer-remove), and any error.
func resolveAllScriptPaths(args []string) (scriptPaths []string, remotePaths []string, err error) {
	seen := make(map[string]bool)

	for _, arg := range args {
		if utils.DirExists(arg) {
			dirScripts, dirErr := findScriptsInDir(arg)
			if dirErr != nil {
				return nil, remotePaths, dirErr
			}
			for _, p := range dirScripts {
				if !seen[p] {
					seen[p] = true
					scriptPaths = append(scriptPaths, p)
				}
			}
			continue
		}

		resolved, isRemote, resolveErr := resolveScriptPath(arg)
		if resolveErr != nil {
			return nil, remotePaths, resolveErr
		}
		if isRemote {
			remotePaths = append(remotePaths, resolved)
		}
		if !seen[resolved] {
			seen[resolved] = true
			scriptPaths = append(scriptPaths, resolved)
		}
	}

	return scriptPaths, remotePaths, nil
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
		if err := os.Chmod(tempPath, utils.PermExec); err != nil {
			os.Remove(tempPath)
			return "", false, fmt.Errorf("failed to make script executable: %w", err)
		}

		utils.PrintMessage("Downloaded build script to %s", utils.StylePath(tempPath))
		return tempPath, true, nil
	}

	return "", false, fmt.Errorf("build script for %s not found", normalized)
}
