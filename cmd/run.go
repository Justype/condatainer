package cmd

import (
	"fmt"
	"os"
	"path/filepath"
	"strings"

	"github.com/Justype/condatainer/internal/apptainer"
	"github.com/Justype/condatainer/internal/build"
	"github.com/Justype/condatainer/internal/config"
	execpkg "github.com/Justype/condatainer/internal/exec"
	"github.com/Justype/condatainer/internal/utils"
	"github.com/spf13/cobra"
)

var (
	runWritableImg  bool
	runBaseImage    string
	runAutoInstall  bool
	runEnvSettings  []string
	runBindPaths    []string
	runFakeroot     bool
)

var runCmd = &cobra.Command{
	Use:   "run [script]",
	Short: "Run a script and auto-solve the dependencies by #DEP tags",
	Long: `Execute a script with automatic dependency resolution.

The script can contain special comment tags:
  #DEP: package/version  - Declares a dependency
  #CNT args              - Additional arguments to pass to condatainer

Supported #CNT arguments:
  -w, --writable-img     Make .img overlays writable
  -b, --base-image PATH  Use custom base image
  --env KEY=VALUE        Set environment variable
  --bind HOST:CONTAINER  Bind mount path
  --fakeroot             Run with fakeroot privileges

Dependencies will be automatically loaded, and missing ones can be auto-installed
with the --auto-install flag.`,
	Example: `  condatainer run script.sh              # Run with dependency check
  condatainer run script.sh -a           # Auto-install missing deps
  condatainer run script.sh -w           # Make .img overlays writable
  condatainer run script.sh -b base.sif  # Use custom base image

  # In script.sh, you can use #CNT to set options:
  # #CNT -w
  # #CNT --env MYVAR=value
  # #CNT --bind /data:/mnt/data`,
	Args: cobra.ExactArgs(1),
	RunE: runScript,
}

func init() {
	rootCmd.AddCommand(runCmd)
	runCmd.Flags().BoolVarP(&runWritableImg, "writable-img", "w", false, "Make .img overlays writable (default: read-only)")
	runCmd.Flags().StringVarP(&runBaseImage, "base-image", "b", "", "Base image to use instead of default")
	runCmd.Flags().BoolVarP(&runAutoInstall, "auto-install", "a", false, "Automatically install missing dependencies")
	runCmd.Flags().BoolP("install", "i", false, "Alias for --auto-install")
}

func runScript(cmd *cobra.Command, args []string) error {
	scriptPath := args[0]

	// Check if -i flag was used (alias for -a)
	if installFlag, _ := cmd.Flags().GetBool("install"); installFlag {
		runAutoInstall = true
	}

	// Ensure base image exists
	if err := apptainer.EnsureBaseImage(); err != nil {
		return err
	}

	// Validate script exists
	if !utils.FileExists(scriptPath) {
		return fmt.Errorf("script file %s not found", utils.StylePath(scriptPath))
	}

	scriptPath, err := filepath.Abs(scriptPath)
	if err != nil {
		return fmt.Errorf("failed to get absolute path: %w", err)
	}

	// Parse additional arguments from script (#CNT tags)
	scriptArgs, err := parseArgsInScript(scriptPath)
	if err != nil {
		return err
	}

	// Apply #CNT arguments
	if len(scriptArgs) > 0 {
		utils.PrintDebug("[RUN] Additional script arguments found: %v", scriptArgs)
		applyScriptArgs(scriptArgs)
	}

	// Apply base image if specified
	if runBaseImage != "" {
		config.Global.BaseImage = runBaseImage
	}

	// Get dependencies from script
	deps, err := utils.GetDependenciesFromScript(scriptPath)
	if err != nil {
		return fmt.Errorf("failed to parse dependencies: %w", err)
	}

	// Check for missing dependencies
	installedOverlays, err := getInstalledOverlaysMap()
	if err != nil {
		return err
	}

	missingDeps := []string{}
	for _, dep := range deps {
		normalized := utils.NormalizeNameVersion(dep)
		if _, ok := installedOverlays[normalized]; !ok {
			missingDeps = append(missingDeps, dep)
		}
	}

	// Handle missing dependencies
	if len(missingDeps) > 0 {
		if !runAutoInstall {
			utils.PrintMessage("Missing dependencies:")
			for _, md := range missingDeps {
				utils.PrintMessage("  - %s", utils.StyleWarning(md))
			}
			utils.PrintMessage("Please run %s to install missing dependencies.",
				utils.StyleAction("condatainer check -a "+scriptPath))
			return nil
		}

		// Auto-install missing dependencies
		utils.PrintMessage("Attempting to auto-install missing dependencies...")

		// Check for custom overlays
		for _, md := range missingDeps {
			if utils.IsOverlay(md) {
				return fmt.Errorf("custom overlay %s is missing - cannot proceed with auto-installation", utils.StyleName(md))
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
	}

	// Resolve overlay paths
	overlays := make([]string, len(deps))
	for i, dep := range deps {
		if utils.IsOverlay(dep) {
			overlays[i] = dep
		} else {
			normalized := utils.NormalizeNameVersion(dep)
			if path, ok := installedOverlays[normalized]; ok {
				overlays[i] = path
			} else {
				overlays[i] = dep
			}
		}
	}

	// Prepare execution command
	// Disable module commands and run the script
	executionScript := `module() { :; }
ml() { :; }
export -f module ml
`

	// Check if script is executable
	fileInfo, err := os.Stat(scriptPath)
	if err != nil {
		return fmt.Errorf("failed to stat script: %w", err)
	}

	if fileInfo.Mode()&0111 != 0 {
		// Script is executable
		executionScript += scriptPath
	} else {
		// Run with bash
		executionScript += "bash " + scriptPath
	}

	// Execute with overlays
	options := execpkg.Options{
		Overlays:     overlays,
		Command:      []string{"/bin/bash", "-c", executionScript},
		WritableImg:  runWritableImg,
		EnvSettings:  runEnvSettings,
		BindPaths:    runBindPaths,
		Fakeroot:     runFakeroot,
		BaseImage:    config.Global.BaseImage,
		ApptainerBin: config.Global.ApptainerBin,
	}

	return execpkg.Run(options)
}

// applyScriptArgs parses #CNT arguments and applies them to run options
func applyScriptArgs(scriptArgs []string) {
	for _, argLine := range scriptArgs {
		parts := strings.Fields(argLine)
		for i := 0; i < len(parts); i++ {
			arg := parts[i]
			switch arg {
			case "-w", "--writable-img":
				runWritableImg = true
			case "--fakeroot":
				runFakeroot = true
			case "-b", "--base-image":
				if i+1 < len(parts) {
					i++
					runBaseImage = parts[i]
				}
			case "--env":
				if i+1 < len(parts) {
					i++
					runEnvSettings = append(runEnvSettings, parts[i])
				}
			case "--bind":
				if i+1 < len(parts) {
					i++
					runBindPaths = append(runBindPaths, parts[i])
				}
			default:
				// Handle --env=VALUE and --bind=VALUE formats
				if strings.HasPrefix(arg, "--env=") {
					runEnvSettings = append(runEnvSettings, strings.TrimPrefix(arg, "--env="))
				} else if strings.HasPrefix(arg, "--bind=") {
					runBindPaths = append(runBindPaths, strings.TrimPrefix(arg, "--bind="))
				} else if strings.HasPrefix(arg, "-b=") {
					runBaseImage = strings.TrimPrefix(arg, "-b=")
				} else if strings.HasPrefix(arg, "--base-image=") {
					runBaseImage = strings.TrimPrefix(arg, "--base-image=")
				}
			}
		}
	}
}

// parseArgsInScript extracts arguments from #CNT comments in the script
func parseArgsInScript(scriptPath string) ([]string, error) {
	data, err := os.ReadFile(scriptPath)
	if err != nil {
		return nil, fmt.Errorf("failed to read script: %w", err)
	}

	args := []string{}
	lines := strings.Split(string(data), "\n")

	for _, line := range lines {
		line = strings.TrimSpace(line)
		if strings.HasPrefix(line, "#CNT") {
			argLine := strings.TrimSpace(strings.TrimPrefix(line, "#CNT"))
			if argLine != "" {
				args = append(args, argLine)
			}
		}
	}

	return args, nil
}
