package cmd

import (
	"fmt"
	"os"
	"os/exec"
	"path/filepath"
	"strings"

	"github.com/Justype/condatainer/internal/apptainer"
	"github.com/Justype/condatainer/internal/config"
	"github.com/Justype/condatainer/internal/utils"
	"github.com/spf13/cobra"
	"github.com/spf13/viper"
)

var (
	showPath     bool
	initLocation string // user, portable, system, or a custom path
)

// configKeys is the list of known configuration keys for shell completion
var configKeys = []string{
	"logs_dir",
	"apptainer_bin",
	"scheduler_bin",
	"submit_job",
	"branch",
	"prefer_remote",
	"extra_base_dirs",
	"scheduler.nodes",
	"scheduler.tasks_per_node",
	"scheduler.ncpus_per_task",
	"scheduler.mem_mb_per_node",
	"scheduler.time",
	"build.ncpus",
	"build.mem_mb",
	"build.time",
	"build.tmp_size_mb",
	"build.compress_args",
	"build.overlay_type",
}

// configKeysCompletion returns config keys for shell completion
func configKeysCompletion(cmd *cobra.Command, args []string, toComplete string) ([]string, cobra.ShellCompDirective) {
	if len(args) == 0 {
		// First arg: complete config keys
		return configKeys, cobra.ShellCompDirectiveNoFileComp
	}
	if len(args) == 1 {
		// Second arg: complete values based on the key
		return configValueCompletion(args[0]), cobra.ShellCompDirectiveNoFileComp
	}
	return nil, cobra.ShellCompDirectiveNoFileComp
}

// configValueCompletion returns suggested values for a config key
func configValueCompletion(key string) []string {
	switch key {
	case "submit_job":
		return []string{"true", "false"}
	case "prefer_remote":
		return []string{"true", "false"}
	case "branch":
		return []string{"main", "dev"}
	case "scheduler.nodes":
		return []string{"1", "2", "4", "8"}
	case "scheduler.tasks_per_node":
		return []string{"1", "2", "4", "8"}
	case "scheduler.ncpus_per_task", "build.ncpus":
		return []string{"4", "8", "16", "32"}
	case "scheduler.mem_mb_per_node", "build.mem_mb":
		return []string{"4096", "8192", "16384", "32768"}
	case "scheduler.time", "build.time":
		return []string{"1h", "2h", "4h", "8h"}
	case "build.tmp_size_mb":
		return []string{"10240", "20480", "40960"}
	case "build.overlay_type":
		return []string{"ext3", "squashfs"}
	default:
		return nil
	}
}

var configCmd = &cobra.Command{
	Use:   "config",
	Short: "Manage condatainer configuration",
	Long: `Manage condatainer configuration settings.

Configuration file priority (highest to lowest):
  1. Command-line flags
  2. Environment variables (CONDATAINER_*)
  3. User config file (~/.config/condatainer/config.yaml)
  4. Portable config (<install-dir>/config.yaml, if in dedicated folder)
  5. System config file (/etc/condatainer/config.yaml)
  6. Defaults

Data directory priority (for read/write operations):
  1. Extra base directories (CONDATAINER_EXTRA_BASE_DIRS or config)
  2. Portable directory (group/shared, from executable location)
  3. User scratch directory ($SCRATCH/condatainer, HPC systems)
  4. User XDG directory (~/.local/share/condatainer)`,
}

var configShowCmd = &cobra.Command{
	Use:   "show",
	Short: "Show current configuration",
	Long: `Display current configuration values and their sources.

Shows:
  - Config file search paths and which one is in use
  - Base directory search paths with writable status
  - All configuration settings (directories, binaries, build config)
  - Environment variable overrides`,
	Run: func(cmd *cobra.Command, args []string) {
		if showPath {
			configPath, err := config.GetUserConfigPath()
			if err != nil {
				ExitWithError("Failed to get config path: %v", err)
			}
			fmt.Println(configPath)
			return
		}

		// Show config file search paths
		fmt.Println(utils.StyleTitle("Config File Search Paths:"))
		searchPaths := config.GetConfigSearchPaths()
		foundActive := false
		for i, sp := range searchPaths {
			status := ""
			if sp.InUse {
				status = " " + utils.StyleSuccess("← in use")
				foundActive = true
			} else if sp.Exists {
				status = " " + utils.StyleInfo("(exists)")
			}
			fmt.Printf("  %d. [%s] %s%s\n", i+1, sp.Type, sp.Path, status)
		}
		if !foundActive {
			fmt.Printf("  %s (use 'condatainer config init' to create)\n", utils.StyleWarning("No config file found"))
		}
		fmt.Println()

		// Show base directory search paths
		fmt.Println(utils.StyleTitle("Base Directory Search Paths:"))

		// Helper function to check if directory is writable
		isWritable := func(dir string) bool {
			testFile := filepath.Join(dir, ".write-test")
			f, err := os.Create(testFile)
			if err != nil {
				return false
			}
			f.Close()
			os.Remove(testFile)
			return true
		}

		pathIndex := 1

		// Extra base directories
		extraBaseDirs := config.GetExtraBaseDirs()
		for _, dir := range extraBaseDirs {
			status := ""
			if _, err := os.Stat(dir); err == nil {
				if isWritable(dir) {
					status = " " + utils.StyleSuccess("(writable)")
				} else {
					status = " " + utils.StyleWarning("(read-only)")
				}
			}
			fmt.Printf("  %d. [extra] %s%s\n", pathIndex, dir, status)
			pathIndex++
		}

		// Portable directory
		if portableDir := config.GetPortableDataDir(); portableDir != "" {
			status := ""
			if _, err := os.Stat(portableDir); err == nil {
				if isWritable(portableDir) {
					status = " " + utils.StyleSuccess("(writable)")
				} else {
					status = " " + utils.StyleWarning("(read-only)")
				}
			}
			fmt.Printf("  %d. [portable] %s%s\n", pathIndex, portableDir, status)
			pathIndex++
		}

		// Scratch directory
		if scratchDir := config.GetScratchDataDir(); scratchDir != "" {
			status := ""
			if _, err := os.Stat(scratchDir); err == nil {
				if isWritable(scratchDir) {
					status = " " + utils.StyleSuccess("(writable)")
				} else {
					status = " " + utils.StyleWarning("(read-only)")
				}
			}
			fmt.Printf("  %d. [scratch] %s%s\n", pathIndex, scratchDir, status)
			pathIndex++
		}

		// User directory
		if userDir := config.GetUserDataDir(); userDir != "" {
			status := ""
			if _, err := os.Stat(userDir); err == nil {
				if isWritable(userDir) {
					status = " " + utils.StyleSuccess("(writable)")
				} else {
					status = " " + utils.StyleWarning("(read-only)")
				}
			}
			fmt.Printf("  %d. [user] %s%s\n", pathIndex, userDir, status)
			pathIndex++
		}

		fmt.Println()

		// Show all settings
		fmt.Println(utils.StyleTitle("Current Configuration:"))
		fmt.Println()

		// Directories
		fmt.Println(utils.StyleTitle("Directories:"))
		logsDir := viper.GetString("logs_dir")
		if logsDir == "" {
			logsDir = config.Global.LogsDir + " (default)"
		}
		fmt.Printf("  logs_dir:       %s\n", logsDir)
		fmt.Println()

		// Binary paths
		fmt.Println(utils.StyleTitle("Binaries:"))
		fmt.Printf("  apptainer_bin:  %s\n", viper.GetString("apptainer_bin"))
		schedulerBin := viper.GetString("scheduler_bin")
		schedulerType := config.GetSchedulerTypeFromBin(schedulerBin)
		if schedulerBin != "" {
			fmt.Printf("  scheduler_bin:  %s (%s)\n", schedulerBin, schedulerType)
		} else {
			fmt.Printf("  scheduler_bin:  %s\n", schedulerBin)
		}
		fmt.Println()

		// Runtime settings
		fmt.Println(utils.StyleTitle("Runtime:"))
		submitJobConfig := viper.GetBool("submit_job")
		submitJobActual := config.Global.SubmitJob
		if submitJobConfig && !submitJobActual {
			fmt.Printf("  submit_job:     %v (disabled: scheduler not accessible)\n", submitJobConfig)
		} else {
			fmt.Printf("  submit_job:     %v\n", submitJobActual)
		}
		fmt.Printf("  branch:         %s\n", viper.GetString("branch"))
		fmt.Printf("  prefer_remote:  %v\n", config.Global.PreferRemote)
		fmt.Println()

		// Scheduler default specs
		fmt.Println(utils.StyleTitle("Scheduler Default Specs:"))
		fmt.Printf("  nodes:           %d\n", viper.GetInt("scheduler.nodes"))
		fmt.Printf("  tasks_per_node:  %d\n", viper.GetInt("scheduler.tasks_per_node"))
		fmt.Printf("  ncpus_per_task:  %d\n", viper.GetInt("scheduler.ncpus_per_task"))
		schedMemMB := viper.GetInt64("scheduler.mem_mb_per_node")
		if schedMemMB > 0 {
			fmt.Printf("  mem_mb_per_node: %d\n", schedMemMB)
		} else {
			fmt.Printf("  mem_mb_per_node: %s\n", utils.StyleInfo("not set"))
		}
		schedTime := viper.GetString("scheduler.time")
		if schedTime != "" {
			fmt.Printf("  time:            %s\n", schedTime)
		} else {
			fmt.Printf("  time:            %s\n", utils.StyleInfo("not set"))
		}
		fmt.Println()

		// Build settings
		fmt.Println(utils.StyleTitle("Build Configuration:"))
		fmt.Printf("  ncpus:           %d\n", viper.GetInt("build.ncpus"))
		fmt.Printf("  mem_mb:          %d\n", viper.GetInt64("build.mem_mb"))
		fmt.Printf("  time:            %s\n", viper.GetString("build.time"))
		fmt.Printf("  tmp_size_mb:     %d\n", viper.GetInt("build.tmp_size_mb"))
		// Show actual compress_args (may be auto-detected based on apptainer version)
		compressArgs := viper.GetString("build.compress_args")
		actualCompressArgs := config.Global.Build.CompressArgs
		if compressArgs != actualCompressArgs {
			fmt.Printf("  compress_args:   %s\n", actualCompressArgs)
		} else {
			fmt.Printf("  compress_args:   %s\n", compressArgs)
		}
		fmt.Printf("  overlay_type:    %s\n", viper.GetString("build.overlay_type"))
		fmt.Println()

		// Extra base directories
		fmt.Println(utils.StyleTitle("Extra Base Directories:"))
		extraDirs := config.GetExtraBaseDirs()
		if len(extraDirs) > 0 {
			fmt.Printf("  extra_base_dirs:\n")
			for _, dir := range extraDirs {
				fmt.Printf("    - %s\n", dir)
			}
		} else {
			fmt.Printf("  extra_base_dirs: %s\n", utils.StyleInfo("none"))
		}
		fmt.Println()

		// Show environment variable overrides
		fmt.Println(utils.StyleTitle("Environment Variable Overrides:"))
		envVars := []string{
			"CONDATAINER_LOGS_DIR",
			"CONDATAINER_APPTAINER_BIN",
			"CONDATAINER_SCHEDULER_BIN",
			"CONDATAINER_SUBMIT_JOB",
			"CONDATAINER_SCHEDULER_NODES",
			"CONDATAINER_SCHEDULER_TASKS_PER_NODE",
			"CONDATAINER_SCHEDULER_NCPUS_PER_TASK",
			"CONDATAINER_SCHEDULER_MEM_MB_PER_NODE",
			"CONDATAINER_SCHEDULER_TIME",
			"CONDATAINER_BUILD_NCPUS",
			"CONDATAINER_BUILD_MEM_MB",
			"CONDATAINER_BUILD_TIME",
		}
		hasEnvOverrides := false
		for _, envVar := range envVars {
			if val := os.Getenv(envVar); val != "" {
				fmt.Printf("  %s=%s\n", envVar, val)
				hasEnvOverrides = true
			}
		}
		if !hasEnvOverrides {
			fmt.Printf("  %s\n", utils.StyleInfo("none"))
		}
	},
}

var configGetCmd = &cobra.Command{
	Use:   "get <key>",
	Short: "Get a configuration value",
	Long:  `Get a specific configuration value.`,
	Example: `  condatainer config get apptainer_bin
  condatainer config get scheduler.ncpus_per_task
  condatainer config get build.ncpus
  condatainer config get submit_job`,
	Args:              cobra.ExactArgs(1),
	ValidArgsFunction: configKeysCompletion,
	Run: func(cmd *cobra.Command, args []string) {
		key := args[0]
		value := viper.Get(key)
		if value == nil {
			ExitWithError("Unknown config key: %s", key)
		}
		fmt.Println(value)
	},
}

var configSetCmd = &cobra.Command{
	Use:   "set <key> <value>",
	Short: "Set a configuration value",
	Long: `Set a configuration value and save to config file.

Time duration format (for build.time):
  Go style:  2h, 30m, 1h30m, 90s
  HPC style: 02:00:00, 2:30:00, 1:30 (HH:MM:SS or HH:MM)`,
	Example: `  condatainer config set apptainer_bin /usr/bin/apptainer
  condatainer config set scheduler.ncpus_per_task 8
  condatainer config set scheduler.time 4h
  condatainer config set build.ncpus 8
  condatainer config set build.time 02:00:00
  condatainer config set submit_job false`,
	Args:              cobra.ExactArgs(2),
	ValidArgsFunction: configKeysCompletion,
	Run: func(cmd *cobra.Command, args []string) {
		key := args[0]
		value := args[1]

		// Validate known keys
		knownKeys := map[string]bool{
			"logs_dir":                  true,
			"apptainer_bin":             true,
			"scheduler_bin":             true,
			"submit_job":                true,
			"branch":                    true,
			"prefer_remote":             true,
			"scheduler.nodes":           true,
			"scheduler.tasks_per_node":  true,
			"scheduler.ncpus_per_task":  true,
			"scheduler.mem_mb_per_node": true,
			"scheduler.time":            true,
			"build.ncpus":               true,
			"build.mem_mb":              true,
			"build.time":                true,
			"build.tmp_size_mb":         true,
			"build.compress_args":       true,
			"build.overlay_type":        true,
			"extra_base_dirs":           true,
		}

		// Array keys that should be edited via config edit or env var
		if key == "extra_base_dirs" {
			utils.PrintError("'%s' is an array setting. Use 'condatainer config edit' or environment variable.", key)
			utils.PrintHint("Config file (YAML array):\n  extra_base_dirs:\n    - /path/to/dir1\n    - /path/to/dir2\n\nEnvironment variable (colon-separated):\n  export CONDATAINER_EXTRA_BASE_DIRS=/path/to/dir1:/path/to/dir2")
			os.Exit(ExitCodeError)
		}

		if !knownKeys[key] {
			utils.PrintWarning("Warning: '%s' is not a standard config key", key)
		}

		// Validate value based on key type
		if key == "scheduler.time" || key == "build.time" {
			if _, err := utils.ParseDuration(value); err != nil {
				utils.PrintError("Invalid duration format: %s", value)
				utils.PrintHint("Use format like: 2h, 30m, 1h30m, or 02:00:00")
				os.Exit(ExitCodeError)
			}
		}

		if key == "branch" {
			lower := strings.ToLower(value)
			if lower != "main" && lower != "dev" {
				utils.PrintWarning("Unknown branch '%s', falling back to 'main'", value)
				value = "main"
			} else {
				value = lower
			}
		}

		// Set the value
		viper.Set(key, value)

		// Get the config path that will be updated
		configPath, err := config.GetActiveOrUserConfigPath()
		if err != nil {
			ExitWithError("Failed to get config path: %v", err)
		}

		// Save to config file
		if err := config.SaveConfig(); err != nil {
			ExitWithError("Failed to save config: %v", err)
		}

		utils.PrintSuccess("Set %s = %s", utils.StyleInfo(key), utils.StyleInfo(value))
		utils.PrintNote("Config saved to: %s", configPath)
	},
}

var configInitCmd = &cobra.Command{
	Use:   "init",
	Short: "Create a config file with defaults",
	Long: `Create a configuration file with default values and auto-detected settings.

Config location options (-l, --location):
  u, user      ~/.config/condatainer/config.yaml (default for standard installations)
  p, portable  <install-dir>/config.yaml (default if installed in a dedicated folder)
  s, system    /etc/condatainer/config.yaml (requires appropriate permissions)

By default, the location is chosen based on the installation:
  - If the executable is in a dedicated folder (e.g., /apps/condatainer/bin/),
    the config is created in the parent folder (portable mode).
  - Otherwise, it uses the user's config directory.`,
	Example: `  condatainer config init              # Create config with smart default location
  condatainer config init -l user      # Force user config location
  condatainer config init -l portable  # Force portable config location`,
	Run: func(cmd *cobra.Command, args []string) {
		var configPath string
		var err error
		var locationType string

		// Determine config path based on --location flag or smart default
		if initLocation != "" {
			// User specified a location
			configPath, err = config.GetConfigPathByLocation(initLocation)
			if err != nil {
				ExitWithError("Invalid location: %v", err)
			}
			// Normalize shortcut to full name for display
			switch initLocation {
			case "u":
				locationType = "user"
			case "p":
				locationType = "portable"
			case "s":
				locationType = "system"
			default:
				locationType = initLocation
			}
		} else {
			// Smart default: portable if in dedicated installation, otherwise user
			if config.IsPortableInstallation() {
				configPath = config.GetPortableConfigPath()
				locationType = "portable"
			} else {
				configPath, err = config.GetUserConfigPath()
				if err != nil {
					ExitWithError("Failed to get config path: %v", err)
				}
				locationType = "user"
			}
		}

		// Check if config already exists
		if _, err := os.Stat(configPath); err == nil {
			utils.PrintWarning("Config file already exists: %s", configPath)
			shouldOverwrite := false
			if utils.ShouldAnswerYes() {
				shouldOverwrite = true
			} else {
				fmt.Print("Overwrite? [y/N]: ")
				var response string
				fmt.Scanln(&response)
				response = strings.ToLower(strings.TrimSpace(response))
				shouldOverwrite = (response == "y" || response == "yes")
			}
			if !shouldOverwrite {
				utils.PrintNote("Cancelled")
				return
			}
		}

		// Warn if inside container
		if config.IsInsideContainer() {
			utils.PrintWarning("Running inside a container - config changes may not persist on the host")
		}

		// Force re-detect binaries from current environment and save to specified path
		updated, err := config.ForceDetectAndSaveTo(configPath)
		if err != nil {
			ExitWithError("Failed to save config: %v", err)
		}

		// Check if apptainer was found
		apptainerBin := viper.GetString("apptainer_bin")

		// Auto-detect compression based on apptainer version and save to config
		detectedCompression := "-comp lz4" // default
		if apptainerBin != "" && config.ValidateBinary(apptainerBin) {
			if err := apptainer.SetBin(apptainerBin); err == nil {
				if version, err := apptainer.GetVersion(); err == nil {
					if apptainer.CheckZstdSupport(version) {
						detectedCompression = "-comp zstd -Xcompression-level 8"
					}
				}
			}
		}
		// Save detected compression to config
		if err := config.SetCompressArgsInConfigTo(detectedCompression, configPath); err != nil {
			utils.PrintWarning("Failed to save compression setting: %v", err)
		}

		if updated {
			utils.PrintSuccess("Config file created with auto-detected settings")
		} else {
			utils.PrintSuccess("Config file created")
		}
		fmt.Printf("  Location: %s (%s)\n", utils.StylePath(configPath), locationType)

		// Show what was detected
		fmt.Println()
		fmt.Println(utils.StyleTitle("Detected settings:"))
		fmt.Printf("  Apptainer: %s\n", viper.GetString("apptainer_bin"))
		if schedulerBin := viper.GetString("scheduler_bin"); schedulerBin != "" {
			fmt.Printf("  Scheduler: %s (%s)\n", schedulerBin, config.GetSchedulerTypeFromBin(schedulerBin))
		} else {
			fmt.Printf("  Scheduler: %s\n", utils.StyleWarning("not found"))
		}
		fmt.Printf("  Compression: %s\n", detectedCompression)
	},
}

var configEditCmd = &cobra.Command{
	Use:   "edit",
	Short: "Edit config file in default editor",
	Long:  "Open the configuration file in your default text editor ($EDITOR)",
	Run: func(cmd *cobra.Command, args []string) {
		configPath, err := config.GetActiveOrUserConfigPath()
		if err != nil {
			ExitWithError("Failed to get config path: %v", err)
		}

		// Create config if it doesn't exist
		if _, err := os.Stat(configPath); os.IsNotExist(err) {
			utils.PrintNote("Config file doesn't exist, creating it first...")
			if err := config.SaveConfig(); err != nil {
				ExitWithError("Failed to create config: %v", err)
			}
		}

		// Get editor from environment
		editor := os.Getenv("EDITOR")
		if editor == "" {
			editor = "vi" // fallback to vi
		}

		// Open editor
		editorCmd := exec.Command(editor, configPath)
		editorCmd.Stdin = os.Stdin
		editorCmd.Stdout = os.Stdout
		editorCmd.Stderr = os.Stderr

		if err := editorCmd.Run(); err != nil {
			ExitWithError("Failed to open editor: %v", err)
		}
	},
}

var configPathsCmd = &cobra.Command{
	Use:   "paths",
	Short: "Show data search paths",
	Long: `Show all search paths for images, build-scripts, and helper-scripts.

Search paths are checked in priority order (first match wins for reads):
  1. Extra base directories (highest priority)
  2. Portable (group/shared directory)
  3. Scratch (user HPC large storage, $SCRATCH/condatainer)
  4. User XDG directory (~/.local/share/condatainer)

Write operations use the first writable directory in the same order.
Portable is preferred for group/shared use, user is the fallback.`,
	Run: func(cmd *cobra.Command, args []string) {
		fmt.Println(utils.StyleTitle("Data Search Paths"))
		fmt.Println()

		// Images
		fmt.Println(utils.StyleTitle("Images:"))
		for i, dir := range config.GetImageSearchPaths() {
			status := ""
			if !config.DirExists(dir) {
				status = " " + utils.StyleWarning("(not found)")
			}
			fmt.Printf("  %d. %s%s\n", i+1, dir, status)
		}
		if writableDir, err := config.GetWritableImagesDir(); err == nil {
			fmt.Printf("  → Writable: %s\n", utils.StyleSuccess(writableDir))
		}
		fmt.Println()

		// Build scripts
		fmt.Println(utils.StyleTitle("Build Scripts:"))
		for i, dir := range config.GetBuildScriptSearchPaths() {
			status := ""
			if !config.DirExists(dir) {
				status = " " + utils.StyleWarning("(not found)")
			}
			fmt.Printf("  %d. %s%s\n", i+1, dir, status)
		}
		if writableDir, err := config.GetWritableBuildScriptsDir(); err == nil {
			fmt.Printf("  → Writable: %s\n", utils.StyleSuccess(writableDir))
		}
		fmt.Println()

		// Helper scripts
		fmt.Println(utils.StyleTitle("Helper Scripts:"))
		for i, dir := range config.GetHelperScriptSearchPaths() {
			status := ""
			if !config.DirExists(dir) {
				status = " " + utils.StyleWarning("(not found)")
			}
			fmt.Printf("  %d. %s%s\n", i+1, dir, status)
		}
		if writableDir, err := config.GetWritableHelperScriptsDir(); err == nil {
			fmt.Printf("  → Writable: %s\n", utils.StyleSuccess(writableDir))
		}
		fmt.Println()

		// Directory info
		fmt.Println(utils.StyleTitle("Directory Locations:"))
		if portableDir := config.GetPortableDataDir(); portableDir != "" {
			fmt.Printf("  Portable: %s\n", portableDir)
		} else {
			fmt.Printf("  Portable: %s\n", utils.StyleWarning("not detected"))
		}
		if scratchDir := config.GetScratchDataDir(); scratchDir != "" {
			fmt.Printf("  Scratch:  %s\n", scratchDir)
		} else {
			fmt.Printf("  Scratch:  %s\n", utils.StyleWarning("$SCRATCH not set"))
		}
		if userDir := config.GetUserDataDir(); userDir != "" {
			fmt.Printf("  User XDG: %s\n", userDir)
		}
	},
}

var configValidateCmd = &cobra.Command{
	Use:   "validate",
	Short: "Validate configuration",
	Long:  "Check if the current configuration is valid and all binaries are accessible",
	Run: func(cmd *cobra.Command, args []string) {
		valid := true

		if !utils.QuietMode {
			fmt.Println(utils.StyleTitle("Validating configuration..."))
			fmt.Println()
		}

		// Check apptainer binary
		apptainerBin := viper.GetString("apptainer_bin")
		if config.ValidateBinary(apptainerBin) {
			if !utils.QuietMode {
				fmt.Printf("%s Apptainer binary: %s\n", utils.StyleSuccess("✓"), apptainerBin)
			}
		} else {
			fmt.Printf("%s Apptainer binary not found: %s\n", utils.StyleError("✗"), apptainerBin)
			valid = false
		}

		// Check scheduler binary
		schedulerBin := viper.GetString("scheduler_bin")
		if schedulerBin != "" {
			if config.ValidateBinary(schedulerBin) {
				if !utils.QuietMode {
					fmt.Printf("%s Scheduler binary: %s\n", utils.StyleSuccess("✓"), schedulerBin)
				}
			} else {
				fmt.Printf("%s Scheduler binary not found: %s\n", utils.StyleError("✗"), schedulerBin)
				valid = false
			}
		} else {
			if !utils.QuietMode {
				fmt.Printf("%s Scheduler binary: %s\n", utils.StyleWarning("⚠"), "not configured")
			}
		}

		// Check build config
		ncpus := viper.GetInt("build.ncpus")
		if ncpus > 0 {
			if !utils.QuietMode {
				fmt.Printf("%s Build CPUs: %d\n", utils.StyleSuccess("✓"), ncpus)
			}
		} else {
			fmt.Printf("%s Build CPUs must be > 0: %d\n", utils.StyleError("✗"), ncpus)
			valid = false
		}

		memMB := viper.GetInt64("build.mem_mb")
		if memMB > 0 {
			if !utils.QuietMode {
				fmt.Printf("%s Build Memory: %d MB\n", utils.StyleSuccess("✓"), memMB)
			}
		} else {
			fmt.Printf("%s Build Memory must be > 0: %d\n", utils.StyleError("✗"), memMB)
			valid = false
		}

		if !utils.QuietMode {
			fmt.Println()
		}
		if valid {
			if !utils.QuietMode {
				utils.PrintSuccess("Configuration is valid")
			}
		} else {
			ExitWithError("Configuration has errors")
		}
	},
}

func init() {
	// Add flags
	configShowCmd.Flags().BoolVar(&showPath, "path", false, "Show only the config file path")
	configInitCmd.Flags().StringVarP(&initLocation, "location", "l", "", "Config location: user (u), portable (p), or system (s)")

	// Add subcommands
	configCmd.AddCommand(configShowCmd)
	configCmd.AddCommand(configGetCmd)
	configCmd.AddCommand(configSetCmd)
	configCmd.AddCommand(configInitCmd)
	configCmd.AddCommand(configEditCmd)
	configCmd.AddCommand(configPathsCmd)
	configCmd.AddCommand(configValidateCmd)

	// Add to root command
	rootCmd.AddCommand(configCmd)
}
