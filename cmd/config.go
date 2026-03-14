package cmd

import (
	"fmt"
	"os"
	"os/exec"
	"sort"
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
	"default_distro",
	"submit_job",
	"scripts_link",
	"extra_scripts_links",
	"prefer_remote",
	"extra_base_dirs",
	"parse_module_load",
	"scheduler_timeout",
	"metadata_cache_ttl",
	"build.ncpus",
	"build.mem_mb",
	"build.time",
	"build.compress_args",
	"build.overlay_type",
	"build.block_size",
	"build.data_block_size",
	"build.use_tmp_overlay",
	"build.tmp_overlay_size_mb",
}

// arrayConfigKeys is the set of config keys that hold string slices.
var arrayConfigKeys = []string{"extra_base_dirs", "extra_scripts_links"}

func isArrayKey(key string) bool {
	for _, k := range arrayConfigKeys {
		if k == key {
			return true
		}
	}
	return false
}

// modifyArrayConfig reads the current slice for key, applies modify, writes back.
func modifyArrayConfig(key string, modify func([]string) []string) error {
	if !isArrayKey(key) {
		return fmt.Errorf("'%s' is not an array config key; array keys: %s",
			key, strings.Join(arrayConfigKeys, ", "))
	}
	current := viper.GetStringSlice(key)
	updated := modify(current)
	viper.Set(key, updated)
	return config.SaveConfig()
}

func arrayKeyCompletion(cmd *cobra.Command, args []string, toComplete string) ([]string, cobra.ShellCompDirective) {
	if len(args) == 0 {
		return arrayConfigKeys, cobra.ShellCompDirectiveNoFileComp
	}
	if len(args) == 1 {
		switch args[0] {
		case "extra_base_dirs":
			return nil, cobra.ShellCompDirectiveFilterDirs
		case "extra_scripts_links":
			return nil, cobra.ShellCompDirectiveDefault
		}
	}
	return nil, cobra.ShellCompDirectiveNoFileComp
}

func arrayRemoveValueCompletion(cmd *cobra.Command, args []string, toComplete string) ([]string, cobra.ShellCompDirective) {
	if len(args) == 1 && isArrayKey(args[0]) {
		return viper.GetStringSlice(args[0]), cobra.ShellCompDirectiveNoFileComp
	}
	return arrayKeyCompletion(cmd, args, toComplete)
}

// getConfigEnvVars returns a list of environment variables corresponding to
// the supported configuration keys.  We base the list on the known config
// keys (configKeys) so that the envelope is automatically updated when new
// settings are added.  The variable names are simply the uppercased key with
// dots replaced by underscores and a CNT_ prefix.
func getConfigEnvVars() []string {
	vars := make([]string, 0, len(configKeys))
	for _, key := range configKeys {
		env := "CNT_" + strings.ToUpper(strings.ReplaceAll(key, ".", "_"))
		vars = append(vars, env)
	}
	// sort for deterministic output
	sort.Strings(vars)
	return vars
}

// configKeysCompletion returns config keys for shell completion
func configKeysCompletion(cmd *cobra.Command, args []string, toComplete string) ([]string, cobra.ShellCompDirective) {
	if len(args) == 0 {
		// First arg: complete config keys
		return configKeys, cobra.ShellCompDirectiveNoFileComp
	}
	if len(args) == 1 {
		// Second arg: complete values based on the key
		//  extra_base_dirs should only complete directories
		if args[0] == "extra_base_dirs" {
			return nil, cobra.ShellCompDirectiveFilterDirs
		}
		// extra_scripts_links is an array setting; allow free-form input
		if args[0] == "extra_scripts_links" {
			return nil, cobra.ShellCompDirectiveDefault
		}
		return configValueCompletion(args[0]), cobra.ShellCompDirectiveNoFileComp
	}
	return nil, cobra.ShellCompDirectiveNoFileComp
}

// configValueCompletion returns suggested values for a config key
func configValueCompletion(key string) []string {
	switch key {
	case "submit_job":
		return []string{"true", "false"}
	case "prefer_remote", "parse_module_load":
		return []string{"true", "false"}
	case "default_distro":
		return config.GetAvailableDistros()
	case "build.ncpus":
		return []string{"4", "8", "16", "32"}
	case "build.mem_mb":
		return []string{"4096", "8192", "16384", "32768"}
	case "build.time":
		return []string{"1h", "2h", "4h", "8h"}
	case "build.overlay_type":
		return []string{"ext3", "squashfs"}
	case "build.compress_args":
		return config.CompressNames()
	case "build.block_size", "build.data_block_size":
		return config.BlockSizeCompletions
	case "build.use_tmp_overlay":
		return []string{"true", "false"}
	case "build.tmp_overlay_size_mb":
		return []string{"10240", "20480", "40960"}
	case "metadata_cache_ttl":
		return []string{"1", "3", "7", "14", "0"}
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
  2. Environment variables (CNT_*)
  3. User config file (~/.config/condatainer/config.yaml)
  4. Portable config (<install-dir>/config.yaml, if in dedicated folder)
  5. System config file (/etc/condatainer/config.yaml)
  6. Defaults

Data directory priority (for read/write operations):
  1. Extra base directories (CNT_EXTRA_BASE_DIRS or config)
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

		pathIndex := 1

		// Extra base directories
		extraBaseDirs := config.GetExtraBaseDirs()
		for _, dir := range extraBaseDirs {
			status := ""
			if _, err := os.Stat(dir); err == nil {
				if utils.IsWritableDir(dir) {
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
				if utils.IsWritableDir(portableDir) {
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
				if utils.IsWritableDir(scratchDir) {
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
				if utils.IsWritableDir(userDir) {
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

		// Distros and default
		fmt.Println(utils.StyleTitle("Distros:"))
		fmt.Printf("  Available:      %s\n", strings.Join(config.GetAvailableDistros(), ", "))
		fmt.Printf("  default_distro: %s\n", config.Global.DefaultDistro)
		fmt.Println()

		// Remote sources
		fmt.Println(utils.StyleTitle("Remote Sources:"))
		fmt.Printf("  scripts_link:  %s\n", config.Global.ScriptsLink)
		extraScriptsLinks := config.GetExtraScriptsLinks()
		if len(extraScriptsLinks) > 0 {
			fmt.Printf("  extra_scripts_links:\n")
			for _, l := range extraScriptsLinks {
				fmt.Printf("    - %s\n", l)
			}
		} else {
			fmt.Printf("  extra_scripts_links: %s\n", utils.StyleInfo("none"))
		}
		fmt.Printf("  prefer_remote: %v\n", config.Global.PreferRemote)
		fmt.Println()

		// Options
		fmt.Println(utils.StyleTitle("Options:"))
		submitJobConfig := viper.GetBool("submit_job")
		submitJobActual := config.Global.SubmitJob
		if submitJobConfig && !submitJobActual {
			fmt.Printf("  submit_job:        %v (disabled: scheduler not accessible)\n", submitJobConfig)
		} else {
			fmt.Printf("  submit_job:        %v\n", submitJobActual)
		}
		if config.Global.SchedulerTimeout == 0 {
			fmt.Printf("  scheduler_timeout: 0 (disabled)\n")
		} else {
			fmt.Printf("  scheduler_timeout: %ds\n", int(config.Global.SchedulerTimeout.Seconds()))
		}
		fmt.Printf("  parse_module_load: %v\n", config.Global.ParseModuleLoad)
		if config.Global.MetadataCacheTTL == 0 {
			fmt.Printf("  metadata_cache_ttl: 0 (disabled)\n")
		} else {
			fmt.Printf("  metadata_cache_ttl: %dd\n", int(config.Global.MetadataCacheTTL.Hours()/24))
		}
		fmt.Println()

		// Scheduler default specs

		// Build settings
		fmt.Printf("%s %s\n", utils.StyleTitle("Build Configuration:"), "build.*")
		fmt.Printf("  ncpus:                %d\n", viper.GetInt("build.ncpus"))
		fmt.Printf("  mem_mb:               %d\n", viper.GetInt64("build.mem_mb"))
		fmt.Printf("  time:                 %s\n", viper.GetString("build.time"))
		// Show actual compress_args (may be auto-detected based on apptainer version)
		compressArgs := viper.GetString("build.compress_args")
		actualCompressArgs := config.Global.Build.CompressArgs
		if compressArgs != actualCompressArgs {
			fmt.Printf("  compress_args:        %s\n", actualCompressArgs)
		} else {
			fmt.Printf("  compress_args:        %s\n", compressArgs)
		}
		fmt.Printf("  overlay_type:         %s\n", viper.GetString("build.overlay_type"))
		fmt.Printf("  block_size:           %s\n", config.Global.Build.BlockSize)
		fmt.Printf("  data_block_size:      %s\n", config.Global.Build.DataBlockSize)
		fmt.Printf("  use_tmp_overlay:      %v\n", config.Global.Build.UseTmpOverlay)
		fmt.Printf("  tmp_overlay_size_mb:  %d\n", viper.GetInt("build.tmp_overlay_size_mb"))
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
		envVars := getConfigEnvVars()
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
  condatainer config get build.ncpus
  condatainer config get submit_job`,
	Args:              cobra.ExactArgs(1),
	ValidArgsFunction: configKeysCompletion,
	Run: func(cmd *cobra.Command, args []string) {
		key := args[0]
		value := viper.Get(key)
		if value == nil {
			ExitWithUsageError("Unknown config key: %s", key)
		}
		if isArrayKey(key) {
			for _, v := range viper.GetStringSlice(key) {
				fmt.Println(v)
			}
		} else {
			fmt.Println(value)
		}
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
  condatainer config set build.ncpus 8
  condatainer config set build.time 02:00:00
  condatainer config set submit_job false`,
	Args:              cobra.ExactArgs(2),
	ValidArgsFunction: configKeysCompletion,
	Run: func(cmd *cobra.Command, args []string) {
		key := args[0]
		value := args[1]

		// Build set of known keys from configKeys (single source of truth)
		knownKeys := make(map[string]bool, len(configKeys))
		for _, k := range configKeys {
			knownKeys[k] = true
		}

		// Validate default_distro against known distros
		if key == "default_distro" {
			if !config.IsValidDistro(value) {
				utils.PrintError("Unknown distro '%s'. Available: %s", value, strings.Join(config.GetAvailableDistros(), ", "))
				os.Exit(ExitCodeUsage)
			}
		}

		// Array keys require append/prepend/remove subcommands
		if isArrayKey(key) {
			utils.PrintError("'%s' is an array setting. Use append/prepend/remove subcommands.", key)
			utils.PrintHint("  condatainer config append  %s <value>\n  condatainer config prepend %s <value>\n  condatainer config remove  %s <value>", key, key, key)
			os.Exit(ExitCodeUsage)
		}

		if !knownKeys[key] {
			utils.PrintWarning("Warning: '%s' is not a standard config key", key)
		}

		// Validate value based on key type
		if key == "scheduler_timeout" {
			var n int
			if _, err := fmt.Sscan(value, &n); err != nil || n < 0 {
				utils.PrintError("Invalid value for scheduler_timeout: %s (must be a non-negative integer; 0 disables the timeout)", value)
				os.Exit(ExitCodeUsage)
			}
		}

		if key == "metadata_cache_ttl" {
			var n int
			if _, err := fmt.Sscan(value, &n); err != nil || n < 0 {
				utils.PrintError("Invalid value for metadata_cache_ttl: %s (must be a non-negative integer in days; 0 disables the cache)", value)
				os.Exit(ExitCodeUsage)
			}
		}

		if key == "build.time" {
			if _, err := utils.ParseWalltime(value); err != nil {
				utils.PrintError("Invalid duration format: %s", value)
				utils.PrintHint("Use format like: 4d12h, 2h30m, 1:30, or 01:30:00")
				os.Exit(ExitCodeUsage)
			}
		}

		if key == "build.compress_args" {
			value = config.NormalizeCompressArgs(value)
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
				response, readErr := utils.ReadLineContext(cmd.Context())
				shouldOverwrite = readErr == nil && (response == "y" || response == "yes")
			}
			if !shouldOverwrite {
				utils.PrintNote("Cancelled")
				return
			}
		}

		// Warn if inside container
		if config.IsInsideContainer() {
			ExitWithError("Cannot initialize config inside a container. Please run this command on the host system.")
		}

		// If neither apptainer nor singularity is in PATH, exit with an error
		if config.DetectApptainerBin() == "" {
			ExitWithError("Neither 'apptainer' nor 'singularity' binary found in PATH.")
		}

		// Force re-detect binaries from current environment and save to specified path
		updated, err := config.ForceDetectAndSaveTo(configPath)
		if err != nil {
			ExitWithError("Failed to save config: %v", err)
		}

		// Auto-detect compression based on binary type and version, then save to config.
		// Clear any previously set viper value so AutoDetectCompression always runs during init.
		detectedCompression := config.Global.Build.CompressArgs // fallback to runtime default
		apptainerBin := viper.GetString("apptainer_bin")
		if apptainerBin != "" && config.ValidateBinary(apptainerBin) {
			if err := apptainer.SetBin(apptainerBin); err == nil {
				if version, err := apptainer.GetVersion(); err == nil {
					viper.Set("build.compress_args", "")
					config.AutoDetectCompression(apptainer.CheckZstdSupport(version), apptainer.IsSingularity())
					detectedCompression = config.Global.Build.CompressArgs
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

		// Check default distro validity
		distro := config.Global.DefaultDistro
		if !config.IsValidDistro(distro) {
			fmt.Printf("%s Default distro invalid: %s\n", utils.StyleError("✗"), distro)
			valid = false
		} else {
			if !utils.QuietMode {
				fmt.Printf("%s Default distro: %s\n", utils.StyleSuccess("✓"), distro)
			}
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

var configAppendCmd = &cobra.Command{
	Use:               "append <key> <value>",
	Short:             "Append a value to an array config key",
	Example:           "  condatainer config append extra_base_dirs /scratch/shared\n  condatainer config append extra_scripts_links https://example.com/scripts",
	Args:              cobra.ExactArgs(2),
	ValidArgsFunction: arrayKeyCompletion,
	SilenceUsage:      true,
	Run: func(cmd *cobra.Command, args []string) {
		key, value := args[0], args[1]
		if err := modifyArrayConfig(key, func(cur []string) []string {
			return append(cur, value)
		}); err != nil {
			ExitWithError("%v", err)
		}
		configPath, _ := config.GetActiveOrUserConfigPath()
		utils.PrintSuccess("Appended %s to %s", utils.StyleInfo(value), utils.StyleInfo(key))
		utils.PrintNote("Config saved to: %s", configPath)
	},
}

var configPrependCmd = &cobra.Command{
	Use:               "prepend <key> <value>",
	Short:             "Prepend a value to an array config key (highest priority)",
	Example:           "  condatainer config prepend extra_base_dirs /project/common\n  condatainer config prepend extra_scripts_links https://myorg.com/scripts",
	Args:              cobra.ExactArgs(2),
	ValidArgsFunction: arrayKeyCompletion,
	SilenceUsage:      true,
	Run: func(cmd *cobra.Command, args []string) {
		key, value := args[0], args[1]
		if err := modifyArrayConfig(key, func(cur []string) []string {
			return append([]string{value}, cur...)
		}); err != nil {
			ExitWithError("%v", err)
		}
		configPath, _ := config.GetActiveOrUserConfigPath()
		utils.PrintSuccess("Prepended %s to %s", utils.StyleInfo(value), utils.StyleInfo(key))
		utils.PrintNote("Config saved to: %s", configPath)
	},
}

var configRemoveCmd = &cobra.Command{
	Use:               "remove <key> <value>",
	Short:             "Remove a value from an array config key",
	Example:           "  condatainer config remove extra_base_dirs /scratch/shared\n  condatainer config remove extra_scripts_links https://example.com/scripts",
	Args:              cobra.ExactArgs(2),
	ValidArgsFunction: arrayRemoveValueCompletion,
	SilenceUsage:      true,
	Run: func(cmd *cobra.Command, args []string) {
		key, value := args[0], args[1]
		removed := false
		if err := modifyArrayConfig(key, func(cur []string) []string {
			var out []string
			for _, v := range cur {
				if v == value {
					removed = true
					continue
				}
				out = append(out, v)
			}
			return out
		}); err != nil {
			ExitWithError("%v", err)
		}
		configPath, _ := config.GetActiveOrUserConfigPath()
		if removed {
			utils.PrintSuccess("Removed %s from %s", utils.StyleInfo(value), utils.StyleInfo(key))
			utils.PrintNote("Config saved to: %s", configPath)
		} else {
			utils.PrintWarning("%s not found in %s", utils.StyleInfo(value), utils.StyleInfo(key))
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
	configCmd.AddCommand(configAppendCmd)
	configCmd.AddCommand(configPrependCmd)
	configCmd.AddCommand(configRemoveCmd)
	configCmd.AddCommand(configInitCmd)
	configCmd.AddCommand(configEditCmd)
	configCmd.AddCommand(configPathsCmd)
	configCmd.AddCommand(configValidateCmd)

	// Add to root command
	rootCmd.AddCommand(configCmd)
}
