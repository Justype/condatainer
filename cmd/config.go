package cmd

import (
	"fmt"
	"os"
	"os/exec"
	"strings"

	"github.com/Justype/condatainer/internal/apptainer"
	"github.com/Justype/condatainer/internal/config"
	"github.com/Justype/condatainer/internal/utils"
	"github.com/spf13/cobra"
	"github.com/spf13/viper"
)

var (
	showPath bool
)

var configCmd = &cobra.Command{
	Use:   "config",
	Short: "Manage condatainer configuration",
	Long: `Manage condatainer configuration settings.

Configuration priority (highest to lowest):
  1. Command-line flags
  2. Environment variables (CONDATAINER_*)
  3. User config file (~/.config/condatainer/config.yaml)
  4. System config file (/etc/condatainer/config.yaml)
  5. Defaults`,
}

var configShowCmd = &cobra.Command{
	Use:   "show",
	Short: "Show current configuration",
	Long:  "Display all current configuration values and their sources",
	Run: func(cmd *cobra.Command, args []string) {
		if showPath {
			configPath, err := config.GetUserConfigPath()
			if err != nil {
				utils.PrintError("Failed to get config path: %v", err)
				os.Exit(1)
			}
			fmt.Println(configPath)
			return
		}

		// Show config file location
		configPath, _ := config.GetUserConfigPath()
		fmt.Printf("Config file: %s\n", utils.StylePath(configPath))
		if _, err := os.Stat(configPath); os.IsNotExist(err) {
			fmt.Printf("  %s (use 'condatainer config init' to create)\n\n", utils.StyleWarning("not found"))
		} else {
			fmt.Printf("  %s\n\n", utils.StyleSuccess("exists"))
		}

		// Show all settings
		fmt.Println(utils.StyleTitle("Current Configuration:"))
		fmt.Println()

		// Directories
		fmt.Println(utils.StyleTitle("Directories:"))
		baseDir := viper.GetString("base_dir")
		if baseDir == "" {
			baseDir = config.Global.BaseDir + " (auto-detected)"
		}
		fmt.Printf("  base_dir:       %s\n", baseDir)
		logsDir := viper.GetString("logs_dir")
		if logsDir == "" {
			logsDir = config.Global.LogsDir + " (default)"
		}
		fmt.Printf("  logs_dir:       %s\n", logsDir)
		fmt.Println()

		// Binary paths
		fmt.Println(utils.StyleTitle("Binaries:"))
		fmt.Printf("  apptainer_bin:  %s\n", viper.GetString("apptainer_bin"))
		fmt.Printf("  scheduler_bin:  %s\n", viper.GetString("scheduler_bin"))
		fmt.Printf("  scheduler_type: %s\n", viper.GetString("scheduler_type"))
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
		fmt.Println()

		// Build settings
		fmt.Println(utils.StyleTitle("Build Configuration:"))
		fmt.Printf("  default_cpus:    %d\n", viper.GetInt("build.default_cpus"))
		fmt.Printf("  default_mem_mb:  %d\n", viper.GetInt64("build.default_mem_mb"))
		fmt.Printf("  default_time:    %s\n", viper.GetString("build.default_time"))
		fmt.Printf("  tmp_size_mb:     %d\n", viper.GetInt("build.tmp_size_mb"))
		// Show actual compress_args (may be auto-detected based on apptainer version)
		compressArgs := viper.GetString("build.compress_args")
		actualCompressArgs := config.Global.Build.CompressArgs
		if compressArgs != actualCompressArgs {
			fmt.Printf("  compress_args:   %s (auto-detected: zstd supported)\n", actualCompressArgs)
		} else {
			fmt.Printf("  compress_args:   %s\n", compressArgs)
		}
		fmt.Printf("  overlay_type:    %s\n", viper.GetString("build.overlay_type"))
		fmt.Println()

		// Show environment variable overrides
		fmt.Println(utils.StyleTitle("Environment Variable Overrides:"))
		envVars := []string{
			"CONDATAINER_BASE_DIR",
			"CONDATAINER_LOGS_DIR",
			"CONDATAINER_APPTAINER_BIN",
			"CONDATAINER_SCHEDULER_BIN",
			"CONDATAINER_SUBMIT_JOB",
			"CONDATAINER_BUILD_DEFAULT_CPUS",
			"CONDATAINER_BUILD_DEFAULT_MEM_MB",
			"CONDATAINER_BUILD_DEFAULT_TIME",
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
	Long: `Get a specific configuration value.

Examples:
  condatainer config get apptainer_bin
  condatainer config get build.default_cpus
  condatainer config get submit_job`,
	Args: cobra.ExactArgs(1),
	Run: func(cmd *cobra.Command, args []string) {
		key := args[0]
		value := viper.Get(key)
		if value == nil {
			utils.PrintError("Unknown config key: %s", key)
			os.Exit(1)
		}
		fmt.Println(value)
	},
}

var configSetCmd = &cobra.Command{
	Use:   "set <key> <value>",
	Short: "Set a configuration value",
	Long: `Set a configuration value and save to config file.

Examples:
  condatainer config set apptainer_bin /usr/bin/apptainer
  condatainer config set build.default_cpus 8
  condatainer config set build.default_time 4h
  condatainer config set build.default_time 02:00:00
  condatainer config set submit_job false

Time duration format (for build.default_time):
  Go style:  2h, 30m, 1h30m, 90s
  HPC style: 02:00:00, 2:30:00, 1:30 (HH:MM:SS or HH:MM)`,
	Args: cobra.ExactArgs(2),
	Run: func(cmd *cobra.Command, args []string) {
		key := args[0]
		value := args[1]

		// Validate known keys
		knownKeys := map[string]bool{
			"base_dir":             true,
			"logs_dir":             true,
			"apptainer_bin":        true,
			"scheduler_bin":        true,
			"scheduler_type":       true,
			"submit_job":           true,
			"build.default_cpus":   true,
			"build.default_mem_mb": true,
			"build.default_time":   true,
			"build.tmp_size_mb":    true,
			"build.compress_args":  true,
			"build.overlay_type":   true,
		}

		if !knownKeys[key] {
			utils.PrintWarning("Warning: '%s' is not a standard config key", key)
		}

		// Validate value based on key type
		if key == "build.default_time" {
			if _, err := utils.ParseDuration(value); err != nil {
				utils.PrintError("Invalid duration format: %s", value)
				utils.PrintHint("Use format like: 2h, 30m, 1h30m, or 02:00:00")
				os.Exit(1)
			}
		}

		// Set the value
		viper.Set(key, value)

		// Save to config file
		if err := config.SaveConfig(); err != nil {
			utils.PrintError("Failed to save config: %v", err)
			os.Exit(1)
		}

		configPath, _ := config.GetUserConfigPath()
		utils.PrintSuccess("Set %s = %s", utils.StyleInfo(key), utils.StyleInfo(value))
		utils.PrintNote("Config saved to: %s", configPath)
	},
}

var configInitCmd = &cobra.Command{
	Use:   "init",
	Short: "Create a config file with defaults",
	Long:  "Create a configuration file with default values and auto-detected settings",
	Run: func(cmd *cobra.Command, args []string) {
		configPath, err := config.GetUserConfigPath()
		if err != nil {
			utils.PrintError("Failed to get config path: %v", err)
			os.Exit(1)
		}

		// Check if config already exists
		if _, err := os.Stat(configPath); err == nil {
			utils.PrintWarning("Config file already exists: %s", configPath)
			fmt.Print("Overwrite? [y/N]: ")
			var response string
			fmt.Scanln(&response)
			response = strings.ToLower(strings.TrimSpace(response))
			if response != "y" && response != "yes" {
				utils.PrintNote("Cancelled")
				return
			}
		}

		// Warn if inside container
		if config.IsInsideContainer() {
			utils.PrintWarning("Running inside a container - config changes may not persist on the host")
		}

		// Force re-detect binaries from current environment
		updated, err := config.ForceDetectAndSave()
		if err != nil {
			utils.PrintError("Failed to save config: %v", err)
			os.Exit(1)
		}

		// Check if apptainer was found
		apptainerBin := viper.GetString("apptainer_bin")
		if apptainerBin == "" || !config.ValidateBinary(apptainerBin) {
			utils.PrintWarning("Cannot find apptainer. Please load apptainer module:")
			utils.PrintHint("ml apptainer")
		}

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
		if err := config.SetCompressArgsInConfig(detectedCompression); err != nil {
			utils.PrintWarning("Failed to save compression setting: %v", err)
		}

		if updated {
			utils.PrintSuccess("Config file created with auto-detected settings: %s", configPath)
		} else {
			utils.PrintSuccess("Config file created: %s", configPath)
		}

		// Show what was detected
		fmt.Println()
		fmt.Println(utils.StyleTitle("Detected settings:"))
		fmt.Printf("  Apptainer: %s\n", viper.GetString("apptainer_bin"))
		if schedulerBin := viper.GetString("scheduler_bin"); schedulerBin != "" {
			fmt.Printf("  Scheduler: %s (%s)\n", schedulerBin, viper.GetString("scheduler_type"))
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
		configPath, err := config.GetUserConfigPath()
		if err != nil {
			utils.PrintError("Failed to get config path: %v", err)
			os.Exit(1)
		}

		// Create config if it doesn't exist
		if _, err := os.Stat(configPath); os.IsNotExist(err) {
			utils.PrintNote("Config file doesn't exist, creating it first...")
			if err := config.SaveConfig(); err != nil {
				utils.PrintError("Failed to create config: %v", err)
				os.Exit(1)
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
			utils.PrintError("Failed to open editor: %v", err)
			os.Exit(1)
		}
	},
}

var configValidateCmd = &cobra.Command{
	Use:   "validate",
	Short: "Validate configuration",
	Long:  "Check if the current configuration is valid and all binaries are accessible",
	Run: func(cmd *cobra.Command, args []string) {
		valid := true

		fmt.Println(utils.StyleTitle("Validating configuration..."))
		fmt.Println()

		// Check apptainer binary
		apptainerBin := viper.GetString("apptainer_bin")
		if config.ValidateBinary(apptainerBin) {
			fmt.Printf("%s Apptainer binary: %s\n", utils.StyleSuccess("✓"), apptainerBin)
		} else {
			fmt.Printf("%s Apptainer binary not found: %s\n", utils.StyleError("✗"), apptainerBin)
			valid = false
		}

		// Check scheduler binary
		schedulerBin := viper.GetString("scheduler_bin")
		if schedulerBin != "" {
			if config.ValidateBinary(schedulerBin) {
				fmt.Printf("%s Scheduler binary: %s\n", utils.StyleSuccess("✓"), schedulerBin)
			} else {
				fmt.Printf("%s Scheduler binary not found: %s\n", utils.StyleError("✗"), schedulerBin)
				valid = false
			}
		} else {
			fmt.Printf("%s Scheduler binary: %s\n", utils.StyleWarning("⚠"), "not configured")
		}

		// Check build config
		defaultCPUs := viper.GetInt("build.default_cpus")
		if defaultCPUs > 0 {
			fmt.Printf("%s Build CPUs: %d\n", utils.StyleSuccess("✓"), defaultCPUs)
		} else {
			fmt.Printf("%s Build CPUs must be > 0: %d\n", utils.StyleError("✗"), defaultCPUs)
			valid = false
		}

		defaultMemMB := viper.GetInt64("build.default_mem_mb")
		if defaultMemMB > 0 {
			fmt.Printf("%s Build Memory: %d MB\n", utils.StyleSuccess("✓"), defaultMemMB)
		} else {
			fmt.Printf("%s Build Memory must be > 0: %d\n", utils.StyleError("✗"), defaultMemMB)
			valid = false
		}

		fmt.Println()
		if valid {
			utils.PrintSuccess("Configuration is valid")
		} else {
			utils.PrintError("Configuration has errors")
			os.Exit(1)
		}
	},
}

func init() {
	// Add flags
	configShowCmd.Flags().BoolVar(&showPath, "path", false, "Show only the config file path")

	// Add subcommands
	configCmd.AddCommand(configShowCmd)
	configCmd.AddCommand(configGetCmd)
	configCmd.AddCommand(configSetCmd)
	configCmd.AddCommand(configInitCmd)
	configCmd.AddCommand(configEditCmd)
	configCmd.AddCommand(configValidateCmd)

	// Add to root command
	rootCmd.AddCommand(configCmd)
}
