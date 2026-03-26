package cmd

import (
	"fmt"
	"os"
	"path/filepath"
	"slices"
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
	showSources  bool
	initLocation string // user, root, extra-root, system, or a custom path
	getLocation  string // target location for config get (single-layer read)
	setLocation  string // target location for config set/append/prepend/remove
)

// configKeyDefs maps every known config key to whether it holds a string slice (array).
// true = array key (use append/prepend/remove); false = scalar key (use set).
var configKeyDefs = map[string]bool{
	"logs_dir":               false,
	"apptainer_bin":          false,
	"scheduler_bin":          false,
	"default_distro":         false,
	"submit_job":             false,
	"scripts_link":           false,
	"extra_scripts_links":    true,
	"prefer_remote":          false,
	"extra_image_dirs":       true,
	"extra_build_dirs":       true,
	"extra_helper_dirs":      true,
	"parse_module_load":      false,
	"scheduler_timeout":      false,
	"metadata_cache_ttl":     false,
	"build.ncpus":            false,
	"build.mem":              false,
	"build.time":             false,
	"build.compress_args":    false,
	"build.block_size":       false,
	"build.data_block_size":  false,
	"build.use_tmp_overlay":  false,
	"build.tmp_overlay_size": false,
	"channels":               true,
}

func isArrayKey(key string) bool { return configKeyDefs[key] }

// modifyArrayConfig reads the current slice for key from the target config file,
// applies modify, and writes the result back. Returns the config path written.
func modifyArrayConfig(key string, modify func([]string) []string) (string, error) {
	if !isArrayKey(key) {
		var arrayKeys []string
		for k, isArr := range configKeyDefs {
			if isArr {
				arrayKeys = append(arrayKeys, k)
			}
		}
		sort.Strings(arrayKeys)
		return "", fmt.Errorf("'%s' is not an array config key; array keys: %s",
			key, strings.Join(arrayKeys, ", "))
	}
	configPath, _, err := config.ResolveWritableConfigPath(setLocation)
	if err != nil {
		return "", err
	}
	current := config.ReadConfigSliceKey(configPath, key)
	updated := modify(current)
	if len(updated) == 0 {
		return configPath, config.DeleteConfigKey(configPath, key)
	}
	return configPath, config.UpdateConfigKey(configPath, key, updated)
}

func arrayKeyCompletion(cmd *cobra.Command, args []string, toComplete string) ([]string, cobra.ShellCompDirective) {
	if len(args) == 0 {
		var keys []string
		for k, isArr := range configKeyDefs {
			if isArr {
				keys = append(keys, k)
			}
		}
		return keys, cobra.ShellCompDirectiveNoFileComp
	}
	if len(args) == 1 {
		switch args[0] {
		case "extra_scripts_links":
			return nil, cobra.ShellCompDirectiveDefault
		}
	}
	return nil, cobra.ShellCompDirectiveNoFileComp
}

func arrayRemoveValueCompletion(cmd *cobra.Command, args []string, toComplete string) ([]string, cobra.ShellCompDirective) {
	loc, _ := cmd.Flags().GetString("location")
	if len(args) == 0 {
		// Suggest only keys that are actually set in the target config file.
		if configPath, _, err := config.ResolveWritableConfigPath(loc); err == nil {
			v := config.ReadConfigAllKeys(configPath)
			if len(v) > 0 {
				sort.Strings(v)
				return v, cobra.ShellCompDirectiveNoFileComp
			}
		}
		// Fallback: all known keys
		var keys []string
		for k := range configKeyDefs {
			keys = append(keys, k)
		}
		sort.Strings(keys)
		return keys, cobra.ShellCompDirectiveNoFileComp
	}
	if len(args) == 1 && isArrayKey(args[0]) {
		// Suggest values present in the target config file for this array key.
		if configPath, _, err := config.ResolveWritableConfigPath(loc); err == nil {
			return config.ReadConfigSliceKey(configPath, args[0]), cobra.ShellCompDirectiveNoFileComp
		}
		return viper.GetStringSlice(args[0]), cobra.ShellCompDirectiveNoFileComp
	}
	return nil, cobra.ShellCompDirectiveNoFileComp
}

// getConfigEnvVars returns a sorted list of CNT_* environment variables for all known config keys.
func getConfigEnvVars() []string {
	vars := make([]string, 0, len(configKeyDefs))
	for key := range configKeyDefs {
		vars = append(vars, "CNT_"+strings.ToUpper(strings.ReplaceAll(key, ".", "_")))
	}
	sort.Strings(vars)
	return vars
}

// scalarConfigKeys returns all non-array config keys.
func scalarConfigKeys() []string {
	var out []string
	for k, isArr := range configKeyDefs {
		if !isArr {
			out = append(out, k)
		}
	}
	return out
}

// configSetKeysCompletion returns only scalar (non-array) config keys for `config set`.
func configSetKeysCompletion(cmd *cobra.Command, args []string, toComplete string) ([]string, cobra.ShellCompDirective) {
	if len(args) == 0 {
		return scalarConfigKeys(), cobra.ShellCompDirectiveNoFileComp
	}
	if len(args) == 1 {
		return configValueCompletion(args[0]), cobra.ShellCompDirectiveNoFileComp
	}
	return nil, cobra.ShellCompDirectiveNoFileComp
}

// configKeysCompletion returns config keys for shell completion
func configKeysCompletion(cmd *cobra.Command, args []string, toComplete string) ([]string, cobra.ShellCompDirective) {
	if len(args) == 0 {
		// First arg: complete all config keys
		var keys []string
		for k := range configKeyDefs {
			keys = append(keys, k)
		}
		return keys, cobra.ShellCompDirectiveNoFileComp
	}
	if len(args) == 1 {
		// Second arg: complete values based on the key
		//  extra_image_dirs / extra_build_dirs / extra_helper_dirs should only complete directories
		if args[0] == "extra_image_dirs" || args[0] == "extra_build_dirs" || args[0] == "extra_helper_dirs" {
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
	case "build.mem":
		return []string{"4g", "8g", "16g", "32g"}
	case "build.time":
		return []string{"1h", "2h", "4h", "8h"}
	case "build.compress_args":
		return config.CompressNames()
	case "build.block_size", "build.data_block_size":
		return config.BlockSizeCompletions
	case "build.use_tmp_overlay":
		return []string{"true", "false"}
	case "build.tmp_overlay_size":
		return []string{"10g", "20g", "40g"}
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
  4. Extra-root config ($CNT_EXTRA_ROOT/config.yaml, group/lab layer)
  5. App-root config (<install-dir>/config.yaml, if in dedicated folder or CNT_ROOT set)
  6. System config file (/etc/condatainer/config.yaml)
  7. Defaults

Data directory priority (for read/write operations):
  1. Extra base directories (CNT_EXTRA_BASE_DIRS or config)
  2. Extra-root directory ($CNT_EXTRA_ROOT, group/lab layer)
  3. App-root directory (auto-detected or CNT_ROOT)
  4. User scratch directory ($SCRATCH/condatainer, HPC systems)
  5. User XDG directory (~/.local/share/condatainer)`,
}

var configShowCmd = &cobra.Command{
	Use:   "show",
	Short: "Show current configuration",
	Long: `Display current configuration values and their sources.

Shows:
  - Config file search paths and which one is in use
  - Base directory search paths with writable status
  - All configuration settings (directories, binaries, build config)
  - Environment variable overrides

Use --sources (-s) to annotate each value with which config layer provides it.`,
	Run: func(cmd *cobra.Command, args []string) {
		if showPath {
			configPath, err := config.GetUserConfigPath()
			if err != nil {
				ExitWithError("Failed to get config path: %v", err)
			}
			fmt.Println(configPath)
			return
		}

		// --sources helpers: inline [layer] / [env] annotations
		layers := config.GetConfigLayerInfos()

		// srcTag returns "  [layer]" (green) or "  [env: VAR]" (yellow) for a scalar key.
		// Always returns "  [default]" (gray) if the key is not set in any layer.
		srcTag := func(key string) string {
			envVar := "CNT_" + strings.ToUpper(strings.ReplaceAll(key, ".", "_"))
			if os.Getenv(envVar) != "" {
				if showSources {
					return "  " + utils.StyleWarning("[env: "+envVar+"]")
				}
				return ""
			}
			for _, layer := range layers {
				if layer.InConfig(key) {
					if showSources {
						return "  " + utils.StyleSuccess("["+layer.Type+"]")
					}
					return ""
				}
			}
			return "  " + utils.StyleDebug("[default]")
		}

		// srcEntryTag returns "  [layer]" (green/yellow) for a value inside an array key.
		srcEntryTag := func(key, val string) string {
			if !showSources {
				return ""
			}
			envVar := "CNT_" + strings.ToUpper(strings.ReplaceAll(key, ".", "_"))
			if os.Getenv(envVar) != "" {
				return "  " + utils.StyleWarning("[env]")
			}
			for _, layer := range layers {
				if !layer.InConfig(key) {
					continue
				}
				for _, v := range layer.GetStringSlice(key) {
					if v == val {
						return "  " + utils.StyleSuccess("["+layer.Type+"]")
					}
				}
			}
			return ""
		}

		// printOverridden prints dimmed lines for lower-priority layers that also set key.
		// indent should pad to the same column as the active value on the line above.
		printOverridden := func(indent, key string) {
			if !showSources {
				return
			}
			envActive := os.Getenv("CNT_"+strings.ToUpper(strings.ReplaceAll(key, ".", "_"))) != ""
			skippedFirst := false
			for _, layer := range layers {
				if !layer.InConfig(key) {
					continue
				}
				if !skippedFirst && !envActive {
					skippedFirst = true
					continue
				}
				suffix := "(overridden)"
				if envActive {
					suffix = "(overridden by env)"
				}
				fmt.Printf("%s%s\n", indent,
					utils.StyleDebug(layer.GetString(key)+" ["+layer.Type+"] "+suffix))
			}
		}

		// Show config file search paths
		fmt.Println(utils.StyleTitle("Config File Search Paths:"))
		searchPaths := config.GetConfigSearchPaths()
		foundActive := false
		for i, sp := range searchPaths {
			status := ""
			if sp.InUse {
				foundActive = true
			} else if !sp.Exists {
				status = " " + utils.StyleWarning("(not found)")
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
		seenBaseDirs := make(map[string]bool)

		addBaseDir := func(label, path string, ro bool) {
			if path == "" || seenBaseDirs[path] {
				return
			}
			seenBaseDirs[path] = true
			status := ""
			if _, err := os.Stat(path); err == nil {
				if ro {
					status = " " + utils.StyleWarning("(search-only)")
				} else if utils.CanWriteToDir(path) {
					status = " " + utils.StyleSuccess("(writable)")
				} else {
					status = " " + utils.StyleWarning("(read-only)")
				}
			}
			fmt.Printf("  %d. [%s] %s%s\n", pathIndex, label, path, status)
			pathIndex++
		}

		// Explicit extra image/build/helper directories
		for _, entry := range config.GetExtraImageDirs() {
			path, ro := config.ParseDirEntry(entry)
			addBaseDir("extra-images", path, ro)
		}
		for _, path := range config.GetExtraBuildDirs() {
			addBaseDir("extra-build", path, false)
		}
		for _, entry := range config.GetExtraHelperDirs() {
			path, ro := config.ParseDirEntry(entry)
			addBaseDir("extra-helper", path, ro)
		}

		// Tier base directories (extra-root → root → scratch → user)
		addBaseDir("extra-root", config.GetExtraRootDir(), false)
		addBaseDir("app-root", config.GetRootDir(), false)
		addBaseDir("scratch", config.GetScratchDataDir(), false)
		addBaseDir("user", config.GetUserDataDir(), false)

		fmt.Println()

		// Show all settings
		fmt.Println(utils.StyleTitle("Current Configuration:"))
		fmt.Println()

		// Paths: binaries then directories
		fmt.Println(utils.StyleTitle("Paths:"))
		fmt.Printf("  apptainer_bin: %s%s\n", viper.GetString("apptainer_bin"), srcTag("apptainer_bin"))
		printOverridden("                 ", "apptainer_bin")
		schedulerBin := viper.GetString("scheduler_bin")
		schedulerType := config.GetSchedulerTypeFromBin(schedulerBin)
		if schedulerBin != "" {
			fmt.Printf("  scheduler_bin: %s (%s)%s\n", schedulerBin, schedulerType, srcTag("scheduler_bin"))
		} else {
			fmt.Printf("  scheduler_bin: %s%s\n", schedulerBin, srcTag("scheduler_bin"))
		}
		printOverridden("                 ", "scheduler_bin")
		logsDir := viper.GetString("logs_dir")
		if logsDir == "" {
			logsDir = config.Global.LogsDir + " (default)"
		}
		fmt.Printf("  logs_dir:      %s%s\n", logsDir, srcTag("logs_dir"))
		printOverridden("                 ", "logs_dir")
		extraImageDirs := config.GetExtraImageDirs()
		if len(extraImageDirs) > 0 {
			fmt.Printf("  extra_image_dirs:\n")
			for _, entry := range extraImageDirs {
				path, ro := config.ParseDirEntry(entry)
				roSuffix := ""
				if ro {
					roSuffix = " " + utils.StyleWarning("(search-only)")
				}
				fmt.Printf("    - %s%s%s\n", path, roSuffix, srcEntryTag("extra_image_dirs", entry))
			}
		} else {
			fmt.Printf("  extra_image_dirs: %s\n", utils.StyleInfo("none"))
		}
		extraBuildDirs := config.GetExtraBuildDirs()
		if len(extraBuildDirs) > 0 {
			fmt.Printf("  extra_build_dirs:\n")
			for _, path := range extraBuildDirs {
				fmt.Printf("    - %s%s\n", path, srcEntryTag("extra_build_dirs", path))
			}
		} else {
			fmt.Printf("  extra_build_dirs: %s\n", utils.StyleInfo("none"))
		}
		extraHelperDirs := config.GetExtraHelperDirs()
		if len(extraHelperDirs) > 0 {
			fmt.Printf("  extra_helper_dirs:\n")
			for _, path := range extraHelperDirs {
				fmt.Printf("    - %s%s\n", path, srcEntryTag("extra_helper_dirs", path))
			}
		} else {
			fmt.Printf("  extra_helper_dirs: %s\n", utils.StyleInfo("none"))
		}
		fmt.Println()

		// Remote sources
		fmt.Println(utils.StyleTitle("Remote Sources:"))
		fmt.Printf("  scripts_link:  %s%s\n", config.Global.ScriptsLink, srcTag("scripts_link"))
		printOverridden("                 ", "scripts_link")
		extraScriptsLinks := config.GetExtraScriptsLinks()
		if len(extraScriptsLinks) > 0 {
			fmt.Printf("  extra_scripts_links:\n")
			for _, l := range extraScriptsLinks {
				fmt.Printf("    - %s%s\n", l, srcEntryTag("extra_scripts_links", l))
			}
		} else {
			fmt.Printf("  extra_scripts_links: %s\n", utils.StyleInfo("none"))
		}
		fmt.Printf("  prefer_remote: %v%s\n", config.Global.PreferRemote, srcTag("prefer_remote"))
		printOverridden("                 ", "prefer_remote")
		fmt.Println()

		// Options (longest key: metadata_cache_ttl = 18 chars)
		fmt.Println(utils.StyleTitle("Options:"))
		fmt.Printf("  %-19s %s%s\n", "default_distro:", config.Global.DefaultDistro, srcTag("default_distro"))
		printOverridden("                      ", "default_distro")
		submitJobConfig := viper.GetBool("submit_job")
		submitJobActual := config.Global.SubmitJob
		if submitJobConfig && !submitJobActual {
			fmt.Printf("  %-19s %v (disabled: scheduler not accessible)%s\n", "submit_job:", submitJobConfig, srcTag("submit_job"))
		} else {
			fmt.Printf("  %-19s %v%s\n", "submit_job:", submitJobActual, srcTag("submit_job"))
		}
		printOverridden("                      ", "submit_job")
		if config.Global.SchedulerTimeout == 0 {
			fmt.Printf("  %-19s 0 (disabled)%s\n", "scheduler_timeout:", srcTag("scheduler_timeout"))
		} else {
			fmt.Printf("  %-19s %s%s\n", "scheduler_timeout:", utils.FormatDuration(config.Global.SchedulerTimeout), srcTag("scheduler_timeout"))
		}
		printOverridden("                      ", "scheduler_timeout")
		fmt.Printf("  %-19s %v%s\n", "parse_module_load:", config.Global.ParseModuleLoad, srcTag("parse_module_load"))
		printOverridden("                      ", "parse_module_load")
		if config.Global.MetadataCacheTTL == 0 {
			fmt.Printf("  %-19s 0 (disabled)%s\n", "metadata_cache_ttl:", srcTag("metadata_cache_ttl"))
		} else {
			fmt.Printf("  %-19s %dd%s\n", "metadata_cache_ttl:", int(config.Global.MetadataCacheTTL.Hours()/24), srcTag("metadata_cache_ttl"))
		}
		printOverridden("                      ", "metadata_cache_ttl")
		channels := config.Global.Build.Channels
		if len(channels) > 0 {
			fmt.Printf("  %-19s\n", "channels:")
			for _, ch := range channels {
				fmt.Printf("    - %s%s\n", ch, srcEntryTag("channels", ch))
			}
		} else {
			fmt.Printf("  %-19s %s\n", "channels:", utils.StyleInfo("none"))
		}
		fmt.Println()

		// Build settings (longest key: tmp_overlay_size = 16 chars)
		fmt.Printf("%s %s\n", utils.StyleTitle("Build Configuration:"), "build.*")
		fmt.Printf("  %-17s %d%s\n", "ncpus:", config.Global.Build.Defaults.CpusPerTask, srcTag("build.ncpus"))
		printOverridden("                    ", "build.ncpus")
		fmt.Printf("  %-17s %s%s\n", "mem:", utils.FormatMemoryMB(config.Global.Build.Defaults.MemPerNodeMB), srcTag("build.mem"))
		printOverridden("                    ", "build.mem")
		fmt.Printf("  %-17s %s%s\n", "time:", utils.FormatDuration(config.Global.Build.Defaults.Time), srcTag("build.time"))
		printOverridden("                    ", "build.time")
		// Show actual compress_args (may be auto-detected based on apptainer version)
		compressArgs := viper.GetString("build.compress_args")
		actualCompressArgs := config.Global.Build.CompressArgs
		if compressArgs != actualCompressArgs {
			fmt.Printf("  %-17s %s%s\n", "compress_args:", actualCompressArgs, srcTag("build.compress_args"))
		} else {
			fmt.Printf("  %-17s %s%s\n", "compress_args:", compressArgs, srcTag("build.compress_args"))
		}
		printOverridden("                    ", "build.compress_args")
		fmt.Printf("  %-17s %s%s\n", "block_size:", config.Global.Build.BlockSize, srcTag("build.block_size"))
		printOverridden("                    ", "build.block_size")
		fmt.Printf("  %-17s %s%s\n", "data_block_size:", config.Global.Build.DataBlockSize, srcTag("build.data_block_size"))
		printOverridden("                    ", "build.data_block_size")
		fmt.Printf("  %-17s %v%s\n", "use_tmp_overlay:", config.Global.Build.UseTmpOverlay, srcTag("build.use_tmp_overlay"))
		printOverridden("                    ", "build.use_tmp_overlay")
		fmt.Printf("  %-17s %s%s\n", "tmp_overlay_size:", utils.FormatMemoryMB(int64(config.Global.Build.TmpSizeMB)), srcTag("build.tmp_overlay_size"))
		printOverridden("                    ", "build.tmp_overlay_size")
		fmt.Println()

		// Show environment variable overrides
		fmt.Println(utils.StyleTitle("Environment Variable Overrides:"))
		envVars := getConfigEnvVars()
		hasEnvOverrides := false
		// Special vars not derived from config keys (system → group → personal)
		for _, envVar := range []string{"CNT_ROOT", "CNT_EXTRA_ROOT", "SCRATCH"} {
			if val := os.Getenv(envVar); val != "" {
				fmt.Printf("  %s=%s\n", envVar, val)
				hasEnvOverrides = true
			}
		}
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

Without -l, returns the effective merged value (env overrides + all config layers).
With -l, reads only the specified config layer file.`,
	Example: `  condatainer config get apptainer_bin
  condatainer config get build.ncpus
  condatainer config get extra_image_dirs
  condatainer config get extra_image_dirs -l app-root   # app-root layer only`,
	Args:              cobra.ExactArgs(1),
	ValidArgsFunction: configKeysCompletion,
	Run: func(cmd *cobra.Command, args []string) {
		key := args[0]

		// Layer-specific read
		if getLocation != "" {
			configPath, _, err := config.ResolveReadableConfigPath(getLocation)
			if err != nil {
				ExitWithError("Invalid location: %v", err)
			}
			if isArrayKey(key) {
				for _, v := range config.ReadConfigSliceKey(configPath, key) {
					fmt.Println(v)
				}
			} else {
				fmt.Println(config.ReadConfigKey(configPath, key))
			}
			return
		}

		// Merged effective value: env > layers in priority order > viper default
		if _, known := configKeyDefs[key]; !known {
			ExitWithUsageError("Unknown config key: %s", key)
		}
		if isArrayKey(key) {
			for _, v := range viper.GetStringSlice(key) {
				fmt.Println(v)
			}
		} else {
			envVar := "CNT_" + strings.ToUpper(strings.ReplaceAll(key, ".", "_"))
			if ev := os.Getenv(envVar); ev != "" {
				fmt.Println(ev)
			} else {
				found := false
				for _, layer := range config.GetConfigLayerInfos() {
					if layer.InConfig(key) {
						fmt.Println(layer.GetString(key))
						found = true
						break
					}
				}
				if !found {
					fmt.Println(viper.GetString(key)) // default
				}
			}
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
	ValidArgsFunction: configSetKeysCompletion,
	Run: func(cmd *cobra.Command, args []string) {
		key := args[0]
		value := args[1]
		if strings.TrimSpace(value) == "" {
			ExitWithError("value cannot be empty")
		}

		// configKeyDefs is the single source of truth for known keys
		knownKeys := configKeyDefs

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

		if _, known := knownKeys[key]; !known {
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

		if key == "build.mem" {
			if mb, err := utils.ParseMemoryMB(value); err != nil || mb <= 0 {
				utils.PrintError("Invalid memory format: %s", value)
				utils.PrintHint("Use format like: 8GB, 16384MB, 8192")
				os.Exit(ExitCodeUsage)
			} else {
				value = fmt.Sprintf("%d", mb)
			}
		}

		if key == "build.tmp_overlay_size" {
			if mb, err := utils.ParseMemoryMB(value); err != nil || mb <= 0 {
				utils.PrintError("Invalid memory format: %s", value)
				utils.PrintHint("Use format like: 10g, 20480m, 20480, 1t")
				os.Exit(ExitCodeUsage)
			} else {
				value = fmt.Sprintf("%d", mb)
			}
		}

		if key == "build.compress_args" {
			value = config.NormalizeCompressArgs(value)
		}

		configPath, locationType, err := config.ResolveWritableConfigPath(setLocation)
		if err != nil {
			ExitWithError("%v", err)
		}
		if err := config.UpdateConfigKey(configPath, key, value); err != nil {
			ExitWithError("Failed to save config: %v", err)
		}

		utils.PrintSuccess("Set %s = %s", utils.StyleInfo(key), utils.StyleInfo(value))
		utils.PrintNote("Config saved to: %s (%s)", configPath, locationType)
	},
}

var configInitCmd = &cobra.Command{
	Use:   "init",
	Short: "Create a config file with defaults",
	Long: `Create a configuration file with default values and auto-detected settings.

Config location options (-l, --location):
  u, user        ~/.config/condatainer/config.yaml (default for standard installations)
  r, app-root    <install-dir>/config.yaml (default if installed in a dedicated folder, or CNT_ROOT set)
  e, extra-root  $CNT_EXTRA_ROOT/config.yaml (group/lab layer, requires CNT_EXTRA_ROOT)
  s, system      /etc/condatainer/config.yaml (requires appropriate permissions)

By default, the location is chosen based on the installation:
  - If the executable is in a dedicated folder (e.g., /apps/condatainer/bin/),
    the config is created in the parent folder (standalone layout, -l app-root).
  - Otherwise, it uses the user's config directory.`,
	Example: `  condatainer config init              # Create config with smart default location
  condatainer config init -l user      # Force user config location
  condatainer config init -l app-root  # Force app-root config location
  condatainer config init -l r         # Same as -l app-root
  CNT_EXTRA_ROOT=/lab condatainer config init -l extra-root  # Group/lab layer`,
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
			locationType = config.NormalizeConfigLocation(initLocation)
			// Fail early with an actionable error if the root dir is read-only
			if locationType == "app-root" {
				if rootDir := filepath.Dir(configPath); !utils.CanWriteToDir(rootDir) {
					ExitWithError(
						"App-root config directory is read-only: %s\n"+
							"Run as a privileged user, or use '-l user' to create a personal config instead.",
						rootDir,
					)
				}
			}
		} else {
			// Smart default: root if in dedicated installation AND writable, otherwise user
			rootPath := config.GetRootConfigPath()
			if rootPath != "" && utils.CanWriteToDir(filepath.Dir(rootPath)) {
				configPath = rootPath
				locationType = "app-root"
			} else {
				if rootPath != "" {
					// Standalone layout detected but dir is read-only — fall back gracefully
					utils.PrintNote("Standalone layout detected but directory is read-only; using user config instead.")
				}
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

		// lowerLayersFor returns loaded config layers that are lower priority than loc.
		// Keys already set in these layers will be skipped when saving.
		layerOrder := []string{"user", "extra-root", "app-root", "system"}
		lowerLayersFor := func(loc string) []config.ConfigLayerInfo {
			cutIdx := slices.Index(layerOrder, loc) + 1
			var result []config.ConfigLayerInfo
			for _, l := range config.GetConfigLayerInfos() {
				if slices.Index(layerOrder, l.Type) >= cutIdx {
					result = append(result, l)
				}
			}
			return result
		}

		// Detect the three minimal keys: apptainer_bin, scheduler_bin, build.compress_args
		detectedApptainerBin := config.DetectApptainerBin()
		if detectedApptainerBin == "" {
			ExitWithError("Neither 'apptainer' nor 'singularity' binary found (checked PATH and 'module avail').")
		}
		viper.Set("apptainer_bin", detectedApptainerBin)

		detectedSchedulerBin := config.DetectSchedulerBin()
		if detectedSchedulerBin != "" {
			viper.Set("scheduler_bin", detectedSchedulerBin)
		}

		detectedCompression := config.Global.Build.CompressArgs // fallback to runtime default
		if err := apptainer.SetBin(detectedApptainerBin); err == nil {
			if version, err := apptainer.GetVersion(); err == nil {
				config.AutoDetectCompression(apptainer.CheckZstdSupport(version), apptainer.IsSingularity())
				detectedCompression = config.Global.Build.CompressArgs
			}
		}
		viper.Set("build.compress_args", detectedCompression)

		if err := config.SaveMinimalConfigTo(configPath, detectedApptainerBin, detectedSchedulerBin, detectedCompression, lowerLayersFor(locationType)); err != nil {
			ExitWithError("Failed to save config: %v", err)
		}

		utils.PrintSuccess("Config file created")
		fmt.Printf("  Location: %s (%s)\n", utils.StylePath(configPath), locationType)

		// Show what was detected
		fmt.Println()
		fmt.Println(utils.StyleTitle("Detected settings:"))
		fmt.Printf("  Apptainer: %s\n", detectedApptainerBin)
		if detectedSchedulerBin != "" {
			fmt.Printf("  Scheduler: %s (%s)\n", detectedSchedulerBin, config.GetSchedulerTypeFromBin(detectedSchedulerBin))
		} else {
			fmt.Printf("  Scheduler: %s\n", utils.StyleWarning("not found"))
		}
		fmt.Printf("  Compression: %s\n", detectedCompression)
	},
}

var configPathsCmd = &cobra.Command{
	Use:   "paths",
	Short: "Show data search paths",
	Long: `Show all search paths for images, build-scripts, and helper-scripts.

Search paths are checked in priority order (first match wins for reads):
  1. Extra base directories (highest priority)
  2. Extra-root ($CNT_EXTRA_ROOT, group/lab layer)
  3. App-root (auto-detected or $CNT_ROOT)
  4. Scratch (user HPC large storage, $SCRATCH/condatainer)
  5. User XDG directory (~/.local/share/condatainer)

Write operations use the first writable directory in the same order.
App-root is preferred for group/shared use, user is the fallback.`,
	Run: func(cmd *cobra.Command, args []string) {
		fmt.Println(utils.StyleTitle("Data Search Paths"))
		fmt.Println()

		// pathStatus returns inline status tags for a search-path entry.
		// writeTarget is the resolved writable directory for this section (empty = not applicable).
		pathStatus := func(dir, writeTarget string) string {
			if !config.DirExists(dir) {
				return " " + utils.StyleWarning("(not found)")
			}
			var tags string
			if utils.CanWriteToDir(dir) {
				if writeTarget != "" && dir == writeTarget {
					tags = " " + utils.StyleSuccess("(writable, target)")
				} else {
					tags = " " + utils.StyleSuccess("(writable)")
				}
			} else {
				tags = " " + utils.StyleWarning("(read-only)")
			}
			return tags
		}

		imagesWritable, _ := config.GetWritableImagesDir()
		helperWritable, _ := config.GetWritableHelperScriptsDir()
		cacheWritable, _ := config.GetWritableCacheDir()

		// Images
		fmt.Println(utils.StyleTitle("Images:"))
		for i, dir := range config.GetImageSearchPaths() {
			fmt.Printf("  %d. %s%s\n", i+1, dir, pathStatus(dir, imagesWritable))
		}
		fmt.Println()

		// Build scripts (read-only search — nothing writes here)
		fmt.Println(utils.StyleTitle("Build Scripts:"))
		for i, dir := range config.GetBuildScriptSearchPaths() {
			fmt.Printf("  %d. %s%s\n", i+1, dir, pathStatus(dir, ""))
		}
		fmt.Println()

		// Helper scripts
		fmt.Println(utils.StyleTitle("Helper Scripts:"))
		for i, dir := range config.GetHelperScriptSearchPaths() {
			fmt.Printf("  %d. %s%s\n", i+1, dir, pathStatus(dir, helperWritable))
		}
		fmt.Println()

		// Cache
		fmt.Println(utils.StyleTitle("Cache:"))
		for i, dir := range config.GetCacheSearchPaths() {
			fmt.Printf("  %d. %s%s\n", i+1, dir, pathStatus(dir, cacheWritable))
		}
		fmt.Println()

		// Config file search paths
		fmt.Println(utils.StyleTitle("Config Files:"))
		for _, cp := range config.GetConfigSearchPaths() {
			var status string
			if !cp.Exists {
				status = " " + utils.StyleWarning("(not found)")
			}
			fmt.Printf("  %-12s %s%s\n", cp.Type+":", cp.Path, status)
		}
		fmt.Println()

		// Directory info
		fmt.Println(utils.StyleTitle("Directory Locations:"))
		if rootDir := config.GetRootDir(); rootDir != "" {
			fmt.Printf("  App-root: %s\n", rootDir)
		} else {
			fmt.Printf("  App-root: %s\n", utils.StyleWarning("not detected (set CNT_ROOT to override)"))
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
		ncpus := config.Global.Build.Defaults.CpusPerTask
		if ncpus > 0 {
			if !utils.QuietMode {
				fmt.Printf("%s Build CPUs: %d\n", utils.StyleSuccess("✓"), ncpus)
			}
		} else {
			fmt.Printf("%s Build CPUs must be > 0: %d\n", utils.StyleError("✗"), ncpus)
			valid = false
		}

		memMB := config.Global.Build.Defaults.MemPerNodeMB
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
	Example:           "  condatainer config append extra_image_dirs /shared/lab/images:ro\n  condatainer config append extra_scripts_links https://example.com/scripts",
	Args:              cobra.ExactArgs(2),
	ValidArgsFunction: arrayKeyCompletion,
	SilenceUsage:      true,
	Run: func(cmd *cobra.Command, args []string) {
		key, value := args[0], args[1]
		if strings.TrimSpace(value) == "" {
			ExitWithError("value cannot be empty")
		}
		moved := false
		configPath, err := modifyArrayConfig(key, func(cur []string) []string {
			out := make([]string, 0, len(cur))
			for _, v := range cur {
				if v == value {
					moved = true
					continue
				}
				out = append(out, v)
			}
			return append(out, value)
		})
		if err != nil {
			ExitWithError("%v", err)
		}
		if moved {
			utils.PrintSuccess("Moved %s to end of %s", utils.StyleInfo(value), utils.StyleInfo(key))
		} else {
			utils.PrintSuccess("Appended %s to %s", utils.StyleInfo(value), utils.StyleInfo(key))
		}
		utils.PrintNote("Config saved to: %s", configPath)
	},
}

var configPrependCmd = &cobra.Command{
	Use:               "prepend <key> <value>",
	Short:             "Prepend a value to an array config key (highest priority)",
	Example:           "  condatainer config prepend extra_image_dirs /fast/images\n  condatainer config prepend extra_scripts_links https://myorg.com/scripts",
	Args:              cobra.ExactArgs(2),
	ValidArgsFunction: arrayKeyCompletion,
	SilenceUsage:      true,
	Run: func(cmd *cobra.Command, args []string) {
		key, value := args[0], args[1]
		if strings.TrimSpace(value) == "" {
			ExitWithError("value cannot be empty")
		}
		moved := false
		configPath, err := modifyArrayConfig(key, func(cur []string) []string {
			out := make([]string, 0, len(cur))
			for _, v := range cur {
				if v == value {
					moved = true
					continue
				}
				out = append(out, v)
			}
			return append([]string{value}, out...)
		})
		if err != nil {
			ExitWithError("%v", err)
		}
		if moved {
			utils.PrintSuccess("Moved %s to front of %s", utils.StyleInfo(value), utils.StyleInfo(key))
		} else {
			utils.PrintSuccess("Prepended %s to %s", utils.StyleInfo(value), utils.StyleInfo(key))
		}
		utils.PrintNote("Config saved to: %s", configPath)
	},
}

var configRemoveCmd = &cobra.Command{
	Use:               "remove <key> [value]",
	Short:             "Remove a config key or a value from an array key",
	Example:           "  condatainer config remove apptainer_bin\n  condatainer config remove extra_image_dirs /shared/lab/images:ro",
	Args:              cobra.RangeArgs(1, 2),
	ValidArgsFunction: arrayRemoveValueCompletion,
	SilenceUsage:      true,
	Run: func(cmd *cobra.Command, args []string) {
		key := args[0]
		if len(args) == 1 {
			// Scalar removal: unset the key entirely
			if isArrayKey(key) {
				ExitWithError("'%s' is an array key; specify a value to remove: config remove %s <value>", key, key)
			}
			configPath, _, err := config.ResolveWritableConfigPath(setLocation)
			if err != nil {
				ExitWithError("%v", err)
			}
			if config.ReadConfigKey(configPath, key) == "" {
				utils.PrintWarning("%s is not set in %s", utils.StyleInfo(key), utils.StylePath(configPath))
				return
			}
			if err := config.DeleteConfigKey(configPath, key); err != nil {
				ExitWithError("%v", err)
			}
			utils.PrintSuccess("Removed %s from config", utils.StyleInfo(key))
			utils.PrintNote("Config saved to: %s", configPath)
			return
		}
		// Array removal
		value := args[1]
		removed := false
		configPath, err := modifyArrayConfig(key, func(cur []string) []string {
			var out []string
			for _, v := range cur {
				if v == value {
					removed = true
					continue
				}
				out = append(out, v)
			}
			return out
		})
		if err != nil {
			ExitWithError("%v", err)
		}
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
	configShowCmd.Flags().BoolVarP(&showSources, "sources", "s", false, "Show which config layer provides each value")
	configInitCmd.Flags().StringVarP(&initLocation, "location", "l", "", "Config location: user (u), app-root (r), extra-root (e), or system (s)")
	configGetCmd.Flags().StringVarP(&getLocation, "location", "l", "", "Read from a specific config layer: user (u), app-root (r), extra-root (e), or system (s)")
	configSetCmd.Flags().StringVarP(&setLocation, "location", "l", "", "Config location to write: user (u), app-root (r), extra-root (e), or system (s)")
	configAppendCmd.Flags().StringVarP(&setLocation, "location", "l", "", "Config location to write: user (u), app-root (r), extra-root (e), or system (s)")
	configPrependCmd.Flags().StringVarP(&setLocation, "location", "l", "", "Config location to write: user (u), app-root (r), extra-root (e), or system (s)")
	configRemoveCmd.Flags().StringVarP(&setLocation, "location", "l", "", "Config location to write: user (u), app-root (r), extra-root (e), or system (s)")

	// Add subcommands
	configCmd.AddCommand(configShowCmd)
	configCmd.AddCommand(configGetCmd)
	configCmd.AddCommand(configSetCmd)
	configCmd.AddCommand(configAppendCmd)
	configCmd.AddCommand(configPrependCmd)
	configCmd.AddCommand(configRemoveCmd)
	configCmd.AddCommand(configInitCmd)
	configCmd.AddCommand(configPathsCmd)
	configCmd.AddCommand(configValidateCmd)

	// Add to root command
	rootCmd.AddCommand(configCmd)
}
