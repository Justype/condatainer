package cmd

import (
	"context"
	"fmt"
	"os"
	"path/filepath"
	"slices"
	"sort"
	"strconv"
	"strings"

	"github.com/Justype/condatainer/internal/config"
	"github.com/Justype/condatainer/internal/utils"
	"github.com/spf13/cobra"
	"github.com/spf13/viper"
)

var (
	showPath    bool
	showSources bool
	initLayer   string // user, app-root, extra-root, system, or a custom path
	getLayer    string // target layer for config get (single-layer read)
	setLayer    string // target layer for config set/append/prepend/remove
)

// configKeyDefs maps every known config key to whether it holds a string slice (array).
// true = array key (use append/prepend/remove); false = scalar key (use set).
var configKeyDefs = map[string]bool{
	"logs_dir":                   false,
	"default_distro":             false,
	"submit_job":                 false,
	"sources":                    true,
	"autoload_gpu":               false,
	"notification":               false,
	"metadata_cache_ttl":         false,
	"store_gc_grace":             false,
	"proxy_perjob":               false,
	"helper_bind_all":            false,
	"build.system_apptainer":     false,
	"build.ncpus":                false,
	"build.mem":                  false,
	"build.time":                 false,
	"build.compress_args":        false,
	"build.block_size":           false,
	"build.data_block_size":      false,
	"build.app_tmp_overlay":      false,
	"build.always_submit":        false,
	"build.app_tmp_overlay_size": false,
	"scheduler.bin":              false,
	"scheduler.timeout":          false,
	"scheduler.account":          false,
	"scheduler.partition":        false,
	"scheduler.ncpus":            false,
	"scheduler.mem":              false,
	"scheduler.time":             false,
	"channels":                   true,
}

func isArrayKey(key string) bool { return configKeyDefs[key] }

func isBoolKey(key string) bool {
	switch key {
	case "submit_job", "proxy_perjob", "helper_bind_all", "autoload_gpu",
		"build.app_tmp_overlay", "build.always_submit":
		return true
	}
	return false
}

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
	configPath, _, fellBackFrom, err := config.ResolveWritableConfigPathVerbose(setLayer)
	if err != nil {
		return "", err
	}
	if fellBackFrom != "" {
		utils.PrintWarning("%s config is read-only; saving to your user config instead (applies only to you).", fellBackFrom)
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
		case "sources":
			return nil, cobra.ShellCompDirectiveDefault
		}
	}
	return nil, cobra.ShellCompDirectiveNoFileComp
}

func arrayRemoveValueCompletion(cmd *cobra.Command, args []string, toComplete string) ([]string, cobra.ShellCompDirective) {
	loc, _ := cmd.Flags().GetString("layer")
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
		// sources entries are name=base pairs; allow free-form input
		if args[0] == "sources" {
			return nil, cobra.ShellCompDirectiveDefault
		}
		return configValueCompletion(args[0]), cobra.ShellCompDirectiveNoFileComp
	}
	return nil, cobra.ShellCompDirectiveNoFileComp
}

// recipeExists reports whether a module name resolves in any configured source.
func recipeExists(ctx context.Context, name string) bool {
	cat, err := config.OpenCatalog(ctx)
	if err != nil {
		return false
	}
	_, found, err := cat.Lookup(ctx, name)
	return err == nil && found
}

// configValueCompletion returns suggested values for a config key
func configValueCompletion(key string) []string {
	switch key {
	case "submit_job", "proxy_perjob", "helper_bind_all":
		return []string{"true", "false"}
	case "autoload_gpu":
		return []string{"true", "false"}
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
	case "build.app_tmp_overlay", "build.always_submit":
		return []string{"true", "false"}
	case "build.app_tmp_overlay_size":
		return []string{"10g", "20g", "40g"}
	case "scheduler.ncpus":
		return []string{"1", "2", "4", "8"}
	case "scheduler.mem":
		return []string{"2g", "4g", "8g", "16g"}
	case "scheduler.time":
		return []string{"1h", "2h", "4h", "8h"}
	case "notification":
		return []string{"none", "terminal", "web", "both"}
	case "metadata_cache_ttl":
		return []string{"1", "3", "7", "14", "0"}
	case "store_gc_grace":
		return []string{"7", "30", "90", "180"}
	default:
		return nil
	}
}

// configLayersHelp documents the -l/--layer values. Shared by every config
// subcommand that accepts --layer so the list stays consistent, not copy-pasted.
const configLayersHelp = `Config file layers (-l, --layer):
  u, user        ~/.config/condatainer/config.yaml (standard/home install)
  r, app-root    <install-dir>/config.yaml         (dedicated install, or CNT_ROOT)
  e, extra-root  $CNT_EXTRA_ROOT/config.yaml       (group, needs CNT_EXTRA_ROOT)
  s, system      /etc/condatainer/config.yaml`

var configCmd = &cobra.Command{
	Use:   "config",
	Short: "Manage condatainer configuration",
	Long: `Manage condatainer configuration settings.

Configuration file priority (highest to lowest):
  1. Command-line flags
  2. Environment variables (CNT_*)
  3. User config file (~/.config/condatainer/config.yaml)
  4. Extra-root config ($CNT_EXTRA_ROOT/config.yaml, group layer)
  5. App-root config (<install-dir>/config.yaml, in dedicated folder or CNT_ROOT set)
  6. System config file (/etc/condatainer/config.yaml)
  7. Defaults

Data directory priority — reads go nearest-first, writes furthest-first:

  read (first match wins)          write (first writable)
  1. Scratch ($SCRATCH)            1. Extra-root ($CNT_EXTRA_ROOT, group layer)
  2. User XDG (~/.local/share)     2. App-root (auto-detected or CNT_ROOT)
  3. Extra-root ($CNT_EXTRA_ROOT)  3. Scratch ($SCRATCH/condatainer)
  4. App-root (auto-detected)      4. User XDG (~/.local/share/condatainer)

A build lands as far out as permissions allow, so one copy serves the group;
build into your own directory to have your version win for you.`,
}

var configShowCmd = &cobra.Command{
	Use:   "show",
	Args:  cobra.NoArgs,
	Short: "Show current configuration",
	Long: `Show current settings and where each comes from:
  search paths, all values, and environment-variable overrides.
`,
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

		// Tier base directories (extra-root → root → scratch → user)
		addBaseDir("extra-root", config.GetExtraRootDir(), false)
		addBaseDir("app-root", config.GetRootDir(), false)
		addBaseDir("scratch", config.GetScratchDataDir(), false)
		addBaseDir("user", config.GetUserDataDir(), false)

		fmt.Println()

		// Show all settings
		fmt.Println(utils.StyleTitle("Paths:"))
		fmt.Printf("  logs_dir: %s%s\n", config.Global.LogsDir, srcTag("logs_dir"))
		printOverridden("            ", "logs_dir")
		fmt.Println()

		// Remote sources
		fmt.Println(utils.StyleTitle("Recipe Sources:"))
		if len(config.Global.Sources) > 0 {
			fmt.Printf("  sources:\n")
			for _, src := range config.Global.Sources {
				fmt.Printf("    - %s: %s\n", src.Name, src.Base)
			}
		} else {
			fmt.Printf("  sources: %s\n", utils.StyleInfo("none configured"))
		}
		fmt.Println()

		// Options (longest key: metadata_cache_ttl = 18 chars)
		fmt.Println(utils.StyleTitle("Options:"))
		fmt.Printf("  %-19s %s%s\n", "default_distro:", config.ResolvedDefaultDistro(), srcTag("default_distro"))
		printOverridden("                      ", "default_distro")
		if cat, err := config.OpenCatalog(cmd.Context()); err == nil {
			if def := config.SourceDefaultDistro(cat); def != "" && def != config.ResolvedDefaultDistro() {
				fmt.Printf("                      %s\n",
					utils.StyleHint(fmt.Sprintf("sources now recommend %s — `config set default_distro %s` to switch", def, def)))
			}
		}
		submitJobConfig := viper.GetBool("submit_job")
		submitJobActual := config.Global.SubmitJob
		if submitJobConfig && !submitJobActual {
			fmt.Printf("  %-19s %v (disabled: scheduler not accessible)%s\n", "submit_job:", submitJobConfig, srcTag("submit_job"))
		} else {
			fmt.Printf("  %-19s %v%s\n", "submit_job:", submitJobActual, srcTag("submit_job"))
		}
		printOverridden("                      ", "submit_job")
		fmt.Printf("  %-19s %v%s\n", "autoload_gpu:", config.Global.AutoloadGPU, srcTag("autoload_gpu"))
		printOverridden("                      ", "autoload_gpu")
		if config.Global.MetadataCacheTTL == 0 {
			fmt.Printf("  %-19s 0 (disabled)%s\n", "metadata_cache_ttl:", srcTag("metadata_cache_ttl"))
		} else {
			fmt.Printf("  %-19s %dd%s\n", "metadata_cache_ttl:", int(config.Global.MetadataCacheTTL.Hours()/24), srcTag("metadata_cache_ttl"))
		}
		printOverridden("                      ", "metadata_cache_ttl")
		fmt.Printf("  %-19s %dd%s\n", "store_gc_grace:", int(config.Global.StoreGCGrace.Hours()/24), srcTag("store_gc_grace"))
		printOverridden("                      ", "store_gc_grace")
		fmt.Printf("  %-19s %v%s\n", "proxy_perjob:", config.Global.ProxyPerJob, srcTag("proxy_perjob"))
		printOverridden("                      ", "proxy_perjob")
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

		fmt.Println(utils.StyleTitle("Helper:"))
		fmt.Printf("  %-19s %v%s\n", "helper_bind_all:", config.Global.HelperBindAll, srcTag("helper_bind_all"))
		printOverridden("                       ", "helper_bind_all")
		notif := config.Global.Notification
		if notif == "" {
			notif = "none"
		}
		fmt.Printf("  %-19s %s%s\n", "notification:", notif, srcTag("notification"))
		printOverridden("                       ", "notification")
		fmt.Println()

		// Build settings (longest key: app_tmp_overlay_size = 20 chars)
		fmt.Printf("%s %s\n", utils.StyleTitle("Build Configuration:"), "build.*")
		fmt.Printf("  %-21s %s%s\n", "system_apptainer:", config.Global.Build.SystemApptainer, srcTag("build.system_apptainer"))
		printOverridden("                        ", "build.system_apptainer")
		fmt.Printf("  %-21s %v%s\n", "always_submit:", config.Global.Build.AlwaysSubmit, srcTag("build.always_submit"))
		printOverridden("                        ", "build.always_submit")
		fmt.Printf("  %-21s %d%s\n", "ncpus:", config.Global.Build.Defaults.CpusPerTask, srcTag("build.ncpus"))
		printOverridden("                        ", "build.ncpus")
		fmt.Printf("  %-21s %s%s\n", "mem:", utils.FormatMemoryMB(config.Global.Build.Defaults.MemPerNodeMB), srcTag("build.mem"))
		printOverridden("                        ", "build.mem")
		fmt.Printf("  %-21s %s%s\n", "time:", utils.FormatDuration(config.Global.Build.Defaults.Time), srcTag("build.time"))
		printOverridden("                        ", "build.time")
		// Show the expanded args when a shorthand name (e.g. "zstd-medium")
		// was configured, since the raw value alone would not be mksquashfs-ready.
		compressArgs := viper.GetString("build.compress_args")
		actualCompressArgs := config.Global.Build.CompressArgs
		if compressArgs != actualCompressArgs {
			fmt.Printf("  %-21s %s%s\n", "compress_args:", actualCompressArgs, srcTag("build.compress_args"))
		} else {
			fmt.Printf("  %-21s %s%s\n", "compress_args:", compressArgs, srcTag("build.compress_args"))
		}
		printOverridden("                        ", "build.compress_args")
		fmt.Printf("  %-21s %s%s\n", "block_size:", config.Global.Build.BlockSize, srcTag("build.block_size"))
		printOverridden("                        ", "build.block_size")
		fmt.Printf("  %-21s %s%s\n", "data_block_size:", config.Global.Build.DataBlockSize, srcTag("build.data_block_size"))
		printOverridden("                        ", "build.data_block_size")
		fmt.Printf("  %-21s %v%s\n", "app_tmp_overlay:", config.Global.Build.AppTmpOverlay, srcTag("build.app_tmp_overlay"))
		printOverridden("                        ", "build.app_tmp_overlay")
		fmt.Printf("  %-21s %s%s\n", "app_tmp_overlay_size:", utils.FormatMemoryMB(int64(config.Global.Build.AppTmpOverlaySizeMB)), srcTag("build.app_tmp_overlay_size"))
		printOverridden("                        ", "build.app_tmp_overlay_size")
		fmt.Println()

		// Scheduler settings (longest key: partition = 10 chars)
		fmt.Printf("%s %s\n", utils.StyleTitle("Scheduler Configuration:"), "scheduler.*")
		schedulerBin := config.Global.Scheduler.Bin
		schedulerType := config.GetSchedulerTypeFromBin(schedulerBin)
		if schedulerBin != "" {
			fmt.Printf("  %-12s %s (%s)%s\n", "bin:", schedulerBin, schedulerType, srcTag("scheduler.bin"))
		} else {
			fmt.Printf("  %-12s %s%s\n", "bin:", schedulerBin, srcTag("scheduler.bin"))
		}
		printOverridden("             ", "scheduler.bin")
		if config.Global.Scheduler.Timeout == 0 {
			fmt.Printf("  %-12s 0 (disabled)%s\n", "timeout:", srcTag("scheduler.timeout"))
		} else {
			fmt.Printf("  %-12s %s%s\n", "timeout:", utils.FormatDuration(config.Global.Scheduler.Timeout), srcTag("scheduler.timeout"))
		}
		printOverridden("             ", "scheduler.timeout")
		account := config.Global.Scheduler.Account
		if account == "" {
			account = utils.StyleInfo("(scheduler default)")
		}
		fmt.Printf("  %-12s %s%s\n", "account:", account, srcTag("scheduler.account"))
		printOverridden("             ", "scheduler.account")
		partition := config.Global.Scheduler.Partition
		if partition == "" {
			partition = utils.StyleInfo("(scheduler default)")
		}
		fmt.Printf("  %-12s %s%s\n", "partition:", partition, srcTag("scheduler.partition"))
		printOverridden("             ", "scheduler.partition")
		fmt.Printf("  %-12s %d%s\n", "ncpus:", config.Global.Scheduler.Defaults.CpusPerTask, srcTag("scheduler.ncpus"))
		printOverridden("             ", "scheduler.ncpus")
		fmt.Printf("  %-12s %s%s\n", "mem:", utils.FormatMemoryMB(config.Global.Scheduler.Defaults.MemPerNodeMB), srcTag("scheduler.mem"))
		printOverridden("             ", "scheduler.mem")
		fmt.Printf("  %-12s %s%s\n", "time:", utils.FormatDuration(config.Global.Scheduler.Defaults.Time), srcTag("scheduler.time"))
		printOverridden("             ", "scheduler.time")
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
	Use:   "get [flags] <key>",
	Short: "Get a configuration value",
	Long: `Get a specific configuration value.

Without -l, returns the effective merged value (env + all config layers).
With -l, reads only that config layer file.

` + configLayersHelp,
	Example: `  condatainer config get build.system_apptainer
  condatainer config get build.ncpus
  condatainer config get sources
  condatainer config get sources -l app-root`,
	Args:              cobra.ExactArgs(1),
	ValidArgsFunction: configKeysCompletion,
	Run: func(cmd *cobra.Command, args []string) {
		key := args[0]

		// Layer-specific read
		if getLayer != "" {
			configPath, _, err := config.ResolveReadableConfigPath(getLayer)
			if err != nil {
				ExitWithError("Invalid layer: %v", err)
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
			ExitWithError("Unknown config key: %s", key)
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
						if isBoolKey(key) {
							fmt.Println(layer.GetBool(key))
						} else {
							fmt.Println(layer.GetString(key))
						}
						found = true
						break
					}
				}
				if !found {
					if isBoolKey(key) {
						fmt.Println(viper.GetBool(key))
					} else {
						fmt.Println(viper.GetString(key)) // default
					}
				}
			}
		}
	},
}

var configSetCmd = &cobra.Command{
	Use:   "set [flags] <key> <value>",
	Short: "Set a configuration value",
	Long: `Set a configuration value and save to config file.

Time duration format (for build.time):
  Go style:  2h, 30m, 1h30m, 90s
  HPC style: 02:00:00, 2:30:00, 1:30 (HH:MM:SS or HH:MM)

` + configLayersHelp,
	Example: `  condatainer config set build.system_apptainer /usr/bin/apptainer
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

		// A default distro must name a recipe that some source actually provides.
		if key == "default_distro" {
			if !recipeExists(cmd.Context(), value+"/base") {
				utils.PrintError("No recipe %s/base in any configured source", value)
				os.Exit(ExitCodeError)
			}
		}

		// Array keys require append/prepend/remove subcommands
		if isArrayKey(key) {
			utils.PrintError("'%s' is an array setting. Use append/prepend/remove subcommands.", key)
			utils.PrintHint("  condatainer config append  %s <value>\n  condatainer config prepend %s <value>\n  condatainer config remove  %s <value>", key, key, key)
			os.Exit(ExitCodeError)
		}

		if _, known := knownKeys[key]; !known {
			utils.PrintWarning("Warning: '%s' is not a standard config key", key)
		}

		// Validate value based on key type
		if key == "scheduler.timeout" {
			var n int
			if _, err := fmt.Sscan(value, &n); err != nil || n < 0 {
				utils.PrintError("Invalid value for scheduler.timeout: %s (must be a non-negative integer; 0 disables the timeout)", value)
				os.Exit(ExitCodeError)
			}
		}

		if key == "notification" {
			v := strings.ToLower(value)
			if v != "" && v != "none" && v != "terminal" && v != "web" && v != "both" {
				utils.PrintError("Invalid value for notification: %q (valid values: none, terminal, web, both)", value)
				os.Exit(ExitCodeError)
			}
		}

		if key == "metadata_cache_ttl" {
			var n int
			if _, err := fmt.Sscan(value, &n); err != nil || n < 0 {
				utils.PrintError("Invalid value for metadata_cache_ttl: %s (must be a non-negative integer in days; 0 disables the cache)", value)
				os.Exit(ExitCodeError)
			}
		}

		// No zero: a grace of nothing would make everything collectable the
		// moment it is installed, which is not a policy anyone means to set.
		if key == "store_gc_grace" {
			var n int
			if _, err := fmt.Sscan(value, &n); err != nil || n < 1 {
				utils.PrintError("Invalid value for store_gc_grace: %s (must be a positive integer in days)", value)
				os.Exit(ExitCodeError)
			}
		}

		if key == "build.time" {
			if _, err := utils.ParseWalltime(value); err != nil {
				utils.PrintError("Invalid duration format: %s", value)
				utils.PrintHint("Use format like: 4d12h, 2h30m, 1:30, or 01:30:00")
				os.Exit(ExitCodeError)
			}
		}

		if key == "build.mem" {
			if mb, err := utils.ParseMemoryMB(value); err != nil || mb <= 0 {
				utils.PrintError("Invalid memory format: %s", value)
				utils.PrintHint("Use format like: 8GB, 16384MB, 8192")
				os.Exit(ExitCodeError)
			} else {
				value = fmt.Sprintf("%d", mb)
			}
		}

		if key == "scheduler.time" {
			if _, err := utils.ParseWalltime(value); err != nil {
				utils.PrintError("Invalid duration format: %s", value)
				utils.PrintHint("Use format like: 4d12h, 2h30m, 1:30, or 01:30:00")
				os.Exit(ExitCodeError)
			}
		}

		if key == "scheduler.mem" {
			if mb, err := utils.ParseMemoryMB(value); err != nil || mb <= 0 {
				utils.PrintError("Invalid memory format: %s", value)
				utils.PrintHint("Use format like: 8GB, 16384MB, 8192")
				os.Exit(ExitCodeError)
			} else {
				value = fmt.Sprintf("%d", mb)
			}
		}

		if key == "build.app_tmp_overlay_size" {
			if mb, err := utils.ParseMemoryMB(value); err != nil || mb <= 0 {
				utils.PrintError("Invalid memory format: %s", value)
				utils.PrintHint("Use format like: 10g, 20480m, 20480, 1t")
				os.Exit(ExitCodeError)
			} else {
				value = fmt.Sprintf("%d", mb)
			}
		}

		if key == "build.compress_args" {
			value = config.NormalizeCompressArgs(value)
		}

		if isBoolKey(key) {
			b, err := strconv.ParseBool(value)
			if err != nil {
				utils.PrintError("Invalid value for %s: %q (use true or false)", key, value)
				os.Exit(ExitCodeError)
			}
			value = strconv.FormatBool(b)
		}

		configPath, layerType, fellBackFrom, err := config.ResolveWritableConfigPathVerbose(setLayer)
		if err != nil {
			ExitWithError("%v", err)
		}
		if fellBackFrom != "" {
			utils.PrintWarning("%s config is read-only; saving to your user config instead (applies only to you).", fellBackFrom)
		}
		if err := config.UpdateConfigKey(configPath, key, value); err != nil {
			ExitWithError("Failed to save config: %v", err)
		}

		utils.PrintSuccess("Set %s = %s", utils.StyleInfo(key), utils.StyleInfo(value))
		utils.PrintNote("Config saved to: %s (%s)", configPath, layerType)
	},
}

var configInitCmd = &cobra.Command{
	Use:   "init",
	Args:  cobra.NoArgs,
	Short: "Create a config file with defaults",
	Long: `Create a configuration file with default values and auto-detected settings.

` + configLayersHelp + `

Without -l, the layer is chosen from the install layout:
  - Executable in a dedicated folder (e.g. /apps/condatainer/bin/) -> parent folder.
  - Otherwise -> the user config directory.`,
	Example: `  condatainer config init -l app-root  # write to a specific layer
  condatainer config init -l r         # via shortcut (r = app-root)`,
	Run: func(cmd *cobra.Command, args []string) {
		var configPath string
		var err error
		var layerType string

		// Determine config path based on --layer flag or smart default
		if initLayer != "" {
			// User specified a layer
			configPath, err = config.GetConfigPathByLayer(initLayer)
			if err != nil {
				ExitWithError("Invalid layer: %v", err)
			}
			layerType = config.NormalizeConfigLayer(initLayer)
			// Fail early with an actionable error if the root dir is read-only
			if layerType == "app-root" {
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
				layerType = "app-root"
			} else {
				if rootPath != "" {
					// Standalone layout detected but dir is read-only — fall back gracefully
					utils.PrintNote("Standalone layout detected but directory is read-only; using user config instead.")
				}
				configPath, err = config.GetUserConfigPath()
				if err != nil {
					ExitWithError("Failed to get config path: %v", err)
				}
				layerType = "user"
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

		// Detect the two host-specific keys: build.system_apptainer, scheduler.bin.
		// build.compress_args is not detected — every artifact is always
		// zstd-medium, unconditionally (config.LoadDefaults) — but is still
		// written explicitly so it shows plainly in the generated file.
		detectedApptainerBin := config.DetectApptainerBin()
		if detectedApptainerBin == "" {
			ExitWithError("Neither 'apptainer' nor 'singularity' binary found (checked PATH and 'module avail').")
		}
		viper.Set("build.system_apptainer", detectedApptainerBin)

		detectedSchedulerBin := config.DetectSchedulerBin()
		if detectedSchedulerBin != "" {
			viper.Set("scheduler.bin", detectedSchedulerBin)
		}

		detectedCompression := config.Global.Build.CompressArgs
		viper.Set("build.compress_args", detectedCompression)

		if err := config.SaveMinimalConfigTo(configPath, detectedApptainerBin, detectedSchedulerBin, detectedCompression, lowerLayersFor(layerType)); err != nil {
			ExitWithError("Failed to save config: %v", err)
		}

		utils.PrintSuccess("Config file created")
		fmt.Printf("  Location: %s (%s)\n", utils.StylePath(configPath), layerType)

		// Record the default distro now, from whichever source recommends one.
		// It is written once and never revised, so the container root stays put.
		if cat, err := config.OpenCatalog(cmd.Context()); err == nil {
			if distro := config.EnsureDefaultDistro(cat); distro != "" {
				fmt.Printf("  Distro:   %s (from %s)\n", utils.StyleInfo(distro), utils.StyleInfo("default_distro"))
			}
		}

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
	Args:  cobra.NoArgs,
	Short: "Show data search paths",
	Long: `Show all search paths for images and helper-scripts.

Search paths are checked nearest-first (first match wins for reads):
  1. Scratch (user HPC large storage, $SCRATCH/condatainer)
  2. User XDG directory (~/.local/share/condatainer)
  3. Extra-root ($CNT_EXTRA_ROOT, group layer)
  4. App-root (auto-detected or $CNT_ROOT)

Write operations take the first writable directory in the reverse order —
extra-root, app-root, scratch, user — so a build lands as far out as
permissions allow and one copy serves the whole group.`,
	Run: func(cmd *cobra.Command, args []string) {
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

		// withLayer appends the data layer a directory belongs to, matching the
		// tags in `list` output and the -l/--layer values commands accept.
		withLayer := func(dir string) string {
			return dir + " " + utils.StyleDebug("("+string(config.ClassifyDataDir(dir))+")")
		}

		// Recipe sources, in resolution order (first match wins)
		fmt.Println(utils.StyleTitle("Sources:"))
		if len(config.Global.Sources) == 0 {
			fmt.Println("  (none configured)")
		}
		for i, src := range config.Global.Sources {
			// A source is read-only, so writability says nothing; only whether a
			// local one is actually there.
			status := ""
			if !strings.HasPrefix(src.Base, "http") && !config.DirExists(src.Base) {
				status = " " + utils.StyleWarning("(not found)")
			}
			fmt.Printf("  %d. %s %s%s\n", i+1, src.Name, src.Base, status)
		}
		fmt.Println()

		// Images
		fmt.Println(utils.StyleTitle("Images:"))
		for i, dir := range config.GetImageSearchPaths() {
			fmt.Printf("  %d. %s%s\n", i+1, withLayer(dir), pathStatus(dir, imagesWritable))
		}
		fmt.Println()

		// Helper scripts
		fmt.Println(utils.StyleTitle("Helper Scripts:"))
		for i, dir := range config.GetHelperScriptSearchPaths() {
			fmt.Printf("  %d. %s%s\n", i+1, withLayer(dir), pathStatus(dir, helperWritable))
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
	Args:  cobra.NoArgs,
	Short: "Validate configuration",
	Long:  "Check if the current configuration is valid and all binaries are accessible",
	Run: func(cmd *cobra.Command, args []string) {
		valid := true

		if !utils.QuietMode {
			fmt.Println(utils.StyleTitle("Validating configuration..."))
			fmt.Println()
		}

		// Check apptainer binary
		apptainerBin := viper.GetString("build.system_apptainer")
		if config.ValidateBinary(apptainerBin) {
			if !utils.QuietMode {
				fmt.Printf("%s Apptainer binary: %s\n", utils.StyleSuccess("✓"), apptainerBin)
			}
		} else {
			fmt.Printf("%s Apptainer binary not found: %s\n", utils.StyleError("✗"), apptainerBin)
			valid = false
		}

		// Check scheduler binary
		schedulerBin := viper.GetString("scheduler.bin")
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

		// Check that a base recipe resolves. Unlike the name-expansion paths,
		// validate may open the catalog, so a source's default_distro counts.
		cat, catErr := config.OpenCatalog(cmd.Context())
		baseRecipe := ""
		if catErr == nil {
			baseRecipe = config.BaseRecipeNameFrom(cat)
		}
		switch {
		case baseRecipe == "":
			fmt.Printf("%s No default distro configured and no source declares default_distro\n", utils.StyleError("✗"))
			valid = false
		case !recipeExists(cmd.Context(), baseRecipe):
			fmt.Printf("%s No %s recipe in any configured source\n", utils.StyleError("✗"), baseRecipe)
			valid = false
		default:
			if !utils.QuietMode {
				fmt.Printf("%s Base recipe: %s\n", utils.StyleSuccess("✓"), baseRecipe)
				if def := config.SourceDefaultDistro(cat); def != "" && def+"/base" != baseRecipe {
					fmt.Printf("  %s\n", utils.StyleHint(
						fmt.Sprintf("sources now recommend %s; run `condatainer config set default_distro %s` to switch", def, def)))
				}
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
	Use:   "append [flags] <key> <value>",
	Short: "Append a value to an array config key",
	Long: `Append a value to an array config key (lowest search priority).

` + configLayersHelp,
	Example:           `  condatainer config append sources lab=/shared/lab/recipes`,
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
	Use:   "prepend [flags] <key> <value>",
	Short: "Prepend a value to an array config key (highest priority)",
	Long: `Prepend a value to an array config key (highest search priority).

` + configLayersHelp,
	Example:           `  condatainer config prepend sources lab=/shared/lab/recipes`,
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
	Use:   "remove [flags] <key> [value]",
	Short: "Remove a config key or a value from an array key",
	Long: `Remove a config key entirely, or one value from an array config key.

` + configLayersHelp,
	Example: `  condatainer config remove build.system_apptainer
  condatainer config remove sources lab=/shared/lab/recipes`,
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
			configPath, _, err := config.ResolveWritableConfigPath(setLayer)
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
	configInitCmd.Flags().StringVarP(&initLayer, "layer", "l", "", "Config layer to create in (u/r/e/s; see help)")
	configGetCmd.Flags().StringVarP(&getLayer, "layer", "l", "", "Read from one config layer (u/r/e/s; see help)")
	configSetCmd.Flags().StringVarP(&setLayer, "layer", "l", "", "Config layer to write (u/r/e/s; see help)")
	configAppendCmd.Flags().StringVarP(&setLayer, "layer", "l", "", "Config layer to write (u/r/e/s; see help)")
	configPrependCmd.Flags().StringVarP(&setLayer, "layer", "l", "", "Config layer to write (u/r/e/s; see help)")
	configRemoveCmd.Flags().StringVarP(&setLayer, "layer", "l", "", "Config layer to write (u/r/e/s; see help)")

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
