package config

import (
	"fmt"
	"log/slog"
	"os"
	"os/exec"
	"path/filepath"
	"regexp"
	"strconv"
	"strings"
	"time"

	"github.com/Justype/condatainer/internal/utils"
	"github.com/spf13/viper"
)

// ConfigFilename is the name of the config file
const ConfigFilename = "config"

// ConfigType is the type of config file (yaml, json, toml)
const ConfigType = "yaml"

// configLayers holds per-file viper instances in priority order (user → extra-root → root → system).
// Populated by InitViper(). Used by layerSources() to merge the `sources` key.
var configLayers []*viper.Viper

// ConfigLayerInfo exposes one loaded config layer for per-layer display.
type ConfigLayerInfo struct {
	Path string
	Type string // "user", "extra-root", "root", "system"
	v    *viper.Viper
}

// InConfig reports whether key is explicitly set in this config file (not just a default).
func (l ConfigLayerInfo) InConfig(key string) bool { return l.v.InConfig(key) }

// GetString returns the string value of key from this config layer.
func (l ConfigLayerInfo) GetString(key string) string { return l.v.GetString(key) }

// GetBool returns the bool value of key from this config layer.
func (l ConfigLayerInfo) GetBool(key string) bool { return l.v.GetBool(key) }

// GetStringSlice returns the string slice value of key from this config layer.
func (l ConfigLayerInfo) GetStringSlice(key string) []string { return l.v.GetStringSlice(key) }

// loadedLayers mirrors configLayers but with labels, for per-layer display in config show --sources.
var loadedLayers []ConfigLayerInfo

// GetConfigLayerInfos returns all loaded config layers in priority order.
// Used by config show --sources to display per-layer value annotations.
func GetConfigLayerInfos() []ConfigLayerInfo { return loadedLayers }

// InitViper initializes Viper with the search paths and defaults. A scalar key
// is won by the highest-priority layer that sets it; array keys merge across all
// layers, channels excepted. See the README's Configuration Hierarchy.
func InitViper() error {
	viper.SetConfigName(ConfigFilename)
	viper.SetConfigType(ConfigType)
	viper.SetEnvPrefix("CNT")
	viper.AutomaticEnv()
	setDefaults()

	// Collect config paths in priority order: user > extra-root > root > system
	type configSource struct{ path, label string }
	var sources []configSource
	if userPath, err := GetUserConfigPath(); err == nil {
		sources = append(sources, configSource{userPath, "user"})
	}
	if extraRoot := GetExtraRootDir(); extraRoot != "" {
		sources = append(sources, configSource{filepath.Join(extraRoot, ConfigFilename+"."+ConfigType), "extra-root"})
	}
	if rootPath := GetRootConfigPath(); rootPath != "" {
		sources = append(sources, configSource{rootPath, "app-root"})
	}
	sources = append(sources, configSource{GetSystemConfigPath(), "system"})

	// Load each existing config file.
	// The first found becomes the global viper primary (for scalar keys + env var compat).
	// configLayers[0] reuses the global viper instance to avoid reading the primary file twice.
	// seenPaths prevents loading the same file twice (e.g. CNT_EXTRA_ROOT == CNT_ROOT).
	primaryLoaded := false
	configLayers = nil
	loadedLayers = nil
	seenPaths := make(map[string]bool)
	for _, src := range sources {
		if !fileExists(src.path) {
			continue
		}
		if seenPaths[src.path] {
			continue
		}
		seenPaths[src.path] = true
		if !primaryLoaded {
			viper.SetConfigFile(src.path)
			if err := viper.ReadInConfig(); err != nil {
				return fmt.Errorf("error reading config %s: %w", src.path, err)
			}
			configLayers = append(configLayers, viper.GetViper()) // reuse, no re-read
			loadedLayers = append(loadedLayers, ConfigLayerInfo{src.path, src.label, viper.GetViper()})
			primaryLoaded = true
		} else {
			v := viper.New()
			v.SetConfigFile(src.path)
			if err := v.ReadInConfig(); err != nil {
				return fmt.Errorf("error reading config %s: %w", src.path, err)
			}
			configLayers = append(configLayers, v)
			loadedLayers = append(loadedLayers, ConfigLayerInfo{src.path, src.label, v})
		}
	}
	return nil
}

// setDefaults sets default values for all config keys
func setDefaults() {
	viper.SetDefault("submit_job", true)
	viper.SetDefault("logs_dir", DefaultLogsDir())

	// Build config defaults
	viper.SetDefault("build.system_apptainer", "apptainer")
	viper.SetDefault("build.ncpus", DefaultNcpus)
	viper.SetDefault("build.mem", DefaultMemMB)
	viper.SetDefault("build.time", DefaultBuildTime)
	viper.SetDefault("build.compress_args", ArgsForCompress("zstd-medium"))
	viper.SetDefault("build.block_size", DefaultBlockSize)
	viper.SetDefault("build.data_block_size", DefaultDataBlockSize)
	viper.SetDefault("build.always_submit", false)

	// Scheduler config defaults
	viper.SetDefault("scheduler.bin", "")
	viper.SetDefault("scheduler.timeout", 0) // seconds; 0 = no timeout
	viper.SetDefault("scheduler.account", "")
	viper.SetDefault("scheduler.partition", "")
	viper.SetDefault("scheduler.ncpus", DefaultSchedulerNcpus)
	viper.SetDefault("scheduler.mem", DefaultSchedulerMemMB)
	viper.SetDefault("scheduler.time", DefaultSchedulerTime)

	viper.SetDefault("channels", DefaultChannels())
	viper.SetDefault("autoload_gpu", true)
	viper.SetDefault("nested_run", DefaultNestedRun)
	// "web" = browser notification via dashboard; "terminal" = bell; "both" = terminal + web; "" or "none" = silent
	viper.SetDefault("notification", DefaultNotification)
	viper.SetDefault("metadata_cache_ttl", DefaultCacheTTLDay) // days
	viper.SetDefault("store_gc_grace", DefaultGCGraceDay)      // days
	viper.SetDefault("proxy_perjob", false)
	viper.SetDefault("helper_bind_all", false)
}

// GetUserConfigPath returns the path to the user config file
func GetUserConfigPath() (string, error) {
	userConfigDir, err := os.UserConfigDir()
	if err != nil {
		// Fallback to home directory
		home, err := os.UserHomeDir()
		if err != nil {
			return "", err
		}
		return filepath.Join(home, ".condatainer", ConfigFilename+"."+ConfigType), nil
	}

	return filepath.Join(userConfigDir, "condatainer", ConfigFilename+"."+ConfigType), nil
}

// GetRootConfigPath returns the root config path (CNT_ROOT or standalone layout).
// Returns empty string if the root dir is not set.
func GetRootConfigPath() string {
	if rootDir := GetRootDir(); rootDir != "" {
		return filepath.Join(rootDir, ConfigFilename+"."+ConfigType)
	}
	return ""
}

// GetExtraRootConfigPath returns the extra-root config path ($CNT_EXTRA_ROOT/config.yaml).
// Returns empty string if CNT_EXTRA_ROOT is not set.
func GetExtraRootConfigPath() string {
	if extraRoot := GetExtraRootDir(); extraRoot != "" {
		return filepath.Join(extraRoot, ConfigFilename+"."+ConfigType)
	}
	return ""
}

// GetSystemConfigPath returns the system-wide config path
func GetSystemConfigPath() string {
	return filepath.Join("/etc", "condatainer", ConfigFilename+"."+ConfigType)
}

// ConfigSearchPath represents a config file location with metadata
type ConfigSearchPath struct {
	Path   string // Full path to config file
	Type   string // Type: "root", "extra-root", "user", "system"
	Exists bool   // Whether the file exists
	InUse  bool   // Whether this is the active config file
}

// fileExists checks if a path exists and is a regular file (not a directory)
func fileExists(path string) bool {
	info, err := os.Stat(path)
	if err != nil {
		return false
	}
	return info.Mode().IsRegular()
}

// GetLoadedConfigPaths returns the file paths of all config files actually loaded by
// InitViper, in priority order. These are the files that actively contribute to the
// effective configuration (via scalar precedence + array merging).
func GetLoadedConfigPaths() []string {
	paths := make([]string, 0, len(configLayers))
	for _, v := range configLayers {
		if p := v.ConfigFileUsed(); p != "" {
			paths = append(paths, p)
		}
	}
	return paths
}

// GetConfigSearchPaths returns all config search paths in priority order.
// This matches the order used by InitViper.
// InUse is true for every config file that was loaded (all contribute via merging),
// not just the primary. Duplicate paths (e.g. CNT_EXTRA_ROOT == CNT_ROOT) are skipped.
func GetConfigSearchPaths() []ConfigSearchPath {
	// Build loaded-path set for InUse: all loaded layers contribute, not just primary.
	loadedPaths := make(map[string]bool)
	for _, p := range GetLoadedConfigPaths() {
		loadedPaths[p] = true
	}

	var paths []ConfigSearchPath
	seenPaths := make(map[string]bool)

	add := func(path, typ string) {
		if path == "" || seenPaths[path] {
			return
		}
		seenPaths[path] = true
		paths = append(paths, ConfigSearchPath{
			Path:   path,
			Type:   typ,
			Exists: fileExists(path),
			InUse:  loadedPaths[path],
		})
	}

	// Priority order matches InitViper: user > extra-root > root > system
	if userConfigDir, err := os.UserConfigDir(); err == nil {
		add(filepath.Join(userConfigDir, "condatainer", ConfigFilename+"."+ConfigType), "user")
	}
	if extraRoot := GetExtraRootDir(); extraRoot != "" {
		add(filepath.Join(extraRoot, ConfigFilename+"."+ConfigType), "extra-root")
	}
	if rootPath := GetRootConfigPath(); rootPath != "" {
		add(rootPath, "app-root")
	}
	add(GetSystemConfigPath(), "system")

	return paths
}

// NormalizeConfigLayer expands a config layer shorthand to its full name.
// Returns the full name ("user", "app-root", "extra-root", "system") or the input unchanged.
// "root" and "r" are accepted as aliases for "app-root".
func NormalizeConfigLayer(layer string) string {
	switch layer {
	case "u":
		return "user"
	case "r", "root":
		return "app-root"
	case "e":
		return "extra-root"
	case "s":
		return "system"
	default:
		return layer
	}
}

// GetConfigPathByLayer returns the config path for the specified config layer.
// Supported layers: "user"/"u", "app-root"/"root"/"r", "extra-root"/"e", "system"/"s"
func GetConfigPathByLayer(layer string) (string, error) {
	switch layer {
	case "user", "u":
		return GetUserConfigPath()
	case "app-root", "root", "r":
		if path := GetRootConfigPath(); path != "" {
			return path, nil
		}
		return "", fmt.Errorf("app-root dir not available (not a standalone layout and CNT_ROOT not set)")
	case "extra-root", "e":
		if path := GetExtraRootConfigPath(); path != "" {
			return path, nil
		}
		return "", fmt.Errorf("CNT_EXTRA_ROOT is not set")
	case "system", "s":
		return GetSystemConfigPath(), nil
	default:
		return "", fmt.Errorf("invalid layer '%s': use 'user' (u), 'app-root' (r), 'extra-root' (e), or 'system' (s)", layer)
	}
}

// inferConfigLayer returns "app-root", "extra-root", "system", or "user" by comparing path
// against known config file paths. Returns "user" for any unrecognised path.
func inferConfigLayer(path string) string {
	if p := GetExtraRootConfigPath(); p != "" && path == p {
		return "extra-root"
	}
	if p := GetRootConfigPath(); p != "" && path == p {
		return "app-root"
	}
	if path == GetSystemConfigPath() {
		return "system"
	}
	return "user"
}

// ResolveWritableConfigPath determines where a config-mutating command should write.
// If layer is non-empty it is resolved and checked for writability.
// If layer is empty the active config file is used when writable; otherwise falls back
// to the user config path, mirroring the config init smart-default behaviour.
//
// Writability is checked on the config file itself, not just its directory: a
// read-only file in a writable directory (a frozen shared config) must not be
// treated as a valid write target.
func ResolveWritableConfigPath(layer string) (path, layerType string, err error) {
	path, layerType, _, err = resolveWritableConfigPath(layer)
	return path, layerType, err
}

// ResolveWritableConfigPathVerbose is ResolveWritableConfigPath plus the layer it
// fell back from. fellBackFrom is empty unless an active config existed but was
// read-only, forcing the write down to the user layer — callers report that, since
// the change then applies only to the current user.
func ResolveWritableConfigPathVerbose(layer string) (path, layerType, fellBackFrom string, err error) {
	return resolveWritableConfigPath(layer)
}

func resolveWritableConfigPath(layer string) (path, layerType, fellBackFrom string, err error) {
	if layer != "" {
		path, err = GetConfigPathByLayer(layer)
		if err != nil {
			return "", "", "", err
		}
		layerType = NormalizeConfigLayer(layer)
		if !utils.CanWriteToFile(path) {
			return "", "", "", fmt.Errorf(
				"config layer '%s' is read-only: %s\nuse a different layer or run with appropriate permissions",
				layerType, path,
			)
		}
		return path, layerType, "", nil
	}

	// Auto-detect: active config file when writable, else user config
	if active := viper.ConfigFileUsed(); active != "" {
		if utils.CanWriteToFile(active) {
			return active, inferConfigLayer(active), "", nil
		}
		fellBackFrom = inferConfigLayer(active)
	}
	path, err = GetUserConfigPath()
	if err != nil {
		return "", "", "", fmt.Errorf("failed to determine user config path: %w", err)
	}
	if fellBackFrom == "user" {
		// The user's own config is unwritable; falling back to itself would loop.
		return "", "", "", fmt.Errorf(
			"user config is read-only: %s\nfix its permissions to change settings", path,
		)
	}
	return path, "user", fellBackFrom, nil
}

// ResolveReadableConfigPath resolves a config layer name to a config file path for reading.
// Unlike ResolveWritableConfigPath it does not check write permissions.
func ResolveReadableConfigPath(layer string) (path, layerType string, err error) {
	if layer == "" {
		path, err = GetUserConfigPath()
		if err != nil {
			return "", "", fmt.Errorf("failed to determine user config path: %w", err)
		}
		return path, "user", nil
	}
	path, err = GetConfigPathByLayer(layer)
	if err != nil {
		return "", "", err
	}
	return path, NormalizeConfigLayer(layer), nil
}

// ReadConfigSliceKey reads a string-slice key from a config file without touching
// global viper state. Returns nil if the file doesn't exist or the key is not set.
func ReadConfigSliceKey(configPath, key string) []string {
	if configPath == "" {
		return nil
	}
	v := viper.New()
	v.SetConfigFile(configPath)
	if err := v.ReadInConfig(); err != nil {
		return nil
	}
	return v.GetStringSlice(key)
}

// UpdateConfigKey reads path into a fresh viper instance, updates key to value,
// and writes the result back to path. No global viper state is affected.
// Creates the directory and file if needed (first-use creation).
func UpdateConfigKey(path, key string, value any) error {
	dir := filepath.Dir(path)
	if err := utils.MkdirAllShared(dir); err != nil {
		return fmt.Errorf("failed to create config directory: %w", err)
	}
	v := viper.New()
	v.SetConfigType(ConfigType)
	v.SetConfigFile(path)
	_ = v.ReadInConfig() // ignore error: file may not exist yet
	v.Set(key, value)
	if err := v.WriteConfigAs(path); err != nil {
		return fmt.Errorf("failed to write config to %s: %w", path, err)
	}
	utils.ShareWithParentGroup(path)
	return nil
}

// ReadConfigAllKeys returns all keys explicitly set in the config file at path.
func ReadConfigAllKeys(path string) []string {
	v := viper.New()
	v.SetConfigType(ConfigType)
	v.SetConfigFile(path)
	if err := v.ReadInConfig(); err != nil {
		return nil
	}
	return v.AllKeys()
}

// ReadConfigKey reads a single key from a config file without touching global viper state.
// Returns "" if the file doesn't exist, can't be read, or the key is not set.
func ReadConfigKey(configPath, key string) string {
	if configPath == "" {
		return ""
	}
	v := viper.New()
	v.SetConfigFile(configPath)
	if err := v.ReadInConfig(); err != nil {
		return ""
	}
	return v.GetString(key)
}

// SaveMinimalConfigTo writes only detected keys to path using a fresh viper instance,
// so no default values bleed in. Keys already provided by lowerLayers are skipped.
func SaveMinimalConfigTo(path, apptainerBin, schedulerBin, compressArgs string, lowerLayers []ConfigLayerInfo) error {
	dir := filepath.Dir(path)
	if err := utils.MkdirAllShared(dir); err != nil {
		return fmt.Errorf("failed to create config directory: %w", err)
	}
	alreadySet := func(key string) bool {
		for _, l := range lowerLayers {
			if l.InConfig(key) {
				return true
			}
		}
		return false
	}
	v := viper.New()
	v.SetConfigType(ConfigType)
	if apptainerBin != "" && !alreadySet("build.system_apptainer") {
		v.Set("build.system_apptainer", apptainerBin)
	}
	if schedulerBin != "" && !alreadySet("scheduler.bin") {
		v.Set("scheduler.bin", schedulerBin)
	}
	if compressArgs != "" && !alreadySet("build.compress_args") {
		v.Set("build.compress_args", compressArgs)
	}
	if err := v.WriteConfigAs(path); err != nil {
		return fmt.Errorf("failed to write config to %s: %w", path, err)
	}
	utils.ShareWithParentGroup(path)
	return nil
}

// SetConfigKey sets a single key in the config file at path, creating it if needed.
// Existing keys are preserved. No-op if the file already has the key set to the same value.
func SetConfigKey(path, key, value string) error {
	v := viper.New()
	v.SetConfigType(ConfigType)
	v.SetConfigFile(path)
	if err := v.ReadInConfig(); err != nil && !os.IsNotExist(err) {
		return fmt.Errorf("failed to read config %s: %w", path, err)
	}
	if v.GetString(key) == value {
		return nil // already set to this value
	}
	v.Set(key, value)
	dir := filepath.Dir(path)
	if err := utils.MkdirAllShared(dir); err != nil {
		return fmt.Errorf("failed to create config directory: %w", err)
	}
	if err := v.WriteConfigAs(path); err != nil {
		return fmt.Errorf("failed to write config %s: %w", path, err)
	}
	utils.ShareWithParentGroup(path)
	return nil
}

// DeleteConfigKey removes a key from the config file at path.
// Nested keys (e.g. "build.compress_args") are handled by navigating the settings map.
// If removing a nested key leaves the parent empty, the parent is removed too.
func DeleteConfigKey(path, key string) error {
	v := viper.New()
	v.SetConfigType(ConfigType)
	v.SetConfigFile(path)
	if err := v.ReadInConfig(); err != nil {
		if os.IsNotExist(err) {
			return nil
		}
		return fmt.Errorf("failed to read config %s: %w", path, err)
	}
	all := v.AllSettings()
	deleteNestedKey(all, strings.Split(key, "."))
	v2 := viper.New()
	v2.SetConfigType(ConfigType)
	for k, val := range all {
		v2.Set(k, val)
	}
	if err := v2.WriteConfigAs(path); err != nil {
		return fmt.Errorf("failed to write config %s: %w", path, err)
	}
	utils.ShareWithParentGroup(path)
	return nil
}

func deleteNestedKey(m map[string]any, parts []string) {
	if len(parts) == 0 {
		return
	}
	if len(parts) == 1 {
		delete(m, parts[0])
		return
	}
	if sub, ok := m[parts[0]].(map[string]any); ok {
		deleteNestedKey(sub, parts[1:])
		if len(sub) == 0 {
			delete(m, parts[0])
		}
	}
}

// ValidateBinary checks if a binary exists and is executable
func ValidateBinary(binPath string) bool {
	if binPath == "" {
		return false
	}

	// If it's a full path, check directly
	if filepath.IsAbs(binPath) {
		info, err := os.Stat(binPath)
		if err != nil {
			return false
		}
		// Check if it's executable (unix-style check)
		return info.Mode()&0111 != 0
	}

	// Otherwise, try to find it in PATH
	_, err := exec.LookPath(binPath)
	return err == nil
}

// FindApptainerBin returns the apptainer (or singularity) on PATH, else the
// newest one a module provides. "" when there is none.
func FindApptainerBin() string {
	if path := apptainerOnPath(); path != "" {
		return path
	}
	return detectApptainerFromModules()
}

// InvalidSystemApptainer returns the configured build.system_apptainer when it
// is set but not a usable binary, "" otherwise.
func InvalidSystemApptainer() string {
	if bin := layerString("build.system_apptainer"); bin != "" && !ValidateBinary(bin) {
		return bin
	}
	return ""
}

func detectApptainerFromModules() string {
	slog.Default().Info("searching modules for apptainer/singularity via 'module avail'")

	bestModule := ""
	bestVersion := ""

	for _, moduleName := range []string{"apptainer", "singularity"} {
		module, version := detectLatestModule(moduleName)
		if module == "" {
			continue
		}

		if bestModule == "" || compareModuleVersion(version, bestVersion) > 0 {
			bestModule = module
			bestVersion = version
		}
	}

	if bestModule == "" {
		slog.Default().Debug("no apptainer/singularity modules found via module avail")
		return ""
	}

	slog.Default().Info("using module candidate", "module", bestModule)

	// Resolve the actual binary path after loading the selected module.
	// Use both names as fallback because module name and binary name can differ.
	cmdStr := fmt.Sprintf("module load %q >/dev/null 2>&1 && (command -v apptainer || command -v singularity)", bestModule)
	out, err := exec.Command("bash", "-lc", cmdStr).Output()
	if err != nil {
		return ""
	}

	resolved := strings.TrimSpace(string(out))
	if resolved == "" {
		return ""
	}

	line := strings.Split(resolved, "\n")[0]
	line = strings.TrimSpace(line)
	if line == "" {
		return ""
	}

	if filepath.IsAbs(line) {
		if info, statErr := os.Stat(line); statErr == nil && info.Mode()&0111 != 0 {
			return line
		}
	}

	if fullPath, lookErr := exec.LookPath(line); lookErr == nil {
		return fullPath
	}

	return ""
}

func detectLatestModule(moduleName string) (module string, version string) {
	cmdStr := fmt.Sprintf("module -t avail %s 2>&1 || true", moduleName)
	out, err := exec.Command("bash", "-lc", cmdStr).Output()
	if err != nil {
		return "", ""
	}

	module, version = parseLatestModuleFromAvailOutput(moduleName, string(out))
	return module, version
}

func parseLatestModuleFromAvailOutput(moduleName, output string) (module string, version string) {
	lineRe := regexp.MustCompile(`^` + regexp.QuoteMeta(moduleName) + `(?:/([^\s()]+))?(?:\([^)]*\))?$`)

	bestModule := ""
	bestVersion := ""

	for _, rawLine := range strings.Split(output, "\n") {
		line := strings.TrimSpace(rawLine)
		if line == "" {
			continue
		}

		if strings.Contains(line, ":") || strings.HasPrefix(line, "-") || strings.HasPrefix(line, "Lmod") {
			continue
		}

		fields := strings.Fields(line)
		if len(fields) == 0 {
			continue
		}

		for _, token := range fields {
			token = strings.TrimSpace(token)
			token = strings.TrimSuffix(token, "*")

			match := lineRe.FindStringSubmatch(token)
			if len(match) == 0 {
				continue
			}

			moduleVersion := ""
			if len(match) > 1 {
				moduleVersion = match[1]
			}

			cleanModule := moduleName
			if moduleVersion != "" {
				cleanModule = moduleName + "/" + moduleVersion
			}

			if bestModule == "" || compareModuleVersion(moduleVersion, bestVersion) > 0 {
				bestModule = cleanModule
				bestVersion = moduleVersion
			}
		}
	}

	return bestModule, bestVersion
}

func compareModuleVersion(a, b string) int {
	if a == b {
		return 0
	}
	if a == "" {
		return -1
	}
	if b == "" {
		return 1
	}

	partsA := splitVersionParts(a)
	partsB := splitVersionParts(b)
	maxLen := len(partsA)
	if len(partsB) > maxLen {
		maxLen = len(partsB)
	}

	for i := 0; i < maxLen; i++ {
		va := ""
		vb := ""
		if i < len(partsA) {
			va = partsA[i]
		}
		if i < len(partsB) {
			vb = partsB[i]
		}

		na, errA := strconv.Atoi(va)
		nb, errB := strconv.Atoi(vb)

		switch {
		case errA == nil && errB == nil:
			if na > nb {
				return 1
			}
			if na < nb {
				return -1
			}
		case errA == nil && errB != nil:
			return 1
		case errA != nil && errB == nil:
			return -1
		default:
			if va > vb {
				return 1
			}
			if va < vb {
				return -1
			}
		}
	}

	return 0
}

func splitVersionParts(v string) []string {
	parts := strings.FieldsFunc(v, func(r rune) bool {
		switch r {
		case '.', '-', '_', '+':
			return true
		default:
			return false
		}
	})

	if len(parts) == 0 {
		return []string{v}
	}

	return parts
}

// DetectSchedulerBin attempts to find scheduler binary
// Returns the binary path if found, empty string otherwise
func DetectSchedulerBin() string {
	// Try SLURM first (most common in HPC)
	if path, err := exec.LookPath("sbatch"); err == nil {
		return path
	}

	// Try PBS
	if path, err := exec.LookPath("qsub"); err == nil {
		return path
	}

	// Try LSF
	if path, err := exec.LookPath("bsub"); err == nil {
		return path
	}

	return ""
}

// GetSchedulerTypeFromBin derives the scheduler type from the binary path/name.
// The type is always inferred from the binary - it cannot be set independently.
func GetSchedulerTypeFromBin(binPath string) string {
	if binPath == "" {
		return ""
	}

	baseName := filepath.Base(binPath)

	switch baseName {
	case "sbatch", "srun", "salloc", "scancel", "squeue":
		return "SLURM"
	case "qsub", "qdel", "qstat":
		// SGE also uses qsub, check for SGE-specific env
		if _, exists := os.LookupEnv("SGE_ROOT"); exists {
			return "SGE"
		}
		return "PBS"
	case "bsub", "bjobs", "bkill":
		return "LSF"
	default:
		return ""
	}
}

// layerString returns the value of a scalar string key from the first source that
// explicitly sets it. Priority: env var > user > extra-root > root > system.
// Returns "" if no source explicitly sets the key.
func layerString(key string) string {
	envKey := "CNT_" + strings.ToUpper(strings.ReplaceAll(key, ".", "_"))
	if ev := os.Getenv(envKey); ev != "" {
		return ev
	}
	for _, v := range configLayers {
		if v.InConfig(key) {
			return v.GetString(key)
		}
	}
	return ""
}

// layerBool returns the value of a scalar bool key and whether it was explicitly set.
// Priority: env var > user > extra-root > root > system.
func layerBool(key string) (bool, bool) {
	envKey := "CNT_" + strings.ToUpper(strings.ReplaceAll(key, ".", "_"))
	if ev := os.Getenv(envKey); ev != "" {
		if val, err := strconv.ParseBool(ev); err == nil {
			return val, true
		}
	}
	for _, v := range configLayers {
		if v.InConfig(key) {
			return v.GetBool(key), true
		}
	}
	return false, false
}

// layerStringSet returns the value of a scalar string key and whether it was
// explicitly set. Priority: env var > user > extra-root > root > system.
//
// Needed where empty is itself a legal value — notification, where it means
// silent — since layerString returns "" for unset and for explicitly empty alike.
func layerStringSet(key string) (string, bool) {
	envKey := "CNT_" + strings.ToUpper(strings.ReplaceAll(key, ".", "_"))
	if ev := os.Getenv(envKey); ev != "" {
		return ev, true
	}
	for _, v := range configLayers {
		if v.InConfig(key) {
			return v.GetString(key), true
		}
	}
	return "", false
}

// layerInt returns the value of a scalar int key and whether it was explicitly set.
// Priority: env var > user > extra-root > root > system.
func layerInt(key string) (int, bool) {
	envKey := "CNT_" + strings.ToUpper(strings.ReplaceAll(key, ".", "_"))
	if ev := os.Getenv(envKey); ev != "" {
		if n, err := strconv.Atoi(ev); err == nil {
			return n, true
		}
	}
	for _, v := range configLayers {
		if v.InConfig(key) {
			return v.GetInt(key), true
		}
	}
	return 0, false
}

// LoadFromViper loads config from Viper into Global struct
func LoadFromViper() {
	// A valid configured apptainer overrides the one LoadDefaults found on PATH.
	if bin := layerString("build.system_apptainer"); bin != "" && ValidateBinary(bin) {
		Global.Build.SystemApptainer = bin
	}

	if bin := layerString("scheduler.bin"); bin != "" {
		Global.Scheduler.Bin = bin
	}

	// Handle submit_job: disable if scheduler is not accessible
	if submitJob, ok := layerBool("submit_job"); ok && !submitJob {
		Global.SubmitJob = false
	} else {
		// Auto-disable if no scheduler binary is available
		if Global.Scheduler.Bin == "" {
			Global.Scheduler.Bin = DetectSchedulerBin()
		}
		if Global.Scheduler.Bin == "" || !ValidateBinary(Global.Scheduler.Bin) {
			Global.SubmitJob = false
		}
	}

	if timeout, ok := layerInt("scheduler.timeout"); ok {
		Global.Scheduler.Timeout = time.Duration(timeout) * time.Second
	}

	if account := layerString("scheduler.account"); account != "" {
		Global.Scheduler.Account = account
	}

	if partition := layerString("scheduler.partition"); partition != "" {
		Global.Scheduler.Partition = partition
	}

	if ncpus, ok := layerInt("scheduler.ncpus"); ok && ncpus > 0 {
		Global.Scheduler.Defaults.CpusPerTask = ncpus
	}

	if memStr := layerString("scheduler.mem"); memStr != "" {
		if memMB, err := utils.ParseMemoryMB(memStr); err == nil && memMB > 0 {
			Global.Scheduler.Defaults.MemPerNodeMB = memMB
		}
	}

	if schedTime := layerString("scheduler.time"); schedTime != "" {
		if dur, err := utils.ParseWalltime(schedTime); err == nil {
			Global.Scheduler.Defaults.Time = dur
		}
	}

	// Load logs_dir from config (overrides default $HOME/logs)
	if logsDir := layerString("logs_dir"); logsDir != "" {
		logsDir = os.ExpandEnv(logsDir)
		if absLogsDir, err := filepath.Abs(logsDir); err == nil {
			logsDir = absLogsDir
		}
		Global.LogsDir = logsDir
	}

	// Recipe collections, in order; earlier entries shadow later ones.
	Global.Sources = layerSources()
	Global.DefaultDistro = layerString("default_distro")

	// Load build config from Viper
	if ncpus, ok := layerInt("build.ncpus"); ok && ncpus > 0 {
		Global.Build.Defaults.CpusPerTask = ncpus
	}

	if memStr := layerString("build.mem"); memStr != "" {
		if memMB, err := utils.ParseMemoryMB(memStr); err == nil && memMB > 0 {
			Global.Build.Defaults.MemPerNodeMB = memMB
		}
	}

	if buildTime := layerString("build.time"); buildTime != "" {
		if dur, err := utils.ParseWalltime(buildTime); err == nil {
			Global.Build.Defaults.Time = dur
		}
	}

	// Only override compress_args if explicitly set in config; otherwise
	// LoadDefaults' unconditional zstd-medium stands.
	if compressArgs := layerString("build.compress_args"); compressArgs != "" {
		Global.Build.CompressArgs = NormalizeCompressArgs(compressArgs)
	}

	if v := layerString("build.block_size"); v != "" {
		if IsValidBlockSize(v) {
			Global.Build.BlockSize = v
		} else {
			slog.Default().Warn("invalid build.block_size, using default", "value", v, "default", DefaultBlockSize)
			Global.Build.BlockSize = DefaultBlockSize
		}
	}
	if v := layerString("build.data_block_size"); v != "" {
		if IsValidBlockSize(v) {
			Global.Build.DataBlockSize = v
		} else {
			slog.Default().Warn("invalid build.data_block_size, using default", "value", v, "default", DefaultDataBlockSize)
			Global.Build.DataBlockSize = DefaultDataBlockSize
		}
	}

	if alwaysSubmit, ok := layerBool("build.always_submit"); ok {
		Global.Build.AlwaysSubmit = alwaysSubmit
	}

	if ch := GetChannels(); len(ch) > 0 {
		Global.Build.Channels = ch
	}

	if autoloadGPU, ok := layerBool("autoload_gpu"); ok {
		Global.AutoloadGPU = autoloadGPU
	}

	if v, ok := layerStringSet("nested_run"); ok {
		if normalized, valid := ParseNestedRun(v); valid {
			Global.NestedRun = normalized
		} else {
			fmt.Fprintf(os.Stderr, "[WARN] Unknown nested_run value %q. Valid values: auto, true, false. Using auto.\n", v)
			Global.NestedRun = DefaultNestedRun
		}
	}

	// Only when set: "" means silent, so an unset key must keep DefaultNotification.
	if v, ok := layerStringSet("notification"); ok {
		Global.Notification = v
		switch Global.Notification {
		case "", "none", "terminal", "web", "both":
			// valid
		default:
			fmt.Fprintf(os.Stderr, "[WARN] Unknown notification value %q. Valid values: terminal, web, both, none. Treating as none.\n", Global.Notification)
			Global.Notification = ""
		}
	}

	if ttl, ok := layerInt("metadata_cache_ttl"); ok {
		Global.MetadataCacheTTL = time.Duration(ttl) * 24 * time.Hour
	}

	if grace, ok := layerInt("store_gc_grace"); ok {
		Global.StoreGCGrace = time.Duration(grace) * 24 * time.Hour
	}

	if proxyPerJob, ok := layerBool("proxy_perjob"); ok {
		Global.ProxyPerJob = proxyPerJob
	}

	if helperBindAll, ok := layerBool("helper_bind_all"); ok {
		Global.HelperBindAll = helperBindAll
	}

}

// NormalizeCompressArgs is a thin wrapper around ArgsForCompress and exists
// for historical compatibility with earlier versions of the code.
func NormalizeCompressArgs(val string) string {
	return ArgsForCompress(val)
}
