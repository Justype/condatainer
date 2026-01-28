package config

import (
	"fmt"
	"os"
	"os/exec"
	"path/filepath"

	"github.com/Justype/condatainer/internal/utils"
	"github.com/spf13/viper"
)

// ConfigFilename is the name of the config file
const ConfigFilename = "config"

// ConfigType is the type of config file (yaml, json, toml)
const ConfigType = "yaml"

// InitViper initializes Viper with proper search paths and defaults
// Priority (highest to lowest):
// 1. Command-line flags (handled by cobra)
// 2. Environment variables (CONDATAINER_*)
// 3. Portable config (parent of bin/ folder where executable lives)
// 4. User config file (~/.config/condatainer/config.yaml)
// 5. User config fallback (~/.condatainer/config.yaml)
// 6. System config file (/etc/condatainer/config.yaml)
// 7. Current directory (for development)
// 8. Defaults
func InitViper() error {
	viper.SetConfigName(ConfigFilename)
	viper.SetConfigType(ConfigType)

	// Set config search paths (order matters - first found wins)

	// Portable config: if executable is in a "bin" folder, check parent folder
	// e.g., /path/to/condatainer/bin/condatainer -> /path/to/condatainer/config.yaml
	// Excludes common user bin directories ($HOME/bin, $HOME/.local/bin, /usr/bin, etc.)
	if exePath, err := os.Executable(); err == nil {
		exeDir := filepath.Dir(exePath)
		if filepath.Base(exeDir) == "bin" {
			parentDir := filepath.Dir(exeDir)
			// Exclude common bin directories that aren't dedicated installations
			home, _ := os.UserHomeDir()
			excludedParents := []string{
				home,                          // $HOME/bin
				filepath.Join(home, ".local"), // $HOME/.local/bin
				"/usr",                         // /usr/bin
				"/usr/local",                   // /usr/local/bin
				"/opt",                         // /opt/bin
			}
			isExcluded := false
			for _, excluded := range excludedParents {
				if parentDir == excluded {
					isExcluded = true
					break
				}
			}
			if !isExcluded {
				viper.AddConfigPath(parentDir)
			}
		}
	}

	// User config (highest priority for user-specific settings)
	if userConfigDir, err := os.UserConfigDir(); err == nil {
		viper.AddConfigPath(filepath.Join(userConfigDir, "condatainer"))
	}

	// Home directory fallback
	if home, err := os.UserHomeDir(); err == nil {
		viper.AddConfigPath(filepath.Join(home, ".condatainer"))
	}

	// System-wide config (lower priority)
	viper.AddConfigPath("/etc/condatainer")

	// Current directory (for development)
	viper.AddConfigPath(".")

	// Environment variables
	viper.SetEnvPrefix("CONDATAINER")
	viper.AutomaticEnv()

	// Set defaults (lowest priority)
	setDefaults()

	// Read config file (non-fatal if not found)
	if err := viper.ReadInConfig(); err != nil {
		if _, ok := err.(viper.ConfigFileNotFoundError); ok {
			// Config file not found; will use defaults and auto-detect
			return nil
		}
		return fmt.Errorf("error reading config file: %w", err)
	}

	return nil
}

// setDefaults sets default values for all config keys
func setDefaults() {
	// Base directory: empty means use auto-detection ($SCRATCH/condatainer or $HOME/condatainer)
	viper.SetDefault("base_dir", "")

	viper.SetDefault("apptainer_bin", "apptainer")
	viper.SetDefault("scheduler_bin", "")
	viper.SetDefault("scheduler_type", "")
	viper.SetDefault("submit_job", true)
	viper.SetDefault("logs_dir", "")

	// Build config defaults
	viper.SetDefault("build.default_cpus", 4)
	viper.SetDefault("build.default_mem_mb", 8192)
	viper.SetDefault("build.default_time", "2h")
	viper.SetDefault("build.tmp_size_mb", 20480)
	viper.SetDefault("build.compress_args", "") // Empty means auto-detect based on apptainer version
	viper.SetDefault("build.overlay_type", "ext3")
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

// GetPortableConfigPath returns the portable config path if the executable is in a
// dedicated installation (bin/ folder that's not a common user/system bin).
// Returns empty string if not in a portable installation.
func GetPortableConfigPath() string {
	exePath, err := os.Executable()
	if err != nil {
		return ""
	}

	exeDir := filepath.Dir(exePath)
	if filepath.Base(exeDir) != "bin" {
		return ""
	}

	parentDir := filepath.Dir(exeDir)

	// Exclude common bin directories that aren't dedicated installations
	home, _ := os.UserHomeDir()
	excludedParents := []string{
		home,                          // $HOME/bin
		filepath.Join(home, ".local"), // $HOME/.local/bin
		"/usr",                        // /usr/bin
		"/usr/local",                  // /usr/local/bin
		"/opt",                        // /opt/bin
	}

	for _, excluded := range excludedParents {
		if parentDir == excluded {
			return ""
		}
	}

	return filepath.Join(parentDir, ConfigFilename+"."+ConfigType)
}

// GetDefaultConfigPath returns the recommended config path based on installation type.
// - If in a portable installation (bin/ folder), returns the portable path
// - Otherwise, returns the user config path (~/.config/condatainer/config.yaml)
func GetDefaultConfigPath() (string, error) {
	// Check for portable installation first
	if portablePath := GetPortableConfigPath(); portablePath != "" {
		return portablePath, nil
	}

	// Fall back to user config
	return GetUserConfigPath()
}

// IsPortableInstallation returns true if the executable is in a portable installation
func IsPortableInstallation() bool {
	return GetPortableConfigPath() != ""
}

// GetSystemConfigPath returns the system-wide config path
func GetSystemConfigPath() string {
	return filepath.Join("/etc", "condatainer", ConfigFilename+"."+ConfigType)
}

// GetConfigPathByLocation returns the config path for the specified location type.
// Supported locations: "user"/"u", "portable"/"p", "system"/"s"
func GetConfigPathByLocation(location string) (string, error) {
	switch location {
	case "user", "u":
		return GetUserConfigPath()
	case "portable", "p":
		if path := GetPortableConfigPath(); path != "" {
			return path, nil
		}
		return "", fmt.Errorf("not in a portable installation (executable not in a dedicated bin/ folder)")
	case "system", "s":
		return GetSystemConfigPath(), nil
	default:
		return "", fmt.Errorf("invalid location '%s': use 'user' (u), 'portable' (p), or 'system' (s)", location)
	}
}

// SaveConfig saves current Viper config to user config file
func SaveConfig() error {
	configPath, err := GetUserConfigPath()
	if err != nil {
		return fmt.Errorf("failed to get config path: %w", err)
	}
	return SaveConfigTo(configPath)
}

// SaveConfigTo saves current Viper config to the specified path
func SaveConfigTo(configPath string) error {
	// Create directory if it doesn't exist
	configDir := filepath.Dir(configPath)
	if err := os.MkdirAll(configDir, 0755); err != nil {
		return fmt.Errorf("failed to create config directory: %w", err)
	}

	// Write config file
	if err := viper.WriteConfigAs(configPath); err != nil {
		return fmt.Errorf("failed to write config file: %w", err)
	}

	return nil
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

// DetectApptainerBin attempts to find apptainer binary
// Returns the full absolute path if found, empty string otherwise
func DetectApptainerBin() string {
	// Try common names in PATH (returns full path)
	candidates := []string{"apptainer", "singularity"}
	for _, candidate := range candidates {
		if path, err := exec.LookPath(candidate); err == nil {
			// exec.LookPath already returns the full path
			return path
		}
	}

	return ""
}

// DetectSchedulerBin attempts to find scheduler binary
// Returns (binary_path, scheduler_type) if found
func DetectSchedulerBin() (string, string) {
	// Try SLURM first (most common in HPC)
	if path, err := exec.LookPath("sbatch"); err == nil {
		return path, "SLURM"
	}

	// Try PBS/Torque
	if path, err := exec.LookPath("qsub"); err == nil {
		return path, "PBS"
	}

	// Try LSF
	if path, err := exec.LookPath("bsub"); err == nil {
		return path, "LSF"
	}

	// Try SGE (Sun Grid Engine)
	if path, err := exec.LookPath("qsub"); err == nil {
		// SGE also uses qsub, check for SGE-specific env
		if _, exists := os.LookupEnv("SGE_ROOT"); exists {
			return path, "SGE"
		}
	}

	return "", ""
}

// AutoDetectAndSave auto-detects binaries and saves to config if needed
// Returns true if config was updated
func AutoDetectAndSave() (bool, error) {
	updated := false

	if !IsInsideContainer() {
		apptainerBin := viper.GetString("apptainer_bin")
		if !ValidateBinary(apptainerBin) {
			detected := DetectApptainerBin()
			if detected != "" {
				viper.Set("apptainer_bin", detected)
				updated = true
			}
		}
	}

	// Check and detect scheduler binary
	schedulerBin := viper.GetString("scheduler_bin")
	if !ValidateBinary(schedulerBin) {
		detectedBin, detectedType := DetectSchedulerBin()
		if detectedBin != "" {
			viper.Set("scheduler_bin", detectedBin)
			viper.Set("scheduler_type", detectedType)
			updated = true
		}
	}

	// Save if anything was updated
	if updated {
		if err := SaveConfig(); err != nil {
			return false, err
		}
	}

	return updated, nil
}

// ForceDetectAndSave always re-detects binaries from current environment and saves
// This is useful for config init to capture the exact paths from current PATH
// Returns true if config was updated
func ForceDetectAndSave() (bool, error) {
	configPath, err := GetUserConfigPath()
	if err != nil {
		return false, err
	}
	return ForceDetectAndSaveTo(configPath)
}

// ForceDetectAndSaveTo always re-detects binaries and saves to the specified path
// Returns true if config was updated
func ForceDetectAndSaveTo(configPath string) (bool, error) {
	updated := false

	if !IsInsideContainer() {
		// Always re-detect apptainer binary
		detected := DetectApptainerBin()
		if detected != "" {
			currentBin := viper.GetString("apptainer_bin")
			if currentBin != detected {
				viper.Set("apptainer_bin", detected)
				updated = true
			}
		}
	}

	// Always re-detect scheduler binary
	detectedBin, detectedType := DetectSchedulerBin()
	if detectedBin != "" {
		currentBin := viper.GetString("scheduler_bin")
		currentType := viper.GetString("scheduler_type")
		if currentBin != detectedBin || currentType != detectedType {
			viper.Set("scheduler_bin", detectedBin)
			viper.Set("scheduler_type", detectedType)
			updated = true
		}
	}

	// Always save (even if nothing changed, to create the file)
	if err := SaveConfigTo(configPath); err != nil {
		return false, err
	}

	return updated, nil
}

// SetCompressArgsInConfig saves the compress_args to config file
// This is called during config init to persist the auto-detected compression
func SetCompressArgsInConfig(compressArgs string) error {
	configPath, err := GetUserConfigPath()
	if err != nil {
		return err
	}
	return SetCompressArgsInConfigTo(compressArgs, configPath)
}

// SetCompressArgsInConfigTo saves the compress_args to the specified config file
func SetCompressArgsInConfigTo(compressArgs, configPath string) error {
	viper.Set("build.compress_args", compressArgs)
	return SaveConfigTo(configPath)
}

// LoadFromViper loads config from Viper into Global struct
func LoadFromViper() {
	// Load base_dir from config if set (overrides auto-detection)
	if baseDir := viper.GetString("base_dir"); baseDir != "" {
		// Expand environment variables in the path
		baseDir = os.ExpandEnv(baseDir)
		if absBaseDir, err := filepath.Abs(baseDir); err == nil {
			baseDir = absBaseDir
		}
		// Update all derived paths
		Global.BaseDir = baseDir
		Global.ImagesDir = filepath.Join(baseDir, "images")
		Global.BuildScriptsDir = filepath.Join(baseDir, "build-scripts")
		Global.HelperScriptsDir = filepath.Join(baseDir, "helper-scripts")
		Global.TmpDir = filepath.Join(baseDir, "tmp")
		Global.BaseImage = filepath.Join(baseDir, "images", "base_image.sif")
	}

	// Update binary paths from Viper, with fallback to detection
	if bin := viper.GetString("apptainer_bin"); bin != "" && ValidateBinary(bin) {
		Global.ApptainerBin = bin
	} else if bin == "" || !ValidateBinary(bin) {
		// Fallback to detection if config value is empty or invalid
		detected := detectApptainerBin()
		if detected != "" {
			Global.ApptainerBin = detected
		}
	}

	if bin := viper.GetString("scheduler_bin"); bin != "" {
		Global.SchedulerBin = bin
	}

	// Handle submit_job: disable if scheduler is not accessible
	if submitJob := viper.GetBool("submit_job"); !submitJob {
		Global.SubmitJob = false
	} else {
		// Auto-disable if no scheduler binary is available
		schedulerBin := Global.SchedulerBin
		if schedulerBin == "" {
			schedulerBin, _ = DetectSchedulerBin()
		}
		if schedulerBin == "" || !ValidateBinary(schedulerBin) {
			Global.SubmitJob = false
		}
	}

	// Load logs_dir from config (overrides default $HOME/logs)
	if logsDir := viper.GetString("logs_dir"); logsDir != "" {
		logsDir = os.ExpandEnv(logsDir)
		if absLogsDir, err := filepath.Abs(logsDir); err == nil {
			logsDir = absLogsDir
		}
		Global.LogsDir = logsDir
	}

	// Load build config from Viper
	if defaultCPUs := viper.GetInt("build.default_cpus"); defaultCPUs > 0 {
		Global.Build.DefaultCPUs = defaultCPUs
	}

	if defaultMemMB := viper.GetInt64("build.default_mem_mb"); defaultMemMB > 0 {
		Global.Build.DefaultMemMB = defaultMemMB
	}

	if defaultTime := viper.GetString("build.default_time"); defaultTime != "" {
		// Parse time duration from string (e.g., "2h", "30m", "1h30m", or "02:00:00")
		if dur, err := utils.ParseDuration(defaultTime); err == nil {
			Global.Build.DefaultTime = dur
		}
	}

	if tmpSizeMB := viper.GetInt("build.tmp_size_mb"); tmpSizeMB > 0 {
		Global.Build.TmpSizeMB = tmpSizeMB
	}

	// Only override compress_args if explicitly set in config (non-empty)
	// Empty means auto-detect based on apptainer version (handled by AutoDetectCompression)
	if compressArgs := viper.GetString("build.compress_args"); compressArgs != "" {
		Global.Build.CompressArgs = compressArgs
	}

	if overlayType := viper.GetString("build.overlay_type"); overlayType != "" {
		Global.Build.OverlayType = overlayType
	}
}

// AutoDetectCompression sets compression to zstd if supported and user hasn't explicitly set it.
// This should be called after apptainer version is known.
// supportsZstd: whether the current apptainer version supports zstd (>= 1.4)
func AutoDetectCompression(supportsZstd bool) {
	// Only auto-detect if user hasn't explicitly set compress_args in config
	// Empty string in config means "auto-detect"
	if viper.GetString("build.compress_args") != "" {
		// User explicitly set compress_args, respect it
		return
	}

	// Auto-detect: use zstd if supported, otherwise lz4
	if supportsZstd {
		Global.Build.CompressArgs = "-comp zstd -Xcompression-level 8"
		utils.PrintDebug("Auto-detected zstd support, using zstd compression")
	} else {
		Global.Build.CompressArgs = "-comp lz4"
		utils.PrintDebug("Using lz4 compression (zstd not supported)")
	}
}
