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
// 3. User config file (~/.config/condatainer/config.yaml)
// 4. System config file (/etc/condatainer/config.yaml)
// 5. Defaults
func InitViper() error {
	viper.SetConfigName(ConfigFilename)
	viper.SetConfigType(ConfigType)

	// Set config search paths (order matters)
	// User config (highest priority)
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
	viper.SetDefault("apptainer_bin", "apptainer")
	viper.SetDefault("scheduler_bin", "")
	viper.SetDefault("scheduler_type", "")
	viper.SetDefault("submit_job", true)

	// Build config defaults
	viper.SetDefault("build.default_cpus", 4)
	viper.SetDefault("build.default_mem_mb", 8192)
	viper.SetDefault("build.default_time", "2h")
	viper.SetDefault("build.tmp_size_mb", 20480)
	viper.SetDefault("build.compress_args", "-comp lz4")
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

// SaveConfig saves current Viper config to user config file
func SaveConfig() error {
	configPath, err := GetUserConfigPath()
	if err != nil {
		return fmt.Errorf("failed to get config path: %w", err)
	}

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

	// Check and detect apptainer binary
	apptainerBin := viper.GetString("apptainer_bin")
	if !ValidateBinary(apptainerBin) {
		detected := DetectApptainerBin()
		if detected != "" {
			viper.Set("apptainer_bin", detected)
			updated = true
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
	updated := false

	// Always re-detect apptainer binary
	detected := DetectApptainerBin()
	if detected != "" {
		currentBin := viper.GetString("apptainer_bin")
		if currentBin != detected {
			viper.Set("apptainer_bin", detected)
			updated = true
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
	if err := SaveConfig(); err != nil {
		return false, err
	}

	return updated, nil
}

// LoadFromViper loads config from Viper into Global struct
func LoadFromViper() {
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

	if submitJob := viper.GetBool("submit_job"); !submitJob {
		Global.SubmitJob = submitJob
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

	if compressArgs := viper.GetString("build.compress_args"); compressArgs != "" {
		Global.Build.CompressArgs = compressArgs
	}

	if overlayType := viper.GetString("build.overlay_type"); overlayType != "" {
		Global.Build.OverlayType = overlayType
	}
}
