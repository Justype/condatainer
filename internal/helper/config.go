package helper

import (
	"bufio"
	"fmt"
	"os"
	"path/filepath"
	"strings"

	"github.com/Justype/condatainer/internal/config"
	"github.com/Justype/condatainer/internal/utils"
)

// HelperConfigPath returns the path to the per-helper saved config file.
// Path: ~/.config/condatainer/helper-config/{name}.conf
func HelperConfigPath(name string) string {
	configDir := config.GetUserConfigDir()
	if configDir == "" {
		return ""
	}
	return filepath.Join(configDir, "helper-config", name+".conf")
}

func helperConfigPath(name string) string { return HelperConfigPath(name) }

// GetHelperConfigKey reads a single key from the per-helper config.
// Keys are case-insensitive. Returns ("", false) if the file or key is absent.
func GetHelperConfigKey(name, key string) (string, bool) {
	cfg, err := LoadHelperConfig(name)
	if err != nil {
		return "", false
	}
	v, ok := cfg[strings.ToLower(key)]
	return v, ok
}

// SetHelperConfigKey updates or inserts a single key in the per-helper config,
// preserving all other existing keys. Keys are normalized to lowercase.
func SetHelperConfigKey(name, key, value string) error {
	cfg, err := LoadHelperConfig(name)
	if err != nil {
		cfg = map[string]string{}
	}
	cfg[strings.ToLower(key)] = value
	return SaveHelperConfig(name, cfg)
}

// LoadHelperConfig reads the saved KEY=value config for a helper.
// Returns an empty map if the file doesn't exist.
func LoadHelperConfig(name string) (map[string]string, error) {
	path := helperConfigPath(name)
	if path == "" {
		return map[string]string{}, nil
	}
	f, err := os.Open(path)
	if err != nil {
		if os.IsNotExist(err) {
			return map[string]string{}, nil
		}
		return nil, fmt.Errorf("opening helper config: %w", err)
	}
	defer f.Close()

	cfg := make(map[string]string)
	scanner := bufio.NewScanner(f)
	for scanner.Scan() {
		line := strings.TrimSpace(scanner.Text())
		if line == "" || strings.HasPrefix(line, "#") {
			continue
		}
		k, v, _ := strings.Cut(line, "=")
		if k != "" {
			cfg[strings.ToLower(strings.TrimSpace(k))] = strings.TrimSpace(v)
		}
	}
	return cfg, scanner.Err()
}

// SaveHelperConfig writes KEY=value pairs to the per-helper config file.
func SaveHelperConfig(name string, cfg map[string]string) error {
	path := helperConfigPath(name)
	if path == "" {
		return fmt.Errorf("cannot determine helper config path")
	}
	dir := filepath.Dir(path)
	if err := utils.MkdirAllShared(dir); err != nil {
		return err
	}
	f, err := utils.CreateFileWritable(path)
	if err != nil {
		return err
	}
	defer f.Close()
	for k, v := range cfg {
		if _, err := fmt.Fprintf(f, "%s=%s\n", k, v); err != nil {
			return err
		}
	}
	return nil
}
