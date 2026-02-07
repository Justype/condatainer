package container

import (
	"fmt"
	"os"
	"path/filepath"
	"sort"
	"strconv"
	"strings"

	"github.com/Justype/condatainer/internal/config"
	"github.com/Justype/condatainer/internal/overlay"
	"github.com/Justype/condatainer/internal/utils"
)

var (
	commonEnvVars = []string{
		"LC_ALL=C.UTF-8",
		"LANG=C.UTF-8",
		"CURL_CA_BUNDLE=", // Unset CA bundle to avoid host ENV interference
		"SSL_CERT_FILE=",  // Unset SSL cert file for the same reason
	}
)

// SetupConfig holds all configuration needed to set up a container
type SetupConfig struct {
	Overlays       []string // Overlay paths (will be resolved)
	WritableImg    bool     // Whether .img overlays should be writable
	EnvSettings    []string // User-specified environment variables (KEY=VALUE format)
	BindPaths      []string // User-specified bind paths
	Fakeroot       bool     // Whether to use fakeroot
	ApptainerFlags []string // Additional apptainer flags to pass through
}

// SetupResult contains all the processed configuration ready for container execution
type SetupResult struct {
	Overlays       []string          // Resolved and ordered overlay paths
	OverlayArgs    []string          // Overlay paths with :ro/:rw suffixes
	EnvList        []string          // Complete environment variable list
	EnvNotes       map[string]string // Environment variable notes for display
	BindPaths      []string          // Deduplicated bind paths
	Fakeroot       bool              // Final fakeroot setting (may be auto-enabled)
	ApptainerFlags []string          // Apptainer flags including GPU flags
	LastImg        string            // Path to .img overlay if present
}

// Setup processes all container configuration and returns a ready-to-use result
func Setup(cfg SetupConfig) (*SetupResult, error) {
	// Resolve overlay paths
	overlays, err := ResolveOverlayPaths(cfg.Overlays)
	if err != nil {
		return nil, err
	}

	// Ensure at most one .img overlay
	if err := ensureSingleImage(overlays); err != nil {
		return nil, err
	}

	// Put .img overlay last if present
	overlays = putImgToLast(overlays)

	// Process overlays and check availability
	overlayArgs := make([]string, 0, len(overlays))
	var lastImg string
	for _, ol := range overlays {
		isImg := utils.IsImg(ol)
		if isImg {
			lastImg = ol

			// Only lock .img files as requested
			if utils.FileExists(ol) && !utils.DirExists(ol) {
				// If it's the principal image and writableImg is true, we need an exclusive lock.
				// In putImgToLast, the principal image is always the last one.
				isPrincipalImg := (ol == overlays[len(overlays)-1])
				writeLock := isPrincipalImg && cfg.WritableImg

				if err := overlay.CheckAvailable(ol, writeLock); err != nil {
					return nil, err
				}
			}
		}

		overlayArgs = append(overlayArgs, FormatOverlayMount(ol, cfg.WritableImg))
	}

	// Build environment variables
	envList, envNotes := buildEnvironment(overlays, lastImg, cfg)

	// Build bind paths
	bindPaths := BindPaths()
	if len(cfg.BindPaths) > 0 {
		bindPaths = append(bindPaths, cfg.BindPaths...)
	}
	bindPaths = DeduplicateBindPaths(bindPaths)

	// Detect GPU flags
	apptainerFlags := append([]string{}, DetectGPUFlags()...)
	apptainerFlags = append(apptainerFlags, cfg.ApptainerFlags...)

	return &SetupResult{
		Overlays:       overlays,
		OverlayArgs:    overlayArgs,
		EnvList:        envList,
		EnvNotes:       envNotes,
		BindPaths:      bindPaths,
		Fakeroot:       cfg.Fakeroot,
		ApptainerFlags: apptainerFlags,
		LastImg:        lastImg,
	}, nil
}

// buildEnvironment constructs the complete environment variable list
func buildEnvironment(overlays []string, lastImg string, cfg SetupConfig) ([]string, map[string]string) {
	// Collect overlay environment variables (from .env files)
	configs, notes := CollectOverlayEnv(overlays)
	envKeys := make([]string, 0, len(configs))
	for key := range configs {
		envKeys = append(envKeys, key)
	}
	sort.Strings(envKeys)

	// Start with overlay-specific vars
	envList := make([]string, 0, len(envKeys)+20)
	for _, key := range envKeys {
		envList = append(envList, fmt.Sprintf("%s=%s", key, configs[key]))
	}

	// Add layer tracking
	layer := 0
	if raw := os.Getenv("IN_CONDATAINER"); raw != "" {
		if parsed, err := strconv.Atoi(raw); err == nil && parsed >= 1 {
			layer = parsed
		}
	}
	layer++
	ps1Prefix := "CNT"
	if layer > 1 {
		ps1Prefix = fmt.Sprintf("CNT_%d", layer)
	}

	// Build PATH
	pathEnv := BuildPathEnv(overlays)

	// For nested usage: if this is a portable installation, prepend its bin directory to PATH
	if config.IsPortable() {
		if portableBin := config.GetPortableBinDir(); portableBin != "" {
			pathEnv = portableBin + ":" + pathEnv
		}
	}

	// Add standard environment variables
	envList = append(envList,
		fmt.Sprintf("PATH=%s", pathEnv),
		fmt.Sprintf("PS1=%s \\[\\e[0;34m\\]\\w\\[\\e[0m\\]> ", ps1Prefix),
		fmt.Sprintf("IN_CONDATAINER=%d", layer),
	)
	envList = append(envList, commonEnvVars...)

	// Add .img-specific environment variables
	if lastImg != "" {
		if os.Getenv("IN_CONDATAINER") != "" {
			utils.PrintWarning("You are trying to mount an .img overlay inside an existing CondaTainer environment. This may lead to unexpected behavior.")
		}

		envList = append(envList,
			"MAMBA_ROOT_PREFIX=/ext3/env",
			"CONDA_PREFIX=/ext3/env",
			"CONDA_DEFINE_ENV=env",
			"RETICULATE_PYTHON=/ext3/env/bin/python",
		)

		if cfg.WritableImg {
			envList = append(envList, "CNT_CONDA_PREFIX=/ext3/env")
		}
	}

	// Add user-specified environment variables (with validation)
	for _, setting := range cfg.EnvSettings {
		setting = strings.TrimSpace(setting)
		if setting == "" {
			continue
		}
		if !strings.Contains(setting, "=") {
			utils.PrintWarning("Invalid env setting %s. It should be in KEY=VALUE format. Skipping.", utils.StyleName(setting))
			continue
		}
		envList = append(envList, setting)
	}

	// Prepare environment notes for display
	envNotes := make(map[string]string)
	for key, value := range configs {
		note := notes[key]
		if note == "" {
			note = value
		}
		envNotes[key] = note
	}
	if lastImg != "" && cfg.WritableImg {
		envNotes["CNT_CONDA_PREFIX"] = "/ext3/env"
	}

	return envList, envNotes
}

// AutoEnableFakeroot checks if fakeroot should be auto-enabled for writable .img overlays
// Returns the updated fakeroot setting
func AutoEnableFakeroot(lastImg string, writable bool, currentFakeroot bool) bool {
	if lastImg == "" || !writable || currentFakeroot {
		return currentFakeroot
	}

	// Check UID status and auto-enable fakeroot if needed
	if status := overlay.InspectImageUIDStatus(lastImg); status == overlay.UIDStatusRoot {
		utils.PrintNote("Root overlay %s detected. --fakeroot enabled automatically.", utils.StylePath(filepath.Base(lastImg)))
		return true
	} else if status == overlay.UIDStatusDifferentUser {
		utils.PrintWarning("%s's inner UID differs from current user. --fakeroot enabled automatically.", utils.StylePath(filepath.Base(lastImg)))
		return true
	}

	return currentFakeroot
}

// ensureSingleImage checks that at most one .img overlay is specified
func ensureSingleImage(overlays []string) error {
	imgCount := 0
	for _, overlay := range overlays {
		if utils.IsImg(overlay) {
			imgCount++
		}
	}
	if imgCount > 1 {
		return fmt.Errorf("only one .img overlay is allowed, found %d", imgCount)
	}
	return nil
}

// putImgToLast moves any .img overlay to the end of the overlay list
func putImgToLast(overlays []string) []string {
	var imgOverlay string
	var others []string

	for _, overlay := range overlays {
		if utils.IsImg(overlay) {
			imgOverlay = overlay
		} else {
			others = append(others, overlay)
		}
	}

	if imgOverlay != "" {
		return append(others, imgOverlay)
	}
	return overlays
}
