package config

import (
	"os"
	"os/exec"
	"path/filepath"
	"time"
)

const VERSION = "1.0.6"
const GitHubRepo = "Justype/condatainer"

// BuildConfig holds default settings for build operations
type BuildConfig struct {
	DefaultCPUs  int           // Default CPUs for builds (if not specified in script)
	DefaultMemMB int64         // Default memory for builds in MB
	DefaultTime  time.Duration // Default time limit
	TmpSizeMB    int           // Size of temporary overlay in MB
	CompressArgs string        // mksquashfs compression arguments
	OverlayType  string        // Overlay filesystem type: "ext3" or "squashfs"
}

// Config holds global application settings
type Config struct {
	// Runtime settings
	Debug     bool
	SubmitJob bool
	Version   string

	// Directory paths
	ProgramDir       string
	BaseDir          string
	ImagesDir        string
	RefImagesDir     string
	BuildScriptsDir  string
	HelperScriptsDir string
	TmpDir           string
	LogsDir          string

	// Binary paths
	ApptainerBin string
	BaseImage    string
	SchedulerBin string // Optional: path to sbatch/scheduler binary (auto-detected if empty)

	// Build configuration
	Build BuildConfig
}

// Global holds the singleton configuration instance
var Global Config

func LoadDefaults(executablePath string) {
	programDir := filepath.Dir(executablePath)
	if absProgDir, err := filepath.Abs(programDir); err == nil {
		programDir = absProgDir
	}

	// Determine base directory with priority:
	// 1. $SCRATCH/condatainer (if $SCRATCH is set)
	// 2. $HOME/condatainer (fallback)
	baseDir := detectBaseDir()

	if absBaseDir, err := filepath.Abs(baseDir); err == nil {
		baseDir = absBaseDir
	}

	Global = Config{
		Debug:     false,
		SubmitJob: true,
		Version:   VERSION,

		ProgramDir:       programDir,
		BaseDir:          baseDir,
		ImagesDir:        filepath.Join(baseDir, "images"),
		BuildScriptsDir:  filepath.Join(baseDir, "build-scripts"),
		HelperScriptsDir: filepath.Join(baseDir, "helper-scripts"),
		TmpDir:           filepath.Join(baseDir, "tmp"),
		LogsDir:          filepath.Join(os.Getenv("HOME"), "logs"),

		ApptainerBin: detectApptainerBin(),
		BaseImage:    filepath.Join(baseDir, "images", "base_image.sif"),
		SchedulerBin: "", // Auto-detect scheduler binary (empty = search PATH)

		Build: BuildConfig{
			DefaultCPUs:  4,             // 4 CPUs default
			DefaultMemMB: 8192,          // 8GB default memory
			DefaultTime:  2 * time.Hour, // 2 hour default time limit
			TmpSizeMB:    20480,         // 20GB temporary overlay
			CompressArgs: "-comp lz4",   // zstd only compatible with apptainer version > 1.4
			OverlayType:  "ext3",        // ext3 for temporary overlays
		},
	}
}

// IsInsideContainer checks if we're currently running inside a container
func IsInsideContainer() bool {
	// Check for IN_CONDATAINER environment variable (our own containers)
	if os.Getenv("IN_CONDATAINER") != "" {
		return true
	}

	// Check for standard Apptainer/Singularity environment variables
	if os.Getenv("APPTAINER_NAME") != "" || os.Getenv("SINGULARITY_NAME") != "" {
		return true
	}

	// Check for Apptainer/Singularity filesystem markers
	if _, err := os.Stat("/.singularity.d"); err == nil {
		return true
	}
	if _, err := os.Stat("/.apptainer.d"); err == nil {
		return true
	}

	return false
}

// detectBaseDir determines the base directory for condatainer data.
// Priority:
//  1. $SCRATCH/condatainer (if $SCRATCH is set, common on HPC systems)
//  2. $HOME/condatainer (fallback)
func detectBaseDir() string {
	// Try $SCRATCH first (common on HPC systems)
	if scratch := os.Getenv("SCRATCH"); scratch != "" {
		return filepath.Join(scratch, "condatainer")
	}

	// Fallback to $HOME/condatainer
	if home := os.Getenv("HOME"); home != "" {
		return filepath.Join(home, "condatainer")
	}

	// Last resort: use current working directory
	if cwd, err := os.Getwd(); err == nil {
		return filepath.Join(cwd, "condatainer")
	}

	return "condatainer"
}

// detectApptainerBin tries to find the apptainer binary, with special handling for containers
func detectApptainerBin() string {
	// If we're inside a container, apptainer might be in a different location
	// or might not be available at all
	if IsInsideContainer() {
		// Try common container paths first
		containerPaths := []string{
			"/usr/local/bin/apptainer",
			"/usr/bin/apptainer",
			"/opt/apptainer/bin/apptainer",
		}

		for _, path := range containerPaths {
			if _, err := os.Stat(path); err == nil {
				return path
			}
		}
	}

	// Fall back to PATH lookup (works both inside and outside containers)
	if binPath, err := exec.LookPath("apptainer"); err == nil {
		return binPath
	}

	// Default to just "apptainer" and let the user's PATH resolve it
	return "apptainer"
}
