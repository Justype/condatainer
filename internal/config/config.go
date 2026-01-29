package config

import (
	"os"
	"os/exec"
	"path/filepath"
	"time"
)

const VERSION = "1.0.6"
const GitHubRepo = "Justype/condatainer"
const GITHUB_REPO = GitHubRepo // Exported constant for compatibility

// PrebuiltBaseURL is the base URL for downloading prebuilt images and overlays
const PrebuiltBaseURL = "https://github.com/Justype/cnt-prebuilt/releases/download/prebuilt"

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
	ProgramDir string
	BaseDir    string
	LogsDir    string

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

		ProgramDir: programDir,
		BaseDir:    baseDir,
		LogsDir:    filepath.Join(os.Getenv("HOME"), "logs"),

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

// GetImagesDir returns the images directory derived from BaseDir
func GetImagesDir() string {
	return filepath.Join(Global.BaseDir, "images")
}

// GetBuildScriptsDir returns the build-scripts directory derived from BaseDir
func GetBuildScriptsDir() string {
	return filepath.Join(Global.BaseDir, "build-scripts")
}

// GetHelperScriptsDir returns the helper-scripts directory derived from BaseDir
func GetHelperScriptsDir() string {
	return filepath.Join(Global.BaseDir, "helper-scripts")
}

// GetTmpDir returns the tmp directory derived from BaseDir
func GetTmpDir() string {
	return filepath.Join(Global.BaseDir, "tmp")
}

// GetBaseImage returns the path to base_image.sif, searching all image directories.
// If user has explicitly set Global.BaseImage (via --base-image flag), returns that.
// Otherwise searches all directories and returns the first found, or falls back to default path.
func GetBaseImage() string {
	// If explicitly set by user (e.g., via --base-image flag), use that
	if Global.BaseImage != "" && Global.BaseImage != filepath.Join(Global.BaseDir, "images", "base_image.sif") {
		return Global.BaseImage
	}

	// Search all image directories
	if found := FindBaseImage(); found != "" {
		return found
	}

	// Fall back to user's images directory (default location)
	return filepath.Join(Global.BaseDir, "images", "base_image.sif")
}

// isWritableDir checks if a directory is writable by trying to create a test file
func isWritableDir(dir string) bool {
	// Try to create directory if it doesn't exist
	if err := os.MkdirAll(dir, 0755); err != nil {
		return false
	}

	// Test write permission
	testFile := filepath.Join(dir, ".write-test")
	f, err := os.Create(testFile)
	if err != nil {
		return false
	}
	f.Close()
	os.Remove(testFile)
	return true
}

// GetWritableTmpDir returns the first writable tmp directory
func GetWritableTmpDir() string {
	// Check extra base dirs
	for _, baseDir := range GetExtraBaseDirs() {
		if baseDir == "" {
			continue
		}
		tmpDir := filepath.Join(baseDir, "tmp")
		if isWritableDir(tmpDir) {
			return tmpDir
		}
	}

	// Check portable data dir
	if portableDir := GetPortableDataDir(); portableDir != "" {
		tmpDir := filepath.Join(portableDir, "tmp")
		if isWritableDir(tmpDir) {
			return tmpDir
		}
	}

	// Check scratch data dir
	if scratchDir := GetScratchDataDir(); scratchDir != "" {
		tmpDir := filepath.Join(scratchDir, "tmp")
		if isWritableDir(tmpDir) {
			return tmpDir
		}
	}

	// Check user data dir
	if userDir := GetUserDataDir(); userDir != "" {
		tmpDir := filepath.Join(userDir, "tmp")
		if isWritableDir(tmpDir) {
			return tmpDir
		}
	}

	// Fall back to base dir (should always be writable)
	return filepath.Join(Global.BaseDir, "tmp")
}
