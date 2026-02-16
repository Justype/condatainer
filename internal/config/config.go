package config

import (
	"os"
	"os/exec"
	"path/filepath"
	"time"

	"github.com/Justype/condatainer/internal/utils"
)

const VERSION = "1.1.1"
const GitHubRepo = "Justype/condatainer"
const GITHUB_REPO = GitHubRepo // Exported constant for compatibility

// PrebuiltBaseURL is the base URL for downloading prebuilt images and overlays
const PrebuiltBaseURL = "https://github.com/Justype/cnt-prebuilt/releases/download/prebuilt"

// SchedulerConfig holds default resource specs for scheduler script parsing
type SchedulerConfig struct {
	Ncpus  int           // Default CPUs per task (default: 4)
	MemMB  int64         // Default memory in MB (default: 0 = unset)
	Time   time.Duration // Default time limit (default: 0 = unset)
	Nodes  int           // Default number of nodes (default: 1)
	Ntasks int           // Default number of tasks (default: 1)
}

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
	LogsDir    string

	// Binary paths
	ApptainerBin string
	SchedulerBin string // Optional: path to sbatch/scheduler binary (auto-detected if empty)

	// Remote repository settings
	Branch       string // Git branch for fetching remote build scripts and metadata (default: "main")
	PreferRemote bool   // Remote build scripts take precedence over local

	// Scheduler default specs
	Scheduler SchedulerConfig

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

	Global = Config{
		Debug:     false,
		SubmitJob: true,
		Version:   VERSION,

		ProgramDir: programDir,
		LogsDir:    filepath.Join(os.Getenv("HOME"), "logs"),

		ApptainerBin: detectApptainerBin(),
		SchedulerBin: "", // Auto-detect scheduler binary (empty = search PATH)

		Branch: "main", // Default branch for remote build scripts

		Scheduler: SchedulerConfig{
			Ncpus:  2,
			MemMB:  8192,          // 8GB
			Time:   4 * time.Hour, // 4 hours
			Nodes:  1,
			Ntasks: 1,
		},

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

// GetBaseImage returns the path to base_image.sif, searching all image directories.
// Searches all directories and returns the first found, or falls back to first writable path.
func GetBaseImage() string {
	// Search all image directories
	if found := FindBaseImage(); found != "" {
		return found
	}

	// Fall back to first writable path (where base image would be downloaded/created)
	if writePath, err := GetBaseImageWritePath(); err == nil {
		return writePath
	}

	// Last resort: use user data dir (backward compatibility)
	if userDir := GetUserDataDir(); userDir != "" {
		return filepath.Join(userDir, "images", "base_image.sif")
	}

	return "base_image.sif"
}

// isWritableDir checks if a directory is writable by trying to create a test file
func isWritableDir(dir string) bool {
	// Try to create directory if it doesn't exist
	if err := os.MkdirAll(dir, utils.PermDir); err != nil {
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

	// Last resort: current directory
	return "tmp"
}
