package config

import (
	"os"
	"os/exec"
	"path/filepath"
	"slices"
	"strconv"
	"strings"
	"time"

	"github.com/Justype/condatainer/internal/scheduler"
	"github.com/Justype/condatainer/internal/utils"
)

const VERSION = "1.4.1"
const GitHubRepo = "Justype/condatainer"
const GITHUB_REPO = GitHubRepo // Exported constant for compatibility
const DEFAULT_DISTRO = "ubuntu24"

func GetAvailableDistros() []string {
	return []string{
		"ubuntu20",
		"ubuntu22",
		"ubuntu24",
	}
}

// IsValidDistro returns true if the provided name matches one of the known distro slugs.
func IsValidDistro(distro string) bool {
	return slices.Contains(GetAvailableDistros(), distro)
}

// DefaultScriptsLink is the base URL for remote build scripts and helpers
const DefaultScriptsLink = "https://raw.githubusercontent.com/Justype/cnt-scripts/main"

// BuildConfig holds default settings for build operations
type BuildConfig struct {
	Defaults      scheduler.ResourceSpec // Default resource spec for build job submissions
	TmpSizeMB     int                    // Size of temporary overlay in MB
	CompressArgs  string                 // mksquashfs compression arguments
	BlockSize     string                 // mksquashfs block size for app/env/external overlays (default: 128k)
	DataBlockSize string                 // mksquashfs block size for data/ref overlays (default: 1m)
	UseTmpOverlay bool                   // Use a temporary overlay instead of a temp directory (default: false)
	AlwaysSubmit  bool                   // Always submit builds as scheduler jobs even without script directives (default: false)
	Channels      []string               // conda channels in priority order (default: [conda-forge, bioconda])
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

	// Base OS overlay slug (e.g. "ubuntu24"). Determines the base image filename:
	// "ubuntu24" → "ubuntu24--base_image.sif". Defaults to "ubuntu24".
	DefaultDistro string

	// Remote repository settings
	ScriptsLink  string   // Base remote URL (scripts_link config key; lowest priority)
	ScriptsLinks []string // Effective ordered list: [extra_scripts_links..., scripts_link]
	PreferRemote bool     // Remote build scripts take precedence over local

	// Dependency parsing
	ParseModuleLoad bool // Parse "module load" / "ml" lines as dependencies (default: false)

	// Scheduler command timeout in seconds (default: 5)
	SchedulerTimeout time.Duration

	// Notification method when a helper job starts running (default: "" = none).
	// Values: "bell" (terminal bell), "email" (scheduler email), ≥5-char string (ntfy.sh topic), "" or "none" (silent).
	Notification string

	// Max age of the on-disk remote build script metadata cache (default: 1 week)
	MetadataCacheTTL time.Duration

	// Build configuration
	Build BuildConfig
}

// CompressOption defines name, mksquashfs arguments, and description
type CompressOption struct {
	Name        string // e.g. "lz4" or "zstd-fast"
	Args        string // full mksquashfs arguments
	Description string // help text for CLI
}

// CompressOptions lists all supported compression shorthand names.
var CompressOptions = []CompressOption{
	{"gzip", "-comp gzip", "Use gzip compression"},
	{"lz4", "-comp lz4", "Use lz4 compression"},
	{"zstd", "-comp zstd -Xcompression-level 14", "Use zstd compression level 14"},
	{"zstd-fast", "-comp zstd -Xcompression-level 3", "Use zstd compression level 3"},
	{"zstd-medium", "-comp zstd -Xcompression-level 8", "Use zstd compression level 8"},
	{"zstd-high", "-comp zstd -Xcompression-level 19", "Use zstd compression level 19"},
}

// ArgsForCompress returns the full mksquashfs arguments corresponding to a
// recognised shortcut name. Return itself if no match found.
func ArgsForCompress(name string) string {
	for _, o := range CompressOptions {
		if o.Name == name {
			return o.Args
		}
	}
	return name
}

// CompressNames returns a list of the shortcut names (used for completion)
func CompressNames() []string {
	names := make([]string, len(CompressOptions))
	for i, o := range CompressOptions {
		names[i] = o.Name
	}
	return names
}

// BlockSizeCompletions lists common mksquashfs block sizes for shell completion
var BlockSizeCompletions = []string{"64k", "128k", "256k", "512k", "1m"}

// IsValidBlockSize validates a mksquashfs -b value.
// Must be a power of two between 4096 and 1048576 (1M), with optional k/K or m/M suffix.
func IsValidBlockSize(size string) bool {
	if size == "" {
		return false
	}
	s := strings.ToLower(size)
	var multiplier int64 = 1
	if strings.HasSuffix(s, "k") {
		multiplier = 1024
		s = s[:len(s)-1]
	} else if strings.HasSuffix(s, "m") {
		multiplier = 1024 * 1024
		s = s[:len(s)-1]
	}
	n, err := strconv.ParseInt(s, 10, 64)
	if err != nil || n <= 0 {
		return false
	}
	bytes := n * multiplier
	// mksquashfs requires power of two, between 4096 and 1048576 (1M)
	return bytes >= 4096 && bytes <= 1048576 && (bytes&(bytes-1)) == 0
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

		ApptainerBin:  detectApptainerBin(),
		SchedulerBin:  "", // Auto-detect scheduler binary (empty = search PATH)
		DefaultDistro: DEFAULT_DISTRO,

		ScriptsLink:      DefaultScriptsLink,
		ScriptsLinks:     []string{DefaultScriptsLink}, // overwritten in LoadFromViper
		SchedulerTimeout: 5 * time.Second,
		Notification:     "",
		MetadataCacheTTL: 7 * 24 * time.Hour, // 1 week

		Build: BuildConfig{
			Defaults: scheduler.ResourceSpec{
				CpusPerTask:  4,             // 4 CPUs default
				MemPerNodeMB: 8192,          // 8GB default memory
				Time:         2 * time.Hour, // 2 hour default time limit
			},
			TmpSizeMB:     20480,                               // 20GB temporary overlay
			CompressArgs:  "-comp lz4",                         // zstd only compatible with apptainer version > 1.4
			BlockSize:     "128k",                              // mksquashfs block size for app/env/external overlays
			DataBlockSize: "512k",                              // mksquashfs block size for data/ref overlays
			Channels:      []string{"conda-forge", "bioconda"}, // default conda channels
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

// detectApptainerBin tries to find the apptainer (or singularity) binary, with special handling
// for containers. Returns the full path when found, or "" when neither binary is available.
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

	// Try both names in PATH order: apptainer first, singularity as fallback
	for _, name := range []string{"apptainer", "singularity"} {
		if binPath, err := exec.LookPath(name); err == nil {
			return binPath
		}
	}

	// Neither found; callers treat "" as "no binary available"
	return ""
}

// BaseImageSifName returns the expected .sif filename for the configured DefaultDistro.
// e.g. "ubuntu24" → "ubuntu24--base_image.sif"
func BaseImageSifName() string {
	return Global.DefaultDistro + "--base_image.sif"
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

// GetWritableTmpDir returns the first writable tmp directory.
// Uses the same SearchDir + firstWritableDir logic as the other writable-dir resolvers:
// shared dirs (extra-root, root) create tmp only when the parent already exists;
// personal dirs (scratch, user) always create on first use.
func GetWritableTmpDir() string {
	// Global override for all build temp paths
	if os.Getenv("CNT_TMPDIR") != "" {
		return utils.GetTmpDir()
	}

	var dirs []SearchDir
	if extraRoot := GetExtraRootDir(); extraRoot != "" {
		dirs = append(dirs, SearchDir{Path: filepath.Join(extraRoot, "tmp")})
	}
	if rootDir := GetRootDir(); rootDir != "" {
		dirs = append(dirs, SearchDir{Path: filepath.Join(rootDir, "tmp")})
	}
	if scratchDir := GetScratchDataDir(); scratchDir != "" {
		dirs = append(dirs, SearchDir{Path: filepath.Join(scratchDir, "tmp"), Personal: true})
	}
	if userDir := GetUserDataDir(); userDir != "" {
		dirs = append(dirs, SearchDir{Path: filepath.Join(userDir, "tmp"), Personal: true})
	}

	if dir := firstWritableDir(deduplicateWriteDirs(dirs)); dir != "" {
		return dir
	}

	// Last resort: current directory
	return "tmp"
}
