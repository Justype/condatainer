package config

import (
	"fmt"
	"os"
	"os/exec"
	"path/filepath"
	"strconv"
	"strings"
	"time"

	"github.com/Justype/condatainer/catalog"
	"github.com/Justype/condatainer/internal/scheduler"
	"github.com/Justype/condatainer/internal/utils"
)

// ExitCodeJobsSubmitted is returned by condatainer commands when build jobs are
// submitted to the scheduler rather than run immediately. Callers (e.g. helper
// overlay checks) can detect this exit code to show a more helpful message.
const ExitCodeJobsSubmitted = 3

// VERSION is a variable so tagged release builds can set it with -ldflags -X
// without modifying tracked source before compilation.
var VERSION = "1.4.1"

const GitHubRepo = "Justype/condatainer"

// BuildConfig holds default settings for build operations
type BuildConfig struct {
	Defaults        scheduler.ResourceSpec // Default resource spec for build job submissions
	CompressArgs    string                 // mksquashfs compression arguments
	BlockSize       string                 // mksquashfs block size for app/env/external overlays (DefaultBlockSize)
	DataBlockSize   string                 // mksquashfs block size for data/ref overlays (DefaultDataBlockSize)
	AlwaysSubmit    bool                   // Always submit builds as scheduler jobs even without script directives (default: false)
	Channels        []string               // conda channels in priority order (default: [conda-forge, bioconda])
	SystemApptainer string                 // Path to the system/module apptainer or singularity binary (auto-detected if empty)
}

// SchedulerConfig holds scheduler binary/submission settings shared by every
// caller that submits a job without parsing an existing script's own
// directives (e.g. helper) — Account/Partition/Defaults are that caller's
// fallback, not a build-specific one.
type SchedulerConfig struct {
	Bin       string                 // Path to sbatch/scheduler binary (auto-detected if empty)
	Timeout   time.Duration          // Scheduler command timeout (default: 0 = no timeout)
	Account   string                 // Default billing/allocation account (empty = scheduler's own default)
	Partition string                 // Default partition/queue (empty = scheduler's own default)
	Defaults  scheduler.ResourceSpec // Default resource spec for jobs with no script directives
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

	// Recipe collections, in order. Earlier entries shadow later ones.
	Sources []catalog.Spec

	// Default distro for the container root, e.g. "ubuntu24" -> recipes/ubuntu24/base.def.
	// Empty falls back to the first source declaring a default_distro.
	DefaultDistro string

	// Pass --nv / --rocm when the host has the matching device nodes (default: true).
	// Turn off on a node whose driver is present but unusable: the device nodes
	// still exist, so detection fires and the container then fails to start.
	AutoloadGPU bool

	// NestedRun says whether an apptainer is provided inside a container so
	// containers can be launched from within one: "auto" (default) uses what
	// exists and builds nothing, "true" also builds the apptainer overlay when it
	// is missing and fails if it cannot be provided, "false" turns it off.
	NestedRun string

	// Notification method when a helper job starts running (default: "web").
	// Values: "web" (browser notification via dashboard), "terminal" (bell ×2, 1.1 s apart),
	// "both" (terminal + web), "" or "none" (silent).
	Notification string

	// Max age of the on-disk remote build script metadata cache (default: 1 week)
	MetadataCacheTTL time.Duration

	// Age below which `store gc` never reports an entry collectable (default: 30
	// days), overridden per invocation by --grace.
	//
	// Configurable because a store in a group root serves everyone using that
	// install, and its owner is who knows its turnover. It resolves per
	// invocation rather than per store, which is what --layer is for.
	StoreGCGrace time.Duration

	// ProxyPerJob: if true, inject "condatainer proxy start --via <login-node>" at the
	// top of every generated scheduler script so the compute node starts a per-job proxy.
	ProxyPerJob bool

	// HelperBindAll: if true, helper services bind to 0.0.0.0 and the server connects
	// via direct TCP instead of SSH tunnel (for clusters where inter-node SSH is blocked).
	HelperBindAll bool

	// Build configuration
	Build BuildConfig

	// Scheduler configuration
	Scheduler SchedulerConfig
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

// Compiled-in defaults for the scalar config keys.
//
// Named because each is needed both in LoadDefaults, which sets what the program
// runs on, and in setDefaults, which is what `config get` reports for an unset
// key. Two literals drift and the two commands then disagree.
const (
	// Data is large and read sequentially, so it takes the bigger block; an app
	// is many small files where a big block wastes read bandwidth.
	DefaultBlockSize     = "128k"
	DefaultDataBlockSize = "512k"

	DefaultNcpus        = 4    // CPUs for a build job
	DefaultMemMB        = 8192 // memory for a build job
	DefaultBuildTime    = "2h" // walltime for a build job
	DefaultCacheTTLDay  = 7    // remote recipe metadata cache, 1 week
	DefaultGCGraceDay   = 30   // store gc: age below which an entry is never collectable
	DefaultNotification = "web"
	DefaultNestedRun    = NestedRunAuto

	DefaultSchedulerNcpus = 1    // CPUs for a job with no script directives
	DefaultSchedulerMemMB = 2048 // memory for a job with no script directives
	DefaultSchedulerTime  = "2h" // walltime for a job with no script directives
)

// DefaultBuildDuration is DefaultBuildTime as a duration, so the two cannot disagree.
var DefaultBuildDuration = 2 * time.Hour

// DefaultSchedulerDuration is DefaultSchedulerTime as a duration, so the two cannot disagree.
var DefaultSchedulerDuration = 2 * time.Hour

// DefaultChannels are the conda channels micromamba gets, highest priority first.
func DefaultChannels() []string { return []string{"conda-forge", "bioconda"} }

// DefaultLogsDir is where job logs land when nothing configures it.
func DefaultLogsDir() string { return filepath.Join(os.Getenv("HOME"), "logs") }

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
		Debug:       false,
		SubmitJob:   true,
		AutoloadGPU: true,
		NestedRun:   DefaultNestedRun,
		Version:     VERSION,

		ProgramDir: programDir,
		LogsDir:    DefaultLogsDir(),

		Notification:     DefaultNotification,
		MetadataCacheTTL: DefaultCacheTTLDay * 24 * time.Hour,
		StoreGCGrace:     DefaultGCGraceDay * 24 * time.Hour,

		Build: BuildConfig{
			Defaults: scheduler.ResourceSpec{
				CpusPerTask:  DefaultNcpus,
				MemPerNodeMB: DefaultMemMB,
				Time:         DefaultBuildDuration,
			},
			// Every reader is >= 1.4: libexec's apptainer (verified at provision
			// time) or a version-checked system apptainer (apptainer.Normal).
			// A registry consumer outside condatainer's own exec/run is
			// responsible for having a compatible apptainer themselves.
			CompressArgs:    ArgsForCompress("zstd-medium"),
			BlockSize:       DefaultBlockSize,
			DataBlockSize:   DefaultDataBlockSize,
			Channels:        DefaultChannels(),
			SystemApptainer: apptainerOnPath(),
		},

		Scheduler: SchedulerConfig{
			Bin:     "", // Auto-detect scheduler binary (empty = search PATH)
			Timeout: 0,  // no timeout by default
			Defaults: scheduler.ResourceSpec{
				CpusPerTask:  DefaultSchedulerNcpus,
				MemPerNodeMB: DefaultSchedulerMemMB,
				Time:         DefaultSchedulerDuration,
			},
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

// apptainerOnPath returns the apptainer (or singularity) found on PATH, "" when neither is.
func apptainerOnPath() string {
	for _, name := range []string{"apptainer", "singularity"} {
		if binPath, err := exec.LookPath(name); err == nil {
			return binPath
		}
	}
	return ""
}

// BaseImageFileName returns the expected filename for the configured base, e.g.
// "ubuntu24" → "ubuntu24--base.sqf" — the flat name any artifact of that name
// gets. Empty when no base is configured, so callers do not go looking for a
// file called ".sqf".
func BaseImageFileName() string {
	name := BaseRecipeName()
	if name == "" {
		return ""
	}
	return strings.ReplaceAll(name, "/", "--") + ".sqf"
}

// GetBaseImage returns the installed base image, searching every image directory.
//
// It only ever returns a file that exists — conflating this with where a base
// could be written used to hand callers a path to a nonexistent image and let
// Apptainer report it.
func GetBaseImage() (string, error) {
	if found := FindBaseImage(); found != "" {
		return found, nil
	}
	name := BaseRecipeName()
	if name == "" {
		return "", fmt.Errorf("no default distro configured: set `default_distro`, or configure a source declaring default_distro")
	}
	return "", fmt.Errorf("base image %s is not installed (searched %s)",
		name, strings.Join(GetImageSearchPaths(), ", "))
}

// GetWritableTmpDir returns the first writable tmp directory under the data
// dirs — the stable root, for work too large or too long-lived for node-local
// scratch. Shared dirs (extra-root, root) create tmp only when the parent
// already exists; personal dirs (scratch, user) always create on first use.
//
// CNT_TMPDIR does not redirect this: it selects the fast root, and collapsing
// the two would put a data build's payload on scratch that a job wipes.
func GetWritableTmpDir() string {
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

	// No writable data dir at all: fall back to the fast root rather than a
	// relative path, which would put the workspace wherever the caller stood.
	return utils.GetTmpDir()
}

// NestedRun values.
const (
	NestedRunAuto  = "auto"
	NestedRunTrue  = "true"
	NestedRunFalse = "false"
)

// NestedRunValues lists the accepted values of the nested_run key.
var NestedRunValues = []string{NestedRunAuto, NestedRunTrue, NestedRunFalse}

// ParseNestedRun normalizes a nested_run value and reports whether it is valid.
func ParseNestedRun(v string) (string, bool) {
	v = strings.ToLower(strings.TrimSpace(v))
	for _, ok := range NestedRunValues {
		if v == ok {
			return v, true
		}
	}
	return "", false
}
