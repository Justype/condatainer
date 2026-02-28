// Package scheduler provides a unified interface for HPC job schedulers
package scheduler

import (
	"fmt"
	"os"
	"os/exec"
	"path/filepath"
	"strconv"
	"strings"
	"time"
)

// SchedulerType represents the type of job scheduler
type SchedulerType string

const (
	SchedulerUnknown  SchedulerType = ""
	SchedulerSLURM    SchedulerType = "SLURM"
	SchedulerPBS      SchedulerType = "PBS"
	SchedulerLSF      SchedulerType = "LSF"
	SchedulerHTCondor SchedulerType = "HTCondor"
)

// SchedulerInfo holds information about the detected scheduler
type SchedulerInfo struct {
	Type      string // Scheduler type (e.g., "SLURM", "PBS", "LSF")
	Binary    string // Path to scheduler binary (e.g., "/usr/bin/sbatch")
	Version   string // Scheduler version (if available)
	InJob     bool   // Whether we're currently inside a scheduled job
	Available bool   // Whether scheduler is available for job submission
}

// GpuSpec holds GPU requirements parsed from scheduler directives.
// Count is per node (not total). Total GPUs = Count * ResourceSpec.Nodes.
type GpuSpec struct {
	Type  string // GPU type/model (e.g., "h100", "a100", "v100")
	Count int    // Number of GPUs per node
	Raw   string // Raw GPU specification string from scheduler
}

// GpuInfo holds information about available GPU types in the cluster
type GpuInfo struct {
	Type      string // GPU type/model
	Total     int    // Total number of this GPU type
	Available int    // Currently available
	Partition string // Partition/queue name (optional)
}

// ResourceLimits holds scheduler resource limits per partition/queue.
// CPU and memory limits are per-node (matching how schedulers report them).
type ResourceLimits struct {
	MaxNodes        int           // Maximum nodes per job
	MaxCpusPerNode  int           // Maximum CPUs per node in this partition
	MaxMemMBPerNode int64         // Maximum memory per node in MB
	MaxGpus         int           // Maximum GPUs per job (total across all nodes)
	DefaultTime     time.Duration // Default walltime if not specified
	MaxTime         time.Duration // Maximum walltime (e.g., "7-00:00:00")
	Partition       string        // Partition/queue name these limits apply to
}

// ResourceSpec holds compute geometry for a job.
// A nil pointer means resource parsing failed — the job runs in passthrough mode.
type ResourceSpec struct {
	Nodes        int           // Number of nodes
	TasksPerNode int           // MPI ranks / tasks per node
	CpusPerTask  int           // CPU threads per task
	MemPerNodeMB int64         // Total RAM per node in MB
	Gpu          *GpuSpec      // GPU requirements (nil = no GPU; Count = per node)
	Time         time.Duration // Job walltime limit
	Exclusive    bool          // Request exclusive node access (no other jobs on the same node)
}

// RuntimeConfig holds job-level control settings (name, I/O paths, notifications).
type RuntimeConfig struct {
	JobName      string // Identifier for the job
	Stdout       string // Standard output file path
	Stderr       string // Standard error file path
	EmailOnBegin bool   // Notification on job start
	EmailOnEnd   bool   // Notification on job end
	EmailOnFail  bool   // Notification on job failure/abort
	MailUser     string // Target email/user for notifications (empty = submitting user)
	Partition    string // Partition/queue to submit to (cleared on cross-scheduler translation)
}

// ClusterInfo holds cluster configuration information
type ClusterInfo struct {
	AvailableGpus   []GpuInfo        // Available GPU types
	Limits          []ResourceLimits // Resource limits per partition
	MaxCpusPerNode  int              // Maximum CPUs available per node (from node info)
	MaxMemMBPerNode int64            // Maximum memory available per node in MB (from node info)
}

// ScriptSpecs holds the full result of parsing a scheduler script.
type ScriptSpecs struct {
	ScriptPath     string        // Absolute path of the parsed script (for HTCondor .sub: the executable)
	Spec           *ResourceSpec // Compute geometry (nil = resource parse failed → passthrough mode)
	Control        RuntimeConfig // Job control settings
	HasDirectives  bool          // True if any scheduler directive (#SBATCH, #PBS, etc.) was found
	RawFlags       []string      // ALL original directives — immutable audit log
	RemainingFlags []string      // Directives not absorbed by Spec or Control
	ScriptType     SchedulerType // Scheduler type detected from script directives
}

// HasSchedulerSpecs returns true if ScriptSpecs contains meaningful scheduler directives.
// This is used to determine if a script should be submitted to a scheduler,
// rather than relying on RawFlags which may be cleared during cross-scheduler translation.
func HasSchedulerSpecs(specs *ScriptSpecs) bool {
	if specs == nil {
		return false
	}
	// Check if any scheduler directive was found during parsing
	return specs.HasDirectives
}

// IsPassthrough reports having directives but no valid resource spec.
func IsPassthrough(specs *ScriptSpecs) bool {
	return specs != nil && specs.HasDirectives && specs.Spec == nil
}

// specDefaults is the package-level defaults used by all schedulers.
var specDefaults = ResourceSpec{
	CpusPerTask:  2,
	MemPerNodeMB: 8192,
	Time:         4 * time.Hour,
	Nodes:        1,
	TasksPerNode: 1,
}

// SetSpecDefaults overrides the default values used by ReadScriptSpecs.
func SetSpecDefaults(d ResourceSpec) {
	specDefaults = d
}

// GetSpecDefaults returns the current default values for ScriptSpecs.
func GetSpecDefaults() ResourceSpec {
	return specDefaults
}

// Override updates rs in-place with non-zero/non-nil fields from other.
// Fields in other that are zero/nil are skipped (rs retains its value).
// Safe to call with a nil other.
func (rs *ResourceSpec) Override(other *ResourceSpec) {
	if other == nil {
		return
	}
	if other.Nodes > 0 {
		rs.Nodes = other.Nodes
	}
	if other.TasksPerNode > 0 {
		rs.TasksPerNode = other.TasksPerNode
	}
	if other.CpusPerTask > 0 {
		rs.CpusPerTask = other.CpusPerTask
	}
	if other.MemPerNodeMB > 0 {
		rs.MemPerNodeMB = other.MemPerNodeMB
	}
	if other.Gpu != nil {
		rs.Gpu = other.Gpu
	}
	if other.Time > 0 {
		rs.Time = other.Time
	}
	if other.Exclusive {
		rs.Exclusive = other.Exclusive
	}
}

// ArraySpec describes an input-file-driven array job.
// Each task processes one line from InputFile using the scheduler's
// native array mechanism. Limit=0 means no concurrency cap.
type ArraySpec struct {
	InputFile  string   // Absolute path to input list (one entry per line)
	Count      int      // Non-empty line count; computed by cmd/run.go
	Limit      int      // Max concurrently running tasks (0 = unlimited)
	ArgCount   int      // Number of shell-split tokens per line (all lines must match)
	SampleArgs []string // Tokens from the first line, for dry-run display
	BlankLines []int    // 1-based line numbers of all blank lines (nil = none)
}

// JobSpec represents specifications for submitting a batch job
type JobSpec struct {
	Name           string            // Job name (for temp and log file naming)
	Command        string            // Command to execute
	Specs          *ScriptSpecs      // Job specifications
	DepJobIDs      []string          // Job IDs this job depends on
	Metadata       map[string]string // Additional metadata: ScriptPath, BuildSource, etc.
	OverrideOutput bool              // If true, always set Stdout/Stderr from Name (ignores script directives)
	Array          *ArraySpec        // Non-nil → emit array job directives
}

// Scheduler defines the interface for job schedulers
type Scheduler interface {
	// IsAvailable checks if the scheduler is available and we're not already in a job
	IsAvailable() bool

	// ReadScriptSpecs parses scheduler-specific directives from a build script
	// Returns ScriptSpecs with parsed scheduler directives (excludes #DEP and module directives)
	ReadScriptSpecs(scriptPath string) (*ScriptSpecs, error)

	// CreateScriptWithSpec generates a batch script with the given specifications
	// Returns the path to the created script
	CreateScriptWithSpec(spec *JobSpec, outputDir string) (string, error)

	// Submit submits a job script with optional dependency chain
	// Returns the job ID assigned by the scheduler
	Submit(scriptPath string, dependencyJobIDs []string) (string, error)

	// GetClusterInfo retrieves cluster configuration (GPUs, limits)
	// Returns nil if information is not available
	GetClusterInfo() (*ClusterInfo, error)

	// GetInfo returns information about the scheduler
	GetInfo() *SchedulerInfo

	// GetJobResources reads allocated resources from scheduler environment variables.
	// Fields with value 0 were not exposed by the scheduler.
	// Returns nil if not running inside a job of this scheduler type.
	GetJobResources() *ResourceSpec
}

// ResolveResourceSpecFrom merges resources using the priority chain:
//
//	defaults → scriptSpecs.Spec (when HasDirectives=true) → jobResources (non-zero fields)
//
// Checks HasDirectives rather than Spec != nil: when a script has no directives,
// the parser fills Spec with GetSpecDefaults(), so checking only Spec != nil would
// silently replace the caller's defaults with scheduler defaults.
// Always returns a non-nil *ResourceSpec.
func ResolveResourceSpecFrom(defaults ResourceSpec, jobResources *ResourceSpec, scriptSpecs *ScriptSpecs) *ResourceSpec {
	base := defaults
	if scriptSpecs != nil && scriptSpecs.HasDirectives && scriptSpecs.Spec != nil {
		base = *scriptSpecs.Spec
	}
	base.Override(jobResources)
	return &base
}

// ResolveResourceSpec merges resources using the priority chain:
//
//	GetSpecDefaults() → scriptSpecs.Spec (when HasDirectives=true) → jobResources
//
// Always returns a non-nil *ResourceSpec.
func ResolveResourceSpec(jobResources *ResourceSpec, scriptSpecs *ScriptSpecs) *ResourceSpec {
	return ResolveResourceSpecFrom(GetSpecDefaults(), jobResources, scriptSpecs)
}

// ValidateSpecs validates job specs against cluster limits.
// Skips validation when specs.Spec is nil (passthrough mode).
func ValidateSpecs(specs *ScriptSpecs, limits *ResourceLimits) error {
	if limits == nil || specs == nil || specs.Spec == nil {
		return nil // No limits or passthrough mode — skip validation
	}
	s := specs.Spec

	// CPUs per node = CpusPerTask * TasksPerNode
	cpusPerNode := s.CpusPerTask * s.TasksPerNode
	if cpusPerNode <= 0 {
		cpusPerNode = s.CpusPerTask
	}
	if limits.MaxCpusPerNode > 0 && cpusPerNode > limits.MaxCpusPerNode {
		return &ValidationError{
			Field:     "CpusPerNode",
			Requested: cpusPerNode,
			Limit:     limits.MaxCpusPerNode,
			Partition: limits.Partition,
		}
	}

	if limits.MaxMemMBPerNode > 0 && s.MemPerNodeMB > limits.MaxMemMBPerNode {
		return &ValidationError{
			Field:     "MemPerNodeMB",
			Requested: int(s.MemPerNodeMB),
			Limit:     int(limits.MaxMemMBPerNode),
			Partition: limits.Partition,
		}
	}

	if limits.MaxTime > 0 && s.Time > limits.MaxTime {
		return &ValidationError{
			Field:     "Time",
			Requested: int(s.Time.Hours()),
			Limit:     int(limits.MaxTime.Hours()),
			Partition: limits.Partition,
		}
	}

	// GPU count validation: Gpu.Count is per-node; compare total against MaxGpus
	if s.Gpu != nil && limits.MaxGpus > 0 {
		totalGpus := s.Gpu.Count * s.Nodes
		if totalGpus > limits.MaxGpus {
			return &ValidationError{
				Field:     "Gpu",
				Requested: totalGpus,
				Limit:     limits.MaxGpus,
				Partition: limits.Partition,
			}
		}
	}

	return nil
}

// ValidateAndConvertSpecs validates job specs against cluster limits.
// Returns a descriptive error if any resource exceeds the partition limit.
// Never modifies specs in-place.
// Returns nil if cluster info is unavailable (validation skipped).
func ValidateAndConvertSpecs(specs *ScriptSpecs) error {
	sched, err := DetectScheduler()
	if err != nil {
		return nil
	}
	clusterInfo, err := sched.GetClusterInfo()
	if err != nil {
		return nil
	}
	// Skip all resource validation in passthrough mode
	if specs.Spec == nil {
		return nil
	}

	// Filter limits to the requested partition when specified.
	// When Control.Partition is set, only validate against that partition's limits.
	// Falls back to all limits if the requested partition is not found.
	activeLimits := clusterInfo.Limits
	if specs.Control.Partition != "" {
		filtered := make([]ResourceLimits, 0, 1)
		for _, limit := range clusterInfo.Limits {
			if limit.Partition == specs.Control.Partition {
				filtered = append(filtered, limit)
				break
			}
		}
		if len(filtered) > 0 {
			activeLimits = filtered
		}
	}

	// Validate all resources (CPU, memory, time, GPU count) against active partition limits.
	// Specs are valid if they fit at least one partition.
	if len(activeLimits) > 0 {
		validForSomePartition := false
		var firstValidationErr error

		for _, limit := range activeLimits {
			if err := ValidateSpecs(specs, &limit); err == nil {
				validForSomePartition = true
				break
			} else if firstValidationErr == nil {
				firstValidationErr = err
			}
		}

		if !validForSomePartition {
			return firstValidationErr
		}
	}

	// Validate GPU availability — no conversion, just report error with suggestions.
	nodes := specs.Spec.Nodes
	if nodes <= 0 {
		nodes = 1
	}
	if specs.Spec.Gpu != nil && len(clusterInfo.AvailableGpus) > 0 {
		if err := ValidateGpuAvailability(specs.Spec.Gpu, nodes, clusterInfo); err != nil {
			return err // GpuValidationError includes list of available alternatives
		}
	}

	return nil
}

// SubmitWithDependencies submits multiple jobs with dependency chain
func SubmitWithDependencies(scheduler Scheduler, jobs []*JobSpec, outputDir string) (map[string]string, error) {
	jobIDs := make(map[string]string)

	for _, jobSpec := range jobs {
		// Use explicit job ID dependencies only
		depIDs := jobSpec.DepJobIDs

		// Create batch script
		scriptPath, err := scheduler.CreateScriptWithSpec(jobSpec, outputDir)
		if err != nil {
			return nil, fmt.Errorf("failed to create script for %s: %w", jobSpec.Name, err)
		}

		// Submit job
		jobID, err := scheduler.Submit(scriptPath, depIDs)
		if err != nil {
			return nil, fmt.Errorf("failed to submit job %s: %w", jobSpec.Name, err)
		}

		// Store job ID
		jobIDs[jobSpec.Name] = jobID

		fmt.Printf("Submitted job %s with ID %s", jobSpec.Name, jobID)
		if len(depIDs) > 0 {
			fmt.Printf(" (depends on: %s)", strings.Join(depIDs, ", "))
		}
		fmt.Println()
	}

	return jobIDs, nil
}

// DetectScheduler attempts to detect and return an available scheduler.
// Returns the scheduler instance if available, otherwise returns ErrSchedulerNotAvailable or ErrSchedulerNotFound.
func DetectScheduler() (Scheduler, error) {
	sched, err := DetectSchedulerWithBinary("")
	if err != nil {
		return nil, err
	}
	if !sched.IsAvailable() {
		return nil, ErrSchedulerNotAvailable
	}
	return sched, nil
}

// DetectSchedulerWithBinary attempts to initialize a scheduler using a preferred binary path.
// If preferredBin is empty, detection falls back to the default discovery path.
// This function returns a Scheduler instance if the scheduler binary is present, regardless of availability.
// Use DetectScheduler to require availability (not inside a job and submission enabled).
func DetectSchedulerWithBinary(preferredBin string) (Scheduler, error) {
	// If a preferred binary is specified, infer scheduler type from the binary name
	if preferredBin != "" {
		baseName := filepath.Base(preferredBin)
		switch baseName {
		case "qsub", "qdel", "qstat":
			return NewPbsSchedulerWithBinary(preferredBin)
		case "bsub", "bjobs", "bkill":
			return NewLsfSchedulerWithBinary(preferredBin)
		case "condor_submit", "condor_q", "condor_status":
			return NewHTCondorSchedulerWithBinary(preferredBin)
		default:
			// Default to SLURM for sbatch and any other binary
			return NewSlurmSchedulerWithBinary(preferredBin)
		}
	}

	// Try SLURM via PATH (most common)
	slurm, err := NewSlurmScheduler()
	if err == nil {
		return slurm, nil
	}

	// Try PBS via PATH
	pbs, pbsErr := NewPbsScheduler()
	if pbsErr == nil {
		return pbs, nil
	}

	// Try LSF via PATH
	lsf, lsfErr := NewLsfScheduler()
	if lsfErr == nil {
		return lsf, nil
	}

	// Try HTCondor via PATH
	htcondor, htcondorErr := NewHTCondorScheduler()
	if htcondorErr == nil {
		return htcondor, nil
	}

	return nil, ErrSchedulerNotFound
}

// CheckAvailability checks if a scheduler is available on the system
// Returns true if a scheduler is available and we're not inside a job
func CheckAvailability() bool {
	scheduler, err := DetectScheduler()
	if err != nil {
		return false
	}
	return scheduler.IsAvailable()
}

// ReadMaxSpecs reads cluster information and returns the maximum resource limits
// Returns the most permissive limits across all partitions
func ReadMaxSpecs(scheduler Scheduler) (*ResourceLimits, error) {
	info, err := scheduler.GetClusterInfo()
	if err != nil {
		return nil, err
	}

	if len(info.Limits) == 0 {
		return nil, ErrClusterInfoUnavailable
	}

	// Find the most permissive limits across all partitions
	maxLimits := &ResourceLimits{
		Partition: "all",
	}

	for _, limit := range info.Limits {
		if limit.MaxCpusPerNode > maxLimits.MaxCpusPerNode {
			maxLimits.MaxCpusPerNode = limit.MaxCpusPerNode
		}
		if limit.MaxMemMBPerNode > maxLimits.MaxMemMBPerNode {
			maxLimits.MaxMemMBPerNode = limit.MaxMemMBPerNode
		}
		if limit.MaxGpus > maxLimits.MaxGpus {
			maxLimits.MaxGpus = limit.MaxGpus
		}
		if limit.MaxNodes > maxLimits.MaxNodes {
			maxLimits.MaxNodes = limit.MaxNodes
		}
		if limit.MaxTime > maxLimits.MaxTime {
			maxLimits.MaxTime = limit.MaxTime
		}
		if limit.DefaultTime > maxLimits.DefaultTime {
			maxLimits.DefaultTime = limit.DefaultTime
		}
	}

	return maxLimits, nil
}

// ReadScript reads and parses scheduler directives from a build script
// This is an alias for ReadScriptSpecs for convenience
func ReadScript(scheduler Scheduler, scriptPath string) (*ScriptSpecs, error) {
	return scheduler.ReadScriptSpecs(scriptPath)
}

// ValidateScript validates a script's specs against cluster limits
func ValidateScript(scheduler Scheduler, scriptPath string) error {
	// Read script specs
	specs, err := scheduler.ReadScriptSpecs(scriptPath)
	if err != nil {
		return err
	}

	// Get cluster info
	info, err := scheduler.GetClusterInfo()
	if err != nil {
		// If we can't get cluster info, we can't validate
		// This is not necessarily an error (cluster might not expose this info)
		return nil
	}

	// Validate against the requested partition (or all if none specified)
	if len(info.Limits) == 0 {
		return nil // No limits to validate against
	}

	limits := info.Limits
	if specs != nil && specs.Control.Partition != "" {
		for _, limit := range info.Limits {
			if limit.Partition == specs.Control.Partition {
				return ValidateSpecs(specs, &limit)
			}
		}
		// Requested partition not found — fall through to all-partition check
	}

	// Validate against each partition's limits; pass if any accepts the job
	for _, limit := range limits {
		if err := ValidateSpecs(specs, &limit); err != nil {
			continue
		}
		return nil
	}

	return ValidateSpecs(specs, &limits[0])
}

// WriteScript creates a batch script with the given specifications
// This is an alias for CreateScriptWithSpec for convenience
func WriteScript(scheduler Scheduler, jobSpec *JobSpec, outputDir string) (string, error) {
	return scheduler.CreateScriptWithSpec(jobSpec, outputDir)
}

// SubmitJob submits a single job and returns its job ID
func SubmitJob(scheduler Scheduler, scriptPath string, dependencyJobIDs []string) (string, error) {
	return scheduler.Submit(scriptPath, dependencyJobIDs)
}

// SubmitJobs submits multiple jobs with dependency chains and returns their job IDs
func SubmitJobs(scheduler Scheduler, jobs []*JobSpec, outputDir string) (map[string]string, error) {
	return SubmitWithDependencies(scheduler, jobs, outputDir)
}

// Init auto-detects and initializes the active scheduler.
// If preferredBin is provided, it will be used instead of auto-detection.
// Returns the detected scheduler type and any error.
func Init(preferredBin string) (SchedulerType, error) {
	sched, err := DetectSchedulerWithBinary(preferredBin)
	if err != nil {
		ClearActiveScheduler()
		return SchedulerUnknown, err
	}

	SetActiveScheduler(sched)

	// Determine scheduler type from info
	info := sched.GetInfo()
	return SchedulerType(info.Type), nil
}

// InitIfAvailable attempts to initialize a scheduler if one is available.
// Unlike Init, this does not return an error if no scheduler is found.
// Returns the detected scheduler type (or SchedulerUnknown if none).
func InitIfAvailable(preferredBin string) SchedulerType {
	schedType, _ := Init(preferredBin)
	return schedType
}

// DetectType returns the type of scheduler available on the system without initializing it.
// This is useful for checking what scheduler is available before deciding to use it.
func DetectType() SchedulerType {
	// Check for SLURM (sbatch)
	if _, err := exec.LookPath("sbatch"); err == nil {
		return SchedulerSLURM
	}

	// Check for PBS (qsub with PBS-specific behavior)
	if _, err := exec.LookPath("qsub"); err == nil {
		// TODO: Distinguish PBS from other qsub implementations (SGE, etc.)
		return SchedulerPBS
	}

	// Check for LSF (bsub)
	if _, err := exec.LookPath("bsub"); err == nil {
		return SchedulerLSF
	}

	// Check for HTCondor (condor_submit)
	if _, err := exec.LookPath("condor_submit"); err == nil {
		return SchedulerHTCondor
	}

	return SchedulerUnknown
}

// IsInsideJob checks if we're currently running inside a scheduler job.
// This is useful to avoid nested job submission.
func IsInsideJob() bool {
	// Check SLURM
	if _, ok := os.LookupEnv("SLURM_JOB_ID"); ok {
		return true
	}
	// Check PBS
	if _, ok := os.LookupEnv("PBS_JOBID"); ok {
		return true
	}
	// Check LSF
	if _, ok := os.LookupEnv("LSB_JOBID"); ok {
		return true
	}
	// Check HTCondor
	if _, ok := os.LookupEnv("_CONDOR_JOB_AD"); ok {
		return true
	}
	return false
}

// ParsedScript contains the normalized specs and detected script type
type ParsedScript struct {
	Specs      *ScriptSpecs  // Normalized scheduler specifications
	ScriptType SchedulerType // Detected scheduler type from script directives
}

// ParseScriptAny parses scheduler directives from a script using the appropriate parser.
// It tries the current scheduler first, then falls back to other schedulers.
// This enables cross-scheduler compatibility (e.g., PBS can run SLURM scripts).
//
// Returns:
//   - ParsedScript with specs and detected type if directives found
//   - nil ParsedScript if no scheduler directives found (not an error)
//   - error if parsing fails
func ParseScriptAny(scriptPath string) (*ParsedScript, error) {
	// Determine current scheduler type
	currentType := DetectType()

	// Try current scheduler first, then others
	tryOrder := []SchedulerType{}

	// Add current scheduler first (if known)
	if currentType != SchedulerUnknown {
		tryOrder = append(tryOrder, currentType)
	}

	// Add other schedulers
	allTypes := []SchedulerType{SchedulerSLURM, SchedulerPBS, SchedulerLSF, SchedulerHTCondor}
	for _, st := range allTypes {
		if st != currentType {
			tryOrder = append(tryOrder, st)
		}
	}

	// Try each scheduler parser
	for _, schedType := range tryOrder {
		var specs *ScriptSpecs
		var err error

		switch schedType {
		case SchedulerSLURM:
			specs, err = TryParseSlurmScript(scriptPath)
		case SchedulerPBS:
			specs, err = TryParsePbsScript(scriptPath)
		case SchedulerLSF:
			specs, err = TryParseLsfScript(scriptPath)
		case SchedulerHTCondor:
			// Only parse as HTCondor if the file has the native .sub extension.
			if strings.HasSuffix(scriptPath, ".sub") {
				specs, err = TryParseHTCondorScript(scriptPath)
			}
		default:
			continue
		}

		if err != nil {
			return nil, err
		}

		// If we found directives, return the result
		if specs != nil && HasSchedulerSpecs(specs) {
			return &ParsedScript{
				Specs:      specs,
				ScriptType: schedType,
			}, nil
		}
	}

	// No scheduler directives found in any format
	return nil, nil
}

// getEnvInt reads an environment variable and parses it as a positive int.
// Returns nil if unset, empty, or not a valid positive integer.
func getEnvInt(key string) *int {
	val := os.Getenv(key)
	if val == "" {
		return nil
	}
	n, err := strconv.Atoi(val)
	if err != nil || n <= 0 {
		return nil
	}
	return &n
}

// getEnvInt64 reads an environment variable and parses it as a positive int64.
// Returns nil if unset, empty, or not a valid positive integer.
func getEnvInt64(key string) *int64 {
	val := os.Getenv(key)
	if val == "" {
		return nil
	}
	n, err := strconv.ParseInt(val, 10, 64)
	if err != nil || n <= 0 {
		return nil
	}
	return &n
}

// getCudaDeviceCount parses CUDA_VISIBLE_DEVICES and returns the number of devices.
// Returns nil if the variable is unset or empty.
func getCudaDeviceCount() *int {
	val := os.Getenv("CUDA_VISIBLE_DEVICES")
	if val == "" {
		return nil
	}
	// Count comma-separated items (e.g., "0,1,2" → 3)
	count := len(strings.Split(val, ","))
	if count <= 0 {
		return nil
	}
	return &count
}

// ReadScriptSpecsFromPath reads scheduler specs from a script file.
// It uses cross-scheduler translation: parses any scheduler format and returns
// normalized specs that can be used with the active scheduler.
//
// If the script's scheduler type differs from the host scheduler, a warning
// message is returned (but not an error) to allow the build to proceed.
//
// Always returns a non-nil *ScriptSpecs. When no directives are found,
// HasDirectives is false and Spec is filled with GetSpecDefaults() values (Nodes=1).
// Use HasDirectives and IsPassthrough to distinguish the three states.
func ReadScriptSpecsFromPath(scriptPath string) (*ScriptSpecs, error) {
	// Parse script using any available parser
	parsed, err := ParseScriptAny(scriptPath)
	if err != nil {
		return nil, err
	}

	// No scheduler directives found — return a default spec with HasDirectives=false.
	// Callers must check HasDirectives (not nil) to distinguish from normal/passthrough mode.
	if parsed == nil {
		d := GetSpecDefaults()
		d.Nodes = 1 // single-node default for local scripts
		return &ScriptSpecs{
			ScriptPath:    scriptPath,
			HasDirectives: false,
			Spec:          &d,
		}, nil
	}

	// Carry the original scheduler type forward so callers can inspect it.
	parsed.Specs.ScriptType = parsed.ScriptType

	// Check for scheduler mismatch
	hostType := DetectType()
	if hostType != SchedulerUnknown && parsed.ScriptType != hostType {
		// Clear RemainingFlags on cross-scheduler translation:
		// scheduler-specific unrecognized flags cannot be translated.
		// RawFlags is the immutable audit log — never cleared.
		parsed.Specs.RemainingFlags = nil
		// Clear Partition: the original partition name is scheduler-specific
		// and cannot be translated to the host scheduler's partition names.
		parsed.Specs.Control.Partition = ""
	}

	return parsed.Specs, nil
}
