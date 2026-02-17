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

	"github.com/Justype/condatainer/internal/utils"
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

// GpuSpec holds GPU requirements parsed from scheduler directives
type GpuSpec struct {
	Type  string // GPU type/model (e.g., "h100", "a100", "v100")
	Count int    // Number of GPUs requested
	Raw   string // Raw GPU specification string from scheduler
}

// GpuInfo holds information about available GPU types in the cluster
type GpuInfo struct {
	Type      string // GPU type/model
	Total     int    // Total number of this GPU type
	Available int    // Currently available
	Partition string // Partition/queue name (optional)
}

// ResourceLimits holds scheduler resource limits
type ResourceLimits struct {
	MaxCpus     int           // Maximum CPUs per job
	MaxMemMB    int64         // Maximum memory in MB per job
	MaxTime     time.Duration // Maximum walltime (e.g., "7-00:00:00")
	MaxGpus     int           // Maximum GPUs per job
	MaxNodes    int           // Maximum nodes per job
	DefaultTime time.Duration // Default walltime if not specified
	Partition   string        // Partition/queue name these limits apply to
}

// ClusterInfo holds cluster configuration information
type ClusterInfo struct {
	AvailableGpus   []GpuInfo        // Available GPU types
	Limits          []ResourceLimits // Resource limits per partition
	MaxCpusPerNode  int              // Maximum CPUs available per node (from node info)
	MaxMemMBPerNode int64            // Maximum memory available per node in MB (from node info)
}

// ScriptSpecs holds the specifications parsed from a job script
type ScriptSpecs struct {
	JobName       string        // Job name
	Ncpus         int           // Number of CPUs per task
	Ntasks        int           // Number of tasks (default 1)
	Nodes         int           // Number of nodes (default 1)
	MemMB         int64         // Memory in MB
	Time          time.Duration // Time limit
	Stdout        string        // Standard output file path
	Stderr        string        // Standard error file path
	Gpu           *GpuSpec      // GPU requirements (nil if no GPU)
	EmailOnBegin  bool          // Send email when job begins (SLURM: BEGIN, PBS: b, LSF: -B)
	EmailOnEnd    bool          // Send email when job ends (SLURM: END, PBS: e, LSF: -N)
	EmailOnFail   bool          // Send email when job fails/aborts (SLURM: FAIL, PBS: a)
	MailUser      string        // Username or email address for notifications (empty = submitting user)
	HasDirectives bool          // True if any scheduler directive was found during parsing
	RawFlags      []string      // Raw scheduler-specific flags (unrecognized only)
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

// SpecDefaults holds configurable default values for ScriptSpecs.
// These are used by ReadScriptSpecs when a script does not specify a resource.
// Set via SetSpecDefaults() during CLI initialization from config.
type SpecDefaults struct {
	Ncpus  int
	MemMB  int64
	Time   time.Duration
	Nodes  int
	Ntasks int
}

// specDefaults is the package-level defaults used by all schedulers.
var specDefaults = SpecDefaults{
	Ncpus:  2,
	MemMB:  8192,
	Time:   4 * time.Hour,
	Nodes:  1,
	Ntasks: 1,
}

// SetSpecDefaults overrides the default values used by ReadScriptSpecs.
func SetSpecDefaults(d SpecDefaults) {
	specDefaults = d
}

// GetSpecDefaults returns the current default values for ScriptSpecs.
func GetSpecDefaults() SpecDefaults {
	return specDefaults
}

// JobSpec represents specifications for submitting a batch job
type JobSpec struct {
	Name      string            // Job name
	Command   string            // Command to execute
	Specs     *ScriptSpecs      // Job specifications
	DepJobIDs []string          // Job IDs this job depends on
	Metadata  map[string]string // Additional metadata: ScriptPath, BuildSource, etc.
}

// JobResources holds resource allocations for the currently running scheduler job.
// A nil pointer field means the scheduler did not expose that resource via environment variables.
type JobResources struct {
	Ncpus  *int   // Number of allocated CPUs per task
	Ntasks *int   // Number of allocated tasks
	Nodes  *int   // Number of allocated nodes
	MemMB  *int64 // Allocated memory in MB
	Ngpus  *int   // Number of allocated GPUs
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
	// Returns nil if not running inside a job of this scheduler type.
	GetJobResources() *JobResources
}

// ValidateSpecs validates job specs against cluster limits
func ValidateSpecs(specs *ScriptSpecs, limits *ResourceLimits) error {
	if limits == nil {
		return nil // No limits to validate against
	}

	if limits.MaxCpus > 0 && specs.Ncpus > limits.MaxCpus {
		return &ValidationError{
			Field:     "Ncpus",
			Requested: specs.Ncpus,
			Limit:     limits.MaxCpus,
			Partition: limits.Partition,
		}
	}

	if limits.MaxMemMB > 0 && specs.MemMB > limits.MaxMemMB {
		return &ValidationError{
			Field:     "MemMB",
			Requested: int(specs.MemMB),
			Limit:     int(limits.MaxMemMB),
			Partition: limits.Partition,
		}
	}

	if limits.MaxTime > 0 && specs.Time > limits.MaxTime {
		return &ValidationError{
			Field:     "Time",
			Requested: int(specs.Time.Hours()),
			Limit:     int(limits.MaxTime.Hours()),
			Partition: limits.Partition,
		}
	}

	if specs.Gpu != nil && limits.MaxGpus > 0 && specs.Gpu.Count > limits.MaxGpus {
		return &ValidationError{
			Field:     "Gpu",
			Requested: specs.Gpu.Count,
			Limit:     limits.MaxGpus,
			Partition: limits.Partition,
		}
	}

	return nil
}

// ValidateAndConvertSpecs validates job specs against cluster limits and attempts conversions if needed.
// - CPUs: Automatically reduced to fit within limits (prints warning)
//   - Time is scaled proportionally based on CPU reduction
//   - Fails if scaled time exceeds cluster time limit
//
// - Memory: Critical - fails if exceeded
// - Time: Critical - fails if exceeded
// - GPUs: Attempts conversion to compatible type
//   - Upgrades: Allowed with warning
//   - Downgrades: Blocked if time limit specified (runtime unpredictable)
//   - No time limit: Downgrades allowed with warning
//
// Modifies specs in-place if conversions succeed.
// Returns validation errors and conversion details.
// If cluster info is unavailable, returns nil error (validation skipped).
func ValidateAndConvertSpecs(specs *ScriptSpecs) (validationErr error, cpuAdjusted bool, cpuMsg string, gpuConverted bool, gpuMsg string) {
	// Detect scheduler
	sched, err := DetectScheduler()
	if err != nil {
		// No scheduler available - skip validation
		return nil, false, "", false, ""
	}

	// Get cluster info
	clusterInfo, err := sched.GetClusterInfo()
	if err != nil {
		// Can't get cluster info - skip validation
		return nil, false, "", false, ""
	}

	// Auto-adjust CPU if it exceeds limits in all partitions
	// When reducing CPUs, also adjust time estimate proportionally
	if len(clusterInfo.Limits) > 0 && specs.Ncpus > 0 {
		maxAllowedCpu := 0
		maxAllowedTime := time.Duration(0)

		for _, limit := range clusterInfo.Limits {
			if limit.MaxCpus > maxAllowedCpu {
				maxAllowedCpu = limit.MaxCpus
			}
			if limit.MaxTime > maxAllowedTime {
				maxAllowedTime = limit.MaxTime
			}
		}

		if maxAllowedCpu > 0 && specs.Ncpus > maxAllowedCpu {
			originalCpu := specs.Ncpus
			originalTime := specs.Time

			// Reduce CPUs to fit limit
			specs.Ncpus = maxAllowedCpu

			// Adjust time proportionally (assuming linear scaling)
			// If we reduce CPUs by half, job might take twice as long
			if originalTime > 0 {
				cpuReductionRatio := float64(originalCpu) / float64(maxAllowedCpu)
				adjustedTime := time.Duration(float64(originalTime) * cpuReductionRatio)

				// Check if adjusted time exceeds cluster limit
				if maxAllowedTime > 0 && adjustedTime > maxAllowedTime {
					// Cannot reduce CPUs because it would make the job exceed time limit
					specs.Ncpus = originalCpu // Restore original value
					return &ValidationError{
						Field:     "Ncpus/Time",
						Requested: originalCpu,
						Limit:     maxAllowedCpu,
						Partition: fmt.Sprintf("reducing CPUs to %d would require %s, exceeding time limit %s",
							maxAllowedCpu, adjustedTime.Round(time.Minute), maxAllowedTime.Round(time.Minute)),
					}, false, "", false, ""
				}

				// Update time estimate
				specs.Time = adjustedTime
				cpuAdjusted = true
				cpuMsg = fmt.Sprintf("CPUs: %d → %d (reduced to fit limit); Time: %s → %s (adjusted proportionally)",
					originalCpu, maxAllowedCpu,
					originalTime.Round(time.Minute), adjustedTime.Round(time.Minute))
			} else {
				// No time specified, just reduce CPUs
				cpuAdjusted = true
				cpuMsg = fmt.Sprintf("CPUs: %d → %d (reduced to fit cluster limit)", originalCpu, maxAllowedCpu)
			}
		}
	}

	// Validate resource limits (memory and other critical resources) against all partitions
	if len(clusterInfo.Limits) > 0 {
		validForSomePartition := false
		var firstValidationErr error

		for _, limit := range clusterInfo.Limits {
			if err := ValidateSpecs(specs, &limit); err == nil {
				validForSomePartition = true
				break
			} else if firstValidationErr == nil {
				firstValidationErr = err
			}
		}

		if !validForSomePartition {
			// Job exceeds limits in all partitions (likely memory or GPU count)
			return firstValidationErr, cpuAdjusted, cpuMsg, false, ""
		}
	}

	// Validate and convert GPU specs if needed
	if specs.Gpu != nil && len(clusterInfo.AvailableGpus) > 0 {
		originalGpuType := specs.Gpu.Type
		originalGpuCount := specs.Gpu.Count

		// Check if requested GPU is available
		err := ValidateGpuAvailability(specs.Gpu, clusterInfo)
		if err != nil {
			// GPU not available - attempt conversion
			convertedSpec, convErr := ConvertGpuSpec(specs.Gpu, clusterInfo)
			if convErr != nil {
				// Conversion failed - return the validation error
				return err, cpuAdjusted, cpuMsg, false, ""
			}

			// Conversion succeeded - update specs and build message
			specs.Gpu = convertedSpec

			// Determine if upgrade or downgrade
			conversionType := "alternative"
			isDowngrade := false
			options, _ := FindCompatibleGpu(specs.Gpu, clusterInfo)
			if len(options) > 0 && options[0].SuggestedSpec.Type == convertedSpec.Type {
				if options[0].IsUpgrade {
					conversionType = "upgrade"
				} else if options[0].IsDowngrade {
					conversionType = "downgrade"
					isDowngrade = true
				}
			}

			// GPU downgrades are risky when time limits exist because we can't estimate
			// how much longer the job will take on a slower GPU (performance varies by workload)
			if isDowngrade && specs.Time > 0 {
				// Restore original GPU spec
				specs.Gpu = &GpuSpec{Type: originalGpuType, Count: originalGpuCount}

				// Return error with explanation
				return &ValidationError{
					Field:     "Gpu",
					Requested: 0, // Not a count issue
					Limit:     0,
					Partition: fmt.Sprintf("GPU downgrade %s → %s not allowed with time limit %s (runtime unpredictable; consider removing time limit or requesting available GPU)",
						originalGpuType, convertedSpec.Type, specs.Time.Round(time.Minute)),
				}, cpuAdjusted, cpuMsg, false, ""
			}

			// Build conversion message
			if originalGpuCount > 1 {
				gpuMsg = fmt.Sprintf("GPUs: %dx %s → %dx %s (%s)",
					originalGpuCount, originalGpuType, convertedSpec.Count, convertedSpec.Type, conversionType)
			} else {
				gpuMsg = fmt.Sprintf("GPU: %s → %s (%s)",
					originalGpuType, convertedSpec.Type, conversionType)
			}

			return nil, cpuAdjusted, cpuMsg, true, gpuMsg
		}
	}

	return nil, cpuAdjusted, cpuMsg, false, ""
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
		if limit.MaxCpus > maxLimits.MaxCpus {
			maxLimits.MaxCpus = limit.MaxCpus
		}
		if limit.MaxMemMB > maxLimits.MaxMemMB {
			maxLimits.MaxMemMB = limit.MaxMemMB
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

	// Validate against all partition limits
	// If the script doesn't specify a partition, validate against all of them
	if len(info.Limits) == 0 {
		return nil // No limits to validate against
	}

	// Validate against each partition's limits
	for _, limit := range info.Limits {
		if err := ValidateSpecs(specs, &limit); err != nil {
			// If it exceeds limits in one partition, that's ok as long as
			// there's another partition that can handle it
			continue
		}
		// Found a partition that can handle this job
		return nil
	}

	// If we get here, the job exceeds limits in all partitions
	// Return validation error for the first partition
	return ValidateSpecs(specs, &info.Limits[0])
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
	// Check PBS/Torque
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
			specs, err = TryParseHTCondorScript(scriptPath)
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
// Returns nil specs (not an error) if no scheduler directives are found.
func ReadScriptSpecsFromPath(scriptPath string) (*ScriptSpecs, error) {
	// Parse script using any available parser
	parsed, err := ParseScriptAny(scriptPath)
	if err != nil {
		return nil, err
	}

	// No scheduler directives found
	if parsed == nil {
		return nil, nil
	}

	// Check for scheduler mismatch and log warning
	hostType := DetectType()
	if hostType != SchedulerUnknown && parsed.ScriptType != hostType {
		// Log warning about mismatch (the specs will still be used)
		utils.PrintWarning("Script contains %s directives but host has %s scheduler. Specs will be translated.",
			parsed.ScriptType, hostType)

		// Clear RawFlags when doing cross-scheduler translation
		// RawFlags contain scheduler-specific syntax that won't work on different schedulers
		// All necessary info is preserved in the parsed ScriptSpecs fields
		parsed.Specs.RawFlags = nil
	}

	return parsed.Specs, nil
}
