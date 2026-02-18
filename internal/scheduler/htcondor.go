package scheduler

import (
	"bufio"
	"fmt"
	"os"
	"os/exec"
	"path/filepath"
	"regexp"
	"strconv"
	"strings"
	"time"

	"github.com/Justype/condatainer/internal/config"
	"github.com/Justype/condatainer/internal/utils"
)

// HTCondorScheduler implements the Scheduler interface for HTCondor
// EXPERIMENTAL: HTCondor is not tested on real clusters and may have edge cases. Feedback welcome.
type HTCondorScheduler struct {
	condorSubmitBin    string
	condorQBin         string
	condorStatusBin    string
	condorConfigValBin string
	jobIDRe            *regexp.Regexp
}

// NewHTCondorScheduler creates a new HTCondor scheduler instance using condor_submit from PATH
func NewHTCondorScheduler() (*HTCondorScheduler, error) {
	return newHTCondorSchedulerWithBinary("")
}

// NewHTCondorSchedulerWithBinary creates an HTCondor scheduler using an explicit condor_submit path
func NewHTCondorSchedulerWithBinary(condorSubmitBin string) (*HTCondorScheduler, error) {
	return newHTCondorSchedulerWithBinary(condorSubmitBin)
}

func newHTCondorSchedulerWithBinary(condorSubmitBin string) (*HTCondorScheduler, error) {
	binPath := condorSubmitBin
	if binPath == "" {
		var err error
		binPath, err = exec.LookPath("condor_submit")
		if err != nil {
			return nil, fmt.Errorf("%w: %v", ErrSchedulerNotFound, err)
		}
	} else {
		if absPath, err := filepath.Abs(binPath); err == nil {
			binPath = absPath
		}
		info, err := os.Stat(binPath)
		if err != nil {
			return nil, fmt.Errorf("%w: %v", ErrSchedulerNotFound, err)
		}
		if info.IsDir() {
			return nil, fmt.Errorf("%w: %s is a directory", ErrSchedulerNotFound, binPath)
		}
	}

	condorQBin, _ := exec.LookPath("condor_q")
	condorStatusBin, _ := exec.LookPath("condor_status")
	condorConfigValBin, _ := exec.LookPath("condor_config_val")

	return &HTCondorScheduler{
		condorSubmitBin:    binPath,
		condorQBin:         condorQBin,
		condorStatusBin:    condorStatusBin,
		condorConfigValBin: condorConfigValBin,
		jobIDRe:            regexp.MustCompile(`submitted to cluster (\d+)`),
	}, nil
}

// IsAvailable checks if HTCondor is available and we're not inside an HTCondor job
func (h *HTCondorScheduler) IsAvailable() bool {
	if h.condorSubmitBin == "" {
		return false
	}

	// Check if we're already inside an HTCondor job
	if _, inJob := os.LookupEnv("_CONDOR_JOB_AD"); inJob {
		return false
	}

	return true
}

// GetInfo returns information about the HTCondor scheduler
func (h *HTCondorScheduler) GetInfo() *SchedulerInfo {
	_, inJob := os.LookupEnv("_CONDOR_JOB_AD")
	available := h.IsAvailable()

	info := &SchedulerInfo{
		Type:      "HTCondor",
		Binary:    h.condorSubmitBin,
		InJob:     inJob,
		Available: available,
	}

	// Try to get HTCondor version
	if h.condorSubmitBin != "" {
		if version, err := h.getHTCondorVersion(); err == nil {
			info.Version = version
		}
	}

	return info
}

// getHTCondorVersion attempts to get the HTCondor version
func (h *HTCondorScheduler) getHTCondorVersion() (string, error) {
	cmd := exec.Command(h.condorSubmitBin, "-version")
	output, err := cmd.Output()
	if err != nil {
		return "", err
	}

	// Parse version from output like "$CondorVersion: 10.0.0 ..."
	versionStr := strings.TrimSpace(string(output))
	lines := strings.Split(versionStr, "\n")
	if len(lines) > 0 {
		line := lines[0]
		// Try to extract version from $CondorVersion: X.Y.Z ...
		if strings.Contains(line, "$CondorVersion:") {
			parts := strings.Fields(line)
			for i, p := range parts {
				if p == "$CondorVersion:" && i+1 < len(parts) {
					return parts[i+1], nil
				}
			}
		}
		// Fallback: return first line
		return strings.TrimSpace(line), nil
	}

	return versionStr, nil
}

// ReadScriptSpecs parses an HTCondor submit file (.sub) in native key=value format.
// The executable line is extracted as the ScriptPath.
func (h *HTCondorScheduler) ReadScriptSpecs(scriptPath string) (*ScriptSpecs, error) {
	lines, err := readFileLines(scriptPath)
	if err != nil {
		return nil, err
	}
	executable := extractHTCondorExecutable(lines)
	return parseScript(executable, lines, h.extractDirectives, h.parseRuntimeConfig, h.parseResourceSpec)
}

// htcondorStructuralKeys are submit file keywords used by extractHTCondorExecutable
// to locate the executable line. Not used for directive extraction.
var htcondorStructuralKeys = map[string]bool{
	"universe":            true,
	"executable":          true,
	"transfer_executable": true,
	"queue":               true,
	"arguments":           true,
}

// htcondorKnownKeys is the whitelist of recognized HTCondor submit-file directive keys
// that carry resource or control information. Structural keys (universe, executable,
// queue, etc.) are intentionally excluded: they are handled separately by
// extractHTCondorExecutable and must not appear in the directive list consumed by
// parseRuntimeConfig / parseResourceSpec.
var htcondorKnownKeys = map[string]bool{
	// Runtime config (parseRuntimeConfig)
	"output": true, "log": true, "error": true,
	"notify_user": true, "accounting_group": true, "notification": true,
	// Resource spec (parseResourceSpec)
	"request_cpus": true, "request_memory": true, "request_gpus": true,
	"+maxruntime": true,
}

// extractDirectives parses native HTCondor submit file lines.
// Comments (#), empty lines, and lines whose key is not in htcondorKnownKeys are skipped.
// Returns key=value directive strings for parsing by parseRuntimeConfig and parseResourceSpec.
func (h *HTCondorScheduler) extractDirectives(lines []string) []string {
	var out []string
	for _, line := range lines {
		trimmed := strings.TrimSpace(line)
		// Skip comments, empty lines
		if trimmed == "" || strings.HasPrefix(trimmed, "#") {
			continue
		}
		// Parse key from key=value or bare keyword
		key := trimmed
		if idx := strings.IndexByte(trimmed, '='); idx >= 0 {
			key = strings.TrimSpace(trimmed[:idx])
		}
		// Only collect lines with a recognized HTCondor key.
		if !htcondorKnownKeys[strings.ToLower(key)] {
			continue
		}
		out = append(out, trimmed)
	}
	return out
}

// extractHTCondorExecutable finds the executable path from HTCondor submit file lines.
func extractHTCondorExecutable(lines []string) string {
	for _, line := range lines {
		trimmed := strings.TrimSpace(line)
		if strings.HasPrefix(trimmed, "#") || trimmed == "" {
			continue
		}
		parts := strings.SplitN(trimmed, "=", 2)
		if len(parts) == 2 && strings.TrimSpace(strings.ToLower(parts[0])) == "executable" {
			return strings.TrimSpace(parts[1])
		}
	}
	return ""
}

// parseRuntimeConfig consumes job control directives (name, I/O, email) from the directive list.
// Returns the populated RuntimeConfig, unconsumed directives, and any critical error.
func (h *HTCondorScheduler) parseRuntimeConfig(directives []string) (RuntimeConfig, []string, error) {
	var rc RuntimeConfig
	var unconsumed []string

	for _, flag := range directives {
		parts := strings.SplitN(flag, "=", 2)
		if len(parts) != 2 {
			unconsumed = append(unconsumed, flag)
			continue
		}

		key := strings.TrimSpace(strings.ToLower(parts[0]))
		value := strings.TrimSpace(parts[1])
		consumed := true

		switch key {
		case "output", "log":
			rc.Stdout = absPath(value)
		case "error":
			rc.Stderr = absPath(value)
		case "notify_user":
			rc.MailUser = value
		case "accounting_group":
			rc.Partition = value
		case "notification":
			switch strings.ToLower(value) {
			case "always":
				rc.EmailOnBegin = true
				rc.EmailOnEnd = true
				rc.EmailOnFail = true
			case "complete":
				rc.EmailOnEnd = true
			case "error":
				rc.EmailOnFail = true
			case "never":
				rc.EmailOnBegin = false
				rc.EmailOnEnd = false
				rc.EmailOnFail = false
			}
		default:
			consumed = false
		}

		if !consumed {
			unconsumed = append(unconsumed, flag)
		}
	}

	return rc, unconsumed, nil
}

// parseResourceSpec consumes compute resource directives from the directive list.
// Returns a populated ResourceSpec and any remaining (unrecognized) directives.
// Returns nil ResourceSpec on parse error (passthrough mode).
func (h *HTCondorScheduler) parseResourceSpec(directives []string) (*ResourceSpec, []string) {
	defaults := GetSpecDefaults()
	rs := &ResourceSpec{
		Nodes:        1, // HTCondor is inherently single-node
		TasksPerNode: 1,
		CpusPerTask:  defaults.CpusPerTask,
		MemPerNodeMB: defaults.MemPerNodeMB,
		Time:         defaults.Time,
	}

	var remaining []string

	for _, flag := range directives {
		parts := strings.SplitN(flag, "=", 2)
		if len(parts) != 2 {
			remaining = append(remaining, flag)
			continue
		}

		key := strings.TrimSpace(strings.ToLower(parts[0]))
		value := strings.TrimSpace(parts[1])
		recognized := true

		switch key {
		case "request_cpus":
			n, err := strconv.Atoi(value)
			if err != nil {
				utils.PrintWarning("HTCondor: invalid request_cpus value %q: %v", value, err)
				return nil, directives
			}
			rs.CpusPerTask = n

		case "request_memory":
			mem, err := parseHTCondorMemory(value)
			if err != nil {
				utils.PrintWarning("HTCondor: invalid request_memory value %q: %v", value, err)
				return nil, directives
			}
			rs.MemPerNodeMB = mem

		case "request_gpus":
			n, err := strconv.Atoi(value)
			if err != nil {
				utils.PrintWarning("HTCondor: invalid request_gpus value %q: %v", value, err)
				return nil, directives
			}
			rs.Gpu = &GpuSpec{
				Type:  "gpu",
				Count: n,
				Raw:   value,
			}

		case "+maxruntime":
			secs, err := strconv.ParseInt(value, 10, 64)
			if err != nil {
				utils.PrintWarning("HTCondor: invalid +MaxRuntime value %q: %v", value, err)
				return nil, directives
			}
			rs.Time = time.Duration(secs) * time.Second

		default:
			recognized = false
		}

		if !recognized {
			remaining = append(remaining, flag)
		}
	}

	return rs, remaining
}

// CreateScriptWithSpec generates an HTCondor submit description file and wrapper script
func (h *HTCondorScheduler) CreateScriptWithSpec(jobSpec *JobSpec, outputDir string) (string, error) {
	specs := jobSpec.Specs

	// Create output directory if specified
	if outputDir != "" {
		if err := os.MkdirAll(outputDir, utils.PermDir); err != nil {
			return "", NewScriptCreationError(jobSpec.Name, outputDir, err)
		}
	}

	// Generate safe name for files
	name := "job"
	if jobSpec.Name != "" {
		name = safeJobName(jobSpec.Name)
	}

	// Set log path based on job name; only override if caller requests it or script has no output set
	if jobSpec.Name != "" && (jobSpec.OverrideOutput || specs.Control.Stdout == "") {
		specs.Control.Stdout = filepath.Join(outputDir, fmt.Sprintf("%s.log", name))
	}
	if specs.Control.Stderr == "" && specs.Control.Stdout != "" {
		specs.Control.Stderr = specs.Control.Stdout
	}

	// === Create wrapper bash script ===
	shPath := filepath.Join(outputDir, fmt.Sprintf("%s.sh", name))
	shFile, err := os.Create(shPath)
	if err != nil {
		return "", NewScriptCreationError(jobSpec.Name, shPath, err)
	}

	shWriter := bufio.NewWriter(shFile)
	fmt.Fprintln(shWriter, "#!/bin/bash")
	fmt.Fprintln(shWriter, "")

	// HTCondor is always single-node; synthesize a ResourceSpec with Nodes=1,
	// TasksPerNode=1, and copy CpusPerTask/MemPerNodeMB from the actual spec if available.
	htcRS := &ResourceSpec{Nodes: 1, TasksPerNode: 1}
	if specs.Spec != nil {
		htcRS.CpusPerTask = specs.Spec.CpusPerTask
		htcRS.MemPerNodeMB = specs.Spec.MemPerNodeMB
	}

	// Export resource variables for use in build scripts
	writeEnvVars(shWriter, htcRS)
	fmt.Fprintln(shWriter, "")

	// Print job information at start
	writeJobHeader(shWriter, "$_CONDOR_CLUSTER_ID.$_CONDOR_PROC_ID", specs, formatHTCondorTime, jobSpec.Metadata)
	fmt.Fprintln(shWriter, "")

	// Write the command
	fmt.Fprintln(shWriter, jobSpec.Command)

	// Print completion info
	fmt.Fprintln(shWriter, "")
	writeJobFooter(shWriter, "$_CONDOR_CLUSTER_ID.$_CONDOR_PROC_ID")

	shWriter.Flush()
	shFile.Close()

	if err := os.Chmod(shPath, utils.PermExec); err != nil {
		return "", NewScriptCreationError(jobSpec.Name, shPath, err)
	}

	// === Create HTCondor submit file ===
	subPath := filepath.Join(outputDir, fmt.Sprintf("%s.sub", name))
	subFile, err := os.Create(subPath)
	if err != nil {
		return "", NewScriptCreationError(jobSpec.Name, subPath, err)
	}
	defer subFile.Close()

	subWriter := bufio.NewWriter(subFile)
	defer subWriter.Flush()

	// Write submit file header
	fmt.Fprintln(subWriter, "# HTCondor Submit File")
	fmt.Fprintln(subWriter, "universe = vanilla")
	fmt.Fprintf(subWriter, "executable = %s\n", shPath)
	fmt.Fprintln(subWriter, "transfer_executable = false")
	fmt.Fprintln(subWriter, "")

	// Write unrecognized flags (RemainingFlags only contains flags not parsed into typed fields)
	for _, flag := range specs.RemainingFlags {
		fmt.Fprintf(subWriter, "%s\n", flag)
	}

	// Write resource specifications
	if specs.Control.JobName != "" {
		fmt.Fprintf(subWriter, "+JobName = \"%s\"\n", specs.Control.JobName)
	}
	if specs.Control.Partition != "" {
		fmt.Fprintf(subWriter, "accounting_group = %s\n", specs.Control.Partition)
	}

	// Output/error/log
	if specs.Control.Stdout != "" {
		fmt.Fprintf(subWriter, "output = %s\n", specs.Control.Stdout)
	}
	if specs.Control.Stderr != "" {
		fmt.Fprintf(subWriter, "error = %s\n", specs.Control.Stderr)
	}
	condorLogPath := filepath.Join(outputDir, fmt.Sprintf("%s.condor.log", name))
	fmt.Fprintf(subWriter, "log = %s\n", condorLogPath)

	// Resource directives -- only when Spec is available
	if specs.Spec != nil {
		// Time limit
		if specs.Spec.Time > 0 {
			secs := int64(specs.Spec.Time.Seconds())
			fmt.Fprintf(subWriter, "+MaxRuntime = %d\n", secs)
		}
	}

	// Email notifications
	if specs.Control.EmailOnBegin || specs.Control.EmailOnEnd || specs.Control.EmailOnFail {
		if specs.Control.EmailOnBegin && specs.Control.EmailOnEnd && specs.Control.EmailOnFail {
			fmt.Fprintln(subWriter, "notification = Always")
		} else if specs.Control.EmailOnEnd {
			fmt.Fprintln(subWriter, "notification = Complete")
		} else if specs.Control.EmailOnFail {
			fmt.Fprintln(subWriter, "notification = Error")
		}
	} else {
		fmt.Fprintln(subWriter, "notification = Never")
	}
	if specs.Control.MailUser != "" {
		fmt.Fprintf(subWriter, "notify_user = %s\n", specs.Control.MailUser)
	}

	// Queue directive
	fmt.Fprintln(subWriter, "")
	fmt.Fprintln(subWriter, "queue")

	// Self-dispose: the wrapper script removes itself (unless in debug mode)
	if !config.Global.Debug {
		shAppend, err := os.OpenFile(shPath, os.O_APPEND|os.O_WRONLY, utils.PermExec)
		if err == nil {
			fmt.Fprintf(shAppend, "\n# Self-dispose\nrm -f %s %s\n", shPath, subPath)
			shAppend.Close()
		}
	}

	return subPath, nil
}

// Submit submits an HTCondor job with optional dependency chain
func (h *HTCondorScheduler) Submit(scriptPath string, dependencyJobIDs []string) (string, error) {
	// HTCondor does not support simple dependency flags like SLURM/PBS/LSF.
	// Job dependencies require DAGMan, which is out of scope here.
	if len(dependencyJobIDs) > 0 {
		utils.PrintWarning("HTCondor does not support simple job dependencies. " +
			"Use DAGMan for dependency workflows. Dependencies will be ignored.")
	}

	// Execute condor_submit
	cmd := exec.Command(h.condorSubmitBin, scriptPath)
	output, err := cmd.CombinedOutput()
	if err != nil {
		return "", NewSubmissionError("HTCondor", filepath.Base(scriptPath), string(output), err)
	}

	// Parse job ID (cluster ID) from output
	// Example: "1 job(s) submitted to cluster 12345."
	matches := h.jobIDRe.FindStringSubmatch(string(output))
	if len(matches) < 2 {
		return "", fmt.Errorf("%w: %s", ErrJobIDParseFailed, string(output))
	}

	jobID := matches[1]
	return jobID, nil
}

// GetClusterInfo retrieves HTCondor cluster configuration
func (h *HTCondorScheduler) GetClusterInfo() (*ClusterInfo, error) {
	info := &ClusterInfo{
		AvailableGpus: make([]GpuInfo, 0),
		Limits:        make([]ResourceLimits, 0),
	}

	if h.condorStatusBin == "" {
		return info, nil
	}

	// Get node info from condor_status
	maxCpus, maxMem, err := h.getMaxNodeResources()
	if err == nil {
		info.MaxCpusPerNode = maxCpus
		info.MaxMemMBPerNode = maxMem
	}

	// Get GPU info
	gpus, err := h.getGpuInfo()
	if err == nil {
		info.AvailableGpus = gpus
	}

	// Get cluster limits (HTCondor uses global/accounting group limits rather than queues)
	if h.condorConfigValBin != "" || h.condorStatusBin != "" {
		limits, err := h.getClusterLimits(info.AvailableGpus, maxCpus, maxMem)
		if err == nil {
			info.Limits = limits
		}
	}

	return info, nil
}

// getMaxNodeResources queries HTCondor for max CPU and memory per node
func (h *HTCondorScheduler) getMaxNodeResources() (int, int64, error) {
	// condor_status -af TotalSlotCpus TotalSlotMemory
	cmd := exec.Command(h.condorStatusBin, "-compact", "-af", "TotalSlotCpus", "TotalSlotMemory")
	output, err := cmd.CombinedOutput()
	if err != nil {
		return 0, 0, NewClusterError("HTCondor", "query node resources", err)
	}

	var maxCpus int
	var maxMemMB int64

	lines := strings.Split(strings.TrimSpace(string(output)), "\n")
	for _, line := range lines {
		fields := strings.Fields(line)
		if len(fields) < 2 {
			continue
		}

		var cpus int
		var memMB int64
		fmt.Sscanf(fields[0], "%d", &cpus)
		fmt.Sscanf(fields[1], "%d", &memMB)

		if cpus > maxCpus {
			maxCpus = cpus
		}
		if memMB > maxMemMB {
			maxMemMB = memMB
		}
	}

	return maxCpus, maxMemMB, nil
}

// getGpuInfo queries HTCondor for available GPU info
func (h *HTCondorScheduler) getGpuInfo() ([]GpuInfo, error) {
	cmd := exec.Command(h.condorStatusBin, "-compact", "-constraint", "TotalGpus > 0",
		"-af", "TotalGpus", "Machine")
	output, err := cmd.CombinedOutput()
	if err != nil {
		return nil, NewClusterError("HTCondor", "query GPUs", err)
	}

	gpuMap := make(map[string]*GpuInfo)
	lines := strings.Split(strings.TrimSpace(string(output)), "\n")

	for _, line := range lines {
		fields := strings.Fields(line)
		if len(fields) < 2 {
			continue
		}

		var gpuCount int
		fmt.Sscanf(fields[0], "%d", &gpuCount)
		machine := fields[1]

		if gpuCount <= 0 {
			continue
		}

		key := fmt.Sprintf("gpu:%s", machine)
		if existing, ok := gpuMap[key]; ok {
			existing.Total += gpuCount
			existing.Available += gpuCount
		} else {
			gpuMap[key] = &GpuInfo{
				Type:      "gpu",
				Total:     gpuCount,
				Available: gpuCount,
				Partition: machine,
			}
		}
	}

	gpus := make([]GpuInfo, 0, len(gpuMap))
	for _, info := range gpuMap {
		gpus = append(gpus, *info)
	}

	return gpus, nil
}

// getClusterLimits queries HTCondor for cluster-wide resource limits
// HTCondor doesn't have partitions/queues like other schedulers, but uses accounting groups
// and global limits. This function creates a "default" limit representing cluster resources.
func (h *HTCondorScheduler) getClusterLimits(gpuInfo []GpuInfo, maxCpus int, maxMemMB int64) ([]ResourceLimits, error) {
	limit := &ResourceLimits{
		Partition: "default",
	}

	// Use the max node resources as baseline
	if maxCpus > 0 {
		limit.MaxCpusPerNode = maxCpus
	}
	if maxMemMB > 0 {
		limit.MaxMemMBPerNode = maxMemMB
	}

	// Try to get max walltime from config if condor_config_val is available
	if h.condorConfigValBin != "" {
		if maxTime, err := h.getConfigValue("MAX_JOB_RUNTIME"); err == nil && maxTime != "" {
			// MAX_JOB_RUNTIME is in seconds
			if seconds, err := strconv.ParseInt(maxTime, 10, 64); err == nil && seconds > 0 {
				limit.MaxTime = time.Duration(seconds) * time.Second
			}
		}
	}

	// Calculate total GPU count
	totalGpus := 0
	for _, gpu := range gpuInfo {
		totalGpus += gpu.Total
	}
	if totalGpus > 0 {
		limit.MaxGpus = totalGpus
	}

	// If we have any limits set, return them
	if limit.MaxCpusPerNode > 0 || limit.MaxMemMBPerNode > 0 || limit.MaxTime > 0 || limit.MaxGpus > 0 {
		return []ResourceLimits{*limit}, nil
	}

	return []ResourceLimits{}, nil
}

// getConfigValue queries HTCondor configuration for a specific parameter
func (h *HTCondorScheduler) getConfigValue(param string) (string, error) {
	if h.condorConfigValBin == "" {
		return "", fmt.Errorf("condor_config_val not available")
	}

	cmd := exec.Command(h.condorConfigValBin, param)
	output, err := cmd.CombinedOutput()
	if err != nil {
		return "", err
	}

	return strings.TrimSpace(string(output)), nil
}

// GetJobResources reads allocated resources from HTCondor environment variables.
func (h *HTCondorScheduler) GetJobResources() *JobResources {
	if _, ok := os.LookupEnv("_CONDOR_JOB_AD"); !ok {
		return nil
	}
	res := &JobResources{}
	res.Ncpus = getEnvInt("_CONDOR_REQUEST_CPUS")
	// _CONDOR_REQUEST_MEMORY is in MB
	res.MemMB = getEnvInt64("_CONDOR_REQUEST_MEMORY")
	res.Ngpus = getCudaDeviceCount()
	return res
}

// parseHTCondorMemory parses HTCondor memory specifications.
// HTCondor default unit is MB if no suffix is provided.
func parseHTCondorMemory(memStr string) (int64, error) {
	memStr = strings.TrimSpace(memStr)
	if memStr == "" {
		return 0, fmt.Errorf("%w: empty memory string", ErrInvalidMemoryFormat)
	}

	// Check if it's a plain number (default MB in HTCondor)
	if val, err := strconv.ParseInt(memStr, 10, 64); err == nil {
		return val, nil
	}

	// Try with unit suffix
	return parseMemory(memStr)
}

// formatHTCondorTime formats a duration for display (HH:MM:SS)
func formatHTCondorTime(d time.Duration) string {
	if d <= 0 {
		return ""
	}
	total := int64(d.Seconds())
	hours := total / 3600
	minutes := (total % 3600) / 60
	seconds := total % 60
	return fmt.Sprintf("%02d:%02d:%02d", hours, minutes, seconds)
}

// TryParseHTCondorScript attempts to parse an HTCondor submit file without requiring HTCondor binaries.
// This is a static parser that can work in any environment.
func TryParseHTCondorScript(scriptPath string) (*ScriptSpecs, error) {
	parser := &HTCondorScheduler{}
	return parser.ReadScriptSpecs(scriptPath)
}
