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

	"github.com/Justype/condatainer/internal/utils"
)

// HTCondorScheduler implements the Scheduler interface for HTCondor
// EXPERIMENTAL: HTCondor is not tested on real clusters and may have edge cases. Feedback welcome.
type HTCondorScheduler struct {
	condorSubmitBin string
	condorQBin      string
	condorStatusBin string
	directiveRe     *regexp.Regexp
	jobIDRe         *regexp.Regexp
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

	return &HTCondorScheduler{
		condorSubmitBin: binPath,
		condorQBin:      condorQBin,
		condorStatusBin: condorStatusBin,
		directiveRe:     regexp.MustCompile(`^\s*#CONDOR\s+(.+)$`),
		jobIDRe:         regexp.MustCompile(`submitted to cluster (\d+)`),
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

// ReadScriptSpecs parses #CONDOR directives from a build script
func (h *HTCondorScheduler) ReadScriptSpecs(scriptPath string) (*ScriptSpecs, error) {
	file, err := os.Open(scriptPath)
	if err != nil {
		if os.IsNotExist(err) {
			return nil, fmt.Errorf("%w: %s", ErrScriptNotFound, scriptPath)
		}
		return nil, err
	}
	defer file.Close()

	defaults := GetSpecDefaults()
	specs := &ScriptSpecs{
		Ncpus:    defaults.Ncpus,
		Ntasks:   defaults.Ntasks,
		Nodes:    defaults.Nodes,
		MemMB:    defaults.MemMB,
		Time:     defaults.Time,
		RawFlags: make([]string, 0),
	}

	scanner := bufio.NewScanner(file)
	lineNum := 0
	for scanner.Scan() {
		lineNum++
		line := scanner.Text()

		// Parse #CONDOR directives only
		if matches := h.directiveRe.FindStringSubmatch(line); matches != nil {
			flag := strings.TrimSpace(matches[1])
			specs.RawFlags = append(specs.RawFlags, flag)

			// Parse individual HTCondor directives
			if err := h.parseCondorDirective(flag, specs); err != nil {
				return nil, NewParseError("HTCondor", lineNum, line, err.Error())
			}
		}
	}

	if err := scanner.Err(); err != nil {
		return nil, fmt.Errorf("error reading script: %w", err)
	}

	return specs, nil
}

// parseCondorDirective parses individual HTCondor directives and updates specs
func (h *HTCondorScheduler) parseCondorDirective(flag string, specs *ScriptSpecs) error {
	// Split on first '=' to get key and value
	parts := strings.SplitN(flag, "=", 2)
	if len(parts) != 2 {
		// Not a key=value directive, store as raw flag
		return nil
	}

	key := strings.TrimSpace(strings.ToLower(parts[0]))
	value := strings.TrimSpace(parts[1])

	switch key {
	case "request_cpus":
		n, err := strconv.Atoi(value)
		if err != nil {
			return fmt.Errorf("invalid request_cpus value: %w", err)
		}
		specs.Ncpus = n

	case "request_memory":
		mem, err := parseHTCondorMemory(value)
		if err != nil {
			return err
		}
		specs.MemMB = mem

	case "request_gpus":
		n, err := strconv.Atoi(value)
		if err != nil {
			return fmt.Errorf("invalid request_gpus value: %w", err)
		}
		specs.Gpu = &GpuSpec{
			Type:  "gpu",
			Count: n,
			Raw:   value,
		}

	case "+maxruntime":
		secs, err := strconv.ParseInt(value, 10, 64)
		if err != nil {
			return fmt.Errorf("invalid +MaxRuntime value: %w", err)
		}
		specs.Time = time.Duration(secs) * time.Second

	case "output":
		specs.Stdout = value

	case "error":
		specs.Stderr = value

	case "notify_user":
		specs.MailUser = value

	case "notification":
		switch strings.ToLower(value) {
		case "always":
			specs.EmailOnBegin = true
			specs.EmailOnEnd = true
			specs.EmailOnFail = true
		case "complete":
			specs.EmailOnEnd = true
		case "error":
			specs.EmailOnFail = true
		case "never":
			specs.EmailOnBegin = false
			specs.EmailOnEnd = false
			specs.EmailOnFail = false
		}
	}

	return nil
}

// CreateScriptWithSpec generates an HTCondor submit description file and wrapper script
func (h *HTCondorScheduler) CreateScriptWithSpec(jobSpec *JobSpec, outputDir string) (string, error) {
	specs := jobSpec.Specs

	// Normalize node/task defaults (HTCondor is inherently single-node)
	if specs.Nodes <= 0 {
		specs.Nodes = 1
	}
	if specs.Ntasks <= 0 {
		specs.Ntasks = 1
	}

	// Create output directory if specified
	if outputDir != "" {
		if err := os.MkdirAll(outputDir, utils.PermDir); err != nil {
			return "", NewScriptCreationError(jobSpec.Name, outputDir, err)
		}
	}

	// Generate safe name for files
	safeName := "job"
	if jobSpec.Name != "" {
		safeName = strings.ReplaceAll(jobSpec.Name, "/", "--")
	}

	// Set log path based on job name
	if jobSpec.Name != "" {
		specs.Stdout = filepath.Join(outputDir, fmt.Sprintf("%s.log", safeName))
	}

	// === Create wrapper bash script ===
	shPath := filepath.Join(outputDir, fmt.Sprintf("%s.sh", safeName))
	shFile, err := os.Create(shPath)
	if err != nil {
		return "", NewScriptCreationError(jobSpec.Name, shPath, err)
	}

	shWriter := bufio.NewWriter(shFile)
	fmt.Fprintln(shWriter, "#!/bin/bash")
	fmt.Fprintln(shWriter, "")

	// Print job information at start
	fmt.Fprintln(shWriter, "# Print job information")
	fmt.Fprintln(shWriter, "_START_TIME=$SECONDS")
	fmt.Fprintln(shWriter, "_format_time() { local s=$1; printf '%02d:%02d:%02d' $((s/3600)) $((s%3600/60)) $((s%60)); }")
	fmt.Fprintln(shWriter, "echo \"========================================\"")
	fmt.Fprintln(shWriter, "echo \"Job ID:    $_CONDOR_CLUSTER_ID.$_CONDOR_PROC_ID\"")
	fmt.Fprintf(shWriter, "echo \"Job Name:  %s\"\n", specs.JobName)
	if specs.Ncpus > 0 {
		fmt.Fprintf(shWriter, "echo \"CPUs:      %d\"\n", specs.Ncpus)
	}
	if specs.MemMB > 0 {
		fmt.Fprintf(shWriter, "echo \"Memory:    %d MB\"\n", specs.MemMB)
	}
	if specs.Time > 0 {
		fmt.Fprintf(shWriter, "echo \"Time:      %s\"\n", formatHTCondorTime(specs.Time))
	}
	fmt.Fprintln(shWriter, "echo \"PWD:       $(pwd)\"")
	// Print additional metadata if available
	if len(jobSpec.Metadata) > 0 {
		maxLen := 0
		for key := range jobSpec.Metadata {
			if len(key) > maxLen {
				maxLen = len(key)
			}
		}
		for key, value := range jobSpec.Metadata {
			if value != "" {
				padding := maxLen - len(key)
				fmt.Fprintf(shWriter, "echo \"%s:%s %s\"\n", key, strings.Repeat(" ", padding+3), value)
			}
		}
	}
	fmt.Fprintf(shWriter, "%s\n", "echo \"Started:   $(date '+%Y-%m-%d %T')\"")
	fmt.Fprintln(shWriter, "echo \"========================================\"")
	fmt.Fprintln(shWriter, "")

	// Write the command
	fmt.Fprintln(shWriter, jobSpec.Command)

	// Print completion info
	fmt.Fprintln(shWriter, "")
	fmt.Fprintln(shWriter, "echo \"========================================\"")
	fmt.Fprintln(shWriter, "echo \"Job ID:    $_CONDOR_CLUSTER_ID.$_CONDOR_PROC_ID\"")
	fmt.Fprintln(shWriter, "echo \"Elapsed:   $(_format_time $(($SECONDS - $_START_TIME)))\"")
	fmt.Fprintf(shWriter, "%s\n", "echo \"Completed: $(date '+%Y-%m-%d %T')\"")
	fmt.Fprintln(shWriter, "echo \"========================================\"")

	shWriter.Flush()
	shFile.Close()

	if err := os.Chmod(shPath, utils.PermExec); err != nil {
		return "", NewScriptCreationError(jobSpec.Name, shPath, err)
	}

	// === Create HTCondor submit file ===
	subPath := filepath.Join(outputDir, fmt.Sprintf("%s.sub", safeName))
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

	// Write raw flags (excluding ones we'll regenerate)
	for _, flag := range specs.RawFlags {
		parts := strings.SplitN(flag, "=", 2)
		if len(parts) != 2 {
			fmt.Fprintf(subWriter, "# %s\n", flag)
			continue
		}
		key := strings.TrimSpace(strings.ToLower(parts[0]))

		// Skip flags that will be regenerated
		switch key {
		case "output", "error", "notification", "notify_user", "+maxruntime":
			continue
		}
		// Skip job-name flags if we have a name
		if specs.JobName != "" && key == "job_name" {
			continue
		}

		fmt.Fprintf(subWriter, "%s\n", flag)
	}

	// Write resource specifications
	if specs.JobName != "" {
		fmt.Fprintf(subWriter, "+JobName = \"%s\"\n", specs.JobName)
	}

	// Output/error/log
	if specs.Stdout != "" {
		fmt.Fprintf(subWriter, "output = %s\n", specs.Stdout)
	}
	errPath := filepath.Join(outputDir, fmt.Sprintf("%s.err", safeName))
	fmt.Fprintf(subWriter, "error = %s\n", errPath)
	condorLogPath := filepath.Join(outputDir, fmt.Sprintf("%s.condor.log", safeName))
	fmt.Fprintf(subWriter, "log = %s\n", condorLogPath)

	// Time limit
	if specs.Time > 0 {
		secs := int64(specs.Time.Seconds())
		fmt.Fprintf(subWriter, "+MaxRuntime = %d\n", secs)
	}

	// Email notifications
	if specs.EmailOnBegin || specs.EmailOnEnd || specs.EmailOnFail {
		if specs.EmailOnBegin && specs.EmailOnEnd && specs.EmailOnFail {
			fmt.Fprintln(subWriter, "notification = Always")
		} else if specs.EmailOnEnd {
			fmt.Fprintln(subWriter, "notification = Complete")
		} else if specs.EmailOnFail {
			fmt.Fprintln(subWriter, "notification = Error")
		}
	} else {
		fmt.Fprintln(subWriter, "notification = Never")
	}
	if specs.MailUser != "" {
		fmt.Fprintf(subWriter, "notify_user = %s\n", specs.MailUser)
	}

	// Queue directive
	fmt.Fprintln(subWriter, "")
	fmt.Fprintln(subWriter, "queue")

	// Self-dispose: the wrapper script removes itself
	// Re-open the sh file to append self-cleanup
	shAppend, err := os.OpenFile(shPath, os.O_APPEND|os.O_WRONLY, utils.PermExec)
	if err == nil {
		fmt.Fprintf(shAppend, "\n# Self-dispose\nrm -f %s %s\n", shPath, subPath)
		shAppend.Close()
	}

	return subPath, nil
}

// Submit submits an HTCondor job with optional dependency chain
func (h *HTCondorScheduler) Submit(scriptPath string, dependencyJobIDs []string) (string, error) {
	// HTCondor does not support simple dependency flags like SLURM/PBS/LSF.
	// Job dependencies require DAGMan, which is out of scope here.
	if len(dependencyJobIDs) > 0 {
		fmt.Printf("[CNT WARN] HTCondor does not support simple job dependencies. " +
			"Use DAGMan for dependency workflows. Dependencies will be ignored.\n")
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

// TryParseHTCondorScript attempts to parse an HTCondor script without requiring HTCondor binaries.
// This is a static parser that can work in any environment.
func TryParseHTCondorScript(scriptPath string) (*ScriptSpecs, error) {
	parser := &HTCondorScheduler{
		directiveRe: regexp.MustCompile(`^\s*#CONDOR\s+(.+)$`),
	}
	return parser.ReadScriptSpecs(scriptPath)
}
