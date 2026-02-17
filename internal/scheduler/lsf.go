package scheduler

import (
	"bufio"
	"fmt"
	"os"
	"os/exec"
	"path/filepath"
	"regexp"
	"sort"
	"strconv"
	"strings"
	"time"

	"github.com/Justype/condatainer/internal/config"
	"github.com/Justype/condatainer/internal/utils"
)

// LsfScheduler implements the Scheduler interface for IBM Spectrum LSF
// EXPERIMENTAL: LSF is not tested on real clusters and may have edge cases. Feedback welcome.
type LsfScheduler struct {
	bsubBin     string
	bjobsBin    string
	bhostsBin   string
	bqueuesBin  string
	directiveRe *regexp.Regexp
	jobIDRe     *regexp.Regexp
}

// NewLsfScheduler creates a new LSF scheduler instance using bsub from PATH
func NewLsfScheduler() (*LsfScheduler, error) {
	return newLsfSchedulerWithBinary("")
}

// NewLsfSchedulerWithBinary creates an LSF scheduler using an explicit bsub path
func NewLsfSchedulerWithBinary(bsubBin string) (*LsfScheduler, error) {
	return newLsfSchedulerWithBinary(bsubBin)
}

func newLsfSchedulerWithBinary(bsubBin string) (*LsfScheduler, error) {
	binPath := bsubBin
	if binPath == "" {
		var err error
		binPath, err = exec.LookPath("bsub")
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

	bjobsCmd, _ := exec.LookPath("bjobs")
	bhostsCmd, _ := exec.LookPath("bhosts")
	bqueuesCmd, _ := exec.LookPath("bqueues")

	return &LsfScheduler{
		bsubBin:     binPath,
		bjobsBin:    bjobsCmd,
		bhostsBin:   bhostsCmd,
		bqueuesBin:  bqueuesCmd,
		directiveRe: regexp.MustCompile(`^\s*#BSUB\s+(.+)$`),
		jobIDRe:     regexp.MustCompile(`Job <(\d+)> is submitted`),
	}, nil
}

// IsAvailable checks if LSF is available and we're not inside an LSF job
func (l *LsfScheduler) IsAvailable() bool {
	if l.bsubBin == "" {
		return false
	}

	// Check if we're already inside an LSF job
	_, inJob := os.LookupEnv("LSB_JOBID")
	if inJob {
		return false
	}

	return true
}

// GetInfo returns information about the LSF scheduler
func (l *LsfScheduler) GetInfo() *SchedulerInfo {
	_, inJob := os.LookupEnv("LSB_JOBID")
	available := l.IsAvailable()

	info := &SchedulerInfo{
		Type:      "LSF",
		Binary:    l.bsubBin,
		InJob:     inJob,
		Available: available,
	}

	// Try to get LSF version
	if l.bsubBin != "" {
		if version, err := l.getLsfVersion(); err == nil {
			info.Version = version
		}
	}

	return info
}

// getLsfVersion attempts to get the LSF version
func (l *LsfScheduler) getLsfVersion() (string, error) {
	cmd := exec.Command(l.bsubBin, "-V")
	output, err := cmd.CombinedOutput()
	if err != nil {
		return "", err
	}

	// Parse version from output like "IBM Spectrum LSF 10.1.0.0, ..."
	versionStr := strings.TrimSpace(string(output))
	return versionStr, nil
}

// ReadScriptSpecs parses #BSUB directives from a build script using the shared pipeline.
func (l *LsfScheduler) ReadScriptSpecs(scriptPath string) (*ScriptSpecs, error) {
	lines, err := readFileLines(scriptPath)
	if err != nil {
		return nil, err
	}
	return parseScript(scriptPath, lines, l.extractDirectives, l.parseRuntimeConfig, l.parseResourceSpec)
}

// extractDirectives extracts #BSUB directive strings from script lines.
func (l *LsfScheduler) extractDirectives(lines []string) []string {
	directives := make([]string, 0)
	for _, line := range lines {
		if matches := l.directiveRe.FindStringSubmatch(line); matches != nil {
			flag := utils.StripInlineComment(matches[1])
			directives = append(directives, flag)
		}
	}
	return directives
}

// parseRuntimeConfig consumes job control fields from the directive list.
// Returns unconsumed directives for parseResourceSpec.
func (l *LsfScheduler) parseRuntimeConfig(directives []string) (RuntimeConfig, []string, error) {
	rc := RuntimeConfig{}
	remaining := make([]string, 0)
	for _, flag := range directives {
		recognized := true
		switch {
		case flagMatches(flag, "-J"):
			rc.JobName, _ = flagValue(flag, "-J")
		case flagMatches(flag, "-o"):
			rc.Stdout, _ = flagValue(flag, "-o")
		case flagMatches(flag, "-e"):
			rc.Stderr, _ = flagValue(flag, "-e")
		case flag == "-B":
			rc.EmailOnBegin = true
		case flag == "-N":
			rc.EmailOnEnd = true
		case flagMatches(flag, "-u"):
			rc.MailUser, _ = flagValue(flag, "-u")
		case flagMatches(flag, "-q"):
			rc.Partition, _ = flagValue(flag, "-q")
		default:
			recognized = false
		}
		if !recognized {
			remaining = append(remaining, flag)
		}
	}
	return rc, remaining, nil
}

// parseResourceSpec consumes resource fields from the directive list.
// Returns nil + all directives on any parse failure (passthrough mode).
func (l *LsfScheduler) parseResourceSpec(directives []string) (*ResourceSpec, []string) {
	defaults := GetSpecDefaults()
	rs := &ResourceSpec{
		CpusPerTask:  defaults.CpusPerTask,
		TasksPerNode: defaults.TasksPerNode,
		Nodes:        defaults.Nodes,
		MemPerNodeMB: defaults.MemPerNodeMB,
		Time:         defaults.Time,
	}
	remaining := make([]string, 0)
	for _, flag := range directives {
		recognized := true
		var parseErr error

		switch {
		case flagMatches(flag, "-n"):
			_, parseErr = flagScanInt(flag, &rs.CpusPerTask, "-n")
		case flagMatches(flag, "-M"):
			_, parseErr = flagScan(flag, &rs.MemPerNodeMB, parseLsfMemory, "-M")
		case flagMatches(flag, "-W"):
			_, parseErr = flagScan(flag, &rs.Time, parseLsfTime, "-W")
		case flagMatches(flag, "-gpu"):
			gpuStr, _ := flagValue(flag, "-gpu")
			gpuStr = strings.Trim(gpuStr, "\"'")
			if gpu := parseLsfGpuDirective(gpuStr); gpu != nil {
				rs.Gpu = gpu
			}
		case flagMatches(flag, "-R"):
			resStr, _ := flagValue(flag, "-R")
			resStr = strings.Trim(resStr, "\"'")
			parseErr = l.parseLsfResourceIntoSpec(resStr, rs)
		case flag == "-x":
			rs.Exclusive = true
		default:
			recognized = false
		}

		if parseErr != nil {
			utils.PrintWarning("LSF: failed to parse directive %q: %v", flag, parseErr)
			return nil, directives
		}
		if !recognized {
			remaining = append(remaining, flag)
		}
	}
	return rs, remaining
}

// parseLsfResourceIntoSpec parses LSF -R resource requirement strings into a ResourceSpec.
// Supports: rusage[mem=N], rusage[ngpus_physical=N], span[hosts=N], etc.
// Returns an error if a known key has an invalid value.
func (l *LsfScheduler) parseLsfResourceIntoSpec(resStr string, rs *ResourceSpec) error {
	// Look for span[hosts=N] block
	spanIdx := strings.Index(resStr, "span[")
	if spanIdx >= 0 {
		start := spanIdx + len("span[")
		end := strings.Index(resStr[start:], "]")
		if end >= 0 {
			spanContent := resStr[start : start+end]
			pairs := strings.FieldsFunc(spanContent, func(r rune) bool {
				return r == ':' || r == ','
			})
			for _, pair := range pairs {
				kv := strings.SplitN(pair, "=", 2)
				if len(kv) == 2 {
					key := strings.TrimSpace(kv[0])
					value := strings.TrimSpace(kv[1])
					if key == "hosts" {
						n, err := strconv.Atoi(value)
						if err != nil {
							return fmt.Errorf("invalid span[hosts=] value %q: %w", value, err)
						}
						rs.Nodes = n
					}
				}
			}
		}
	}

	// Look for rusage[...] block
	rusageIdx := strings.Index(resStr, "rusage[")
	if rusageIdx < 0 {
		return nil
	}

	// Extract content between brackets
	start := rusageIdx + len("rusage[")
	end := strings.Index(resStr[start:], "]")
	if end < 0 {
		return nil
	}

	rusageContent := resStr[start : start+end]

	// Parse key=value pairs separated by colons or commas
	pairs := strings.FieldsFunc(rusageContent, func(r rune) bool {
		return r == ':' || r == ','
	})

	for _, pair := range pairs {
		kv := strings.SplitN(pair, "=", 2)
		if len(kv) != 2 {
			continue
		}

		key := strings.TrimSpace(kv[0])
		value := strings.TrimSpace(kv[1])

		switch key {
		case "mem":
			mem, err := parseLsfMemory(value)
			if err != nil {
				return fmt.Errorf("invalid rusage[mem=] value %q: %w", value, err)
			}
			rs.MemPerNodeMB = mem
		case "ngpus_physical", "ngpus":
			count, err := strconv.Atoi(value)
			if err != nil {
				return fmt.Errorf("invalid rusage[%s=] value %q: %w", key, value, err)
			}
			rs.Gpu = &GpuSpec{
				Type:  "gpu",
				Count: count,
				Raw:   fmt.Sprintf("%s=%s", key, value),
			}
		}
	}
	return nil
}

// parseLsfGpuDirective parses LSF -gpu directive content
// Format: num=N[:type=T[:...]]
func parseLsfGpuDirective(gpuStr string) *GpuSpec {
	if gpuStr == "" {
		return nil
	}

	spec := &GpuSpec{
		Raw:   gpuStr,
		Count: 1,
		Type:  "gpu",
	}

	// Parse key=value pairs separated by colons
	parts := strings.Split(gpuStr, ":")
	for _, part := range parts {
		kv := strings.SplitN(part, "=", 2)
		if len(kv) != 2 {
			continue
		}

		key := strings.TrimSpace(kv[0])
		value := strings.TrimSpace(kv[1])

		switch key {
		case "num":
			if count, err := strconv.Atoi(value); err == nil {
				spec.Count = count
			}
		case "type":
			spec.Type = value
		}
	}

	return spec
}

// CreateScriptWithSpec generates an LSF batch script
func (l *LsfScheduler) CreateScriptWithSpec(jobSpec *JobSpec, outputDir string) (string, error) {
	specs := jobSpec.Specs

	// Normalize node/task defaults when spec is available
	if specs.Spec != nil {
		if specs.Spec.Nodes <= 0 {
			specs.Spec.Nodes = 1
		}
		if specs.Spec.TasksPerNode <= 0 {
			specs.Spec.TasksPerNode = 1
		}
	}

	// Create output directory if specified
	if outputDir != "" {
		if err := os.MkdirAll(outputDir, utils.PermDir); err != nil {
			return "", NewScriptCreationError(jobSpec.Name, outputDir, err)
		}
	}

	// Set log path based on job name
	if jobSpec.Name != "" {
		safeName := strings.ReplaceAll(jobSpec.Name, "/", "--")
		specs.Control.Stdout = filepath.Join(outputDir, fmt.Sprintf("%s.log", safeName))
	}

	// Generate script filename
	scriptName := "job.lsf"
	if jobSpec.Name != "" {
		safeName := strings.ReplaceAll(jobSpec.Name, "/", "--")
		scriptName = fmt.Sprintf("%s.lsf", safeName)
	}

	scriptPath := filepath.Join(outputDir, scriptName)

	// Create the batch script
	file, err := os.Create(scriptPath)
	if err != nil {
		return "", NewScriptCreationError(jobSpec.Name, scriptPath, err)
	}
	defer file.Close()

	writer := bufio.NewWriter(file)
	defer writer.Flush()

	// Write shebang
	fmt.Fprintln(writer, "#!/bin/bash")

	// Write unrecognized flags (RemainingFlags only contains flags not parsed into typed fields)
	for _, flag := range specs.RemainingFlags {
		fmt.Fprintf(writer, "#BSUB %s\n", flag)
	}

	// Add job name if specified
	if specs.Control.JobName != "" {
		fmt.Fprintf(writer, "#BSUB -J %s\n", specs.Control.JobName)
	}

	// Add custom stdout if specified
	if specs.Control.Stdout != "" {
		fmt.Fprintf(writer, "#BSUB -o %s\n", specs.Control.Stdout)
	}

	// Add email notifications if specified
	if specs.Control.EmailOnBegin {
		fmt.Fprintln(writer, "#BSUB -B")
	}
	if specs.Control.EmailOnEnd {
		fmt.Fprintln(writer, "#BSUB -N")
	}
	if specs.Control.MailUser != "" {
		fmt.Fprintf(writer, "#BSUB -u %s\n", specs.Control.MailUser)
	}
	if specs.Control.Partition != "" {
		fmt.Fprintf(writer, "#BSUB -q %s\n", specs.Control.Partition)
	}

	// Resource directives — only when Spec is available
	if specs.Spec != nil {
		// Generate CPU count
		if specs.Spec.CpusPerTask > 0 {
			fmt.Fprintf(writer, "#BSUB -n %d\n", specs.Spec.CpusPerTask)
		}

		// Force node allocation via span
		fmt.Fprintf(writer, "#BSUB -R \"span[hosts=%d]\"\n", specs.Spec.Nodes)

		// Add walltime if specified
		if specs.Spec.Time > 0 {
			fmt.Fprintf(writer, "#BSUB -W %s\n", formatLsfTime(specs.Spec.Time))
		}
		if specs.Spec.Exclusive {
			fmt.Fprintln(writer, "#BSUB -x")
		}
	}

	// Add blank line
	fmt.Fprintln(writer, "")

	// Export resource variables for use in build scripts
	if specs.Spec != nil {
		fmt.Fprintf(writer, "export NNODES=%d\n", specs.Spec.Nodes)
		fmt.Fprintf(writer, "export NTASKS=%d\n", specs.Spec.TasksPerNode)
		if specs.Spec.CpusPerTask > 0 {
			fmt.Fprintf(writer, "export NCPUS=%d\n", specs.Spec.CpusPerTask)
		} else {
			fmt.Fprintln(writer, "export NCPUS=1")
		}
		if specs.Spec.MemPerNodeMB > 0 {
			fmt.Fprintf(writer, "export MEM=%d\n", specs.Spec.MemPerNodeMB)
			fmt.Fprintf(writer, "export MEM_MB=%d\n", specs.Spec.MemPerNodeMB)
			fmt.Fprintf(writer, "export MEM_GB=%d\n", specs.Spec.MemPerNodeMB/1024)
		}
	} else {
		fmt.Fprintln(writer, "export NNODES=1")
		fmt.Fprintln(writer, "export NTASKS=1")
		fmt.Fprintln(writer, "export NCPUS=1")
	}
	fmt.Fprintln(writer, "")

	// Print job information at start
	fmt.Fprintln(writer, "# Print job information")
	fmt.Fprintln(writer, "_START_TIME=$SECONDS")
	fmt.Fprintln(writer, "_format_time() { local s=$1; printf '%02d:%02d:%02d' $((s/3600)) $((s%3600/60)) $((s%60)); }")
	fmt.Fprintln(writer, "echo \"========================================\"")
	fmt.Fprintln(writer, "echo \"Job ID:    $LSB_JOBID\"")
	fmt.Fprintf(writer, "echo \"Job Name:  %s\"\n", specs.Control.JobName)
	if specs.Spec != nil {
		if specs.Spec.Nodes > 1 {
			fmt.Fprintf(writer, "echo \"Nodes:     %d\"\n", specs.Spec.Nodes)
		}
		if specs.Spec.TasksPerNode > 1 {
			fmt.Fprintf(writer, "echo \"Tasks:     %d\"\n", specs.Spec.TasksPerNode)
		}
		if specs.Spec.CpusPerTask > 0 {
			if specs.Spec.TasksPerNode > 1 || specs.Spec.Nodes > 1 {
				fmt.Fprintf(writer, "echo \"CPUs/Task: %d\"\n", specs.Spec.CpusPerTask)
			} else {
				fmt.Fprintf(writer, "echo \"CPUs:      %d\"\n", specs.Spec.CpusPerTask)
			}
		}
		if specs.Spec.MemPerNodeMB > 0 {
			fmt.Fprintf(writer, "echo \"Memory:    %d MB\"\n", specs.Spec.MemPerNodeMB)
		}
		if specs.Spec.Time > 0 {
			fmt.Fprintf(writer, "echo \"Time:      %s\"\n", formatLsfTime(specs.Spec.Time))
		}
	}
	fmt.Fprintln(writer, "echo \"PWD:       $(pwd)\"")
	// Print additional metadata if available
	if len(jobSpec.Metadata) > 0 {
		// Find max key length for alignment
		maxLen := 0
		for key := range jobSpec.Metadata {
			if len(key) > maxLen {
				maxLen = len(key)
			}
		}
		// Print each metadata field with proper alignment
		for key, value := range jobSpec.Metadata {
			if value != "" {
				padding := maxLen - len(key)
				fmt.Fprintf(writer, "echo \"%s:%s %s\"\n", key, strings.Repeat(" ", padding+3), value)
			}
		}
	}
	fmt.Fprintf(writer, "%s\n", "echo \"Started:   $(date '+%Y-%m-%d %T')\"")
	fmt.Fprintln(writer, "echo \"========================================\"")
	fmt.Fprintln(writer, "")

	// Write the command
	fmt.Fprintln(writer, jobSpec.Command)

	// Print completion info
	fmt.Fprintln(writer, "")
	fmt.Fprintln(writer, "echo \"========================================\"")
	fmt.Fprintln(writer, "echo \"Job ID:    $LSB_JOBID\"")
	fmt.Fprintln(writer, "echo \"Elapsed:   $(_format_time $(($SECONDS - $_START_TIME)))\"")
	fmt.Fprintf(writer, "%s\n", "echo \"Completed: $(date '+%Y-%m-%d %T')\"")
	fmt.Fprintln(writer, "echo \"========================================\"")

	// Self-delete the script after execution (unless in debug mode)
	if !config.Global.Debug {
		fmt.Fprintf(writer, "rm -f %s\n", scriptPath)
	}

	// Make executable
	if err := os.Chmod(scriptPath, utils.PermExec); err != nil {
		return "", NewScriptCreationError(jobSpec.Name, scriptPath, err)
	}

	return scriptPath, nil
}

// Submit submits an LSF job with optional dependency chain
func (l *LsfScheduler) Submit(scriptPath string, dependencyJobIDs []string) (string, error) {
	args := []string{}

	// Add dependency if provided
	// LSF uses: -w "done(id1) && done(id2)"
	if len(dependencyJobIDs) > 0 {
		var conditions []string
		for _, id := range dependencyJobIDs {
			conditions = append(conditions, fmt.Sprintf("done(%s)", id))
		}
		depStr := strings.Join(conditions, " && ")
		args = append(args, "-w", depStr)
	}

	// LSF reads the script from stdin: bsub < script.lsf
	args = append(args, "<", scriptPath)

	// Execute bsub with shell to handle stdin redirection
	shellCmd := fmt.Sprintf("%s %s", l.bsubBin, strings.Join(args, " "))
	cmd := exec.Command("bash", "-c", shellCmd)
	output, err := cmd.CombinedOutput()
	if err != nil {
		return "", NewSubmissionError("LSF", filepath.Base(scriptPath), string(output), err)
	}

	// Parse job ID from output like "Job <12345> is submitted to queue <normal>."
	matches := l.jobIDRe.FindStringSubmatch(string(output))
	if len(matches) < 2 {
		return "", fmt.Errorf("%w: %s", ErrJobIDParseFailed, string(output))
	}

	jobID := matches[1]
	return jobID, nil
}

// GetClusterInfo retrieves cluster configuration (GPUs, limits)
func (l *LsfScheduler) GetClusterInfo() (*ClusterInfo, error) {
	info := &ClusterInfo{
		AvailableGpus: make([]GpuInfo, 0),
		Limits:        make([]ResourceLimits, 0),
	}

	// Get host resources using bhosts
	if l.bhostsBin != "" {
		maxCpus, maxMem, err := l.getHostResources()
		if err == nil {
			info.MaxCpusPerNode = maxCpus
			info.MaxMemMBPerNode = maxMem
		}

		// Get GPU info
		gpus, err := l.getGpuInfo()
		if err == nil {
			info.AvailableGpus = gpus
		}
	}

	// Get queue limits
	if l.bqueuesBin != "" {
		limits, err := l.getQueueLimits(info.AvailableGpus)
		if err == nil {
			info.Limits = limits
		}
	}

	if info.MaxCpusPerNode == 0 && info.MaxMemMBPerNode == 0 && len(info.AvailableGpus) == 0 {
		return nil, ErrClusterInfoUnavailable
	}

	return info, nil
}

// getHostResources queries LSF for host resources using bhosts and lshosts
func (l *LsfScheduler) getHostResources() (int, int64, error) {
	// First get CPU info from bhosts
	cmd := exec.Command(l.bhostsBin, "-w")
	output, err := cmd.CombinedOutput()
	if err != nil {
		return 0, 0, NewClusterError("LSF", "query hosts", err)
	}

	var maxCpus int

	lines := strings.Split(string(output), "\n")
	for i, line := range lines {
		if i == 0 {
			continue // Skip header line
		}

		fields := strings.Fields(line)
		if len(fields) < 4 {
			continue
		}

		// bhosts -w output: HOST_NAME STATUS JL/U MAX NJOBS RUN SSUSP USUSP RSV
		if cpus, err := strconv.Atoi(fields[3]); err == nil && cpus > maxCpus {
			maxCpus = cpus
		}
	}

	// Get memory info from lshosts
	var maxMemMB int64
	lshostsCmd, err := exec.LookPath("lshosts")
	if err == nil {
		cmd := exec.Command(lshostsCmd, "-w")
		output, err := cmd.CombinedOutput()
		if err == nil {
			lines := strings.Split(string(output), "\n")
			for i, line := range lines {
				if i == 0 {
					continue // Skip header line
				}

				fields := strings.Fields(line)
				// lshosts output: HOST_NAME type model cpuf ncpus maxmem maxswp server RESOURCES
				// maxmem is typically in format like "32G" or "65536M"
				if len(fields) >= 6 {
					memStr := fields[5]
					if memMB, err := parseLsfMemoryToMB(memStr); err == nil && memMB > maxMemMB {
						maxMemMB = memMB
					}
				}
			}
		}
	}

	return maxCpus, maxMemMB, nil
}

// getQueueLimits queries LSF for queue resource limits using bqueues
func (l *LsfScheduler) getQueueLimits(gpuInfo []GpuInfo) ([]ResourceLimits, error) {
	cmd := exec.Command(l.bqueuesBin, "-l")
	output, err := cmd.CombinedOutput()
	if err != nil {
		return nil, NewClusterError("LSF", "query queues", err)
	}

	queueLimits := make(map[string]*ResourceLimits)
	var currentQueue string

	lines := strings.Split(string(output), "\n")
	for _, line := range lines {
		line = strings.TrimSpace(line)

		// Queue name line starts with "QUEUE: "
		if strings.HasPrefix(line, "QUEUE: ") {
			currentQueue = strings.TrimSpace(strings.TrimPrefix(line, "QUEUE:"))
			if currentQueue != "" {
				queueLimits[currentQueue] = &ResourceLimits{
					Partition: currentQueue,
				}
			}
			continue
		}

		// Skip if no current queue
		if currentQueue == "" {
			continue
		}

		limit := queueLimits[currentQueue]

		// Parse RUNLIMIT (walltime)
		if strings.Contains(line, "RUNLIMIT") {
			// Format: "RUNLIMIT  600.0 min"
			fields := strings.Fields(line)
			for i, field := range fields {
				if field == "RUNLIMIT" && i+1 < len(fields) {
					// Parse value and unit
					var minutes float64
					if _, err := fmt.Sscanf(fields[i+1], "%f", &minutes); err == nil {
						limit.MaxTime = time.Duration(minutes * 60 * float64(time.Second))
					}
					break
				}
			}
		}

		// Parse MEMLIMIT (per-process memory limit)
		if strings.Contains(line, "MEMLIMIT") {
			fields := strings.Fields(line)
			for i, field := range fields {
				if field == "MEMLIMIT" && i+1 < len(fields) {
					// Parse memory value (in MB)
					var memMB int64
					if _, err := fmt.Sscanf(fields[i+1], "%d", &memMB); err == nil && memMB > 0 {
						limit.MaxMemMBPerNode = memMB
					}
					break
				}
			}
		}

		// Parse CPULIMIT or PROCLIMIT
		if strings.Contains(line, "PROCLIMIT") {
			fields := strings.Fields(line)
			for i, field := range fields {
				if field == "PROCLIMIT" && i+1 < len(fields) {
					var procs int
					if _, err := fmt.Sscanf(fields[i+1], "%d", &procs); err == nil && procs > 0 {
						limit.MaxCpusPerNode = procs
					}
					break
				}
			}
		}
	}

	// Get available resources per queue from actual hosts
	if l.bhostsBin != "" {
		availRes, err := l.getAvailableResourcesByQueue()
		if err == nil {
			// Merge available resources into queue limits
			for queueName, limit := range queueLimits {
				if avail, ok := availRes[queueName]; ok {
					// Use available resources as limits if they're set and larger than config
					if avail.MaxCpusPerNode > 0 {
						limit.MaxCpusPerNode = avail.MaxCpusPerNode
					}
					if avail.MaxMemMBPerNode > 0 {
						limit.MaxMemMBPerNode = avail.MaxMemMBPerNode
					}
				}
			}
		}
	}

	// Convert map to slice
	limits := make([]ResourceLimits, 0, len(queueLimits))

	// Calculate GPU count per queue from GPU info
	gpusByQueue := make(map[string]int)
	for _, gpu := range gpuInfo {
		queue := gpu.Partition
		if queue == "" {
			// If no queue assigned, add to all queues
			for qName := range queueLimits {
				gpusByQueue[qName] += gpu.Total
			}
		} else {
			gpusByQueue[queue] += gpu.Total
		}
	}

	for queueName, limit := range queueLimits {
		// Assign GPU count for this queue
		if maxGpus, ok := gpusByQueue[queueName]; ok {
			limit.MaxGpus = maxGpus
		}
		limits = append(limits, *limit)
	}

	// Sort by queue name for consistent output
	sort.Slice(limits, func(i, j int) bool {
		return limits[i].Partition < limits[j].Partition
	})

	return limits, nil
}

// getAvailableResourcesByQueue queries available CPUs and memory per queue
func (l *LsfScheduler) getAvailableResourcesByQueue() (map[string]ResourceLimits, error) {
	// LSF doesn't directly expose queue-to-host mapping in simple commands
	// We'll use bhosts to get host info and try to match to queues via bqueues
	// This is a best-effort implementation

	cmd := exec.Command(l.bhostsBin, "-w")
	output, err := cmd.CombinedOutput()
	if err != nil {
		return nil, NewClusterError("LSF", "query available resources by queue", err)
	}

	// Track resources per queue
	resources := make(map[string]ResourceLimits)

	// Parse bhosts output to get host capabilities
	type hostInfo struct {
		cpus   int
		memMB  int64
		status string
	}
	hosts := make(map[string]hostInfo)

	lines := strings.Split(string(output), "\n")
	for i, line := range lines {
		if i == 0 {
			continue // Skip header line
		}

		fields := strings.Fields(line)
		if len(fields) < 4 {
			continue
		}

		// bhosts -w output: HOST_NAME STATUS JL/U MAX NJOBS RUN SSUSP USUSP RSV
		hostName := fields[0]
		status := fields[1]
		var cpus int
		fmt.Sscanf(fields[3], "%d", &cpus)

		hosts[hostName] = hostInfo{
			cpus:   cpus,
			status: status,
		}
	}

	// Get memory info from lshosts
	lshostsCmd, err := exec.LookPath("lshosts")
	if err == nil {
		cmd := exec.Command(lshostsCmd, "-w")
		output, err := cmd.CombinedOutput()
		if err == nil {
			lines := strings.Split(string(output), "\n")
			for i, line := range lines {
				if i == 0 {
					continue // Skip header
				}

				fields := strings.Fields(line)
				if len(fields) >= 6 {
					hostName := fields[0]
					if host, ok := hosts[hostName]; ok {
						memStr := fields[5]
						if memMB, err := parseLsfMemoryToMB(memStr); err == nil {
							host.memMB = memMB
							hosts[hostName] = host
						}
					}
				}
			}
		}
	}

	// Query bqueues to get queue-to-host mapping
	// Note: This is simplified - LSF queue-host mapping can be complex
	// For a more accurate implementation, we'd need to parse bqueues -l output
	// or use LSF API if available
	cmd = exec.Command(l.bqueuesBin)
	output, err = cmd.CombinedOutput()
	if err == nil {
		// For now, aggregate all host resources under each queue
		// This is a simplification - in reality, queues may have specific host groups
		lines := strings.Split(string(output), "\n")
		for i, line := range lines {
			if i == 0 {
				continue // Skip header
			}

			fields := strings.Fields(line)
			if len(fields) == 0 {
				continue
			}

			queueName := fields[0]

			// Aggregate all ok/closed hosts for this queue
			// In a real implementation, we'd filter by queue's host groups
			var maxCpus int
			var maxMemMB int64

			for _, host := range hosts {
				if host.status == "ok" || host.status == "closed" {
					if host.cpus > maxCpus {
						maxCpus = host.cpus
					}
					if host.memMB > maxMemMB {
						maxMemMB = host.memMB
					}
				}
			}

			if maxCpus > 0 || maxMemMB > 0 {
				resources[queueName] = ResourceLimits{
					Partition: queueName,
					MaxCpusPerNode:  maxCpus,
					MaxMemMBPerNode: maxMemMB,
				}
			}
		}
	}

	return resources, nil
}

// getGpuInfo queries LSF for GPU information using bhosts
func (l *LsfScheduler) getGpuInfo() ([]GpuInfo, error) {
	// Try bhosts -gpu first (LSF 10.1+)
	cmd := exec.Command(l.bhostsBin, "-gpu")
	output, err := cmd.CombinedOutput()
	if err != nil {
		// If -gpu flag not supported, try parsing bhosts -l output
		return l.getGpuInfoFromHostDetails()
	}

	gpuMap := make(map[string]*GpuInfo)
	lines := strings.Split(string(output), "\n")

	for i, line := range lines {
		if i == 0 {
			continue // Skip header line
		}

		fields := strings.Fields(line)
		if len(fields) < 2 {
			continue
		}

		// bhosts -gpu output format (varies by LSF version):
		// HOST_NAME ngpus gputype ...
		// or HOST_NAME gpu_shared_avg_mut gpu_shared_avg_ut ngpus_physical ...

		var gpuCount int
		var gpuType string

		// Try to find ngpus field (usually contains number of GPUs)
		for _, field := range fields[1:] {
			if count, err := strconv.Atoi(field); err == nil && count > 0 {
				gpuCount = count
				break
			}
		}

		if gpuCount == 0 {
			continue
		}

		// Look for GPU type in the output (often contains "nvidia", "tesla", etc.)
		gpuType = "gpu" // default
		for _, field := range fields {
			lower := strings.ToLower(field)
			if strings.Contains(lower, "nvidia") || strings.Contains(lower, "tesla") ||
				strings.Contains(lower, "a100") || strings.Contains(lower, "v100") ||
				strings.Contains(lower, "h100") || strings.Contains(lower, "gpu") {
				gpuType = field
				break
			}
		}

		// Aggregate by GPU type
		if existing, ok := gpuMap[gpuType]; ok {
			existing.Total += gpuCount
			existing.Available += gpuCount // LSF doesn't easily expose available vs used
		} else {
			gpuMap[gpuType] = &GpuInfo{
				Type:      gpuType,
				Total:     gpuCount,
				Available: gpuCount,
			}
		}
	}

	gpus := make([]GpuInfo, 0, len(gpuMap))
	for _, info := range gpuMap {
		gpus = append(gpus, *info)
	}

	return gpus, nil
}

// getGpuInfoFromHostDetails attempts to extract GPU info from detailed host info
func (l *LsfScheduler) getGpuInfoFromHostDetails() ([]GpuInfo, error) {
	// This is a fallback for older LSF versions
	// Query bhosts for list of hosts, then bhosts -l for each
	cmd := exec.Command(l.bhostsBin)
	output, err := cmd.CombinedOutput()
	if err != nil {
		return nil, NewClusterError("LSF", "query GPU info", err)
	}

	gpuMap := make(map[string]*GpuInfo)
	lines := strings.Split(string(output), "\n")

	// Parse bhosts output to find hosts, then query detailed info
	for i, line := range lines {
		if i == 0 {
			continue
		}
		fields := strings.Fields(line)
		if len(fields) == 0 {
			continue
		}
		hostName := fields[0]

		// Query detailed host info
		detailCmd := exec.Command(l.bhostsBin, "-l", hostName)
		detailOutput, err := detailCmd.CombinedOutput()
		if err != nil {
			continue
		}

		// Parse for GPU-related resources
		// LSF stores GPU info in resources like "ngpus_physical", "gpu", etc.
		detailStr := string(detailOutput)
		if strings.Contains(detailStr, "ngpus") || strings.Contains(detailStr, "gpu") {
			// Found GPU info - parse it
			// This is simplified - real parsing would be more sophisticated
			gpuType := "gpu"
			gpuCount := 1 // Conservative estimate

			if existing, ok := gpuMap[gpuType]; ok {
				existing.Total += gpuCount
				existing.Available += gpuCount
			} else {
				gpuMap[gpuType] = &GpuInfo{
					Type:      gpuType,
					Total:     gpuCount,
					Available: gpuCount,
				}
			}
		}
	}

	gpus := make([]GpuInfo, 0, len(gpuMap))
	for _, info := range gpuMap {
		gpus = append(gpus, *info)
	}

	return gpus, nil
}

// parseLsfMemoryToMB parses LSF memory format from lshosts (e.g., "32G", "65536M") to MB
func parseLsfMemoryToMB(memStr string) (int64, error) {
	memStr = strings.TrimSpace(memStr)
	if memStr == "" || memStr == "-" {
		return 0, fmt.Errorf("no memory info")
	}

	// Check for unit suffix
	memUpper := strings.ToUpper(memStr)
	var value int64
	var unit string

	// Try to parse with unit
	n, err := fmt.Sscanf(memUpper, "%d%s", &value, &unit)
	if err != nil || n < 1 {
		return 0, fmt.Errorf("%w: %s", ErrInvalidMemoryFormat, memStr)
	}

	switch unit {
	case "T", "TB":
		return value * 1024 * 1024, nil
	case "G", "GB":
		return value * 1024, nil
	case "M", "MB", "":
		return value, nil
	case "K", "KB":
		return value / 1024, nil
	default:
		return value, nil
	}
}

// parseLsfTime parses LSF walltime format
// Supports: HH:MM, MM, HH:MM:SS, or just minutes as integer
func parseLsfTime(timeStr string) (time.Duration, error) {
	timeStr = strings.TrimSpace(timeStr)
	if timeStr == "" {
		return 0, nil
	}

	parts := strings.Split(timeStr, ":")
	var hours, minutes, seconds int64

	switch len(parts) {
	case 3:
		hours, _ = strconv.ParseInt(parts[0], 10, 64)
		minutes, _ = strconv.ParseInt(parts[1], 10, 64)
		seconds, _ = strconv.ParseInt(parts[2], 10, 64)
	case 2:
		hours, _ = strconv.ParseInt(parts[0], 10, 64)
		minutes, _ = strconv.ParseInt(parts[1], 10, 64)
	case 1:
		minutes, _ = strconv.ParseInt(parts[0], 10, 64)
	default:
		return 0, fmt.Errorf("%w: %s", ErrInvalidTimeFormat, timeStr)
	}

	totalSeconds := hours*3600 + minutes*60 + seconds
	return time.Duration(totalSeconds) * time.Second, nil
}

// formatLsfTime formats a duration as LSF walltime (HH:MM)
func formatLsfTime(d time.Duration) string {
	if d <= 0 {
		return ""
	}
	total := int64(d.Seconds())
	hours := total / 3600
	minutes := (total % 3600) / 60
	return fmt.Sprintf("%02d:%02d", hours, minutes)
}

// parseLsfMemory parses LSF memory specification and returns MB
// LSF -M flag uses KB by default, but can have suffixes
func parseLsfMemory(memStr string) (int64, error) {
	memStr = strings.TrimSpace(memStr)
	if memStr == "" {
		return 0, fmt.Errorf("empty memory string")
	}

	// Try to parse as plain number (KB by default in LSF)
	if value, err := strconv.ParseInt(memStr, 10, 64); err == nil {
		return value / 1024, nil // Convert KB to MB
	}

	// Try with unit suffix
	memUpper := strings.ToUpper(memStr)
	var value int64
	var unit string

	n, err := fmt.Sscanf(memUpper, "%d%s", &value, &unit)
	if err != nil && n == 0 {
		return 0, fmt.Errorf("%w: %s", ErrInvalidMemoryFormat, memStr)
	}

	switch unit {
	case "G", "GB":
		return value * 1024, nil
	case "M", "MB":
		return value, nil
	case "K", "KB", "":
		return value / 1024, nil
	case "T", "TB":
		return value * 1024 * 1024, nil
	default:
		return value / 1024, nil // Default to KB
	}
}

// GetJobResources reads allocated resources from LSF environment variables.
func (s *LsfScheduler) GetJobResources() *JobResources {
	if _, ok := os.LookupEnv("LSB_JOBID"); !ok {
		return nil
	}
	res := &JobResources{}
	res.Ncpus = getEnvInt("LSB_DJOB_NUMPROC")
	if res.Ncpus == nil {
		res.Ncpus = getEnvInt("LSB_MAX_NUM_PROCESSORS")
	}
	// LSB_MAX_MEM_RUSAGE is in KB
	if memKB := getEnvInt64("LSB_MAX_MEM_RUSAGE"); memKB != nil {
		mb := *memKB / 1024
		if mb > 0 {
			res.MemMB = &mb
		}
	}
	res.Ngpus = getCudaDeviceCount()
	return res
}

// TryParseLsfScript attempts to parse an LSF script without requiring LSF binaries.
// This is a static parser that can work in any environment.
func TryParseLsfScript(scriptPath string) (*ScriptSpecs, error) {
	parser := &LsfScheduler{
		directiveRe: regexp.MustCompile(`^\s*#BSUB\s+(.+)$`),
	}
	return parser.ReadScriptSpecs(scriptPath)
}
