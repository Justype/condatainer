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

// LsfScheduler implements the Scheduler interface for IBM Spectrum LSF
// EXPERIMENTAL: LSF is not tested on real clusters and may have edge cases. Feedback welcome.
type LsfScheduler struct {
	bsubBin     string
	bjobsBin    string
	bhostsBin   string
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

	return &LsfScheduler{
		bsubBin:     binPath,
		bjobsBin:    bjobsCmd,
		bhostsBin:   bhostsCmd,
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

// ReadScriptSpecs parses #BSUB directives from a build script
func (l *LsfScheduler) ReadScriptSpecs(scriptPath string) (*ScriptSpecs, error) {
	file, err := os.Open(scriptPath)
	if err != nil {
		if os.IsNotExist(err) {
			return nil, fmt.Errorf("%w: %s", ErrScriptNotFound, scriptPath)
		}
		return nil, err
	}
	defer file.Close()

	specs := &ScriptSpecs{
		Ncpus:    4, // default
		RawFlags: make([]string, 0),
	}

	scanner := bufio.NewScanner(file)
	lineNum := 0
	for scanner.Scan() {
		lineNum++
		line := scanner.Text()

		// Parse #BSUB directives only
		if matches := l.directiveRe.FindStringSubmatch(line); matches != nil {
			flag := strings.TrimSpace(matches[1])
			specs.RawFlags = append(specs.RawFlags, flag)

			// Parse LSF options
			if err := l.parseBsubFlag(flag, specs); err != nil {
				return nil, NewParseError("LSF", lineNum, line, err.Error())
			}
		}
	}

	if err := scanner.Err(); err != nil {
		return nil, fmt.Errorf("error reading script: %w", err)
	}

	return specs, nil
}

// parseBsubFlag parses individual BSUB flags and updates specs
func (l *LsfScheduler) parseBsubFlag(flag string, specs *ScriptSpecs) error {
	// Job name: -J jobname
	if strings.HasPrefix(flag, "-J ") {
		specs.JobName = strings.TrimSpace(strings.TrimPrefix(flag, "-J"))
		return nil
	}

	// Number of cores: -n cores
	if strings.HasPrefix(flag, "-n ") {
		ncpuStr := strings.TrimSpace(strings.TrimPrefix(flag, "-n"))
		if ncpus, err := strconv.Atoi(ncpuStr); err == nil {
			specs.Ncpus = ncpus
		}
		return nil
	}

	// Memory limit: -M memKB (LSF default is KB per process)
	if strings.HasPrefix(flag, "-M ") {
		memStr := strings.TrimSpace(strings.TrimPrefix(flag, "-M"))
		if mem, err := parseLsfMemory(memStr); err == nil {
			specs.MemMB = mem
		}
		return nil
	}

	// Walltime: -W [HH:]MM or HH:MM:SS
	if strings.HasPrefix(flag, "-W ") {
		timeStr := strings.TrimSpace(strings.TrimPrefix(flag, "-W"))
		if dur, err := parseLsfTime(timeStr); err == nil {
			specs.Time = dur
		}
		return nil
	}

	// Output file: -o path
	if strings.HasPrefix(flag, "-o ") {
		specs.Stdout = strings.TrimSpace(strings.TrimPrefix(flag, "-o"))
		return nil
	}

	// Error file: -e path
	if strings.HasPrefix(flag, "-e ") {
		specs.Stderr = strings.TrimSpace(strings.TrimPrefix(flag, "-e"))
		return nil
	}

	// Email on begin: -B (standalone flag)
	if flag == "-B" {
		specs.EmailOnBegin = true
		return nil
	}

	// Email on end: -N (standalone flag)
	if flag == "-N" {
		specs.EmailOnEnd = true
		return nil
	}

	// Mail user: -u email
	if strings.HasPrefix(flag, "-u ") {
		specs.MailUser = strings.TrimSpace(strings.TrimPrefix(flag, "-u"))
		return nil
	}

	// GPU specification: -gpu "num=N[:type=T]"
	if strings.HasPrefix(flag, "-gpu ") {
		gpuStr := strings.TrimSpace(strings.TrimPrefix(flag, "-gpu"))
		gpuStr = strings.Trim(gpuStr, "\"'")
		gpu := parseLsfGpuDirective(gpuStr)
		if gpu != nil {
			specs.Gpu = gpu
		}
		return nil
	}

	// Resource requirement: -R "rusage[...]"
	if strings.HasPrefix(flag, "-R ") {
		resStr := strings.TrimSpace(strings.TrimPrefix(flag, "-R"))
		resStr = strings.Trim(resStr, "\"'")
		l.parseLsfResource(resStr, specs)
		return nil
	}

	return nil
}

// parseLsfResource parses LSF -R resource requirement strings
// Supports: rusage[mem=N], rusage[ngpus_physical=N], span[hosts=1], etc.
func (l *LsfScheduler) parseLsfResource(resStr string, specs *ScriptSpecs) {
	// Look for rusage[...] block
	rusageIdx := strings.Index(resStr, "rusage[")
	if rusageIdx < 0 {
		return
	}

	// Extract content between brackets
	start := rusageIdx + len("rusage[")
	end := strings.Index(resStr[start:], "]")
	if end < 0 {
		return
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
			if mem, err := parseLsfMemory(value); err == nil {
				specs.MemMB = mem
			}
		case "ngpus_physical", "ngpus":
			if count, err := strconv.Atoi(value); err == nil {
				specs.Gpu = &GpuSpec{
					Type:  "gpu",
					Count: count,
					Raw:   fmt.Sprintf("%s=%s", key, value),
				}
			}
		}
	}
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

	// Create output directory if specified
	if outputDir != "" {
		if err := os.MkdirAll(outputDir, utils.PermDir); err != nil {
			return "", NewScriptCreationError(jobSpec.Name, outputDir, err)
		}
	}

	// Set log path based on job name
	if jobSpec.Name != "" {
		safeName := strings.ReplaceAll(jobSpec.Name, "/", "--")
		specs.Stdout = filepath.Join(outputDir, fmt.Sprintf("%s.log", safeName))
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

	// Write parsed BSUB directives with overrides
	for _, flag := range specs.RawFlags {
		// Skip output/error flags if we're overriding them
		if specs.Stdout != "" && strings.HasPrefix(flag, "-o ") {
			continue
		}
		// Always skip stderr flags - we don't want separate error logs
		if strings.HasPrefix(flag, "-e ") {
			continue
		}
		// Skip job-name flags - we'll use specs.JobName if set
		if specs.JobName != "" && strings.HasPrefix(flag, "-J ") {
			continue
		}
		// Skip email flags - they'll be regenerated from specs
		if flag == "-B" || flag == "-N" || strings.HasPrefix(flag, "-u ") {
			continue
		}
		// Skip walltime flags if we have a parsed time
		if specs.Time > 0 && strings.HasPrefix(flag, "-W ") {
			continue
		}
		fmt.Fprintf(writer, "#BSUB %s\n", flag)
	}

	// Add job name if specified
	if specs.JobName != "" {
		fmt.Fprintf(writer, "#BSUB -J %s\n", specs.JobName)
	}

	// Add custom stdout if specified
	if specs.Stdout != "" {
		fmt.Fprintf(writer, "#BSUB -o %s\n", specs.Stdout)
	}

	// Add email notifications if specified
	if specs.EmailOnBegin {
		fmt.Fprintln(writer, "#BSUB -B")
	}
	if specs.EmailOnEnd {
		fmt.Fprintln(writer, "#BSUB -N")
	}
	if specs.MailUser != "" {
		fmt.Fprintf(writer, "#BSUB -u %s\n", specs.MailUser)
	}

	// Add walltime if specified
	if specs.Time > 0 {
		fmt.Fprintf(writer, "#BSUB -W %s\n", formatLsfTime(specs.Time))
	}

	// Add blank line
	fmt.Fprintln(writer, "")

	// Print job information at start
	fmt.Fprintln(writer, "# Print job information")
	fmt.Fprintln(writer, "_START_TIME=$SECONDS")
	fmt.Fprintln(writer, "_format_time() { local s=$1; printf '%02d:%02d:%02d' $((s/3600)) $((s%3600/60)) $((s%60)); }")
	fmt.Fprintln(writer, "echo \"========================================\"")
	fmt.Fprintln(writer, "echo \"Job ID:    $LSB_JOBID\"")
	fmt.Fprintf(writer, "echo \"Job Name:  %s\"\n", specs.JobName)
	if specs.Ncpus > 0 {
		fmt.Fprintf(writer, "echo \"CPUs:      %d\"\n", specs.Ncpus)
	}
	if specs.MemMB > 0 {
		fmt.Fprintf(writer, "echo \"Memory:    %d MB\"\n", specs.MemMB)
	}
	if specs.Time > 0 {
		fmt.Fprintf(writer, "echo \"Time:      %s\"\n", formatLsfTime(specs.Time))
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

	// Self-delete the script after execution
	fmt.Fprintf(writer, "rm -f %s\n", scriptPath)

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
	}

	if info.MaxCpusPerNode == 0 && info.MaxMemMBPerNode == 0 && len(info.AvailableGpus) == 0 {
		return nil, ErrClusterInfoUnavailable
	}

	return info, nil
}

// getHostResources queries LSF for host resources using bhosts
func (l *LsfScheduler) getHostResources() (int, int64, error) {
	cmd := exec.Command(l.bhostsBin, "-w")
	output, err := cmd.CombinedOutput()
	if err != nil {
		return 0, 0, NewClusterError("LSF", "query hosts", err)
	}

	var maxCpus int
	var maxMemMB int64

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

	return maxCpus, maxMemMB, nil
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

// TryParseLsfScript attempts to parse an LSF script without requiring LSF binaries.
// This is a static parser that can work in any environment.
func TryParseLsfScript(scriptPath string) (*ScriptSpecs, error) {
	parser := &LsfScheduler{
		directiveRe: regexp.MustCompile(`^\s*#BSUB\s+(.+)$`),
	}
	return parser.ReadScriptSpecs(scriptPath)
}
