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
)

// PbsScheduler implements the Scheduler interface for PBS/Torque
type PbsScheduler struct {
	qsubBin     string
	qstatBin    string
	directiveRe *regexp.Regexp
	jobIDRe     *regexp.Regexp
}

// NewPbsScheduler creates a new PBS scheduler instance using qsub from PATH
func NewPbsScheduler() (*PbsScheduler, error) {
	return newPbsSchedulerWithBinary("")
}

// NewPbsSchedulerWithBinary creates a PBS scheduler using an explicit qsub path
func NewPbsSchedulerWithBinary(qsubBin string) (*PbsScheduler, error) {
	return newPbsSchedulerWithBinary(qsubBin)
}

func newPbsSchedulerWithBinary(qsubBin string) (*PbsScheduler, error) {
	binPath := qsubBin
	if binPath == "" {
		var err error
		binPath, err = exec.LookPath("qsub")
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

	qstatCmd, _ := exec.LookPath("qstat")

	return &PbsScheduler{
		qsubBin:     binPath,
		qstatBin:    qstatCmd,
		directiveRe: regexp.MustCompile(`^\s*#PBS\s+(.+)$`),
		jobIDRe:     regexp.MustCompile(`^(\d+\..*|^\d+)$`),
	}, nil
}

// IsAvailable checks if PBS is available and we're not inside a PBS job
func (p *PbsScheduler) IsAvailable() bool {
	if p.qsubBin == "" {
		return false
	}

	// Check if we're already inside a PBS job
	_, inJob := os.LookupEnv("PBS_JOBID")
	if inJob {
		return false
	}

	return true
}

// GetInfo returns information about the PBS scheduler
func (p *PbsScheduler) GetInfo() *SchedulerInfo {
	_, inJob := os.LookupEnv("PBS_JOBID")
	available := p.IsAvailable()

	info := &SchedulerInfo{
		Type:      "PBS",
		Binary:    p.qsubBin,
		InJob:     inJob,
		Available: available,
	}

	// Try to get PBS version
	if p.qsubBin != "" {
		if version, err := p.getPbsVersion(); err == nil {
			info.Version = version
		}
	}

	return info
}

// getPbsVersion attempts to get the PBS version
func (p *PbsScheduler) getPbsVersion() (string, error) {
	cmd := exec.Command(p.qsubBin, "--version")
	output, err := cmd.Output()
	if err != nil {
		return "", err
	}

	versionStr := strings.TrimSpace(string(output))
	return versionStr, nil
}

// ReadScriptSpecs parses #PBS directives from a build script
func (p *PbsScheduler) ReadScriptSpecs(scriptPath string) (*ScriptSpecs, error) {
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

		// Parse #PBS directives only
		if matches := p.directiveRe.FindStringSubmatch(line); matches != nil {
			flag := strings.TrimSpace(matches[1])
			specs.RawFlags = append(specs.RawFlags, flag)

			// Parse PBS options
			if err := p.parsePbsFlag(flag, specs); err != nil {
				return nil, NewParseError("PBS", lineNum, line, err.Error())
			}
		}
	}

	if err := scanner.Err(); err != nil {
		return nil, fmt.Errorf("error reading script: %w", err)
	}

	return specs, nil
}

// parsePbsFlag parses individual PBS flags and updates specs
// PBS uses -l for resources, -N for name, -o/-e for output/error, -m for mail
func (p *PbsScheduler) parsePbsFlag(flag string, specs *ScriptSpecs) error {
	// Job name: -N jobname
	if strings.HasPrefix(flag, "-N ") {
		specs.JobName = strings.TrimSpace(strings.TrimPrefix(flag, "-N"))
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

	// Mail user: -M email
	if strings.HasPrefix(flag, "-M ") {
		specs.MailUser = strings.TrimSpace(strings.TrimPrefix(flag, "-M"))
		return nil
	}

	// Mail options: -m [a|b|e|n] (abort, begin, end, none)
	if strings.HasPrefix(flag, "-m ") {
		mailOpts := strings.TrimSpace(strings.TrimPrefix(flag, "-m"))
		if mailOpts == "n" {
			specs.EmailOnBegin = false
			specs.EmailOnEnd = false
			specs.EmailOnFail = false
		} else {
			for _, c := range mailOpts {
				switch c {
				case 'a': // abort/fail
					specs.EmailOnFail = true
				case 'b': // begin
					specs.EmailOnBegin = true
				case 'e': // end
					specs.EmailOnEnd = true
				}
			}
		}
		return nil
	}

	// Resource list: -l resource=value[,resource=value,...]
	if strings.HasPrefix(flag, "-l ") {
		resourceStr := strings.TrimSpace(strings.TrimPrefix(flag, "-l"))
		return p.parseResourceList(resourceStr, specs)
	}

	return nil
}

// parseResourceList parses PBS -l resource specifications
// Examples:
//   - select=1:ncpus=4:mem=8gb
//   - walltime=01:00:00
//   - ngpus=1
//   - nodes=1:ppn=4
func (p *PbsScheduler) parseResourceList(resourceStr string, specs *ScriptSpecs) error {
	// Split by comma, but be careful with select= which can have colons
	resources := strings.Split(resourceStr, ",")

	for _, res := range resources {
		res = strings.TrimSpace(res)
		if res == "" {
			continue
		}

		// Handle select= which contains colon-separated resources
		if strings.HasPrefix(res, "select=") {
			selectParts := strings.Split(res, ":")
			for _, part := range selectParts {
				if err := p.parseSingleResource(part, specs); err != nil {
					return err
				}
			}
			continue
		}

		// Handle nodes=N:ppn=M format
		if strings.HasPrefix(res, "nodes=") {
			nodeParts := strings.Split(res, ":")
			for _, part := range nodeParts {
				if strings.HasPrefix(part, "ppn=") {
					ppnStr := strings.TrimPrefix(part, "ppn=")
					if ppn, err := strconv.Atoi(ppnStr); err == nil {
						specs.Ncpus = ppn
					}
				}
			}
			continue
		}

		if err := p.parseSingleResource(res, specs); err != nil {
			return err
		}
	}

	return nil
}

// parseSingleResource parses a single resource=value specification
func (p *PbsScheduler) parseSingleResource(res string, specs *ScriptSpecs) error {
	parts := strings.SplitN(res, "=", 2)
	if len(parts) != 2 {
		return nil // Not a key=value, skip
	}

	key := strings.TrimSpace(parts[0])
	value := strings.TrimSpace(parts[1])

	switch key {
	case "ncpus":
		if ncpus, err := strconv.Atoi(value); err == nil {
			specs.Ncpus = ncpus
		}

	case "mem":
		if mem, err := parseMemoryString(value); err == nil {
			specs.MemMB = mem
		}

	case "walltime":
		if dur, err := parsePbsTime(value); err == nil {
			specs.Time = dur
		}

	case "ngpus":
		if count, err := strconv.Atoi(value); err == nil {
			specs.Gpu = &GpuSpec{
				Type:  "gpu",
				Count: count,
				Raw:   res,
			}
		}

	case "gpus":
		gpu, err := parseGpuString(value)
		if err == nil {
			specs.Gpu = gpu
		}
	}

	return nil
}

// CreateScriptWithSpec generates a PBS batch script
func (p *PbsScheduler) CreateScriptWithSpec(jobSpec *JobSpec, outputDir string) (string, error) {
	specs := jobSpec.Specs

	// Create output directory if specified
	if outputDir != "" {
		if err := os.MkdirAll(outputDir, 0775); err != nil {
			return "", NewScriptCreationError(jobSpec.Name, outputDir, err)
		}
	}

	// Create log directory if it doesn't exist
	logDir := filepath.Join(outputDir, "log")
	if err := os.MkdirAll(logDir, 0775); err != nil {
		return "", NewScriptCreationError(jobSpec.Name, logDir, err)
	}

	// Always set log path to our standardized location
	if jobSpec.Name != "" {
		safeName := strings.ReplaceAll(jobSpec.Name, "/", "--")
		specs.Stdout = filepath.Join("log", fmt.Sprintf("%s.log", safeName))
	}

	// Generate script filename
	scriptName := "pbs_job.sh"
	if jobSpec.Name != "" {
		// Replace slashes with -- to avoid creating subdirectories
		safeName := strings.ReplaceAll(jobSpec.Name, "/", "--")
		scriptName = fmt.Sprintf("pbs_%s.sh", safeName)
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

	// Write parsed PBS directives with overrides
	for _, flag := range specs.RawFlags {
		// Skip output/error flags if we're overriding them
		if specs.Stdout != "" && strings.HasPrefix(flag, "-o ") {
			continue
		}
		// Always skip stderr flags - we don't want separate error logs
		if strings.HasPrefix(flag, "-e ") {
			continue
		}
		fmt.Fprintf(writer, "#PBS %s\n", flag)
	}

	// Add custom stdout if specified (stderr is not used - output goes to stdout)
	if specs.Stdout != "" {
		fmt.Fprintf(writer, "#PBS -o %s\n", specs.Stdout)
	}

	// Add email notifications if specified
	if specs.EmailOnBegin || specs.EmailOnEnd || specs.EmailOnFail {
		var mailOpts string
		if specs.EmailOnBegin {
			mailOpts += "b"
		}
		if specs.EmailOnEnd {
			mailOpts += "e"
		}
		if specs.EmailOnFail {
			mailOpts += "a"
		}
		fmt.Fprintf(writer, "#PBS -m %s\n", mailOpts)
	}
	if specs.MailUser != "" {
		fmt.Fprintf(writer, "#PBS -M %s\n", specs.MailUser)
	}

	// Add blank line
	fmt.Fprintln(writer, "")

	// Write the command
	fmt.Fprintln(writer, jobSpec.Command)

	// Add job ID echo for tracking
	fmt.Fprintln(writer, "echo PBS_JOBID $PBS_JOBID")

	// Self-delete the script after execution
	fmt.Fprintf(writer, "rm -f %s\n", scriptPath)

	// Make executable
	if err := os.Chmod(scriptPath, 0775); err != nil {
		return "", NewScriptCreationError(jobSpec.Name, scriptPath, err)
	}

	return scriptPath, nil
}

// Submit submits a PBS job with optional dependency chain
func (p *PbsScheduler) Submit(scriptPath string, dependencyJobIDs []string) (string, error) {
	args := []string{scriptPath}

	// Add dependency if provided
	if len(dependencyJobIDs) > 0 {
		depStr := strings.Join(dependencyJobIDs, ":")
		depArg := fmt.Sprintf("-W depend=afterok:%s", depStr)
		args = append([]string{depArg}, args...)
	}

	// Execute qsub
	cmd := exec.Command(p.qsubBin, args...)
	output, err := cmd.CombinedOutput()
	if err != nil {
		return "", NewSubmissionError("PBS", filepath.Base(scriptPath), string(output), err)
	}

	// Parse job ID from output
	jobID := strings.TrimSpace(string(output))
	if jobID == "" {
		return "", fmt.Errorf("%w: %s", ErrJobIDParseFailed, string(output))
	}

	return jobID, nil
}

// GetClusterInfo retrieves cluster configuration (GPUs, limits)
// Returns nil if information is not available
func (p *PbsScheduler) GetClusterInfo() (*ClusterInfo, error) {
	// PBS cluster info retrieval not yet implemented
	// This would require parsing qstat -Q or pbsnodes output
	return nil, ErrClusterInfoUnavailable
}

// parseMemoryString converts memory strings like "8G", "1024M", "8gb" to MB
func parseMemoryString(memStr string) (int64, error) {
	memStr = strings.ToUpper(strings.TrimSpace(memStr))

	var value int64
	var unit string

	n, err := fmt.Sscanf(memStr, "%d%s", &value, &unit)
	if err != nil && n == 0 {
		return 0, fmt.Errorf("%w: %s", ErrInvalidMemoryFormat, memStr)
	}

	switch unit {
	case "G", "GB":
		return value * 1024, nil
	case "M", "MB", "":
		return value, nil
	case "K", "KB":
		return value / 1024, nil
	case "T", "TB":
		return value * 1024 * 1024, nil
	default:
		return value, nil
	}
}

// parsePbsTime parses PBS walltime format: HH:MM:SS
func parsePbsTime(timeStr string) (time.Duration, error) {
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

// parseGpuString parses GPU specifications from various formats
// Supports: "2", "a100:1", "gpu:a100:2", "nvidia_h100:2"
func parseGpuString(gpuStr string) (*GpuSpec, error) {
	if gpuStr == "" {
		return nil, nil
	}

	spec := &GpuSpec{
		Raw:   gpuStr,
		Count: 1,
	}

	// Check if this is a MIG profile (contains dots and underscores)
	if strings.Contains(gpuStr, ".") && strings.Contains(gpuStr, "_") {
		parts := strings.Split(gpuStr, ":")
		if len(parts) == 1 {
			spec.Type = parts[0]
			spec.Count = 1
		} else if len(parts) == 2 {
			spec.Type = parts[0]
			if count, err := strconv.Atoi(parts[1]); err == nil {
				spec.Count = count
			}
		}
		return spec, nil
	}

	parts := strings.Split(gpuStr, ":")

	switch len(parts) {
	case 1:
		if count, err := strconv.Atoi(parts[0]); err == nil {
			spec.Count = count
			spec.Type = "gpu"
		} else {
			spec.Type = parts[0]
			spec.Count = 1
		}

	case 2:
		spec.Type = parts[0]
		if count, err := strconv.Atoi(parts[1]); err == nil {
			spec.Count = count
		}

	case 3:
		spec.Type = parts[1]
		if count, err := strconv.Atoi(parts[2]); err == nil {
			spec.Count = count
		}

	default:
		lastColon := strings.LastIndex(gpuStr, ":")
		if lastColon > 0 {
			spec.Type = gpuStr[:lastColon]
			if count, err := strconv.Atoi(gpuStr[lastColon+1:]); err == nil {
				spec.Count = count
			}
		} else {
			spec.Type = gpuStr
		}
	}

	return spec, nil
}

// TryParsePbsScript attempts to parse a PBS script without requiring PBS binaries.
// This is a static parser that can work in any environment.
func TryParsePbsScript(scriptPath string) (*ScriptSpecs, error) {
	parser := &PbsScheduler{
		directiveRe: regexp.MustCompile(`^\s*#PBS\s+(.+)$`),
	}
	return parser.ReadScriptSpecs(scriptPath)
}
