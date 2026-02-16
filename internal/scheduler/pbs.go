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

// PbsScheduler implements the Scheduler interface for PBS/Torque
// EXPERIMENTAL: PBS is not tested on real clusters and may have edge cases. Feedback welcome.
type PbsScheduler struct {
	qsubBin     string
	qstatBin    string
	pbsnodesBin string
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
	pbsnodesCmd, _ := exec.LookPath("pbsnodes")

	return &PbsScheduler{
		qsubBin:     binPath,
		qstatBin:    qstatCmd,
		pbsnodesBin: pbsnodesCmd,
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

		// Handle select=N:ncpus=M:mpiprocs=P:mem=X format
		if strings.HasPrefix(res, "select=") {
			selectParts := strings.Split(res, ":")
			var mpiprocs int
			for _, part := range selectParts {
				// Extract chunk count (node count) from select=N
				if strings.HasPrefix(part, "select=") {
					selectVal := strings.TrimPrefix(part, "select=")
					if n, err := strconv.Atoi(selectVal); err == nil {
						specs.Nodes = n
					}
				} else if strings.HasPrefix(part, "mpiprocs=") {
					// mpiprocs is tasks per chunk (per node)
					mpStr := strings.TrimPrefix(part, "mpiprocs=")
					if mp, err := strconv.Atoi(mpStr); err == nil {
						mpiprocs = mp
					}
				} else {
					if err := p.parseSingleResource(part, specs); err != nil {
						return err
					}
				}
			}
			// Compute total tasks from nodes * mpiprocs
			if mpiprocs > 0 {
				specs.Ntasks = specs.Nodes * mpiprocs
			}
			continue
		}

		// Handle nodes=N:ppn=M format
		if strings.HasPrefix(res, "nodes=") {
			nodeParts := strings.Split(res, ":")
			for _, part := range nodeParts {
				if strings.HasPrefix(part, "nodes=") {
					nodesStr := strings.TrimPrefix(part, "nodes=")
					if n, err := strconv.Atoi(nodesStr); err == nil {
						specs.Nodes = n
					}
				} else if strings.HasPrefix(part, "ppn=") {
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

	// Create output directory if specified (scripts still live in outputDir)
	if outputDir != "" {
		if err := os.MkdirAll(outputDir, utils.PermDir); err != nil {
			return "", NewScriptCreationError(jobSpec.Name, outputDir, err)
		}
	}

	// Set log path based on job name (logs go to outputDir, which caller controls)
	if jobSpec.Name != "" {
		safeName := strings.ReplaceAll(jobSpec.Name, "/", "--")
		specs.Stdout = filepath.Join(outputDir, fmt.Sprintf("%s.log", safeName))
	}

	// Generate script filename
	scriptName := "job.pbs"
	if jobSpec.Name != "" {
		// Replace slashes with -- to avoid creating subdirectories
		safeName := strings.ReplaceAll(jobSpec.Name, "/", "--")
		scriptName = fmt.Sprintf("%s.pbs", safeName)
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

	// Normalize node/task defaults
	if specs.Nodes <= 0 {
		specs.Nodes = 1
	}
	if specs.Ntasks <= 0 {
		specs.Ntasks = 1
	}

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
		// Skip job-name flags - we'll use specs.JobName if set
		if specs.JobName != "" && strings.HasPrefix(flag, "-N ") {
			continue
		}
		// Skip email flags - they'll be regenerated from specs
		if strings.HasPrefix(flag, "-m ") || strings.HasPrefix(flag, "-M ") {
			continue
		}
		// Skip resource flags that we'll regenerate (select=, nodes=, walltime=)
		if strings.HasPrefix(flag, "-l ") {
			resourceStr := strings.TrimSpace(strings.TrimPrefix(flag, "-l"))
			if strings.HasPrefix(resourceStr, "select=") || strings.HasPrefix(resourceStr, "nodes=") {
				continue
			}
			if strings.HasPrefix(resourceStr, "walltime=") && specs.Time > 0 {
				continue
			}
		}
		fmt.Fprintf(writer, "#PBS %s\n", flag)
	}

	// Add job name if specified
	if specs.JobName != "" {
		fmt.Fprintf(writer, "#PBS -N %s\n", specs.JobName)
	}

	// Add custom stdout if specified (stderr is not used - output goes to stdout)
	if specs.Stdout != "" {
		fmt.Fprintf(writer, "#PBS -o %s\n", specs.Stdout)
	}

	// Generate resource specification: select=N:ncpus=M[:mpiprocs=P][:mem=Xmb][:ngpus=G]
	selectParts := []string{fmt.Sprintf("select=%d", specs.Nodes)}
	if specs.Ncpus > 0 {
		selectParts = append(selectParts, fmt.Sprintf("ncpus=%d", specs.Ncpus))
	}
	if specs.Ntasks > 1 && specs.Nodes > 0 {
		mpiprocs := specs.Ntasks / specs.Nodes
		if mpiprocs > 1 {
			selectParts = append(selectParts, fmt.Sprintf("mpiprocs=%d", mpiprocs))
		}
	}
	if specs.MemMB > 0 {
		selectParts = append(selectParts, fmt.Sprintf("mem=%dmb", specs.MemMB))
	}
	if specs.Gpu != nil && specs.Gpu.Count > 0 {
		selectParts = append(selectParts, fmt.Sprintf("ngpus=%d", specs.Gpu.Count))
	}
	fmt.Fprintf(writer, "#PBS -l %s\n", strings.Join(selectParts, ":"))

	// Generate walltime
	if specs.Time > 0 {
		hours := int(specs.Time.Hours())
		mins := int(specs.Time.Minutes()) % 60
		secs := int(specs.Time.Seconds()) % 60
		fmt.Fprintf(writer, "#PBS -l walltime=%02d:%02d:%02d\n", hours, mins, secs)
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

	// Print job information at start
	fmt.Fprintln(writer, "# Print job information")
	fmt.Fprintln(writer, "_START_TIME=$SECONDS")
	fmt.Fprintln(writer, "_format_time() { local s=$1; printf '%02d:%02d:%02d' $((s/3600)) $((s%3600/60)) $((s%60)); }")
	fmt.Fprintln(writer, "echo \"========================================\"")
	fmt.Fprintln(writer, "echo \"Job ID:    $PBS_JOBID\"")
	fmt.Fprintf(writer, "echo \"Job Name:  %s\"\n", specs.JobName)
	if specs.Nodes > 1 {
		fmt.Fprintf(writer, "echo \"Nodes:     %d\"\n", specs.Nodes)
	}
	if specs.Ntasks > 1 {
		fmt.Fprintf(writer, "echo \"Tasks:     %d\"\n", specs.Ntasks)
	}
	if specs.Ncpus > 0 {
		fmt.Fprintf(writer, "echo \"CPUs/Task: %d\"\n", specs.Ncpus)
	}
	if specs.MemMB > 0 {
		fmt.Fprintf(writer, "echo \"Memory:    %d MB\"\n", specs.MemMB)
	}
	if specs.Time > 0 {
		hours := int(specs.Time.Hours())
		mins := int(specs.Time.Minutes()) % 60
		secs := int(specs.Time.Seconds()) % 60
		fmt.Fprintf(writer, "echo \"Time:      %02d:%02d:%02d\"\n", hours, mins, secs)
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
	fmt.Fprintln(writer, "echo \"Job ID:    $PBS_JOBID\"")
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
	info := &ClusterInfo{
		AvailableGpus: make([]GpuInfo, 0),
		Limits:        make([]ResourceLimits, 0),
	}

	// Get node resources using pbsnodes
	if p.pbsnodesBin != "" {
		maxCpus, maxMem, gpus, err := p.getNodeResources()
		if err == nil {
			info.MaxCpusPerNode = maxCpus
			info.MaxMemMBPerNode = maxMem
			info.AvailableGpus = gpus
		}
	}

	if info.MaxCpusPerNode == 0 && info.MaxMemMBPerNode == 0 && len(info.AvailableGpus) == 0 {
		return nil, ErrClusterInfoUnavailable
	}

	return info, nil
}

// getNodeResources queries PBS for node resources using pbsnodes
func (p *PbsScheduler) getNodeResources() (int, int64, []GpuInfo, error) {
	cmd := exec.Command(p.pbsnodesBin, "-a")
	output, err := cmd.CombinedOutput()
	if err != nil {
		return 0, 0, nil, NewClusterError("PBS", "query nodes", err)
	}

	var maxCpus int
	var maxMemMB int64
	gpuMap := make(map[string]*GpuInfo)

	// Parse pbsnodes output
	// Format is key = value pairs, with node names as headers
	lines := strings.Split(string(output), "\n")
	var currentNode string
	var nodeCpus int
	var nodeMemKB int64
	var nodeState string

	for _, line := range lines {
		line = strings.TrimSpace(line)

		// Node name line (no leading whitespace in original, but we trimmed)
		if line != "" && !strings.Contains(line, "=") {
			// Save previous node's data
			if currentNode != "" && (nodeState == "free" || nodeState == "job-exclusive" || strings.Contains(nodeState, "busy")) {
				if nodeCpus > maxCpus {
					maxCpus = nodeCpus
				}
				nodeMB := nodeMemKB / 1024
				if nodeMB > maxMemMB {
					maxMemMB = nodeMB
				}
			}
			currentNode = line
			nodeCpus = 0
			nodeMemKB = 0
			nodeState = ""
			continue
		}

		// Parse key = value pairs
		if strings.Contains(line, "=") {
			parts := strings.SplitN(line, "=", 2)
			if len(parts) != 2 {
				continue
			}
			key := strings.TrimSpace(parts[0])
			value := strings.TrimSpace(parts[1])

			switch key {
			case "state":
				nodeState = value
			case "np", "pcpus":
				fmt.Sscanf(value, "%d", &nodeCpus)
			case "resources_available.ncpus":
				fmt.Sscanf(value, "%d", &nodeCpus)
			case "resources_available.mem":
				// Memory can be in kb, mb, gb format
				nodeMemKB, _ = parsePbsMemory(value)
			case "resources_available.ngpus":
				var gpuCount int
				fmt.Sscanf(value, "%d", &gpuCount)
				if gpuCount > 0 {
					if existing, ok := gpuMap["gpu"]; ok {
						existing.Total += gpuCount
						if nodeState == "free" {
							existing.Available += gpuCount
						}
					} else {
						available := 0
						if nodeState == "free" {
							available = gpuCount
						}
						gpuMap["gpu"] = &GpuInfo{
							Type:      "gpu",
							Total:     gpuCount,
							Available: available,
						}
					}
				}
			}
		}
	}

	// Don't forget the last node
	if currentNode != "" && (nodeState == "free" || nodeState == "job-exclusive" || strings.Contains(nodeState, "busy")) {
		if nodeCpus > maxCpus {
			maxCpus = nodeCpus
		}
		nodeMB := nodeMemKB / 1024
		if nodeMB > maxMemMB {
			maxMemMB = nodeMB
		}
	}

	gpus := make([]GpuInfo, 0, len(gpuMap))
	for _, info := range gpuMap {
		gpus = append(gpus, *info)
	}

	return maxCpus, maxMemMB, gpus, nil
}

// parsePbsMemory parses PBS memory format (e.g., "8gb", "1024mb", "1048576kb")
func parsePbsMemory(memStr string) (int64, error) {
	memStr = strings.ToLower(strings.TrimSpace(memStr))

	var value int64
	var unit string

	n, err := fmt.Sscanf(memStr, "%d%s", &value, &unit)
	if err != nil && n == 0 {
		return 0, fmt.Errorf("invalid memory format: %s", memStr)
	}

	switch unit {
	case "tb":
		return value * 1024 * 1024 * 1024, nil // Return in KB
	case "gb":
		return value * 1024 * 1024, nil // Return in KB
	case "mb":
		return value * 1024, nil // Return in KB
	case "kb", "":
		return value, nil // Already in KB
	case "b":
		return value / 1024, nil // Convert to KB
	default:
		return value, nil
	}
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

// GetJobResources reads allocated resources from PBS environment variables.
func (s *PbsScheduler) GetJobResources() *JobResources {
	if _, ok := os.LookupEnv("PBS_JOBID"); !ok {
		return nil
	}
	res := &JobResources{}
	res.Ncpus = getEnvInt("PBS_NCPUS")
	if res.Ncpus == nil {
		res.Ncpus = getEnvInt("NCPUS")
	}
	res.Nodes = getEnvInt("PBS_NUM_NODES")
	res.Ntasks = getEnvInt("PBS_NP")
	if res.Ntasks == nil {
		res.Ntasks = getEnvInt("PBS_TASKNUM")
	}
	// PBS_VMEM is in bytes
	if vmem := getEnvInt64("PBS_VMEM"); vmem != nil {
		mb := *vmem / (1024 * 1024)
		if mb > 0 {
			res.MemMB = &mb
		}
	}
	res.Ngpus = getCudaDeviceCount()
	return res
}

// TryParsePbsScript attempts to parse a PBS script without requiring PBS binaries.
// This is a static parser that can work in any environment.
func TryParsePbsScript(scriptPath string) (*ScriptSpecs, error) {
	parser := &PbsScheduler{
		directiveRe: regexp.MustCompile(`^\s*#PBS\s+(.+)$`),
	}
	return parser.ReadScriptSpecs(scriptPath)
}
