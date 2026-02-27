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
	lines, err := readFileLines(scriptPath)
	if err != nil {
		return nil, err
	}
	return parseScript(scriptPath, lines, p.extractDirectives, p.parseRuntimeConfig, p.parseResourceSpec)
}

// extractDirectives extracts raw directive strings from script lines (strips the #PBS prefix).
func (p *PbsScheduler) extractDirectives(lines []string) []string {
	var out []string
	for _, line := range lines {
		if m := p.directiveRe.FindStringSubmatch(line); m != nil {
			out = append(out, utils.StripInlineComment(m[1]))
		}
	}
	return out
}

// parseRuntimeConfig consumes PBS job control directives (name, I/O, email) from the directive list.
// Returns the populated RuntimeConfig, unconsumed directives, and any critical error.
// PBS control fields never hard-fail — unknown or malformed control flags are passed through.
func (p *PbsScheduler) parseRuntimeConfig(directives []string) (RuntimeConfig, []string, error) {
	var rc RuntimeConfig
	var unconsumed []string

	for _, flag := range directives {
		consumed := true
		switch {
		case flagMatches(flag, "-N"):
			rc.JobName, _ = flagValue(flag, "-N")
		case flagMatches(flag, "-o"):
			v, _ := flagValue(flag, "-o")
			rc.Stdout = absPath(v)
		case flagMatches(flag, "-e"):
			v, _ := flagValue(flag, "-e")
			rc.Stderr = absPath(v)
		case flagMatches(flag, "-q"):
			rc.Partition, _ = flagValue(flag, "-q")
		case flagMatches(flag, "-M"):
			rc.MailUser, _ = flagValue(flag, "-M")
		case flagMatches(flag, "-m"):
			mailOpts, _ := flagValue(flag, "-m")
			if mailOpts == "n" {
				rc.EmailOnBegin = false
				rc.EmailOnEnd = false
				rc.EmailOnFail = false
			} else {
				for _, c := range mailOpts {
					switch c {
					case 'a': // abort/fail
						rc.EmailOnFail = true
					case 'b': // begin
						rc.EmailOnBegin = true
					case 'e': // end
						rc.EmailOnEnd = true
					}
				}
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

// parseResourceSpec consumes PBS resource directives from the directive list.
// Returns a populated ResourceSpec and any unconsumed directives.
// Returns nil ResourceSpec on any parse error (passthrough mode).
func (p *PbsScheduler) parseResourceSpec(directives []string) (*ResourceSpec, []string) {
	defaults := GetSpecDefaults()
	rs := &ResourceSpec{
		CpusPerTask:  defaults.CpusPerTask,
		TasksPerNode: defaults.TasksPerNode,
		Nodes:        defaults.Nodes,
		MemPerNodeMB: defaults.MemPerNodeMB,
		Time:         defaults.Time,
	}

	var unconsumed []string

	for _, flag := range directives {
		// Exclusive node access: -x
		if flag == "-x" {
			rs.Exclusive = true
			continue
		}
		// Only -l flags carry resource specifications
		if !strings.HasPrefix(flag, "-l ") {
			unconsumed = append(unconsumed, flag)
			continue
		}

		resourceStr := strings.TrimSpace(strings.TrimPrefix(flag, "-l"))
		if err := p.parseResourceList(resourceStr, rs); err != nil {
			logParseWarning("PBS: failed to parse resource directive %q: %v", flag, err)
			return nil, directives
		}
	}

	return rs, unconsumed
}

// parseResourceList parses PBS -l resource specifications
// Examples:
//   - select=1:ncpus=4:mem=8gb
//   - walltime=01:00:00
//   - ngpus=1
//   - nodes=1:ppn=4
func (p *PbsScheduler) parseResourceList(resourceStr string, rs *ResourceSpec) error {
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
					n, err := strconv.Atoi(selectVal)
					if err != nil {
						return fmt.Errorf("invalid select value %q: %w", selectVal, err)
					}
					rs.Nodes = n
				} else if strings.HasPrefix(part, "mpiprocs=") {
					// mpiprocs is tasks per chunk (per node)
					mpStr := strings.TrimPrefix(part, "mpiprocs=")
					mp, err := strconv.Atoi(mpStr)
					if err != nil {
						return fmt.Errorf("invalid mpiprocs value %q: %w", mpStr, err)
					}
					mpiprocs = mp
				} else {
					if err := p.parseSingleResource(part, rs); err != nil {
						return err
					}
				}
			}
			// mpiprocs maps to TasksPerNode
			if mpiprocs > 0 {
				rs.TasksPerNode = mpiprocs
			}
			continue
		}

		// Handle nodes=N:ppn=M format
		if strings.HasPrefix(res, "nodes=") {
			nodeParts := strings.Split(res, ":")
			for _, part := range nodeParts {
				if strings.HasPrefix(part, "nodes=") {
					nodesStr := strings.TrimPrefix(part, "nodes=")
					n, err := strconv.Atoi(nodesStr)
					if err != nil {
						return fmt.Errorf("invalid nodes value %q: %w", nodesStr, err)
					}
					rs.Nodes = n
				} else if strings.HasPrefix(part, "ppn=") {
					ppnStr := strings.TrimPrefix(part, "ppn=")
					ppn, err := strconv.Atoi(ppnStr)
					if err != nil {
						return fmt.Errorf("invalid ppn value %q: %w", ppnStr, err)
					}
					rs.CpusPerTask = ppn
				}
			}
			continue
		}

		if err := p.parseSingleResource(res, rs); err != nil {
			return err
		}
	}

	return nil
}

// parseSingleResource parses a single resource=value specification
func (p *PbsScheduler) parseSingleResource(res string, rs *ResourceSpec) error {
	parts := strings.SplitN(res, "=", 2)
	if len(parts) != 2 {
		return nil // Not a key=value, skip
	}

	key := strings.TrimSpace(parts[0])
	value := strings.TrimSpace(parts[1])

	switch key {
	case "ncpus":
		ncpus, err := strconv.Atoi(value)
		if err != nil {
			return fmt.Errorf("invalid ncpus value %q: %w", value, err)
		}
		rs.CpusPerTask = ncpus

	case "mem":
		mem, err := parseMemoryMB(value)
		if err != nil {
			return fmt.Errorf("invalid mem value %q: %w", value, err)
		}
		rs.MemPerNodeMB = mem

	case "walltime":
		dur, err := utils.ParseHMSTime(value)
		if err != nil {
			return fmt.Errorf("invalid walltime value %q: %w", value, err)
		}
		rs.Time = dur

	case "ngpus":
		count, err := strconv.Atoi(value)
		if err != nil {
			return fmt.Errorf("invalid ngpus value %q: %w", value, err)
		}
		rs.Gpu = &GpuSpec{
			Type:  "gpu",
			Count: count,
			Raw:   res,
		}

	case "gpus":
		gpu, err := parseGpuString(value)
		if err != nil {
			return fmt.Errorf("invalid gpus value %q: %w", value, err)
		}
		rs.Gpu = gpu

	case "place":
		// place=excl (or place=excl:...) requests exclusive node access
		if strings.Contains(value, "excl") {
			rs.Exclusive = true
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

	// Set log path based on job name; only override if caller requests it or script has no output set
	// Capture before override: separate output when stderr is explicitly set to a different path
	arraySeparateOutput := jobSpec.Array != nil && specs.Control.Stderr != "" && specs.Control.Stderr != specs.Control.Stdout
	if jobSpec.Array != nil {
		// Array job: silence scheduler output; exec redirect in script body handles per-task logs
		specs.Control.Stdout = "/dev/null"
		specs.Control.Stderr = "/dev/null"
	} else if jobSpec.Name != "" && (jobSpec.OverrideOutput || specs.Control.Stdout == "") {
		specs.Control.Stdout = filepath.Join(outputDir, fmt.Sprintf("%s.log", safeJobName(jobSpec.Name)))
	}
	if specs.Control.Stderr == "" && specs.Control.Stdout != "" {
		specs.Control.Stderr = specs.Control.Stdout
	}

	// Generate script filename
	scriptName := "job.pbs"
	if jobSpec.Name != "" {
		scriptName = fmt.Sprintf("%s.pbs", safeJobName(jobSpec.Name))
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

	// Write passthrough flags (directives not consumed by Spec or Control)
	for _, flag := range specs.RemainingFlags {
		fmt.Fprintf(writer, "#PBS %s\n", flag)
	}

	// Write RuntimeConfig directives
	ctrl := specs.Control
	if ctrl.JobName != "" {
		fmt.Fprintf(writer, "#PBS -N %s\n", ctrl.JobName)
	}
	if ctrl.Stdout != "" {
		fmt.Fprintf(writer, "#PBS -o %s\n", ctrl.Stdout)
	}
	if ctrl.Stderr != "" {
		fmt.Fprintf(writer, "#PBS -e %s\n", ctrl.Stderr)
	}
	if ctrl.Partition != "" {
		fmt.Fprintf(writer, "#PBS -q %s\n", ctrl.Partition)
	}

	// Write ResourceSpec directives (only if not in passthrough mode)
	if specs.Spec != nil {
		rs := specs.Spec
		nodes := rs.Nodes
		if nodes <= 0 {
			nodes = 1
		}
		tasksPerNode := rs.TasksPerNode
		if tasksPerNode <= 0 {
			tasksPerNode = 1
		}

		// Generate resource specification: select=N:ncpus=M[:mpiprocs=P][:mem=Xmb][:ngpus=G]
		selectParts := []string{fmt.Sprintf("select=%d", nodes)}
		if rs.CpusPerTask > 0 {
			selectParts = append(selectParts, fmt.Sprintf("ncpus=%d", rs.CpusPerTask))
		}
		if tasksPerNode > 1 {
			selectParts = append(selectParts, fmt.Sprintf("mpiprocs=%d", tasksPerNode))
		}
		if rs.MemPerNodeMB > 0 {
			selectParts = append(selectParts, fmt.Sprintf("mem=%dmb", rs.MemPerNodeMB))
		}
		if rs.Gpu != nil && rs.Gpu.Count > 0 {
			selectParts = append(selectParts, fmt.Sprintf("ngpus=%d", rs.Gpu.Count))
		}
		fmt.Fprintf(writer, "#PBS -l %s\n", strings.Join(selectParts, ":"))

		if rs.Time > 0 {
			fmt.Fprintf(writer, "#PBS -l walltime=%s\n", formatHMSTime(rs.Time))
		}
		if rs.Exclusive {
			fmt.Fprintln(writer, "#PBS -l place=excl")
		}
	}

	// Add email notifications if specified
	if ctrl.EmailOnBegin || ctrl.EmailOnEnd || ctrl.EmailOnFail {
		var mailOpts string
		if ctrl.EmailOnBegin {
			mailOpts += "b"
		}
		if ctrl.EmailOnEnd {
			mailOpts += "e"
		}
		if ctrl.EmailOnFail {
			mailOpts += "a"
		}
		fmt.Fprintf(writer, "#PBS -m %s\n", mailOpts)
	}
	if ctrl.MailUser != "" {
		fmt.Fprintf(writer, "#PBS -M %s\n", ctrl.MailUser)
	}

	// Array directive
	if jobSpec.Array != nil {
		arr := jobSpec.Array
		r := fmt.Sprintf("1-%d", arr.Count)
		if arr.Limit > 0 {
			r += fmt.Sprintf("%%%d", arr.Limit)
		}
		fmt.Fprintf(writer, "#PBS -J %s\n", r)
	}

	fmt.Fprintln(writer, "")

	// Array job: extract input line, set ARRAY_ARGS, and redirect output
	if jobSpec.Array != nil {
		writeArrayBlock(writer, "$PBS_ARRAY_INDEX",
			jobSpec.Array.InputFile, outputDir, safeJobName(jobSpec.Name),
			jobSpec.Array.Count, arraySeparateOutput)
		jobSpec.Metadata["Array Job ID"] = "$PBS_JOBID"
		jobSpec.Metadata["Array Index"] = "$PBS_ARRAY_INDEX"
		jobSpec.Metadata["Array File"] = jobSpec.Array.InputFile
		jobSpec.Metadata["Array Args"] = "$ARRAY_ARGS"
	}

	// Print job information at start
	writeJobHeader(writer, "$PBS_JOBID", specs, formatHMSTime, jobSpec.Metadata)
	fmt.Fprintln(writer, "")

	// Write the command
	fmt.Fprintln(writer, jobSpec.Command)

	// Print completion info
	fmt.Fprintln(writer, "")
	writeJobFooter(writer, "$PBS_JOBID")

	// Self-delete the script after execution (unless in debug mode)
	if !debugMode {
		fmt.Fprintf(writer, "rm -f %s\n", scriptPath)
	}

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
		maxCpus, maxMem, err := p.getNodeResources()
		if err == nil {
			info.MaxCpusPerNode = maxCpus
			info.MaxMemMBPerNode = maxMem
		}

		// Get GPU info separately
		gpus, err := p.getGpuInfo()
		if err == nil {
			info.AvailableGpus = gpus
		}
	}

	// Get queue limits
	if p.qstatBin != "" {
		limits, err := p.getQueueLimits(info.AvailableGpus)
		if err == nil {
			info.Limits = limits
		}
	}

	if info.MaxCpusPerNode == 0 && info.MaxMemMBPerNode == 0 && len(info.AvailableGpus) == 0 {
		return nil, ErrClusterInfoUnavailable
	}

	return info, nil
}

// getNodeResources queries PBS for node resources using pbsnodes
func (p *PbsScheduler) getNodeResources() (int, int64, error) {
	cmd := exec.Command(p.pbsnodesBin, "-a")
	output, err := cmd.CombinedOutput()
	if err != nil {
		return 0, 0, NewClusterError("PBS", "query nodes", err)
	}

	var maxCpus int
	var maxMemMB int64

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

	return maxCpus, maxMemMB, nil
}

// getGpuInfo queries PBS for GPU information using pbsnodes
func (p *PbsScheduler) getGpuInfo() ([]GpuInfo, error) {
	cmd := exec.Command(p.pbsnodesBin, "-a")
	output, err := cmd.CombinedOutput()
	if err != nil {
		return nil, NewClusterError("PBS", "query GPU info", err)
	}

	gpuMap := make(map[string]*GpuInfo)

	// Parse pbsnodes output for GPU information
	lines := strings.Split(string(output), "\n")
	var currentNode string
	var nodeState string
	var nodeQueue string // Try to detect queue assignment if available

	for _, line := range lines {
		line = strings.TrimSpace(line)

		// Node name line
		if line != "" && !strings.Contains(line, "=") {
			currentNode = line
			nodeState = ""
			nodeQueue = "" // Reset queue for new node
			continue
		}

		if currentNode == "" {
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
			case "queue":
				// Some PBS installations show queue assignment
				nodeQueue = value
			case "resources_available.ngpus":
				var gpuCount int
				fmt.Sscanf(value, "%d", &gpuCount)
				if gpuCount > 0 {
					// Try to get GPU type from other fields
					gpuType := "gpu" // default

					// Build key with queue if available
					key := gpuType
					if nodeQueue != "" {
						key = fmt.Sprintf("%s:%s", nodeQueue, gpuType)
					}

					if existing, ok := gpuMap[key]; ok {
						existing.Total += gpuCount
						if nodeState == "free" {
							existing.Available += gpuCount
						}
					} else {
						available := 0
						if nodeState == "free" {
							available = gpuCount
						}
						gpuMap[key] = &GpuInfo{
							Type:      gpuType,
							Total:     gpuCount,
							Available: available,
							Partition: nodeQueue, // May be empty
						}
					}
				}
			case "resources_available.gpu_model", "gpu_type":
				// Some PBS installations expose GPU model/type
				// This would be used to set gpuType above
				// For now, we keep it simple with "gpu" as the type
			}
		}
	}

	gpus := make([]GpuInfo, 0, len(gpuMap))
	for _, info := range gpuMap {
		gpus = append(gpus, *info)
	}

	return gpus, nil
}

// getQueueLimits queries PBS for queue resource limits using qstat -Qf
func (p *PbsScheduler) getQueueLimits(gpuInfo []GpuInfo) ([]ResourceLimits, error) {
	cmd := exec.Command(p.qstatBin, "-Qf")
	output, err := cmd.CombinedOutput()
	if err != nil {
		return nil, NewClusterError("PBS", "query queues", err)
	}

	queueLimits := make(map[string]*ResourceLimits)
	var currentQueue string

	lines := strings.Split(string(output), "\n")
	for _, line := range lines {
		// Queue name line starts with "Queue: "
		if strings.HasPrefix(line, "Queue: ") {
			currentQueue = strings.TrimSpace(strings.TrimPrefix(line, "Queue:"))
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

		// Parse key = value pairs (indented lines)
		line = strings.TrimSpace(line)
		if !strings.Contains(line, "=") {
			continue
		}

		parts := strings.SplitN(line, "=", 2)
		if len(parts) != 2 {
			continue
		}

		key := strings.TrimSpace(parts[0])
		value := strings.TrimSpace(parts[1])
		limit := queueLimits[currentQueue]

		switch key {
		case "resources_max.walltime", "max_walltime":
			if duration, err := utils.ParseHMSTime(value); err == nil && duration > 0 {
				limit.MaxTime = duration
			}
		case "resources_max.mem", "max_mem":
			if memKB, err := parsePbsMemory(value); err == nil && memKB > 0 {
				memMB := memKB / 1024
				limit.MaxMemMBPerNode = memMB
			}
		case "resources_max.ncpus", "max_ncpus":
			var cpus int
			if _, err := fmt.Sscanf(value, "%d", &cpus); err == nil && cpus > 0 {
				limit.MaxCpusPerNode = cpus
			}
		case "resources_available.nodect", "total_nodes":
			var nodes int
			if _, err := fmt.Sscanf(value, "%d", &nodes); err == nil && nodes > 0 {
				limit.MaxNodes = nodes
			}
		case "enabled":
			// Skip disabled queues
			if value == "False" {
				delete(queueLimits, currentQueue)
			}
		}
	}

	// Get available resources per queue from actual nodes
	if p.pbsnodesBin != "" {
		availRes, err := p.getAvailableResourcesByQueue()
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
func (p *PbsScheduler) getAvailableResourcesByQueue() (map[string]ResourceLimits, error) {
	cmd := exec.Command(p.pbsnodesBin, "-a")
	output, err := cmd.CombinedOutput()
	if err != nil {
		return nil, NewClusterError("PBS", "query available resources by queue", err)
	}

	// Track resources per queue
	resources := make(map[string]ResourceLimits)

	// Parse pbsnodes output
	lines := strings.Split(string(output), "\n")
	var currentNode string
	var nodeCpus int
	var nodeMemKB int64
	var nodeState string
	var nodeQueue string

	for _, line := range lines {
		line = strings.TrimSpace(line)

		// Node name line
		if line != "" && !strings.Contains(line, "=") {
			// Save previous node's data
			if currentNode != "" && nodeQueue != "" && (nodeState == "free" || nodeState == "job-exclusive" || strings.Contains(nodeState, "busy")) {
				nodeMB := nodeMemKB / 1024

				// Update max for this queue
				if existing, ok := resources[nodeQueue]; ok {
					if nodeCpus > existing.MaxCpusPerNode {
						existing.MaxCpusPerNode = nodeCpus
					}
					if nodeMB > existing.MaxMemMBPerNode {
						existing.MaxMemMBPerNode = nodeMB
					}
					resources[nodeQueue] = existing
				} else {
					resources[nodeQueue] = ResourceLimits{
						Partition:       nodeQueue,
						MaxCpusPerNode:  nodeCpus,
						MaxMemMBPerNode: nodeMB,
					}
				}
			}
			currentNode = line
			nodeCpus = 0
			nodeMemKB = 0
			nodeState = ""
			nodeQueue = ""
			continue
		}

		if currentNode == "" {
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
			case "queue":
				nodeQueue = value
			case "np", "pcpus":
				fmt.Sscanf(value, "%d", &nodeCpus)
			case "resources_available.ncpus":
				fmt.Sscanf(value, "%d", &nodeCpus)
			case "resources_available.mem":
				nodeMemKB, _ = parsePbsMemory(value)
			}
		}
	}

	// Don't forget the last node
	if currentNode != "" && nodeQueue != "" && (nodeState == "free" || nodeState == "job-exclusive" || strings.Contains(nodeState, "busy")) {
		nodeMB := nodeMemKB / 1024

		if existing, ok := resources[nodeQueue]; ok {
			if nodeCpus > existing.MaxCpusPerNode {
				existing.MaxCpusPerNode = nodeCpus
			}
			if nodeMB > existing.MaxMemMBPerNode {
				existing.MaxMemMBPerNode = nodeMB
			}
			resources[nodeQueue] = existing
		} else {
			resources[nodeQueue] = ResourceLimits{
				Partition:       nodeQueue,
				MaxCpusPerNode:  nodeCpus,
				MaxMemMBPerNode: nodeMB,
			}
		}
	}

	return resources, nil
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
// GetJobResources reads allocated resources from PBS environment variables.
// Fields with value 0 were not exposed by PBS.
func (s *PbsScheduler) GetJobResources() *ResourceSpec {
	if _, ok := os.LookupEnv("PBS_JOBID"); !ok {
		return nil
	}
	res := &ResourceSpec{}
	if v := getEnvInt("PBS_NUM_NODES"); v != nil {
		res.Nodes = *v
	}
	// PBS_NCPUS: total CPUs per node (some PBS variants)
	if v := getEnvInt("PBS_NCPUS"); v != nil {
		res.CpusPerTask = *v
	} else if v := getEnvInt("NCPUS"); v != nil {
		res.CpusPerTask = *v
	}
	// Derive TasksPerNode from total tasks / nodes
	ntasks := getEnvInt("PBS_NP")
	if ntasks == nil {
		ntasks = getEnvInt("PBS_TASKNUM")
	}
	if ntasks != nil && res.Nodes > 0 {
		res.TasksPerNode = *ntasks / res.Nodes
	}
	// PBS_VMEM is in bytes
	if vmem := getEnvInt64("PBS_VMEM"); vmem != nil {
		if mb := *vmem / (1024 * 1024); mb > 0 {
			res.MemPerNodeMB = mb
		}
	}
	if n := getCudaDeviceCount(); n != nil && *n > 0 {
		res.Gpu = &GpuSpec{Count: *n}
	}
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
