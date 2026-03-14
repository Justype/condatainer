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

// slurmBlacklistedFlags lists SLURM directives that describe CPU topology / task
// distribution that cannot be reliably translated. Any match forces passthrough mode.
var slurmBlacklistedFlags = []string{
	"--gpus-per-socket",
	"--sockets-per-node",
	"--cores-per-socket",
	"--threads-per-core",
	"--ntasks-per-socket",
	"--ntasks-per-core",
	"--distribution",
}

// SlurmScheduler implements the Scheduler interface for SLURM
type SlurmScheduler struct {
	sbatchBin         string
	sinfoCommand      string
	scontrolCommand   string
	directiveRe       *regexp.Regexp
	jobIDRe           *regexp.Regexp
	cachedClusterInfo *ClusterInfo
}

// NewSlurmScheduler creates a new SLURM scheduler instance using sbatch from PATH
func NewSlurmScheduler() (*SlurmScheduler, error) {
	return newSlurmSchedulerWithBinary("")
}

// NewSlurmSchedulerWithBinary creates a SLURM scheduler using an explicit sbatch path
func NewSlurmSchedulerWithBinary(sbatchBin string) (*SlurmScheduler, error) {
	return newSlurmSchedulerWithBinary(sbatchBin)
}

func newSlurmSchedulerWithBinary(sbatchBin string) (*SlurmScheduler, error) {
	binPath := sbatchBin
	if binPath == "" {
		var err error
		binPath, err = exec.LookPath("sbatch")
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

	sinfoCmd := siblingBin(binPath, "sinfo")
	scontrolCmd := siblingBin(binPath, "scontrol")

	return &SlurmScheduler{
		sbatchBin:       binPath,
		sinfoCommand:    sinfoCmd,
		scontrolCommand: scontrolCmd,
		directiveRe:     regexp.MustCompile(`^\s*#SBATCH\s+(.+)$`),
		jobIDRe:         regexp.MustCompile(`Submitted batch job (\d+)`),
	}, nil
}

// GetCurrentJobID returns the SLURM job ID of the currently running job, or "".
func (s *SlurmScheduler) GetCurrentJobID() string { return os.Getenv("SLURM_JOB_ID") }

// IsAvailable checks if the SLURM binary is present on this system.
func (s *SlurmScheduler) IsAvailable() bool {
	return s.sbatchBin != ""
}

// IsInsideJob returns true if the current process is running inside a SLURM job.
func (s *SlurmScheduler) IsInsideJob() bool {
	_, ok := os.LookupEnv("SLURM_JOB_ID")
	return ok
}

// GetInfo returns information about the SLURM scheduler
func (s *SlurmScheduler) GetInfo() *SchedulerInfo {
	info := &SchedulerInfo{
		Type:      "SLURM",
		Binary:    s.sbatchBin,
		InJob:     s.IsInsideJob(),
		Available: s.IsAvailable(),
	}

	// Try to get SLURM version
	if s.sbatchBin != "" {
		if version, err := s.getSlurmVersion(); err == nil {
			info.Version = version
		}
	}

	return info
}

// getSlurmVersion attempts to get the SLURM version
func (s *SlurmScheduler) getSlurmVersion() (string, error) {
	output, err := runCommand("SLURM", "get-version", s.sbatchBin, "--version")
	if err != nil {
		return "", err
	}

	// Parse version from output like "slurm 23.02.6"
	versionStr := strings.TrimSpace(string(output))
	parts := strings.Fields(versionStr)
	if len(parts) >= 2 {
		return parts[1], nil
	}

	return versionStr, nil
}

// ReadScriptSpecs parses #SBATCH directives from a build script
func (s *SlurmScheduler) ReadScriptSpecs(scriptPath string) (*ScriptSpecs, error) {
	lines, err := readFileLines(scriptPath)
	if err != nil {
		return nil, err
	}
	return parseScript(scriptPath, lines, s.extractDirectives, s.parseRuntimeConfig, s.parseResourceSpec)
}

// extractDirectives extracts raw directive strings from script lines (strips the #SBATCH prefix).
func (s *SlurmScheduler) extractDirectives(lines []string) []string {
	var out []string
	for _, line := range lines {
		if m := s.directiveRe.FindStringSubmatch(line); m != nil {
			out = append(out, utils.StripInlineComment(m[1]))
		}
	}
	return out
}

// parseRuntimeConfig consumes job control directives (name, I/O, email) from the directive list.
// Returns the populated RuntimeConfig, unconsumed directives, and any critical error.
func (s *SlurmScheduler) parseRuntimeConfig(directives []string) (RuntimeConfig, []string, error) {
	var rc RuntimeConfig
	var unconsumed []string

	for _, flag := range directives {
		consumed := true
		switch {
		case flagMatches(flag, "--job-name", "-J"):
			rc.JobName, _ = flagValue(flag, "--job-name", "-J")
		case flagMatches(flag, "--chdir", "-D"):
			v, _ := flagValue(flag, "--chdir", "-D")
			rc.WorkDir = absPath(v)
		case flagMatches(flag, "--output", "-o"):
			rc.Stdout, _ = flagValue(flag, "--output", "-o")
		case flagMatches(flag, "--error", "-e"):
			rc.Stderr, _ = flagValue(flag, "--error", "-e")
		case flagMatches(flag, "--partition", "-p"):
			rc.Partition, _ = flagValue(flag, "--partition", "-p")
		case flagMatches(flag, "--mail-user"):
			rc.MailUser, _ = flagValue(flag, "--mail-user")
		case flagMatches(flag, "--mail-type"):
			rawMailType, _ := flagValue(flag, "--mail-type")
			mailType := strings.ToUpper(rawMailType)
			if mailType == "NONE" {
				rc.EmailOnBegin = false
				rc.EmailOnEnd = false
				rc.EmailOnFail = false
			} else {
				for _, t := range strings.Split(mailType, ",") {
					switch strings.TrimSpace(t) {
					case "ALL":
						rc.EmailOnBegin = true
						rc.EmailOnEnd = true
						rc.EmailOnFail = true
					case "BEGIN":
						rc.EmailOnBegin = true
					case "END":
						rc.EmailOnEnd = true
					case "FAIL", "REQUEUE", "INVALID_DEPEND":
						rc.EmailOnFail = true
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

// parseResourceSpec consumes compute resource directives from the directive list.
// Returns a populated ResourceSpec and any unconsumed directives.
// Returns nil ResourceSpec on blacklisted flags or any unresolvable interdependency.
//
// Two-phase design:
//
//	Phase 1: scan all directives into typed temp vars (no writes to ResourceSpec yet).
//	Phase 2: resolve in dependency order — Nodes/Tasks → GPU → CPU → Memory → Time.
func (s *SlurmScheduler) parseResourceSpec(directives []string) (*ResourceSpec, []string) {
	defaults := GetSpecDefaults()

	// ── Phase 1: Scan ────────────────────────────────────────────────────────
	var (
		nodes            int = defaults.Nodes
		ntasksPerNode    int
		totalNtasks      int
		hasNtasksPerNode bool
		hasTotalNtasks   bool

		cpusPerTask    int = defaults.CpusPerTask
		hasCpusPerTask bool
		cpusPerGpu     int
		hasCpusPerGpu  bool

		// GPU: at most one form expected; last one seen wins in Phase 1.
		// Resolution priority in Phase 2: gres > per-node > total > per-task.
		gpuGresStr      string // --gres=gpu:... (per-node by SLURM definition)
		gpuPerNodeStr   string // --gpus-per-node=...
		gpuTotalStr     string // --gpus=... (total across all nodes)
		gpuPerTaskStr   string // --gpus-per-task=...
		ntasksPerGpu    int
		hasNtasksPerGpu bool

		memStr       string // --mem
		memPerCpuStr string // --mem-per-cpu (resolved after CPU)
		memPerGpuStr string // --mem-per-gpu (resolved after GPU)
		timeStr      string

		exclusive bool
	)

	var unconsumed []string

	for _, flag := range directives {
		// Blacklist: topology flags we cannot translate → passthrough immediately.
		for _, bl := range slurmBlacklistedFlags {
			if flag == bl || flagMatches(flag, bl) {
				logParseWarning("SLURM: unsupported topology directive %q; using passthrough mode", flag)
				return nil, directives
			}
		}

		consumed := true
		var parseErr error

		switch {
		case flagMatches(flag, "--nodes", "-N"):
			_, parseErr = flagScanInt(flag, &nodes, "--nodes", "-N")
		case flagMatches(flag, "--ntasks-per-node"):
			_, parseErr = flagScanInt(flag, &ntasksPerNode, "--ntasks-per-node")
			hasNtasksPerNode = true
		case flagMatches(flag, "--ntasks", "-n"):
			_, parseErr = flagScanInt(flag, &totalNtasks, "--ntasks", "-n")
			hasTotalNtasks = true
		case flagMatches(flag, "--cpus-per-task", "-c"):
			_, parseErr = flagScanInt(flag, &cpusPerTask, "--cpus-per-task", "-c")
			hasCpusPerTask = true
		case flagMatches(flag, "--cpus-per-gpu"):
			_, parseErr = flagScanInt(flag, &cpusPerGpu, "--cpus-per-gpu")
			hasCpusPerGpu = true
		case strings.HasPrefix(flag, "--gres=gpu:"):
			gpuGresStr = strings.TrimPrefix(flag, "--gres=")
		case flagMatches(flag, "--gpus"):
			gpuTotalStr, _ = flagValue(flag, "--gpus")
		case flagMatches(flag, "--gpus-per-node"):
			gpuPerNodeStr, _ = flagValue(flag, "--gpus-per-node")
		case flagMatches(flag, "--gpus-per-task"):
			gpuPerTaskStr, _ = flagValue(flag, "--gpus-per-task")
		case flagMatches(flag, "--ntasks-per-gpu"):
			_, parseErr = flagScanInt(flag, &ntasksPerGpu, "--ntasks-per-gpu")
			hasNtasksPerGpu = true
		case flagMatches(flag, "--mem"):
			memStr, _ = flagValue(flag, "--mem")
		case flagMatches(flag, "--mem-per-cpu"):
			memPerCpuStr, _ = flagValue(flag, "--mem-per-cpu")
		case flagMatches(flag, "--mem-per-gpu"):
			memPerGpuStr, _ = flagValue(flag, "--mem-per-gpu")
		case flagMatches(flag, "--time", "-t"):
			timeStr, _ = flagValue(flag, "--time", "-t")
		case flag == "--exclusive", strings.HasPrefix(flag, "--exclusive="):
			exclusive = true
		default:
			consumed = false
		}

		if parseErr != nil {
			logParseWarning("SLURM: failed to parse directive %q: %v", flag, parseErr)
			return nil, directives
		}
		if !consumed {
			unconsumed = append(unconsumed, flag)
		}
	}

	// ── Phase 2: Resolve in dependency order ─────────────────────────────────
	rs := &ResourceSpec{
		CpusPerTask:  defaults.CpusPerTask,
		TasksPerNode: defaults.TasksPerNode,
		Nodes:        defaults.Nodes,
		MemPerNodeMB: defaults.MemPerNodeMB,
		Time:         defaults.Time,
		Exclusive:    exclusive,
	}

	// Conflict Check
	if gpuPerTaskStr != "" && hasNtasksPerGpu {
		logParseWarning("SLURM: --gpus-per-task and --ntasks-per-gpu are mutually exclusive; using passthrough")
		return nil, directives
	}

	// 1. Nodes
	rs.Nodes = nodes

	// 2. GPU Resolution
	switch {
	case gpuGresStr != "":
		gpu, err := parseSlurmGpu(gpuGresStr)
		if err != nil {
			logParseWarning("SLURM: failed to parse --gres value %q: %v; using passthrough", gpuGresStr, err)
			return nil, directives
		}
		if rs.Nodes <= 0 {
			rs.Nodes = 1
		}
		rs.Gpu = gpu
	case gpuPerNodeStr != "":
		gpu, err := parseSlurmGpu(gpuPerNodeStr)
		if err != nil {
			logParseWarning("SLURM: failed to parse --gpus-per-node value %q: %v; using passthrough", gpuPerNodeStr, err)
			return nil, directives
		}
		if rs.Nodes <= 0 {
			rs.Nodes = 1
		}
		rs.Gpu = gpu
	case gpuTotalStr != "":
		gpu, err := parseSlurmGpu(gpuTotalStr)
		if err != nil {
			logParseWarning("SLURM: failed to parse --gpus value %q: %v; using passthrough", gpuTotalStr, err)
			return nil, directives
		}
		if rs.Nodes <= 0 {
			// If total gpus provided but no nodes, default to that many nodes with 1 gpu each
			rs.Nodes = gpu.Count
			gpu.Count = 1
		} else {
			if gpu.Count%rs.Nodes != 0 {
				logParseWarning("SLURM: --gpus=%d not divisible by --nodes=%d; using passthrough", gpu.Count, rs.Nodes)
				return nil, directives
			}
			gpu.Count /= rs.Nodes
		}
		rs.Gpu = gpu
	}

	// 3. Tasks
	// We want to calculate explicitly how many Ntasks we have.
	if hasTotalNtasks {
		rs.Ntasks = totalNtasks
		if hasNtasksPerNode {
			rs.TasksPerNode = ntasksPerNode
		} else if rs.Nodes > 0 {
			if totalNtasks%rs.Nodes == 0 {
				rs.TasksPerNode = totalNtasks / rs.Nodes
			} else if rs.Gpu != nil || gpuPerTaskStr != "" {
				logParseWarning("SLURM: --ntasks=%d not evenly divisible by --nodes=%d with GPU specs; using passthrough mode", totalNtasks, rs.Nodes)
				return nil, directives
			} else {
				// Non-divisible: use ceiling to ensure all tasks fit; this matches memory allocation logic.
				rs.TasksPerNode = (totalNtasks + rs.Nodes - 1) / rs.Nodes
			}
		} else {
			rs.TasksPerNode = 0 // no node constraint: free distribution
		}
	} else if hasNtasksPerGpu {
		if rs.Gpu != nil {
			rs.Ntasks = ntasksPerGpu * (rs.Nodes * rs.Gpu.Count)
		} else {
			logParseWarning("SLURM: --ntasks-per-gpu requires a GPU spec; using passthrough")
			return nil, directives
		}
	} else if hasNtasksPerNode {
		rs.TasksPerNode = ntasksPerNode
	}

	// Back-calculate GPUs if --gpus-per-task was used
	if gpuPerTaskStr != "" {
		gpu, err := parseSlurmGpu(gpuPerTaskStr)
		if err != nil {
			logParseWarning("SLURM: failed to parse --gpus-per-task value %q: %v; using passthrough mode", gpuPerTaskStr, err)
			return nil, directives
		}
		if rs.Nodes > 0 {
			// e.g. nodes=2, total tasks=8 => tasksPerNode=4.
			// gpusPerNode = gpusPerTask * tasksPerNode
			tasksUsed := rs.TasksPerNode
			if tasksUsed == 0 && rs.Ntasks > 0 {
				tasksUsed = rs.Ntasks / rs.Nodes
			}
			if tasksUsed == 0 {
				tasksUsed = 1 // Safe fallback
			}
			gpu.Count *= tasksUsed
		}
		rs.Gpu = gpu
	}

	// 4. CPU (depends on task distribution or GPUs)
	if hasCpusPerTask {
		rs.CpusPerTask = cpusPerTask
	} else if hasCpusPerGpu {
		if rs.Gpu == nil {
			logParseWarning("SLURM: --cpus-per-gpu requires a GPU spec; using passthrough mode")
			return nil, directives
		}
		rs.CpusPerTask = cpusPerGpu * rs.Gpu.Count
		// If MPI tasks spread things out, we need to divide CPUs evenly among tasks
		// on that node. If TasksPerNode is defined, we divide.
		tasksUsed := rs.TasksPerNode
		if tasksUsed == 0 && rs.Nodes > 0 && rs.Ntasks > 0 {
			tasksUsed = rs.Ntasks / rs.Nodes
		}
		if tasksUsed > 0 {
			rs.CpusPerTask /= tasksUsed
		}
		if rs.CpusPerTask < 1 {
			rs.CpusPerTask = 1
		}
	} else {
		// New logic mandates fallback to 1 if not set
		if rs.CpusPerTask <= 0 {
			rs.CpusPerTask = 1
		}
	}

	// 5. Memory — three forms, resolved after CPU and GPU.
	switch {
	case memStr != "":
		mem, err := parseMemoryMB(memStr)
		if err != nil {
			logParseWarning("SLURM: invalid --mem value %q: %v; using passthrough mode", memStr, err)
			return nil, directives
		}
		rs.MemPerNodeMB = mem
	case memPerCpuStr != "":
		memPerCpu, err := parseMemoryMB(memPerCpuStr)
		if err != nil {
			logParseWarning("SLURM: invalid --mem-per-cpu value %q: %v; using passthrough mode", memPerCpuStr, err)
			return nil, directives
		}
		rs.MemPerCpuMB = memPerCpu
		rs.MemPerNodeMB = 0 // clear default so GetMemPerNodeMB() derives correctly
	case memPerGpuStr != "":
		if rs.Gpu == nil {
			logParseWarning("SLURM: --mem-per-gpu requires a GPU spec; using passthrough mode")
			return nil, directives
		}
		memPerGpu, err := parseMemoryMB(memPerGpuStr)
		if err != nil {
			logParseWarning("SLURM: invalid --mem-per-gpu value %q: %v; using passthrough mode", memPerGpuStr, err)
			return nil, directives
		}
		rs.MemPerNodeMB = memPerGpu * int64(rs.Gpu.Count)
	}

	// 6. Time (independent).
	if timeStr != "" {
		t, err := utils.ParseDHMSTime(timeStr)
		if err != nil {
			logParseWarning("SLURM: invalid --time value %q: %v; using passthrough mode", timeStr, err)
			return nil, directives
		}
		rs.Time = t
	}

	return rs, unconsumed
}

// CreateScriptWithSpec generates a SLURM batch script
// Returns the script path or any critical error.
func (s *SlurmScheduler) CreateScriptWithSpec(jobSpec *JobSpec, outputDir string) (string, error) {
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
	scriptName := "job.sbatch"
	if jobSpec.Name != "" {
		scriptName = fmt.Sprintf("%s.sbatch", safeJobName(jobSpec.Name))
	}

	scriptPath := filepath.Join(outputDir, scriptName)

	// Create the batch script
	file, err := utils.CreateFileWritable(scriptPath)
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
		fmt.Fprintf(writer, "#SBATCH %s\n", flag)
	}

	// Write RuntimeConfig directives
	ctrl := specs.Control
	if ctrl.JobName != "" {
		fmt.Fprintf(writer, "#SBATCH --job-name=%s\n", ctrl.JobName)
	}
	if ctrl.WorkDir != "" {
		fmt.Fprintf(writer, "#SBATCH --chdir=%s\n", ctrl.WorkDir)
	}
	if ctrl.Stdout != "" {
		fmt.Fprintf(writer, "#SBATCH --output=%s\n", ctrl.AbsStdout())
	}
	if ctrl.Stderr != "" {
		fmt.Fprintf(writer, "#SBATCH --error=%s\n", ctrl.AbsStderr())
	}
	if ctrl.EmailOnBegin || ctrl.EmailOnEnd || ctrl.EmailOnFail {
		var mailTypes []string
		if ctrl.EmailOnBegin {
			mailTypes = append(mailTypes, "BEGIN")
		}
		if ctrl.EmailOnEnd {
			mailTypes = append(mailTypes, "END")
		}
		if ctrl.EmailOnFail {
			mailTypes = append(mailTypes, "FAIL")
		}
		fmt.Fprintf(writer, "#SBATCH --mail-type=%s\n", strings.Join(mailTypes, ","))
	}
	if ctrl.MailUser != "" {
		fmt.Fprintf(writer, "#SBATCH --mail-user=%s\n", ctrl.MailUser)
	}
	if ctrl.Partition != "" {
		fmt.Fprintf(writer, "#SBATCH --partition=%s\n", ctrl.Partition)
	}

	// Write ResourceSpec directives (only if not in passthrough mode)
	if specs.Spec != nil {
		rs := specs.Spec

		if !rs.IsMPI() {
			// Pure OpenMP: pin to 1 node and 1 task so --mem is unambiguous.
			fmt.Fprintf(writer, "#SBATCH --nodes=1\n")
			fmt.Fprintf(writer, "#SBATCH --ntasks=1\n")
			if rs.CpusPerTask > 0 {
				fmt.Fprintf(writer, "#SBATCH --cpus-per-task=%d\n", rs.CpusPerTask)
			}
		} else {
			fmt.Fprintf(writer, "#SBATCH --ntasks=%d\n", rs.GetNtasks())
			if rs.Nodes > 0 {
				fmt.Fprintf(writer, "#SBATCH --nodes=%d\n", rs.Nodes)
			}
			// Calculate and enforce TasksPerNode when both Nodes and Ntasks are set
			tasksPerNode := rs.TasksPerNode
			if tasksPerNode == 0 && rs.Nodes > 0 && rs.GetNtasks() > 0 {
				// Use ceiling to ensure all tasks fit
				tasksPerNode = (rs.GetNtasks() + rs.Nodes - 1) / rs.Nodes
			}
			if tasksPerNode > 0 {
				fmt.Fprintf(writer, "#SBATCH --ntasks-per-node=%d\n", tasksPerNode)
			}
			if rs.CpusPerTask > 0 {
				fmt.Fprintf(writer, "#SBATCH --cpus-per-task=%d\n", rs.CpusPerTask)
			}
		}

		if rs.MemPerCpuMB > 0 {
			fmt.Fprintf(writer, "#SBATCH --mem-per-cpu=%dmb\n", rs.MemPerCpuMB)
		} else if rs.MemPerNodeMB > 0 {
			fmt.Fprintf(writer, "#SBATCH --mem=%dmb\n", rs.MemPerNodeMB)
		}
		if rs.Time > 0 {
			fmt.Fprintf(writer, "#SBATCH --time=%s\n", formatSlurmTimeSpec(rs.Time))
		}
		if rs.Gpu != nil && rs.Gpu.Count > 0 {
			if rs.Gpu.Type != "" && rs.Gpu.Type != "gpu" {
				fmt.Fprintf(writer, "#SBATCH --gpus-per-node=%s:%d\n", rs.Gpu.Type, rs.Gpu.Count)
			} else {
				fmt.Fprintf(writer, "#SBATCH --gpus-per-node=%d\n", rs.Gpu.Count)
			}
		}
		if rs.Exclusive {
			fmt.Fprintln(writer, "#SBATCH --exclusive")
		}
	}

	// Array directive
	if jobSpec.Array != nil {
		arr := jobSpec.Array
		r := fmt.Sprintf("1-%d", arr.Count)
		if arr.Limit > 0 {
			r += fmt.Sprintf("%%%d", arr.Limit)
		}
		fmt.Fprintf(writer, "#SBATCH --array=%s\n", r)
	}

	fmt.Fprintln(writer, "")

	// Array job: extract input line, set ARRAY_ARGS, and redirect output
	if jobSpec.Array != nil {
		writeArrayBlock(writer, "$SLURM_ARRAY_TASK_ID",
			jobSpec.Array.InputFile, outputDir, safeJobName(jobSpec.Name),
			jobSpec.Array.Count, arraySeparateOutput)
		jobSpec.Metadata["Array Job ID"] = "$SLURM_ARRAY_JOB_ID"
		jobSpec.Metadata["Array Index"] = "$SLURM_ARRAY_TASK_ID"
		jobSpec.Metadata["Array File"] = jobSpec.Array.InputFile
		jobSpec.Metadata["Array Args"] = "$ARRAY_ARGS"
	}

	// Print job information at start
	jobIDVar := "$SLURM_JOB_ID"
	if jobSpec.Array != nil {
		jobIDVar = "$SLURM_ARRAY_JOB_ID"
	}
	writeJobHeader(writer, jobIDVar, specs, formatSlurmTimeSpec, jobSpec.Metadata)
	fmt.Fprintln(writer, "")

	// Write the command and capture exit code
	fmt.Fprintln(writer, jobSpec.Command)
	fmt.Fprintln(writer, "_EXIT_CODE=$?")

	// Print completion info
	fmt.Fprintln(writer, "")
	writeJobFooter(writer, jobIDVar)

	// Self-dispose: remove this script file (unless in debug mode)
	if !debugMode {
		fmt.Fprintf(writer, "rm -f %s\n", scriptPath)
	}

	// Exit with command's exit code
	fmt.Fprintln(writer, "exit $_EXIT_CODE")

	// Make executable
	if err := os.Chmod(scriptPath, utils.PermExec); err != nil {
		return "", NewScriptCreationError(jobSpec.Name, scriptPath, err)
	}

	return scriptPath, nil
}

// Submit submits a SLURM job with optional dependency chain
// buildSlurmDepFlag returns the --dependency flag string for sbatch, or "" if deps is empty.
// Format: --dependency=afterok:ID1:ID2,afternotok:ID3,afterany:ID4
func buildSlurmDepFlag(deps []Dependency) string {
	var parts []string
	for _, dep := range deps {
		if len(dep.JobIDs) > 0 {
			parts = append(parts, dep.Type+":"+strings.Join(dep.JobIDs, ":"))
		}
	}
	if len(parts) == 0 {
		return ""
	}
	return "--dependency=" + strings.Join(parts, ",")
}

// buildSlurmSubmitArgs returns the sbatch argument list for the given deps and script path.
func buildSlurmSubmitArgs(deps []Dependency, scriptPath string) []string {
	args := []string{scriptPath}
	if flag := buildSlurmDepFlag(deps); flag != "" {
		args = append([]string{flag, "--kill-on-invalid-dep=yes"}, args...)
	}
	return args
}

func (s *SlurmScheduler) Submit(scriptPath string, deps []Dependency) (string, error) {
	args := buildSlurmSubmitArgs(deps, scriptPath)

	// Execute sbatch
	output, err := runCommand("SLURM", "submit", s.sbatchBin, args...)
	if err != nil {
		return "", NewSubmissionError("SLURM", filepath.Base(scriptPath), string(output), err)
	}

	// Parse job ID from output
	matches := s.jobIDRe.FindStringSubmatch(string(output))
	if len(matches) < 2 {
		return "", fmt.Errorf("%w: %s", ErrJobIDParseFailed, string(output))
	}

	jobID := matches[1]
	return jobID, nil
}

// GetClusterInfo retrieves SLURM cluster configuration
func (s *SlurmScheduler) GetClusterInfo() (*ClusterInfo, error) {
	if s.cachedClusterInfo != nil {
		return s.cachedClusterInfo, nil
	}

	info := &ClusterInfo{
		AvailableGpus: make([]GpuInfo, 0),
		Limits:        make([]ResourceLimits, 0),
	}

	// Get GPU information
	if s.sinfoCommand != "" {
		gpus, err := s.getGpuInfo()
		if err == nil {
			info.AvailableGpus = gpus
		}

		// Get max node resources (CPUs and memory)
		maxCpus, maxMem, err := s.getMaxNodeResources()
		if err == nil {
			info.MaxCpusPerNode = maxCpus
			info.MaxMemMBPerNode = maxMem
		}
	}

	// Get partition limits (passing GPU info for max GPU calculation)
	if s.scontrolCommand != "" {
		limits, err := s.getPartitionLimits(info.AvailableGpus)
		if err == nil {
			info.Limits = limits
		}
	}

	s.cachedClusterInfo = info
	return info, nil
}

// getGpuInfo queries SLURM for available GPU types
func (s *SlurmScheduler) getGpuInfo() ([]GpuInfo, error) {
	output, err := runCommand("SLURM", "query-gpus", s.sinfoCommand, "-o", "%P|%G|%D|%T", "--noheader")
	if err != nil {
		return nil, NewClusterError("SLURM", "query GPUs", err)
	}

	gpuMap := make(map[string]*GpuInfo)
	lines := strings.Split(strings.TrimSpace(string(output)), "\n")

	for _, line := range lines {
		parts := strings.Split(line, "|")
		if len(parts) < 4 {
			continue
		}

		partition := strings.TrimSpace(strings.TrimSuffix(parts[0], "*"))
		gresStr := strings.TrimSpace(parts[1])
		nodesStr := strings.TrimSpace(parts[2])
		state := strings.TrimSpace(parts[3])

		// Parse multiple GPU types separated by commas
		// Example: gpu:nvidia_h100_80gb_hbm3_3g.40gb:8(S:0-3),gpu:nvidia_h100_80gb_hbm3_2g.20gb:8(S:0-3)
		gpuEntries := strings.Split(gresStr, ",")

		for _, entry := range gpuEntries {
			if !strings.HasPrefix(entry, "gpu:") && !strings.HasPrefix(entry, "gpu(") {
				continue
			}

			// Remove "gpu:" or "gpu(" prefix and any trailing ")"
			entry = strings.TrimPrefix(entry, "gpu:")
			entry = strings.TrimPrefix(entry, "gpu(")

			// Remove socket information like (S:0-3)
			if idx := strings.Index(entry, "("); idx > 0 {
				entry = entry[:idx]
			}
			entry = strings.TrimSuffix(entry, ")")

			gpuType := "gpu"
			gpuCount := 0

			if strings.Contains(entry, ":") {
				entryParts := strings.Split(entry, ":")
				if len(entryParts) >= 1 {
					gpuType = entryParts[0]
				}
				if len(entryParts) >= 2 {
					// The last part should be the count
					fmt.Sscanf(entryParts[len(entryParts)-1], "%d", &gpuCount)
				}
			} else {
				fmt.Sscanf(entry, "%d", &gpuCount)
			}

			if gpuCount == 0 {
				continue
			}

			nodes := 0
			fmt.Sscanf(nodesStr, "%d", &nodes)
			totalGpus := gpuCount * nodes

			key := fmt.Sprintf("%s:%s", partition, gpuType)
			if existing, ok := gpuMap[key]; ok {
				existing.Total += totalGpus
				if state == "idle" || state == "mixed" {
					existing.Available += totalGpus
				}
			} else {
				available := 0
				if state == "idle" || state == "mixed" {
					available = totalGpus
				}
				gpuMap[key] = &GpuInfo{
					Type:      gpuType,
					Total:     totalGpus,
					Available: available,
					Partition: partition,
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

// getPartitionLimits queries SLURM for partition resource limits
func (s *SlurmScheduler) getPartitionLimits(gpuInfo []GpuInfo) ([]ResourceLimits, error) {
	// First, get partition config limits
	output, err := runCommand("SLURM", "query-partition-limits", s.scontrolCommand, "show", "partition", "-o")
	if err != nil {
		return nil, NewClusterError("SLURM", "query partition limits", err)
	}

	limits := make([]ResourceLimits, 0)
	lines := strings.Split(strings.TrimSpace(string(output)), "\n")

	for _, line := range lines {
		limit := s.parsePartitionLine(line)
		if limit != nil {
			limits = append(limits, *limit)
		}
	}

	// Get available resources per partition from actual nodes
	if s.sinfoCommand != "" {
		availRes, err := s.getAvailableResourcesByPartition()
		if err == nil {
			// Merge available resources into limits
			for i := range limits {
				if avail, ok := availRes[limits[i].Partition]; ok {
					// Use available resources as limits if they're set
					if avail.MaxCpusPerNode > 0 {
						limits[i].MaxCpusPerNode = avail.MaxCpusPerNode
					}
					if avail.MaxMemMBPerNode > 0 {
						limits[i].MaxMemMBPerNode = avail.MaxMemMBPerNode
					}
				}
			}
		}
	}

	// Calculate max GPUs per partition from GPU info
	gpusByPartition := make(map[string]int)
	for _, gpu := range gpuInfo {
		partition := gpu.Partition
		if partition == "" {
			partition = "default"
		}
		gpusByPartition[partition] += gpu.Total
	}

	// Add max GPUs to limits
	for i := range limits {
		if maxGpus, ok := gpusByPartition[limits[i].Partition]; ok {
			limits[i].MaxGpus = maxGpus
		}
	}

	return limits, nil
}

// getAvailableResourcesByPartition queries available CPUs and memory per partition
func (s *SlurmScheduler) getAvailableResourcesByPartition() (map[string]ResourceLimits, error) {
	// Query sinfo for partition, CPUs, and memory: %R = partition, %c = CPUs, %m = memory (MB)
	output, err := runCommand("SLURM", "query-available-resources", s.sinfoCommand, "-o", "%R|%c|%m", "--noheader")
	if err != nil {
		return nil, NewClusterError("SLURM", "query available resources", err)
	}

	// Track max resources per partition
	resources := make(map[string]ResourceLimits)

	lines := strings.Split(strings.TrimSpace(string(output)), "\n")
	for _, line := range lines {
		parts := strings.Split(line, "|")
		if len(parts) < 3 {
			continue
		}

		partition := strings.TrimSpace(strings.TrimSuffix(parts[0], "*"))
		var cpus int
		var memMB int64
		fmt.Sscanf(strings.TrimSpace(parts[1]), "%d", &cpus)
		memMB, _ = parseMemoryMB(strings.TrimSpace(parts[2]))

		// Update max for this partition
		if existing, ok := resources[partition]; ok {
			if cpus > existing.MaxCpusPerNode {
				existing.MaxCpusPerNode = cpus
			}
			if memMB > existing.MaxMemMBPerNode {
				existing.MaxMemMBPerNode = memMB
			}
			resources[partition] = existing
		} else {
			resources[partition] = ResourceLimits{
				Partition:       partition,
				MaxCpusPerNode:  cpus,
				MaxMemMBPerNode: memMB,
			}
		}
	}

	return resources, nil
}

// getMaxNodeResources queries SLURM for maximum CPU and memory available per node
func (s *SlurmScheduler) getMaxNodeResources() (int, int64, error) {
	// Query sinfo for CPUs and memory per node: %c = CPUs, %m = memory in MB
	output, err := runCommand("SLURM", "query-node-resources", s.sinfoCommand, "-o", "%c|%m", "--noheader")
	if err != nil {
		return 0, 0, NewClusterError("SLURM", "query node resources", err)
	}

	var maxCpus int
	var maxMemMB int64

	lines := strings.Split(strings.TrimSpace(string(output)), "\n")
	for _, line := range lines {
		parts := strings.Split(line, "|")
		if len(parts) < 2 {
			continue
		}

		var cpus int
		var memMB int64
		fmt.Sscanf(strings.TrimSpace(parts[0]), "%d", &cpus)
		memMB, _ = parseMemoryMB(strings.TrimSpace(parts[1]))

		if cpus > maxCpus {
			maxCpus = cpus
		}
		if memMB > maxMemMB {
			maxMemMB = memMB
		}
	}

	return maxCpus, maxMemMB, nil
}

// parsePartitionLine parses a single partition line from scontrol output
func (s *SlurmScheduler) parsePartitionLine(line string) *ResourceLimits {
	// https://slurm.schedmd.com/slurm.conf.html#SECTION_PARTITION-CONFIGURATION
	limit := &ResourceLimits{}

	fields := strings.Fields(line)
	for _, field := range fields {
		kv := strings.SplitN(field, "=", 2)
		if len(kv) != 2 {
			continue
		}

		key := kv[0]
		value := kv[1]

		switch key {
		case "PartitionName":
			limit.Partition = value
		case "MaxTime":
			if value != "UNLIMITED" {
				if dur, err := utils.ParseDHMSTime(value); err == nil {
					limit.MaxTime = dur
				}
			}
		case "DefaultTime":
			if value != "NONE" && value != "UNLIMITED" {
				if dur, err := utils.ParseDHMSTime(value); err == nil {
					limit.DefaultTime = dur
				}
			}
		case "MaxCPUsPerNode":
			fmt.Sscanf(value, "%d", &limit.MaxCpusPerNode)
		case "MaxMemPerNode":
			if value != "UNLIMITED" {
				limit.MaxMemMBPerNode, _ = parseMemoryMB(value)
			}
		case "MaxNodes":
			if value != "UNLIMITED" {
				fmt.Sscanf(value, "%d", &limit.MaxNodes)
			}
		}
	}

	if limit.Partition == "" {
		return nil
	}

	return limit
}

// parseSlurmGpu parses SLURM GPU specifications
// Supports formats like:
// - gpu:2 (2 generic GPUs)
// - a100:1 (1 A100 GPU)
// - gpu:a100:2 (2 A100 GPUs)
// - nvidia_h100_80gb_hbm3_1g.10gb (MIG profile)
// - nvidia_h100_80gb_hbm3_1g.10gb:2 (2 MIG instances)
func parseSlurmGpu(gpuStr string) (*GpuSpec, error) {
	if gpuStr == "" {
		return nil, nil
	}

	spec := &GpuSpec{
		Raw:   gpuStr,
		Count: 1,
	}

	// Check if this is a MIG profile (contains dots and underscores)
	// MIG format: nvidia_h100_80gb_hbm3_1g.10gb or nvidia_h100_80gb_hbm3_1g.10gb:2
	if strings.Contains(gpuStr, ".") && strings.Contains(gpuStr, "_") {
		// This looks like a MIG profile
		parts := strings.Split(gpuStr, ":")
		if len(parts) == 1 {
			// Just the MIG profile name
			spec.Type = parts[0]
			spec.Count = 1
		} else if len(parts) == 2 {
			// MIG profile with count
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

func formatSlurmTimeSpec(d time.Duration) string {
	if d <= 0 {
		return ""
	}
	days := int64(d.Hours()) / 24
	if days > 0 {
		rem := d - time.Duration(days*24)*time.Hour
		return fmt.Sprintf("%d-%s", days, formatHMSTime(rem))
	}
	return formatHMSTime(d)
}

// GetJobResources reads allocated resources from SLURM environment variables.
// Fields with value 0 were not exposed by SLURM.
func (s *SlurmScheduler) GetJobResources() *ResourceSpec {
	if _, ok := os.LookupEnv("SLURM_JOB_ID"); !ok {
		return nil
	}
	res := &ResourceSpec{}
	if v := getEnvInt("SLURM_JOB_NUM_NODES"); v != nil {
		res.Nodes = *v
	}

	// Prefer MPI library env vars (global across all ranks) over SLURM vars (task-local in MPI jobs).
	if mpiSize := getMpiCommSize(); mpiSize != nil {
		res.Ntasks = *mpiSize
	} else if v := getEnvInt("SLURM_NTASKS"); v != nil {
		res.Ntasks = *v
	}

	// TasksPerNode: prefer MPI local size (actual tasks on this node) over SLURM var
	if mpiLocal := getMpiLocalSize(); mpiLocal != nil {
		res.TasksPerNode = *mpiLocal
	} else if v := getEnvInt("SLURM_NTASKS_PER_NODE"); v != nil {
		res.TasksPerNode = *v
	} else if res.Ntasks > 0 && res.Nodes > 1 && res.Ntasks%res.Nodes == 0 {
		// Derive TasksPerNode when SLURM_NTASKS_PER_NODE is absent but geometry is uniform.
		res.TasksPerNode = res.Ntasks / res.Nodes
	}
	if v := getEnvInt("SLURM_CPUS_PER_TASK"); v != nil {
		res.CpusPerTask = *v
	}
	// SLURM_MEM_PER_NODE is already in MB.
	if v := getEnvInt64("SLURM_MEM_PER_NODE"); v != nil {
		res.MemPerNodeMB = *v
	}
	if n := getCudaDeviceCount(); n != nil && *n > 0 {
		res.Gpu = &GpuSpec{Count: *n}
	}
	return res
}

// GetJobStatus returns the current status of the given SLURM job ID.
// Uses squeue --format=%T; returns JobStatusUnknown conservatively when squeue is unavailable or times out.
// squeue only lists active jobs; empty output means the job has finished or is absent.
func (s *SlurmScheduler) GetJobStatus(jobID string) (JobStatus, error) {
	squeueBin := siblingBin(s.sbatchBin, "squeue")
	if squeueBin == "" {
		return JobStatusUnknown, nil // conservative: can't check
	}
	out, err := runCommand("SLURM", "job-status", squeueBin, "-j", jobID, "--noheader", "--format=%T")
	if err != nil {
		// squeue exits non-zero when the job ID is no longer known (completed/cancelled).
		// "Invalid job id" covers both old IDs and IDs that never existed.
		if strings.Contains(strings.ToLower(string(out)), "invalid job id") {
			return JobStatusDone, nil
		}
		return JobStatusUnknown, nil // conservative: timeout or unavailable
	}
	state := strings.TrimSpace(string(out))
	if state == "" {
		return JobStatusDone, nil // not in queue → finished or absent
	}
	switch strings.ToUpper(state) {
	case "PENDING", "CONFIGURING", "STAGE_IN":
		return JobStatusPending, nil
	default: // RUNNING, COMPLETING, STAGE_OUT, or any other active state
		return JobStatusRunning, nil
	}
}

// TryParseSlurmScript attempts to parse a SLURM script without requiring SLURM binaries.
// This is a static parser that can work in any environment.
func TryParseSlurmScript(scriptPath string) (*ScriptSpecs, error) {
	parser := &SlurmScheduler{
		directiveRe: regexp.MustCompile(`^\s*#SBATCH\s+(.+)$`),
	}
	return parser.ReadScriptSpecs(scriptPath)
}
