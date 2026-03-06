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
		case flagMatches(flag, "-cwd"):
			v, _ := flagValue(flag, "-cwd")
			rc.WorkDir = absPath(v)
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

// lsfSpanInfo accumulates span[...], affinity[cores(...)] and per-slot memory info from -R.
// rawMLimitPerSlotMB is parsed from -M (memory limit/ulimit) and used only as a fallback
// for MemPerNodeMB when no rusage[mem=X] is present. -M itself passes through to RemainingFlags.
type lsfSpanInfo struct {
	singleNode         bool  // span[hosts=1]
	ptile              int   // span[ptile=M], 0 = not set
	cores              int   // affinity[cores(T)], 0 = not set
	rawMemPerSlotMB    int64 // rusage[mem=X] converted to MB, 0 = not set (priority)
	rawMLimitPerSlotMB int64 // -M X converted to MB, 0 = not set (fallback only)
}

// lsfAffinityRe matches affinity[cores(N)] in an LSF -R resource string.
var lsfAffinityRe = regexp.MustCompile(`affinity\[cores\((\d+)\)\]`)

// parseResourceSpec consumes resource fields from the directive list using a two-pass
// approach so that -n semantics can be resolved after all -R span/affinity info is known.
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

	// Pass 1: collect all directives; defer -n interpretation until span context is known.
	var rawN int
	var hasN bool
	var span lsfSpanInfo

	for _, flag := range directives {
		recognized := true
		var parseErr error

		switch {
		case flagMatches(flag, "-n"):
			_, parseErr = flagScanInt(flag, &rawN, "-n")
			hasN = true
		case flagMatches(flag, "-M"):
			_, parseErr = flagScan(flag, &span.rawMLimitPerSlotMB, parseLsfMemory, "-M")
			recognized = false // -M is a limit (ulimit), not allocation; re-emit unchanged
		case flagMatches(flag, "-W"):
			_, parseErr = flagScan(flag, &rs.Time, utils.ParseHMSTime, "-W")
		case flagMatches(flag, "-gpu"):
			gpuStr, _ := flagValue(flag, "-gpu")
			gpuStr = strings.Trim(gpuStr, "\"'")
			if gpu := parseLsfGpuDirective(gpuStr); gpu != nil {
				rs.Gpu = gpu
			}
		case flagMatches(flag, "-R"):
			resStr, _ := flagValue(flag, "-R")
			resStr = strings.Trim(resStr, "\"'")
			parseErr = l.parseLsfResourceIntoSpec(resStr, rs, &span)
		case flag == "-x":
			rs.Exclusive = true
		default:
			recognized = false
		}

		if parseErr != nil {
			logParseWarning("LSF: failed to parse directive %q: %v", flag, parseErr)
			return nil, directives
		}
		if !recognized {
			remaining = append(remaining, flag)
		}
	}

	// Pass 2: resolve -n meaning based on span context, then resolve per-slot memory.
	if hasN {
		switch {
		case span.singleNode:
			rs.Nodes = 1
			if span.cores > 0 {
				// Hybrid single-node: -n is number of processes; affinity[cores(T)] is CPUs per process.
				// slot = one process = T cores; -n M → M processes × T cores on one node.
				rs.Ntasks = rawN
				rs.TasksPerNode = rawN
				rs.CpusPerTask = span.cores
			} else {
				// OpenMP single-node: slot = one core; -n N means N CPU slots (cores).
				// 1 process uses all N cores; rusage[mem=X] is X per core (per slot).
				rs.Ntasks = 1
				rs.TasksPerNode = 1
				rs.CpusPerTask = rawN
			}
		case span.ptile > 0:
			// MPI (±OpenMP): -n is total tasks; ptile is tasks per node.
			// Use ceiling division to calculate nodes when not evenly divisible.
			rs.Ntasks = rawN
			rs.TasksPerNode = span.ptile
			rs.Nodes = (rawN + span.ptile - 1) / span.ptile // ceiling division
			if span.cores > 0 {
				rs.CpusPerTask = span.cores // MPI+OpenMP: affinity sets thread count
			} else {
				rs.CpusPerTask = 1 // pure MPI: single-threaded tasks
			}
		default:
			// No span: -n is total tasks, free distribution (no ptile constraint).
			// Nodes=0 and TasksPerNode=0 signal "unspecified" — generation omits span.
			// affinity[cores(T)] without ptile → hybrid free-dist; otherwise pure MPI.
			rs.Nodes = 0
			rs.Ntasks = rawN
			rs.TasksPerNode = 0
			if span.cores > 0 {
				rs.CpusPerTask = span.cores // affinity[cores(T)] without ptile → hybrid
			} else {
				rs.CpusPerTask = 1 // pure MPI
			}
		}
	}
	// Resolve memory: both rusage[mem=X] and -M are per-slot (per-task in LSF).
	// rusage[mem] takes priority over -M.
	// rusage[mem=N] bare numbers are MB (LSF docs); -M bare numbers are KB (LSF convention).
	// Both respect LSF_UNIT_FOR_LIMITS when set.
	// Slot semantics depend on affinity:
	//   span[hosts=1] no affinity : slot = 1 core  → MemPerCpuMB = X
	//   span[hosts=1] affinity[T] : slot = T cores → MemPerCpuMB = X/T
	//   span[ptile=M]             : slot = 1 task  → MemPerCpuMB = X/1 (pure MPI)
	//                                                MemPerCpuMB = X/T (hybrid, affinity[T])
	//   no span (free-dist)       : slot = 1 task  → same as ptile
	// When -n is absent: treat memory value as per-node directly.
	cpt := span.cores
	if cpt <= 0 {
		cpt = 1
	}
	if hasN {
		if span.rawMemPerSlotMB > 0 {
			rs.MemPerCpuMB = span.rawMemPerSlotMB / int64(cpt)
			rs.MemPerNodeMB = 0 // clear default so GetMemPerNodeMB() derives correctly
		} else if span.rawMLimitPerSlotMB > 0 {
			rs.MemPerCpuMB = span.rawMLimitPerSlotMB / int64(cpt)
			rs.MemPerNodeMB = 0 // clear default so GetMemPerNodeMB() derives correctly
		}
	} else {
		// No -n: can't determine slot count; treat memory value as per-node directly.
		if span.rawMemPerSlotMB > 0 {
			rs.MemPerNodeMB = span.rawMemPerSlotMB
		} else if span.rawMLimitPerSlotMB > 0 {
			rs.MemPerNodeMB = span.rawMLimitPerSlotMB
		}
	}

	return rs, remaining
}

// parseLsfResourceIntoSpec parses LSF -R resource requirement strings into a ResourceSpec
// and accumulates span/affinity info in span.
// Supports: rusage[mem=N], rusage[ngpus_physical=N], span[hosts=1], span[ptile=M],
// affinity[cores(T)].
// Returns an error if a known key has an invalid value.
func (l *LsfScheduler) parseLsfResourceIntoSpec(resStr string, rs *ResourceSpec, span *lsfSpanInfo) error {
	// Look for span[...] block
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
				if len(kv) != 2 {
					continue
				}
				key := strings.TrimSpace(kv[0])
				value := strings.TrimSpace(kv[1])
				switch key {
				case "hosts":
					n, err := strconv.Atoi(value)
					if err != nil {
						return fmt.Errorf("invalid span[hosts=] value %q: %w", value, err)
					}
					if n == 1 {
						span.singleNode = true
					} else {
						logParseWarning("LSF: span[hosts=%d] > 1 is not supported; ignoring (use span[ptile=M] for multi-node)", n)
					}
				case "ptile":
					n, err := strconv.Atoi(value)
					if err != nil {
						return fmt.Errorf("invalid span[ptile=] value %q: %w", value, err)
					}
					span.ptile = n
				}
			}
		}
	}

	// Look for affinity[cores(N)] block
	if m := lsfAffinityRe.FindStringSubmatch(resStr); m != nil {
		n, err := strconv.Atoi(m[1])
		if err == nil {
			span.cores = n
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
			// rusage[mem=N] defaults to MB per LSF docs; LSF_UNIT_FOR_LIMITS overrides.
			mem, err := parseLsfMemoryWithUnits(value)
			if err != nil {
				return fmt.Errorf("invalid rusage[mem=] value %q: %w", value, err)
			}
			span.rawMemPerSlotMB = mem // per-slot; resolved to per-cpu in Pass 2
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
		case "type", "gmodel":
			spec.Type = value
		}
	}

	return spec
}

// CreateScriptWithSpec generates an LSF batch script
func (l *LsfScheduler) CreateScriptWithSpec(jobSpec *JobSpec, outputDir string) (string, error) {
	specs := jobSpec.Specs

	// No TasksPerNode normalization: TasksPerNode=0 means "freely distributed" and
	// routes to the no-span MPI path in generation. span[ptile=0] is never emitted.

	// Create output directory if specified
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
	scriptName := "job.lsf"
	if jobSpec.Name != "" {
		scriptName = fmt.Sprintf("%s.lsf", safeJobName(jobSpec.Name))
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

	// Add job name if specified; for array jobs embed the range inside the name
	if specs.Control.JobName != "" {
		if jobSpec.Array != nil {
			arr := jobSpec.Array
			r := fmt.Sprintf("[1-%d", arr.Count)
			if arr.Limit > 0 {
				r += fmt.Sprintf("%%%d", arr.Limit)
			}
			fmt.Fprintf(writer, "#BSUB -J \"%s%s]\"\n", specs.Control.JobName, r)
		} else {
			fmt.Fprintf(writer, "#BSUB -J %s\n", specs.Control.JobName)
		}
	}

	// Add custom stdout if specified
	if specs.Control.WorkDir != "" {
		fmt.Fprintf(writer, "#BSUB -cwd %s\n", specs.Control.WorkDir)
	}
	if specs.Control.Stdout != "" {
		fmt.Fprintf(writer, "#BSUB -o %s\n", specs.Control.AbsStdout())
	}
	if specs.Control.Stderr != "" {
		fmt.Fprintf(writer, "#BSUB -e %s\n", specs.Control.AbsStderr())
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
		rs := specs.Spec

		// rusage[mem=X] is the per-slot memory allocation request in LSF.
		// Ceiling-divide per-node MB by slotsPerNode to get per-slot MB.
		// MPI jobs: slot = one task (slotsPerNode = TasksPerNode).
		// OpenMP / single-node / no-span: slot = one CPU (slotsPerNode = CpusPerTask).
		// isMPI = GetNtasks() > 1 (uses Ntasks if set, else Nodes×TasksPerNode, min 1).
		// Three generation paths:
		//   !isMPI                     → OpenMP: span[hosts=1]
		//   isMPI && TasksPerNode > 0  → MPI/Hybrid with known ptile: span[ptile=M]
		//   isMPI && TasksPerNode == 0 → free distribution: no span, just -n Ntasks
		isMPI := rs.IsMPI()
		// Calculate TasksPerNode when both Nodes and Ntasks are set but TasksPerNode is not
		tasksPerNode := rs.TasksPerNode
		if tasksPerNode == 0 && rs.Nodes > 0 && rs.GetNtasks() > 0 {
			// Use ceiling to ensure all tasks fit
			tasksPerNode = (rs.GetNtasks() + rs.Nodes - 1) / rs.Nodes
		}

		memStr := ""
		if rs.MemPerCpuMB > 0 {
			// Inverse of parse: rusage_mem (per-slot) = MemPerCpuMB × affinity_cores.
			// Hybrid uses affinity[cores(CpusPerTask)]; all other cases affinity_cores = 1.
			cpuPerSlot := 1
			if isMPI && rs.CpusPerTask > 1 {
				cpuPerSlot = rs.CpusPerTask
			}
			memStr = fmt.Sprintf(" rusage[mem=%dMB]", rs.MemPerCpuMB*int64(cpuPerSlot))
		} else if rs.MemPerNodeMB > 0 {
			// Per-node → per-slot: divide by slotsPerNode.
			// Use calculated tasksPerNode if available; otherwise fall back to slotsPerNode=1.
			slotsPerNode := 1
			if isMPI {
				if tasksPerNode > 0 {
					slotsPerNode = tasksPerNode
				}
			} else {
				if rs.CpusPerTask > 0 {
					slotsPerNode = rs.CpusPerTask
				}
			}
			memPerSlotMB := (rs.MemPerNodeMB + int64(slotsPerNode) - 1) / int64(slotsPerNode)
			memStr = fmt.Sprintf(" rusage[mem=%dMB]", memPerSlotMB)
		}

		switch {
		case !isMPI:
			// OpenMP / single-node: -n = total CPUs; span[hosts=1]
			if rs.CpusPerTask > 0 {
				fmt.Fprintf(writer, "#BSUB -n %d\n", rs.CpusPerTask)
			}
			fmt.Fprintf(writer, "#BSUB -R \"span[hosts=1]%s\"\n", memStr)
		case tasksPerNode > 0:
			// MPI/Hybrid with known tasks-per-node: span[ptile=M]
			fmt.Fprintf(writer, "#BSUB -n %d\n", rs.GetNtasks())
			if rs.CpusPerTask > 1 {
				fmt.Fprintf(writer, "#BSUB -R \"span[ptile=%d] affinity[cores(%d)]%s\"\n",
					tasksPerNode, rs.CpusPerTask, memStr)
			} else {
				fmt.Fprintf(writer, "#BSUB -R \"span[ptile=%d]%s\"\n", tasksPerNode, memStr)
			}
		default:
			// MPI/Hybrid, free distribution: no Nodes specified, only Ntasks
			fmt.Fprintf(writer, "#BSUB -n %d\n", rs.GetNtasks())
			rStr := ""
			if rs.CpusPerTask > 1 {
				rStr = fmt.Sprintf("affinity[cores(%d)]%s", rs.CpusPerTask, memStr)
			} else {
				rStr = strings.TrimSpace(memStr)
			}
			if rStr != "" {
				fmt.Fprintf(writer, "#BSUB -R \"%s\"\n", rStr)
			}
		}
		if rs.Gpu != nil && rs.Gpu.Count > 0 {
			if rs.Gpu.Type != "" && rs.Gpu.Type != "gpu" {
				fmt.Fprintf(writer, "#BSUB -gpu \"num=%d:gmodel=%s\"\n", rs.Gpu.Count, rs.Gpu.Type)
			} else {
				fmt.Fprintf(writer, "#BSUB -gpu \"num=%d\"\n", rs.Gpu.Count)
			}
		}
		if rs.Time > 0 {
			fmt.Fprintf(writer, "#BSUB -W %s\n", formatLsfTime(rs.Time))
		}
		if rs.Exclusive {
			fmt.Fprintln(writer, "#BSUB -x")
		}
	}

	fmt.Fprintln(writer, "")

	// Export normalized resource env vars: LSF's native env vars (LSB_DJOB_NUMPROC, etc.)
	// don't directly provide the job geometry. Export NCPUS, NNODES, etc. so that:
	// 1. User scripts can access them directly (e.g., `salmon -p $NCPUS`)
	// 2. GetJobResources() can read them back reliably
	// 3. Nested "condatainer run" inherits correct resource context
	if specs.Spec != nil {
		for _, kv := range ResourceEnvVars(specs.Spec) {
			fmt.Fprintf(writer, "export %s\n", kv)
		}
		fmt.Fprintln(writer, "")
	}

	// Array job: extract input line, set ARRAY_ARGS, and redirect output
	if jobSpec.Array != nil {
		// LSF is 1-indexed; LSB_JOBINDEX maps directly to the sed line number
		writeArrayBlock(writer, "$LSB_JOBINDEX",
			jobSpec.Array.InputFile, outputDir, safeJobName(jobSpec.Name),
			jobSpec.Array.Count, arraySeparateOutput)
		jobSpec.Metadata["Array Job ID"] = "$LSB_JOBID"
		jobSpec.Metadata["Array Index"] = "$LSB_JOBINDEX"
		jobSpec.Metadata["Array File"] = jobSpec.Array.InputFile
		jobSpec.Metadata["Array Args"] = "$ARRAY_ARGS"
	}

	// Print job information at start
	writeJobHeader(writer, "$LSB_JOBID", specs, formatLsfTime, jobSpec.Metadata)
	fmt.Fprintln(writer, "")

	// Write the command and capture exit code
	fmt.Fprintln(writer, jobSpec.Command)
	fmt.Fprintln(writer, "_EXIT_CODE=$?")

	// Print completion info
	fmt.Fprintln(writer, "")
	writeJobFooter(writer, "$LSB_JOBID")

	// Self-delete the script after execution (unless in debug mode)
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

// Submit submits an LSF job with optional dependency chain
// buildLsfDepCondition returns the condition expression for bsub -w, or "" if deps is empty.
// Format: done(ID1) && exit(ID2) && ended(ID3)
// afterok→done(), afternotok→exit(), afterany→ended()
func buildLsfDepCondition(deps []Dependency) string {
	var conditions []string
	for _, dep := range deps {
		var lsfFn string
		switch dep.Type {
		case DependencyAfterNotOK:
			lsfFn = "exit"
		case DependencyAfterAny:
			lsfFn = "ended"
		default: // DependencyAfterOK and any unknown type
			lsfFn = "done"
		}
		for _, id := range dep.JobIDs {
			conditions = append(conditions, fmt.Sprintf("%s(%s)", lsfFn, id))
		}
	}
	return strings.Join(conditions, " && ")
}

// buildLsfArgs returns the bsub arguments derived from deps (not including "< scriptPath").
func buildLsfArgs(deps []Dependency) []string {
	var args []string
	if cond := buildLsfDepCondition(deps); cond != "" {
		args = append(args, "-w", cond, "-ti")
	}
	return args
}

func (l *LsfScheduler) Submit(scriptPath string, deps []Dependency) (string, error) {
	args := buildLsfArgs(deps)

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
					Partition:       queueName,
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
	return parseMemoryMB(memStr)
}

// parseLsfMemoryWithUnits parses an LSF rusage[mem=] or select[mem=] value to MB.
// Per LSF docs, bare numbers in rusage[] and select[] default to MB.
// LSF_UNIT_FOR_LIMITS overrides the default when set.
func parseLsfMemoryWithUnits(memStr string) (int64, error) {
	memStr = strings.TrimSpace(memStr)
	if memStr == "" {
		return 0, fmt.Errorf("empty memory string")
	}
	defaultUnit := strings.ToUpper(strings.TrimSpace(os.Getenv("LSF_UNIT_FOR_LIMITS")))
	if defaultUnit == "" {
		defaultUnit = "MB" // rusage[mem=] and select[mem=] default to MB per LSF docs
	}
	result, err := utils.ParseMemoryMBWithDefault(memStr, defaultUnit)
	if err != nil {
		return 0, fmt.Errorf("%w: %s", ErrInvalidMemoryFormat, memStr)
	}
	return result, nil
}

// parseLsfMemory parses an LSF -M (hard memory limit/ulimit) value to MB.
// Per LSF convention, bare numbers for -M default to KB.
// LSF_UNIT_FOR_LIMITS overrides the default when set.
func parseLsfMemory(memStr string) (int64, error) {
	memStr = strings.TrimSpace(memStr)
	if memStr == "" {
		return 0, fmt.Errorf("empty memory string")
	}
	defaultUnit := strings.ToUpper(strings.TrimSpace(os.Getenv("LSF_UNIT_FOR_LIMITS")))
	if defaultUnit == "" {
		defaultUnit = "KB" // -M defaults to KB per LSF convention
	}
	result, err := utils.ParseMemoryMBWithDefault(memStr, defaultUnit)
	if err != nil {
		return 0, fmt.Errorf("%w: %s", ErrInvalidMemoryFormat, memStr)
	}
	return result, nil
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

// GetJobResources reads allocated resources from LSF environment variables and
//
// LSF does not expose Nodes or TasksPerNode natively, so geometry is read from
// our injected NNODES / NTASKS_PER_NODE / NCPUS / NTASKS vars first.
//
// LSB_DJOB_NUMPROC / LSB_MAX_NUM_PROCESSORS is the total number of allocated
// slots (= Ntasks × CpusPerTask). It is used as a Ntasks fallback when our
// injected NTASKS is absent; CpusPerTask (from NCPUS) is used to split the
// slot count correctly.
//
// LSB_MAX_MEM_RUSAGE is per-slot (per-task) memory in KB → MemPerCpuMB.
func (s *LsfScheduler) GetJobResources() *ResourceSpec {
	if _, ok := os.LookupEnv("LSB_JOBID"); !ok {
		return nil
	}
	res := &ResourceSpec{}

	// 1. Geometry from normalized env vars or MPI library vars (most reliable sources).
	if v := getEnvInt("NNODES"); v != nil {
		res.Nodes = *v
	}
	// TasksPerNode: prefer MPI local size (actual tasks on this node)
	if mpiLocal := getMpiLocalSize(); mpiLocal != nil {
		res.TasksPerNode = *mpiLocal
	} else if v := getEnvInt("NTASKS_PER_NODE"); v != nil {
		res.TasksPerNode = *v
	}
	if v := getEnvInt("NCPUS"); v != nil {
		res.CpusPerTask = *v
	}
	// Prefer MPI library vars (global) over our normalized NTASKS var
	if mpiSize := getMpiCommSize(); mpiSize != nil {
		res.Ntasks = *mpiSize
	} else if v := getEnvInt("NTASKS"); v != nil {
		res.Ntasks = *v
	}

	// 2. Total slots from LSF native vars; used only when NTASKS was not injected.
	if res.Ntasks == 0 {
		var totalSlots int
		if v := getEnvInt("LSB_DJOB_NUMPROC"); v != nil {
			totalSlots = *v
		} else if v := getEnvInt("LSB_MAX_NUM_PROCESSORS"); v != nil {
			totalSlots = *v
		}
		if totalSlots > 0 {
			if res.CpusPerTask > 1 {
				// Hybrid: total slots = Ntasks × CpusPerTask.
				res.Ntasks = totalSlots / res.CpusPerTask
			} else {
				// Pure MPI or unknown: treat each slot as one task.
				res.Ntasks = totalSlots
				if res.CpusPerTask == 0 {
					res.CpusPerTask = 1
				}
			}
		}
	}

	// 3. Memory: LSB_MAX_MEM_RUSAGE is per-slot (per-task) in KB → MemPerCpuMB.
	if memKB := getEnvInt64("LSB_MAX_MEM_RUSAGE"); memKB != nil {
		if mb := *memKB / 1024; mb > 0 {
			res.MemPerCpuMB = mb
		}
	}

	if n := getCudaDeviceCount(); n != nil && *n > 0 {
		res.Gpu = &GpuSpec{Count: *n}
	}
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
