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

// SlurmScheduler implements the Scheduler interface for SLURM
type SlurmScheduler struct {
	sbatchBin       string
	sinfoCommand    string
	scontrolCommand string
	directiveRe     *regexp.Regexp
	jobIDRe         *regexp.Regexp
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

	sinfoCmd, _ := exec.LookPath("sinfo")
	scontrolCmd, _ := exec.LookPath("scontrol")

	return &SlurmScheduler{
		sbatchBin:       binPath,
		sinfoCommand:    sinfoCmd,
		scontrolCommand: scontrolCmd,
		directiveRe:     regexp.MustCompile(`^\s*#SBATCH\s+(.+)$`),
		jobIDRe:         regexp.MustCompile(`Submitted batch job (\d+)`),
	}, nil
}

// IsAvailable checks if SLURM is available and we're not inside a SLURM job
func (s *SlurmScheduler) IsAvailable() bool {
	if s.sbatchBin == "" {
		return false
	}

	// Check if we're already inside a SLURM job
	_, inJob := os.LookupEnv("SLURM_JOB_ID")
	if inJob {
		return false
	}

	return true
}

// GetInfo returns information about the SLURM scheduler
func (s *SlurmScheduler) GetInfo() *SchedulerInfo {
	_, inJob := os.LookupEnv("SLURM_JOB_ID")
	available := s.IsAvailable()

	info := &SchedulerInfo{
		Type:      "SLURM",
		Binary:    s.sbatchBin,
		InJob:     inJob,
		Available: available,
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
	cmd := exec.Command(s.sbatchBin, "--version")
	output, err := cmd.Output()
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

	var ntasksPerNode int
	var hasExplicitNtasks bool

	scanner := bufio.NewScanner(file)
	lineNum := 0
	for scanner.Scan() {
		lineNum++
		line := scanner.Text()

		// Parse #SBATCH directives only
		if matches := s.directiveRe.FindStringSubmatch(line); matches != nil {
			specs.HasDirectives = true
			flag := utils.StripInlineComment(matches[1])

			// Track --ntasks-per-node separately (needs post-processing with Nodes)
			if strings.HasPrefix(flag, "--ntasks-per-node=") {
				if _, err := fmt.Sscanf(flag, "--ntasks-per-node=%d", &ntasksPerNode); err != nil {
					return nil, NewParseError("SLURM", lineNum, line, "invalid --ntasks-per-node value")
				}
			}

			// Track explicit --ntasks to avoid overwriting with ntasks-per-node
			if strings.HasPrefix(flag, "--ntasks=") || strings.HasPrefix(flag, "-n ") {
				hasExplicitNtasks = true
			}

			// Parse common SBATCH options - returns true if recognized
			recognized, err := s.parseSbatchFlag(flag, specs)
			if err != nil {
				return nil, NewParseError("SLURM", lineNum, line, err.Error())
			}

			// Only keep unrecognized flags in RawFlags
			if !recognized {
				specs.RawFlags = append(specs.RawFlags, flag)
			}
		}
	}

	if err := scanner.Err(); err != nil {
		return nil, fmt.Errorf("error reading script: %w", err)
	}

	// Post-process: compute Ntasks from ntasks-per-node if --ntasks was not explicitly set
	if ntasksPerNode > 0 && !hasExplicitNtasks {
		specs.Ntasks = specs.Nodes * ntasksPerNode
	}

	return specs, nil
}

// parseSbatchFlag parses individual SBATCH flags and updates specs
// Returns true if the flag was recognized and parsed, false otherwise
func (s *SlurmScheduler) parseSbatchFlag(flag string, specs *ScriptSpecs) (bool, error) {
	// Job name
	if strings.HasPrefix(flag, "--job-name=") {
		specs.JobName = strings.TrimPrefix(flag, "--job-name=")
		return true, nil
	} else if strings.HasPrefix(flag, "-J ") {
		specs.JobName = strings.TrimSpace(strings.TrimPrefix(flag, "-J"))
		return true, nil
	}

	// CPUs per task
	if strings.HasPrefix(flag, "--cpus-per-task=") {
		if _, err := fmt.Sscanf(flag, "--cpus-per-task=%d", &specs.Ncpus); err != nil {
			return false, fmt.Errorf("invalid cpus-per-task value: %w", err)
		}
		return true, nil
	} else if strings.HasPrefix(flag, "-c ") {
		if _, err := fmt.Sscanf(flag, "-c %d", &specs.Ncpus); err != nil {
			return false, fmt.Errorf("invalid -c value: %w", err)
		}
		return true, nil
	}

	// Nodes
	if strings.HasPrefix(flag, "--nodes=") {
		if _, err := fmt.Sscanf(flag, "--nodes=%d", &specs.Nodes); err != nil {
			return false, fmt.Errorf("invalid --nodes value: %w", err)
		}
		return true, nil
	} else if strings.HasPrefix(flag, "-N ") {
		if _, err := fmt.Sscanf(flag, "-N %d", &specs.Nodes); err != nil {
			return false, fmt.Errorf("invalid -N value: %w", err)
		}
		return true, nil
	}

	// Ntasks
	if strings.HasPrefix(flag, "--ntasks=") {
		if _, err := fmt.Sscanf(flag, "--ntasks=%d", &specs.Ntasks); err != nil {
			return false, fmt.Errorf("invalid --ntasks value: %w", err)
		}
		return true, nil
	} else if strings.HasPrefix(flag, "-n ") {
		if _, err := fmt.Sscanf(flag, "-n %d", &specs.Ntasks); err != nil {
			return false, fmt.Errorf("invalid -n value: %w", err)
		}
		return true, nil
	}

	// Ntasks per node (also recognized)
	if strings.HasPrefix(flag, "--ntasks-per-node=") {
		return true, nil
	}

	// Memory
	if strings.HasPrefix(flag, "--mem=") {
		memStr := strings.TrimPrefix(flag, "--mem=")
		mem, err := parseMemory(memStr)
		if err != nil {
			return false, err
		}
		specs.MemMB = mem
		return true, nil
	}

	// Time
	if strings.HasPrefix(flag, "--time=") {
		raw := strings.TrimPrefix(flag, "--time=")
		dur, err := parseSlurmTimeSpec(raw)
		if err != nil {
			return false, fmt.Errorf("invalid --time value: %w", err)
		}
		specs.Time = dur
		return true, nil
	} else if strings.HasPrefix(flag, "-t ") {
		raw := strings.TrimSpace(strings.TrimPrefix(flag, "-t"))
		dur, err := parseSlurmTimeSpec(raw)
		if err != nil {
			return false, fmt.Errorf("invalid -t value: %w", err)
		}
		specs.Time = dur
		return true, nil
	}

	// Output
	if strings.HasPrefix(flag, "--output=") {
		specs.Stdout = strings.TrimPrefix(flag, "--output=")
		return true, nil
	} else if strings.HasPrefix(flag, "-o ") {
		specs.Stdout = strings.TrimSpace(strings.TrimPrefix(flag, "-o"))
		return true, nil
	}

	// Error
	if strings.HasPrefix(flag, "--error=") {
		specs.Stderr = strings.TrimPrefix(flag, "--error=")
		return true, nil
	} else if strings.HasPrefix(flag, "-e ") {
		specs.Stderr = strings.TrimSpace(strings.TrimPrefix(flag, "-e"))
		return true, nil
	}

	// GPU
	if strings.HasPrefix(flag, "--gres=gpu:") {
		gpu, err := parseSlurmGpu(strings.TrimPrefix(flag, "--gres="))
		if err != nil {
			return false, err
		}
		specs.Gpu = gpu
		return true, nil
	} else if strings.HasPrefix(flag, "--gpus=") {
		gpu, err := parseSlurmGpu(strings.TrimPrefix(flag, "--gpus="))
		if err != nil {
			return false, err
		}
		specs.Gpu = gpu
		return true, nil
	} else if strings.HasPrefix(flag, "--gpus-per-node=") {
		gpu, err := parseSlurmGpu(strings.TrimPrefix(flag, "--gpus-per-node="))
		if err != nil {
			return false, err
		}
		specs.Gpu = gpu
		return true, nil
	} else if strings.HasPrefix(flag, "--gpus-per-task=") {
		gpu, err := parseSlurmGpu(strings.TrimPrefix(flag, "--gpus-per-task="))
		if err != nil {
			return false, err
		}
		specs.Gpu = gpu
		return true, nil
	}

	// Email notifications
	if strings.HasPrefix(flag, "--mail-type=") {
		mailType := strings.ToUpper(strings.TrimPrefix(flag, "--mail-type="))
		// Parse comma-separated values (e.g., "BEGIN,END,FAIL" or "ALL")
		if mailType == "NONE" {
			// Explicitly disable all notifications
			specs.EmailOnBegin = false
			specs.EmailOnEnd = false
			specs.EmailOnFail = false
		} else {
			types := strings.Split(mailType, ",")
			for _, t := range types {
				t = strings.TrimSpace(t)
				switch t {
				case "ALL":
					specs.EmailOnBegin = true
					specs.EmailOnEnd = true
					specs.EmailOnFail = true
				case "BEGIN":
					specs.EmailOnBegin = true
				case "END":
					specs.EmailOnEnd = true
				case "FAIL", "REQUEUE", "INVALID_DEPEND":
					specs.EmailOnFail = true
					// Advanced types like TIME_LIMIT_*, STAGE_OUT, ARRAY_TASKS are ignored
					// as they don't map to the common begin/end/fail pattern
				}
			}
		}
		return true, nil
	}

	if strings.HasPrefix(flag, "--mail-user=") {
		specs.MailUser = strings.TrimPrefix(flag, "--mail-user=")
		return true, nil
	}

	// Flag not recognized
	return false, nil
}

// CreateScriptWithSpec generates a SLURM batch script
func (s *SlurmScheduler) CreateScriptWithSpec(jobSpec *JobSpec, outputDir string) (string, error) {
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
	scriptName := "job.sbatch"
	if jobSpec.Name != "" {
		// Replace slashes with -- to avoid creating subdirectories
		safeName := strings.ReplaceAll(jobSpec.Name, "/", "--")
		scriptName = fmt.Sprintf("%s.sbatch", safeName)
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

	// Write unrecognized flags (RawFlags only contains flags not parsed into typed fields)
	for _, flag := range specs.RawFlags {
		fmt.Fprintf(writer, "#SBATCH %s\n", flag)
	}

	// Add job name if specified
	if specs.JobName != "" {
		fmt.Fprintf(writer, "#SBATCH --job-name=%s\n", specs.JobName)
	}

	// Add custom stdout if specified (stderr is not used - output goes to stdout)
	if specs.Stdout != "" {
		fmt.Fprintf(writer, "#SBATCH --output=%s\n", specs.Stdout)
	}

	// Add node/task/cpu directives
	fmt.Fprintf(writer, "#SBATCH --nodes=%d\n", specs.Nodes)
	fmt.Fprintf(writer, "#SBATCH --ntasks=%d\n", specs.Ntasks)
	if specs.Ncpus > 0 {
		fmt.Fprintf(writer, "#SBATCH --cpus-per-task=%d\n", specs.Ncpus)
	}

	// Add email notifications if specified
	if specs.EmailOnBegin || specs.EmailOnEnd || specs.EmailOnFail {
		var mailTypes []string
		if specs.EmailOnBegin {
			mailTypes = append(mailTypes, "BEGIN")
		}
		if specs.EmailOnEnd {
			mailTypes = append(mailTypes, "END")
		}
		if specs.EmailOnFail {
			mailTypes = append(mailTypes, "FAIL")
		}
		fmt.Fprintf(writer, "#SBATCH --mail-type=%s\n", strings.Join(mailTypes, ","))
	}
	if specs.MailUser != "" {
		fmt.Fprintf(writer, "#SBATCH --mail-user=%s\n", specs.MailUser)
	}
	if specs.Time > 0 {
		fmt.Fprintf(writer, "#SBATCH --time=%s\n", formatSlurmTimeSpec(specs.Time))
	}

	// Add blank line
	fmt.Fprintln(writer, "")

	// Export resource variables for use in build scripts
	fmt.Fprintf(writer, "export NNODES=%d\n", specs.Nodes)
	fmt.Fprintf(writer, "export NTASKS=%d\n", specs.Ntasks)
	if specs.Ncpus > 0 {
		fmt.Fprintf(writer, "export NCPUS=%d\n", specs.Ncpus)
	} else {
		fmt.Fprintln(writer, "export NCPUS=1")
	}
	if specs.MemMB > 0 {
		fmt.Fprintf(writer, "export MEM=%d\n", specs.MemMB)
		fmt.Fprintf(writer, "export MEM_MB=%d\n", specs.MemMB)
		fmt.Fprintf(writer, "export MEM_GB=%d\n", specs.MemMB/1024)
	}
	fmt.Fprintln(writer, "")

	// Print job information at start
	fmt.Fprintln(writer, "# Print job information")
	fmt.Fprintln(writer, "_START_TIME=$SECONDS")
	fmt.Fprintln(writer, "_format_time() { local s=$1; printf '%02d:%02d:%02d' $((s/3600)) $((s%3600/60)) $((s%60)); }")
	fmt.Fprintln(writer, "echo \"========================================\"")
	fmt.Fprintln(writer, "echo \"Job ID:    $SLURM_JOB_ID\"")
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
		fmt.Fprintf(writer, "echo \"Time:      %s\"\n", formatSlurmTimeSpec(specs.Time))
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
	fmt.Fprintln(writer, "echo \"Job ID:    $SLURM_JOB_ID\"")
	fmt.Fprintln(writer, "echo \"Elapsed:   $(_format_time $(($SECONDS - $_START_TIME)))\"")
	fmt.Fprintf(writer, "%s\n", "echo \"Completed: $(date '+%Y-%m-%d %T')\"")
	fmt.Fprintln(writer, "echo \"========================================\"")

	// Self-dispose: remove this script file (unless in debug mode)
	if !config.Global.Debug {
		fmt.Fprintf(writer, "rm -f %s\n", scriptPath)
	}

	// Make executable
	if err := os.Chmod(scriptPath, utils.PermExec); err != nil {
		return "", NewScriptCreationError(jobSpec.Name, scriptPath, err)
	}

	return scriptPath, nil
}

// Submit submits a SLURM job with optional dependency chain
func (s *SlurmScheduler) Submit(scriptPath string, dependencyJobIDs []string) (string, error) {
	args := []string{scriptPath}

	// Add dependency if provided
	if len(dependencyJobIDs) > 0 {
		// Use comma-separated job IDs, see: https://bioinformaticsworkbook.org/Appendix/HPC/SLURM/submitting-dependency-jobs-using-slurm.html
		depStr := strings.Join(dependencyJobIDs, ",")
		depArg := fmt.Sprintf("--dependency=afterok:%s", depStr)
		args = append([]string{depArg}, args...)
	}

	// Execute sbatch
	cmd := exec.Command(s.sbatchBin, args...)
	output, err := cmd.CombinedOutput()
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

	return info, nil
}

// getGpuInfo queries SLURM for available GPU types
func (s *SlurmScheduler) getGpuInfo() ([]GpuInfo, error) {
	cmd := exec.Command(s.sinfoCommand, "-o", "%P|%G|%D|%T", "--noheader")
	output, err := cmd.CombinedOutput()
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

			// Dynamically add MIG profiles to the GPU database if not already present
			if IsMigProfile(gpuType) {
				EnsureMigProfileInDatabase(gpuType)
			}

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
	cmd := exec.Command(s.scontrolCommand, "show", "partition", "-o")
	output, err := cmd.CombinedOutput()
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
					if avail.MaxCpus > 0 {
						limits[i].MaxCpus = avail.MaxCpus
					}
					if avail.MaxMemMB > 0 {
						limits[i].MaxMemMB = avail.MaxMemMB
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
	cmd := exec.Command(s.sinfoCommand, "-o", "%R|%c|%m", "--noheader")
	output, err := cmd.CombinedOutput()
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
		memMB, _ = parseMemory(strings.TrimSpace(parts[2]))

		// Update max for this partition
		if existing, ok := resources[partition]; ok {
			if cpus > existing.MaxCpus {
				existing.MaxCpus = cpus
			}
			if memMB > existing.MaxMemMB {
				existing.MaxMemMB = memMB
			}
			resources[partition] = existing
		} else {
			resources[partition] = ResourceLimits{
				Partition: partition,
				MaxCpus:   cpus,
				MaxMemMB:  memMB,
			}
		}
	}

	return resources, nil
}

// getMaxNodeResources queries SLURM for maximum CPU and memory available per node
func (s *SlurmScheduler) getMaxNodeResources() (int, int64, error) {
	// Query sinfo for CPUs and memory per node: %c = CPUs, %m = memory in MB
	cmd := exec.Command(s.sinfoCommand, "-o", "%c|%m", "--noheader")
	output, err := cmd.CombinedOutput()
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
		memMB, _ = parseMemory(strings.TrimSpace(parts[1]))

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
				if dur, err := parseSlurmTimeSpec(value); err == nil {
					limit.MaxTime = dur
				}
			}
		case "DefaultTime":
			if value != "NONE" && value != "UNLIMITED" {
				if dur, err := parseSlurmTimeSpec(value); err == nil {
					limit.DefaultTime = dur
				}
			}
		case "MaxCPUsPerNode":
			fmt.Sscanf(value, "%d", &limit.MaxCpus)
		case "MaxMemPerNode":
			if value != "UNLIMITED" {
				limit.MaxMemMB, _ = parseMemory(value)
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

// parseMemory converts memory strings like "8G", "1024M" to MB
func parseMemory(memStr string) (int64, error) {
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

func parseSlurmTimeSpec(timeStr string) (time.Duration, error) {
	timeStr = strings.TrimSpace(timeStr)
	if timeStr == "" {
		return 0, nil
	}

	var days int64
	hms := timeStr
	if idx := strings.Index(hms, "-"); idx >= 0 {
		dayPart := hms[:idx]
		if parsed, err := strconv.ParseInt(dayPart, 10, 64); err == nil {
			days = parsed
		} else {
			return 0, fmt.Errorf("%w: %s", ErrInvalidTimeFormat, timeStr)
		}
		hms = strings.TrimSpace(hms[idx+1:])
	}

	var hours, minutes, seconds int64
	if hms != "" {
		parts := strings.Split(hms, ":")
		switch len(parts) {
		case 3:
			h, err := strconv.ParseInt(parts[0], 10, 64)
			if err != nil {
				return 0, fmt.Errorf("%w: %s", ErrInvalidTimeFormat, timeStr)
			}
			hours = h
			m, err := strconv.ParseInt(parts[1], 10, 64)
			if err != nil {
				return 0, fmt.Errorf("%w: %s", ErrInvalidTimeFormat, timeStr)
			}
			minutes = m
			s, err := strconv.ParseInt(parts[2], 10, 64)
			if err != nil {
				return 0, fmt.Errorf("%w: %s", ErrInvalidTimeFormat, timeStr)
			}
			seconds = s
		case 2:
			h, err := strconv.ParseInt(parts[0], 10, 64)
			if err != nil {
				return 0, fmt.Errorf("%w: %s", ErrInvalidTimeFormat, timeStr)
			}
			hours = h
			m, err := strconv.ParseInt(parts[1], 10, 64)
			if err != nil {
				return 0, fmt.Errorf("%w: %s", ErrInvalidTimeFormat, timeStr)
			}
			minutes = m
		case 1:
			m, err := strconv.ParseInt(parts[0], 10, 64)
			if err != nil {
				return 0, fmt.Errorf("%w: %s", ErrInvalidTimeFormat, timeStr)
			}
			minutes = m
		default:
			return 0, fmt.Errorf("%w: %s", ErrInvalidTimeFormat, timeStr)
		}
	}

	totalSeconds := days*24*3600 + hours*3600 + minutes*60 + seconds
	return time.Duration(totalSeconds) * time.Second, nil
}

func formatSlurmTimeSpec(d time.Duration) string {
	if d <= 0 {
		return ""
	}
	total := int64(d.Seconds())
	days := total / (24 * 3600)
	rem := total % (24 * 3600)
	hours := rem / 3600
	rem %= 3600
	minutes := rem / 60
	seconds := rem % 60
	if days > 0 {
		return fmt.Sprintf("%d-%02d:%02d:%02d", days, hours, minutes, seconds)
	}
	return fmt.Sprintf("%02d:%02d:%02d", hours, minutes, seconds)
}

// GetJobResources reads allocated resources from SLURM environment variables.
func (s *SlurmScheduler) GetJobResources() *JobResources {
	if _, ok := os.LookupEnv("SLURM_JOB_ID"); !ok {
		return nil
	}
	res := &JobResources{}
	res.Ncpus = getEnvInt("SLURM_CPUS_PER_TASK")
	res.Ntasks = getEnvInt("SLURM_NTASKS")
	res.Nodes = getEnvInt("SLURM_JOB_NUM_NODES")
	// SLURM_MEM_PER_NODE is already in MB
	res.MemMB = getEnvInt64("SLURM_MEM_PER_NODE")
	res.Ngpus = getCudaDeviceCount()
	return res
}

// TryParseSlurmScript attempts to parse a SLURM script without requiring SLURM binaries.
// This is a static parser that can work in any environment.
func TryParseSlurmScript(scriptPath string) (*ScriptSpecs, error) {
	parser := &SlurmScheduler{
		directiveRe: regexp.MustCompile(`^\s*#SBATCH\s+(.+)$`),
	}
	return parser.ReadScriptSpecs(scriptPath)
}
