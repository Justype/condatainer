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

	specs := &ScriptSpecs{
		Ncpus:    4, // default
		RawFlags: make([]string, 0),
	}

	scanner := bufio.NewScanner(file)
	lineNum := 0
	for scanner.Scan() {
		lineNum++
		line := scanner.Text()

		// Parse #SBATCH directives only
		if matches := s.directiveRe.FindStringSubmatch(line); matches != nil {
			flag := strings.TrimSpace(matches[1])
			specs.RawFlags = append(specs.RawFlags, flag)

			// Parse common SBATCH options
			if err := s.parseSbatchFlag(flag, specs); err != nil {
				return nil, NewParseError("SLURM", lineNum, line, err.Error())
			}
		}
	}

	if err := scanner.Err(); err != nil {
		return nil, fmt.Errorf("error reading script: %w", err)
	}

	return specs, nil
}

// parseSbatchFlag parses individual SBATCH flags and updates specs
func (s *SlurmScheduler) parseSbatchFlag(flag string, specs *ScriptSpecs) error {
	// Job name
	if strings.HasPrefix(flag, "--job-name=") {
		specs.JobName = strings.TrimPrefix(flag, "--job-name=")
	} else if strings.HasPrefix(flag, "-J ") {
		specs.JobName = strings.TrimSpace(strings.TrimPrefix(flag, "-J"))
	}

	// CPUs
	if strings.HasPrefix(flag, "--cpus-per-task=") {
		if _, err := fmt.Sscanf(flag, "--cpus-per-task=%d", &specs.Ncpus); err != nil {
			return fmt.Errorf("invalid cpus-per-task value: %w", err)
		}
	} else if strings.HasPrefix(flag, "-c ") {
		if _, err := fmt.Sscanf(flag, "-c %d", &specs.Ncpus); err != nil {
			return fmt.Errorf("invalid -c value: %w", err)
		}
	}

	// Memory
	if strings.HasPrefix(flag, "--mem=") {
		memStr := strings.TrimPrefix(flag, "--mem=")
		mem, err := parseMemory(memStr)
		if err != nil {
			return err
		}
		specs.MemMB = mem
	}

	// Time
	if strings.HasPrefix(flag, "--time=") {
		raw := strings.TrimPrefix(flag, "--time=")
		dur, err := parseSlurmTimeSpec(raw)
		if err != nil {
			return fmt.Errorf("invalid --time value: %w", err)
		}
		specs.Time = dur
	} else if strings.HasPrefix(flag, "-t ") {
		raw := strings.TrimSpace(strings.TrimPrefix(flag, "-t"))
		dur, err := parseSlurmTimeSpec(raw)
		if err != nil {
			return fmt.Errorf("invalid -t value: %w", err)
		}
		specs.Time = dur
	}

	// Output
	if strings.HasPrefix(flag, "--output=") {
		specs.Stdout = strings.TrimPrefix(flag, "--output=")
	} else if strings.HasPrefix(flag, "-o ") {
		specs.Stdout = strings.TrimSpace(strings.TrimPrefix(flag, "-o"))
	}

	// Error
	if strings.HasPrefix(flag, "--error=") {
		specs.Stderr = strings.TrimPrefix(flag, "--error=")
	} else if strings.HasPrefix(flag, "-e ") {
		specs.Stderr = strings.TrimSpace(strings.TrimPrefix(flag, "-e"))
	}

	// GPU
	if strings.HasPrefix(flag, "--gres=") {
		gpu, err := parseSlurmGpu(strings.TrimPrefix(flag, "--gres="))
		if err != nil {
			return err
		}
		specs.Gpu = gpu
	} else if strings.HasPrefix(flag, "--gpus=") {
		gpu, err := parseSlurmGpu(strings.TrimPrefix(flag, "--gpus="))
		if err != nil {
			return err
		}
		specs.Gpu = gpu
	} else if strings.HasPrefix(flag, "--gpus-per-node=") {
		gpu, err := parseSlurmGpu(strings.TrimPrefix(flag, "--gpus-per-node="))
		if err != nil {
			return err
		}
		specs.Gpu = gpu
	} else if strings.HasPrefix(flag, "--gpus-per-task=") {
		gpu, err := parseSlurmGpu(strings.TrimPrefix(flag, "--gpus-per-task="))
		if err != nil {
			return err
		}
		specs.Gpu = gpu
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
	}

	if strings.HasPrefix(flag, "--mail-user=") {
		specs.MailUser = strings.TrimPrefix(flag, "--mail-user=")
	}

	return nil
}

// CreateScriptWithSpec generates a SLURM batch script
func (s *SlurmScheduler) CreateScriptWithSpec(jobSpec *JobSpec, outputDir string) (string, error) {
	specs := jobSpec.Specs

	// Create output directory if specified
	if outputDir != "" {
		if err := os.MkdirAll(outputDir, 0775); err != nil {
			return "", NewScriptCreationError(jobSpec.Name, outputDir, err)
		}
	}

	// Generate script filename
	scriptName := "sbatch_job.sh"
	if jobSpec.Name != "" {
		scriptName = fmt.Sprintf("sbatch_%s.sh", jobSpec.Name)
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

	// Write parsed SBATCH directives with overrides
	for _, flag := range specs.RawFlags {
		// Skip output/error flags if we're overriding them
		if specs.Stdout != "" && (strings.HasPrefix(flag, "--output=") || strings.HasPrefix(flag, "-o ")) {
			continue
		}
		if specs.Stderr != "" && (strings.HasPrefix(flag, "--error=") || strings.HasPrefix(flag, "-e ")) {
			continue
		}
		if specs.Time > 0 && (strings.HasPrefix(flag, "--time=") || strings.HasPrefix(flag, "-t ")) {
			continue
		}
		// Skip email flags - we'll generate them from the parsed boolean fields
		if strings.HasPrefix(flag, "--mail-type=") || strings.HasPrefix(flag, "--mail-user=") {
			continue
		}
		fmt.Fprintf(writer, "#SBATCH %s\n", flag)
	}

	// Add custom stdout/stderr if specified
	if specs.Stdout != "" {
		fmt.Fprintf(writer, "#SBATCH --output=%s\n", specs.Stdout)
	}
	if specs.Stderr != "" {
		fmt.Fprintf(writer, "#SBATCH --error=%s\n", specs.Stderr)
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

	// Write the command
	fmt.Fprintln(writer, jobSpec.Command)

	// Add job ID echo for tracking
	fmt.Fprintln(writer, "echo SLURM_JOB_ID $SLURM_JOB_ID")

	// Self-delete the script after execution
	fmt.Fprintf(writer, "rm -f %s\n", scriptPath)

	// Make executable
	if err := os.Chmod(scriptPath, 0775); err != nil {
		return "", NewScriptCreationError(jobSpec.Name, scriptPath, err)
	}

	return scriptPath, nil
}

// Submit submits a SLURM job with optional dependency chain
func (s *SlurmScheduler) Submit(scriptPath string, dependencyJobIDs []string) (string, error) {
	args := []string{scriptPath}

	// Add dependency if provided
	if len(dependencyJobIDs) > 0 {
		depStr := strings.Join(dependencyJobIDs, ":")
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
	}

	// Get partition limits
	if s.scontrolCommand != "" {
		limits, err := s.getPartitionLimits()
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
func (s *SlurmScheduler) getPartitionLimits() ([]ResourceLimits, error) {
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

	return limits, nil
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
