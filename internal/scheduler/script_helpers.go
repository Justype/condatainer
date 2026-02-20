package scheduler

import (
	"bufio"
	"fmt"
	"io"
	"os"
	"path/filepath"
	"strconv"
	"strings"
	"time"

	"github.com/Justype/condatainer/internal/utils"
)

// logParseWarning prints a parser diagnostic warning, except in a job.
func logParseWarning(format string, args ...any) {
	if !IsInsideJob() {
		utils.PrintWarning(format, args...)
	}
}

// readFileLines opens a file and returns all its lines.
// Shared helper used by all scheduler ReadScriptSpecs implementations.
func readFileLines(path string) ([]string, error) {
	file, err := os.Open(path)
	if err != nil {
		if os.IsNotExist(err) {
			return nil, fmt.Errorf("%w: %s", ErrScriptNotFound, path)
		}
		return nil, err
	}
	defer file.Close()

	var lines []string
	scanner := bufio.NewScanner(file)
	for scanner.Scan() {
		lines = append(lines, scanner.Text())
	}
	if err := scanner.Err(); err != nil {
		return nil, fmt.Errorf("error reading script: %w", err)
	}
	return lines, nil
}

// parseScript is the shared two-stage pipeline called by each scheduler's ReadScriptSpecs.
//
// Stage 1 (critical): ParseRuntimeConfig — job control settings.
//   - Returns an error on critical failure (hard stop).
//   - RuntimeConfig fields are always consumed; they never appear in RemainingFlags.
//
// Stage 2 (best-effort): parseResourceSpec — compute geometry.
//   - Returns nil on any parse failure → passthrough mode (warning printed).
//   - Unconsumed directives become RemainingFlags.
//
// RawFlags is the immutable audit log of ALL directives in the script.
func parseScript(
	scriptPath string,
	lines []string,
	extractor func([]string) []string,
	rcParser func([]string) (RuntimeConfig, []string, error),
	rsParser func([]string) (*ResourceSpec, []string),
) (*ScriptSpecs, error) {
	directives := extractor(lines)

	rc, unconsumed, err := rcParser(directives)
	if err != nil {
		return nil, err // Critical: RuntimeConfig parse failure stops everything
	}

	rs, remaining := rsParser(unconsumed)
	if rs == nil && len(unconsumed) > 0 {
		// Duplicate warning, already printed by the scheduler-specific parser; no need to print again here.
		// utils.PrintWarning("Could not parse resource directives; using passthrough mode")
	}

	return &ScriptSpecs{
		ScriptPath:     scriptPath,
		Spec:           rs,
		Control:        rc,
		HasDirectives:  len(directives) > 0,
		RawFlags:       directives,
		RemainingFlags: remaining,
	}, nil
}

// flagValue extracts the value from a CLI flag, trying each prefix in order.
// Handles "prefix=value" and "prefix value" (space-separated) forms.
// Returns ("", false) if no prefix matches.
func flagValue(flag string, prefixes ...string) (string, bool) {
	for _, prefix := range prefixes {
		if v, ok := strings.CutPrefix(flag, prefix+"="); ok {
			return v, true
		}
		if v, ok := strings.CutPrefix(flag, prefix+" "); ok {
			return strings.TrimSpace(v), true
		}
	}
	return "", false
}

// flagMatches reports whether flag matches any prefix in "prefix=…" or "prefix …" form.
// Intended for use in switch case expressions.
func flagMatches(flag string, prefixes ...string) bool {
	_, ok := flagValue(flag, prefixes...)
	return ok
}

// flagScan extracts a value from a CLI flag and writes it into *dest using the provided parser.
// Generic base — works for any type (int, int64, time.Duration, *GpuSpec, …).
// Returns (false, nil) when no prefix matches; (true, err) on parse failure.
func flagScan[T any](flag string, dest *T, parser func(string) (T, error), prefixes ...string) (bool, error) {
	v, ok := flagValue(flag, prefixes...)
	if !ok {
		return false, nil
	}
	result, err := parser(v)
	if err == nil {
		*dest = result
	}
	return true, err
}

// flagScanInt is a convenience wrapper for flagScan using strconv.Atoi.
func flagScanInt(flag string, dest *int, prefixes ...string) (bool, error) {
	return flagScan(flag, dest, strconv.Atoi, prefixes...)
}

// safeJobName converts a job name to a filesystem-safe string by replacing "/" with "--".
func safeJobName(name string) string {
	return strings.ReplaceAll(name, "/", "--")
}

// absPath returns the absolute form of path. If path is already absolute or
// filepath.Abs fails, it returns path unchanged.
func absPath(path string) string {
	if filepath.IsAbs(path) {
		return path
	}
	if abs, err := filepath.Abs(path); err == nil {
		return abs
	}
	return path
}

// ResourceEnvVars returns KEY=VALUE strings for resource-related environment
// variables derived from rs. All counts default to 1 when rs is nil or unset.
//
//   - NNODES          — number of nodes
//   - NTASKS_PER_NODE — tasks (MPI ranks) per node
//   - NTASKS          — total tasks (NNODES × NTASKS_PER_NODE)
//   - NCPUS           — CPUs per node (CpusPerTask × TasksPerNode, PBS-style)
//   - NCPUS_PER_TASK  — CPUs per task (SLURM-style --cpus-per-task)
//   - MEM / MEM_MB / MEM_GB — included only when MemPerNodeMB > 0
func ResourceEnvVars(rs *ResourceSpec) []string {
	nodes, tasksPerNode, cpusPerTask := 1, 1, 1
	if rs != nil {
		if rs.Nodes > 0 {
			nodes = rs.Nodes
		}
		if rs.TasksPerNode > 0 {
			tasksPerNode = rs.TasksPerNode
		}
		if rs.CpusPerTask > 0 {
			cpusPerTask = rs.CpusPerTask
		}
	}
	env := []string{
		fmt.Sprintf("NNODES=%d", nodes),
		fmt.Sprintf("NTASKS_PER_NODE=%d", tasksPerNode),
		fmt.Sprintf("NTASKS=%d", nodes*tasksPerNode),
		fmt.Sprintf("NCPUS=%d", cpusPerTask*tasksPerNode), // CPUs per node (PBS-style)
		fmt.Sprintf("NCPUS_PER_TASK=%d", cpusPerTask),
	}
	if rs != nil && rs.MemPerNodeMB > 0 {
		env = append(env,
			fmt.Sprintf("MEM=%d", rs.MemPerNodeMB),
			fmt.Sprintf("MEM_MB=%d", rs.MemPerNodeMB),
			fmt.Sprintf("MEM_GB=%d", rs.MemPerNodeMB/1024),
		)
	}
	return env
}

// writeJobHeader writes the job info header echo block to w.
//
// - jobIDVar is the shell expression for the job ID (e.g. "$SLURM_JOB_ID").
//
// - *ScriptSpecs; resource lines (Nodes, Tasks, CPUs/Task, Memory, Time)
// are printed only when specs != nil && specs.Spec != nil.
// specs.ScriptPath is printed when available.
//
// - formatTime formats rs.Time for display; pass nil to skip the Time line.
func writeJobHeader(w io.Writer, jobIDVar string, specs *ScriptSpecs, formatTime func(time.Duration) string, metadata map[string]string) {
	fmt.Fprintln(w, "# Print job information")
	fmt.Fprintln(w, "_START_TIME=$SECONDS")
	fmt.Fprintln(w, "_format_time() { local s=$1; printf '%02d:%02d:%02d' $((s/3600)) $((s%3600/60)) $((s%60)); }")
	fmt.Fprintln(w, "echo \"========================================\"")
	fmt.Fprintf(w, "echo \"Job ID:    %s\"\n", jobIDVar)
	fmt.Fprintf(w, "echo \"Job Name:  %s\"\n", specs.Control.JobName)
	if specs != nil && specs.ScriptPath != "" {
		fmt.Fprintf(w, "echo \"Script:    %s\"\n", specs.ScriptPath)
	}

	// Extract ResourceSpec (may be nil → passthrough mode)
	rs := (*ResourceSpec)(nil)
	if specs != nil {
		rs = specs.Spec
	}
	if rs != nil {
		nodes := rs.Nodes
		if nodes <= 0 {
			nodes = 1
		}
		tasksPerNode := rs.TasksPerNode
		if tasksPerNode <= 0 {
			tasksPerNode = 1
		}
		if nodes > 1 {
			fmt.Fprintf(w, "echo \"Nodes:     %d\"\n", nodes)
		}
		if nodes*tasksPerNode > 1 {
			fmt.Fprintf(w, "echo \"Tasks:     %d\"\n", nodes*tasksPerNode)
		}
		if rs.CpusPerTask > 0 {
			fmt.Fprintf(w, "echo \"CPUs/Task: %d\"\n", rs.CpusPerTask)
		}
		if rs.MemPerNodeMB > 0 {
			fmt.Fprintf(w, "echo \"Memory:    %d MB\"\n", rs.MemPerNodeMB)
		}
		if rs.Time > 0 && formatTime != nil {
			fmt.Fprintf(w, "echo \"Time:      %s\"\n", formatTime(rs.Time))
		}
	}
	fmt.Fprintln(w, "echo \"PWD:       $(pwd)\"")
	if len(metadata) > 0 {
		maxLen := 0
		for key := range metadata {
			if len(key) > maxLen {
				maxLen = len(key)
			}
		}
		for key, value := range metadata {
			if value != "" {
				padding := maxLen - len(key)
				fmt.Fprintf(w, "echo \"%s:%s %s\"\n", key, strings.Repeat(" ", padding+3), value)
			}
		}
	}
	fmt.Fprintf(w, "%s\n", "echo \"Started:   $(date '+%Y-%m-%d %T')\"")
	fmt.Fprintln(w, "echo \"========================================\"")
}

// writeJobFooter writes the job completion footer echo block to w.
// jobIDVar is the shell expression for the job ID (e.g. "$SLURM_JOB_ID").
func writeJobFooter(w io.Writer, jobIDVar string) {
	fmt.Fprintln(w, "echo \"========================================\"")
	fmt.Fprintf(w, "echo \"Job ID:    %s\"\n", jobIDVar)
	fmt.Fprintln(w, "echo \"Elapsed:   $(_format_time $(($SECONDS - $_START_TIME)))\"")
	fmt.Fprintf(w, "%s\n", "echo \"Completed: $(date '+%Y-%m-%d %T')\"")
	fmt.Fprintln(w, "echo \"========================================\"")
}
