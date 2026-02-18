package scheduler

import (
	"bufio"
	"fmt"
	"io"
	"os"
	"strconv"
	"strings"
	"time"

	"github.com/Justype/condatainer/internal/utils"
)

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
		utils.PrintWarning("Could not parse resource directives; using passthrough mode")
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

// writeEnvVars writes resource environment variable exports to w.
// Must only be called when rs != nil.
//
// Always exports: NNODES, NTASKS, NTASKS_PER_NODE, NCPUS (defaults to 1 when CpusPerTask == 0)
// Conditionally: MEM, MEM_MB, MEM_GB (only when MemPerNodeMB > 0)
func writeEnvVars(w io.Writer, rs *ResourceSpec) {
	nodes := rs.Nodes
	if nodes <= 0 {
		nodes = 1
	}
	tasksPerNode := rs.TasksPerNode
	if tasksPerNode <= 0 {
		tasksPerNode = 1
	}
	fmt.Fprintf(w, "export NNODES=%d\n", nodes)
	fmt.Fprintf(w, "export NTASKS=%d\n", nodes*tasksPerNode)
	fmt.Fprintf(w, "export NTASKS_PER_NODE=%d\n", tasksPerNode)
	if rs.CpusPerTask > 0 {
		fmt.Fprintf(w, "export NCPUS=%d\n", rs.CpusPerTask)
	} else {
		fmt.Fprintln(w, "export NCPUS=1")
	}
	if rs.MemPerNodeMB > 0 {
		fmt.Fprintf(w, "export MEM=%d\n", rs.MemPerNodeMB)
		fmt.Fprintf(w, "export MEM_MB=%d\n", rs.MemPerNodeMB)
		fmt.Fprintf(w, "export MEM_GB=%d\n", rs.MemPerNodeMB/1024)
	}
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
