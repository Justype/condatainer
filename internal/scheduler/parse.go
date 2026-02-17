package scheduler

import (
	"bufio"
	"fmt"
	"os"
	"strconv"
	"strings"

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
