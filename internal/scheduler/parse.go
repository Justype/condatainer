package scheduler

import (
	"bufio"
	"fmt"
	"os"

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
		Spec:           rs,
		Control:        rc,
		HasDirectives:  len(directives) > 0,
		RawFlags:       directives,
		RemainingFlags: remaining,
	}, nil
}
