package utils

import (
	"bufio"
	"fmt"
	"os"
	"regexp"
	"sort"
	"strconv"
	"strings"
	"time"
)

// StripInlineComment removes everything after the first '#' character (inline comment).
// Returns the trimmed string without the comment.
func StripInlineComment(s string) string {
	if idx := strings.Index(s, "#"); idx >= 0 {
		s = s[:idx]
	}
	return strings.TrimSpace(s)
}

// ParseMemoryMB parses memory strings like "8G", "1024M", "512K", "1T" into MB (int64).
// Default unit is MB when no suffix is given.
func ParseMemoryMB(memStr string) (int64, error) {
	return ParseMemoryMBWithDefault(memStr, "MB")
}

// FormatMemoryMB formats a MB value into a human-readable string.
// Uses "GB" when divisible by 1024, otherwise "MB". Examples: 8192 → "8GB", 1536 → "1536MB".
func FormatMemoryMB(mb int64) string {
	if mb%1024 == 0 {
		return fmt.Sprintf("%dGB", mb/1024)
	}
	return fmt.Sprintf("%dMB", mb)
}

// FormatDuration formats a duration dropping zero trailing components.
// Examples: 2h0m0s → "2h", 48h → "2d", 25h30m → "1d1h30m", 1h30m15s → "1h30m15s".
func FormatDuration(d time.Duration) string {
	total := int(d.Seconds())
	days := total / 86400
	hours := (total % 86400) / 3600
	mins := (total % 3600) / 60
	secs := total % 60

	var parts []string
	if days > 0 {
		parts = append(parts, fmt.Sprintf("%dd", days))
	}
	if hours > 0 {
		parts = append(parts, fmt.Sprintf("%dh", hours))
	}
	if mins > 0 {
		parts = append(parts, fmt.Sprintf("%dm", mins))
	}
	if secs > 0 {
		parts = append(parts, fmt.Sprintf("%ds", secs))
	}
	if len(parts) == 0 {
		return "0s"
	}
	return strings.Join(parts, "")
}

// ParseMemoryMBWithDefault parses a memory string into MB.
// When the string has an explicit unit suffix (G/GB/M/MB/K/KB/T/TB) it is used directly.
// For bare numbers the provided defaultUnit is applied (e.g. "KB", "MB", "GB").
func ParseMemoryMBWithDefault(memStr, defaultUnit string) (int64, error) {
	memStr = strings.ToUpper(strings.TrimSpace(memStr))

	var value int64
	var unit string

	n, err := fmt.Sscanf(memStr, "%d%s", &value, &unit)
	if err != nil && n == 0 {
		return 0, fmt.Errorf("invalid memory format: %s", memStr)
	}

	if unit == "" {
		unit = strings.ToUpper(strings.TrimSpace(defaultUnit))
	}

	switch unit {
	case "G", "GB":
		return value * 1024, nil
	case "M", "MB":
		return value, nil
	case "K", "KB":
		return value / 1024, nil
	case "T", "TB":
		return value * 1024 * 1024, nil
	default:
		return 0, fmt.Errorf("invalid memory unit %q in: %s", unit, memStr)
	}
}

// ParseSizeToMB converts strings like "10G", "500M", "1024" into Megabytes (int).
// Default unit is MB if no suffix is provided. Delegates to ParseMemoryMB.
func ParseSizeToMB(sizeStr string) (int, error) {
	mb, err := ParseMemoryMB(sizeStr)
	if err != nil {
		return 0, fmt.Errorf("invalid size format: %s (expected '10G', '500M', etc.)", sizeStr)
	}
	return int(mb), nil
}

// SplitDepConstraint splits a raw dep string into its nameVersion and optional
// version constraint. Supports ">=" and ">" operators.
// Example: "samtools/1.21>=1.16" → ("samtools/1.21", ">=", "1.16")
// Example: "samtools/1.21" → ("samtools/1.21", "", "")
func SplitDepConstraint(raw string) (nameVersion, op, minVersion string) {
	for _, sep := range []string{">=", ">"} {
		if idx := strings.Index(raw, sep); idx >= 0 {
			return raw[:idx], sep, raw[idx+len(sep):]
		}
	}
	return raw, "", ""
}

// versionParts splits a version string like "1.21.3" into integer components.
// Stops at the first non-numeric segment (e.g. "1.16rc1" → [1, 16]).
func versionParts(v string) []int {
	var parts []int
	for _, s := range strings.Split(v, ".") {
		n, err := strconv.Atoi(s)
		if err != nil {
			break
		}
		parts = append(parts, n)
	}
	return parts
}

// CompareVersions compares two partial version strings component by component.
// Missing trailing components are treated as 0 (so "1.16" == "1.16.0").
// Returns -1 if a < b, 0 if equal, 1 if a > b.
func CompareVersions(a, b string) int {
	ap, bp := versionParts(a), versionParts(b)
	n := len(ap)
	if len(bp) > n {
		n = len(bp)
	}
	for i := 0; i < n; i++ {
		av, bv := 0, 0
		if i < len(ap) {
			av = ap[i]
		}
		if i < len(bp) {
			bv = bp[i]
		}
		if av < bv {
			return -1
		}
		if av > bv {
			return 1
		}
	}
	return 0
}

// DepSatisfiedByVersion returns true if installedVersion satisfies the constraint
// (op + minVersion) and does not exceed preferredVersion (upper bound).
// op must be ">=" or ">". preferredVersion may be empty to skip the upper bound check.
func DepSatisfiedByVersion(installedVersion, op, minVersion, preferredVersion string) bool {
	cmp := CompareVersions(installedVersion, minVersion)
	var lower bool
	switch op {
	case ">=":
		lower = cmp >= 0
	case ">":
		lower = cmp > 0
	default:
		return false
	}
	if !lower {
		return false
	}
	// Upper bound: installed must not exceed the preferred version.
	if preferredVersion != "" && CompareVersions(installedVersion, preferredVersion) > 0 {
		return false
	}
	return true
}

// NormalizeNameVersion normalizes package spec formats so that
// "name/version", "name=version", "name@version" are treated the same.
// Converts = and @ to /, and -- to /, then strips whitespace.
// Any version constraint suffix (e.g. ">=1.10") is preserved as-is.
func NormalizeNameVersion(nameVersion string) string {
	nv, op, minVer := SplitDepConstraint(strings.TrimSpace(nameVersion))
	s := strings.ReplaceAll(nv, "=", "/")
	s = strings.ReplaceAll(s, "@", "/")
	s = strings.ReplaceAll(s, "--", "/")
	return s + op + minVer
}

// ParseHMSTime parses colon-separated walltime "HH:MM:SS", "HH:MM", or "MM" (bare minutes).
func ParseHMSTime(timeStr string) (time.Duration, error) {
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
		return 0, fmt.Errorf("invalid time format: %s", timeStr)
	}

	return time.Duration(hours*3600+minutes*60+seconds) * time.Second, nil
}

// ParseDHMSTime handles colon-separated and D-HH:MM:SS walltime formats.
func ParseDHMSTime(timeStr string) (time.Duration, error) {
	var days int64
	rest := timeStr
	if before, after, found := strings.Cut(timeStr, "-"); found {
		if d, err := strconv.ParseInt(before, 10, 64); err == nil {
			days = d
			rest = after
		}
	}
	dur, err := ParseHMSTime(rest)
	if err != nil {
		return 0, err
	}
	return time.Duration(days)*24*time.Hour + dur, nil
}

// parseCompoundDuration parses Go-style duration strings with an optional integer days prefix.
// Converts "Nd" to "N*24h" and delegates to time.ParseDuration.
// Examples: "4d12h" → 108h, "2h30m", "1.5h", "1d2h30m45s".
func parseCompoundDuration(s string) (time.Duration, error) {
	lower := strings.ToLower(strings.TrimSpace(s))
	if lower == "" {
		return 0, nil
	}
	if idx := strings.IndexByte(lower, 'd'); idx >= 0 {
		n, err := strconv.ParseInt(lower[:idx], 10, 64)
		if err != nil {
			return 0, fmt.Errorf("invalid time format: %s", s)
		}
		lower = fmt.Sprintf("%dh%s", n*24, lower[idx+1:])
	}
	d, err := time.ParseDuration(lower)
	if err != nil {
		return 0, fmt.Errorf("invalid time format: %s", s)
	}
	return d, nil
}

// ParseWalltime parses a walltime string into a duration. Supported formats:
//
// Compound (Go-style, with int day): 4d12h, 2h30m, 3h, 90m, 1d2h30m45s, 1.5h
//
// Colon-separated: D-HH:MM:SS, HH:MM:SS, HH:MM, MM
func ParseWalltime(timeStr string) (time.Duration, error) {
	timeStr = strings.TrimSpace(timeStr)
	if timeStr == "" {
		return 0, nil
	}
	// Compound duration: contains letters (d/h/m/s) with no colon
	if strings.ContainsAny(strings.ToLower(timeStr), "dhms") && !strings.Contains(timeStr, ":") {
		return parseCompoundDuration(timeStr)
	}
	// Colon-separated or D-HH:MM:SS: must contain only digits, colons, and dashes
	for _, c := range timeStr {
		if (c < '0' || c > '9') && c != ':' && c != '-' {
			return 0, fmt.Errorf("invalid time format: %s", timeStr)
		}
	}
	return ParseDHMSTime(timeStr)
}

// GetWhatIsFromScript reads a script and extracts the first #WHATIS: line.
// Returns the trimmed description string, or empty string if not found.
func GetWhatIsFromScript(scriptPath string) string {
	file, err := os.Open(scriptPath)
	if err != nil {
		return ""
	}
	defer file.Close()

	scanner := bufio.NewScanner(file)
	for scanner.Scan() {
		line := scanner.Text()
		if strings.HasPrefix(line, "#WHATIS:") {
			return strings.TrimSpace(line[len("#WHATIS:"):])
		}
	}
	return ""
}

// GetDependenciesFromScript parses a build script and extracts dependencies
// from #DEP: comments and, optionally, module load / ml commands.
// Set parseModuleLoad=true to also parse "module load" and "ml" lines.
// Returns a list of normalized dependency names with duplicates removed.
func GetDependenciesFromScript(scriptPath string, parseModuleLoad bool) ([]string, error) {
	if !FileExists(scriptPath) {
		return nil, fmt.Errorf("build script not found at %s", scriptPath)
	}

	file, err := os.Open(scriptPath)
	if err != nil {
		return nil, fmt.Errorf("failed to open script: %w", err)
	}
	defer file.Close()

	dependencies := []string{}
	seen := make(map[string]bool)

	var moduleLoadRegex *regexp.Regexp
	if parseModuleLoad {
		moduleLoadRegex = regexp.MustCompile(`^\s*(module\s+load)\s+(.+)$`)
	}
	scanner := bufio.NewScanner(file)

	for scanner.Scan() {
		line := scanner.Text()

		// Check for #DEP: comments
		if strings.HasPrefix(line, "#DEP:") {
			depLine := StripInlineComment(line[5:])
			if depLine != "" {
				if IsOverlay(depLine) {
					if !seen[depLine] {
						dependencies = append(dependencies, depLine)
						seen[depLine] = true
					}
				} else {
					// Split constraint before normalizing so "=" in ">=" isn't mangled.
					nv, op, minVer := SplitDepConstraint(depLine)
					normalized := NormalizeNameVersion(nv) + op + minVer
					if !seen[normalized] {
						dependencies = append(dependencies, normalized)
						seen[normalized] = true
					}
				}
			}
			continue
		}

		// Check for module load commands (only when enabled)
		if parseModuleLoad {
			line = strings.TrimSpace(line)
			if matches := moduleLoadRegex.FindStringSubmatch(line); matches != nil {
				parts := strings.Fields(line)
				if len(parts) >= 3 {
					for _, mod := range parts[2:] {
						normalized := NormalizeNameVersion(mod)
						if !seen[normalized] {
							dependencies = append(dependencies, normalized)
							seen[normalized] = true
						}
					}
				}
			} else if strings.HasPrefix(line, "ml") {
				parts := strings.Fields(line)
				if len(parts) >= 2 {
					for _, mod := range parts[1:] {
						// Skip ml subcommands that indicate no package names follow
						if mod == "purge" || mod == "list" || mod == "avail" || mod == "av" {
							break
						}
						// "ml load foo" may include the literal "load"; skip it and continue
						if mod == "load" {
							continue
						}
						normalized := NormalizeNameVersion(mod)
						if !seen[normalized] {
							dependencies = append(dependencies, normalized)
							seen[normalized] = true
						}
					}
				}
			}
		}
	}

	if err := scanner.Err(); err != nil {
		return nil, fmt.Errorf("error reading script: %w", err)
	}

	return dependencies, nil
}

// GetExternalBuildTypeFromScript parses external build script TYPE metadata.
// Supported tag forms:
//   - #TYPE:<value>
//   - TYPE:<value>
//
// If TYPE is not present, it defaults to "app".
// Supported values (case-insensitive):
//   - app aliases: app, env, tool, conda, small
//   - data aliases: data, ref, large
func GetExternalBuildTypeFromScript(scriptPath string) (string, error) {
	if !FileExists(scriptPath) {
		return "", fmt.Errorf("build script not found at %s", scriptPath)
	}

	file, err := os.Open(scriptPath)
	if err != nil {
		return "", fmt.Errorf("failed to open script: %w", err)
	}
	defer file.Close()

	scanner := bufio.NewScanner(file)
	for scanner.Scan() {
		line := strings.TrimSpace(scanner.Text())

		var raw string
		switch {
		case strings.HasPrefix(line, "#TYPE:"):
			raw = strings.TrimSpace(line[len("#TYPE:"):])
		case strings.HasPrefix(line, "TYPE:"):
			raw = strings.TrimSpace(line[len("TYPE:"):])
		default:
			continue
		}

		value := strings.ToLower(StripInlineComment(raw))
		if value == "" {
			continue
		}

		switch value {
		case "app", "env", "tool", "conda", "small":
			return "app", nil
		case "data", "ref", "large":
			return "data", nil
		default:
			return "", fmt.Errorf("invalid TYPE value %q: valid values are app or data", value)
		}
	}

	if err := scanner.Err(); err != nil {
		return "", fmt.Errorf("error reading script: %w", err)
	}

	return "app", nil
}

// GetInteractivePromptsFromScript parses a build script and extracts interactive
// prompt lines beginning with "#INTERACTIVE:". Returns a list of prompt
// strings (without the prefix) or an error if the file cannot be read.
func GetInteractivePromptsFromScript(scriptPath string) ([]string, error) {
	if !FileExists(scriptPath) {
		return nil, fmt.Errorf("build script not found at %s", scriptPath)
	}

	file, err := os.Open(scriptPath)
	if err != nil {
		return nil, fmt.Errorf("failed to open script: %w", err)
	}
	defer file.Close()

	prompts := []string{}
	scanner := bufio.NewScanner(file)
	for scanner.Scan() {
		line := scanner.Text()
		if strings.HasPrefix(line, "#INTERACTIVE:") {
			content := strings.TrimSpace(line[len("#INTERACTIVE:"):])
			if content != "" {
				prompts = append(prompts, content)
			}
		}
	}

	if err := scanner.Err(); err != nil {
		return nil, fmt.Errorf("error reading script: %w", err)
	}

	return prompts, nil
}

// PlaceholderDef holds a single #PL: declaration from a build script.
// Values are sorted descending (natural version sort); if Open is true,
// the last element of Values is "*" indicating any input is accepted.
type PlaceholderDef struct {
	Name   string
	Values []string // sorted desc concrete values; last element is "*" if Open
	Open   bool     // true if * was present in the original list
}

// SortVersionsDescending sorts version strings in descending natural order.
// Segments are split on ".", "-", or "_" and compared numerically when both
// segments are integers, otherwise lexicographically. The highest version
// comes first.
func SortVersionsDescending(values []string) []string {
	result := make([]string, len(values))
	copy(result, values)
	sort.Slice(result, func(i, j int) bool {
		return naturalVersionGreater(result[i], result[j])
	})
	return result
}

// naturalVersionGreater returns true if a > b in natural version order.
func naturalVersionGreater(a, b string) bool {
	segsA := splitVersionSegments(a)
	segsB := splitVersionSegments(b)
	n := len(segsA)
	if len(segsB) < n {
		n = len(segsB)
	}
	for i := 0; i < n; i++ {
		ia, aErr := strconv.ParseInt(segsA[i], 10, 64)
		ib, bErr := strconv.ParseInt(segsB[i], 10, 64)
		if aErr == nil && bErr == nil {
			if ia != ib {
				return ia > ib
			}
		} else {
			if segsA[i] != segsB[i] {
				return segsA[i] > segsB[i]
			}
		}
	}
	return len(segsA) > len(segsB)
}

var versionSplitRe = regexp.MustCompile(`[.\-_]`)

func splitVersionSegments(v string) []string {
	segs := versionSplitRe.Split(v, -1)
	result := make([]string, 0, len(segs))
	for _, s := range segs {
		if s != "" {
			result = append(result, s)
		}
	}
	return result
}

// parsePLValues parses the raw value string from a #PL: line.
// Handles comma lists, integer ranges (a-b), and the "*" open-ended marker.
// Returns sorted-desc concrete values (with "*" appended at end if open).
func parsePLValues(raw string) ([]string, bool, error) {
	tokens := strings.Split(raw, ",")
	var concrete []string
	open := false
	seen := make(map[string]bool)

	rangeRe := regexp.MustCompile(`^(\d+)-(\d+)$`)

	for _, tok := range tokens {
		tok = strings.TrimSpace(tok)
		if tok == "" {
			continue
		}
		if tok == "*" {
			open = true
			continue
		}
		if m := rangeRe.FindStringSubmatch(tok); m != nil {
			a, _ := strconv.ParseInt(m[1], 10, 64)
			b, _ := strconv.ParseInt(m[2], 10, 64)
			lo, hi := a, b
			if lo > hi {
				lo, hi = hi, lo
			}
			for v := lo; v <= hi; v++ {
				s := strconv.FormatInt(v, 10)
				if !seen[s] {
					concrete = append(concrete, s)
					seen[s] = true
				}
			}
		} else {
			if !seen[tok] {
				concrete = append(concrete, tok)
				seen[tok] = true
			}
		}
	}

	if len(concrete) == 0 && !open {
		return nil, false, fmt.Errorf("empty value list")
	}

	concrete = SortVersionsDescending(concrete)
	if open {
		concrete = append(concrete, "*")
	}
	return concrete, open, nil
}

// GetPlaceholdersFromScript parses a build script and returns all #PL:
// declarations in declaration order. Values are sorted descending; if the
// list contained "*" the last element is "*" and Open is true.
func GetPlaceholdersFromScript(scriptPath string) ([]PlaceholderDef, error) {
	if !FileExists(scriptPath) {
		return nil, fmt.Errorf("build script not found at %s", scriptPath)
	}

	file, err := os.Open(scriptPath)
	if err != nil {
		return nil, fmt.Errorf("failed to open script: %w", err)
	}
	defer file.Close()

	var defs []PlaceholderDef
	seenNames := make(map[string]bool)

	scanner := bufio.NewScanner(file)
	for scanner.Scan() {
		line := scanner.Text()
		if !strings.HasPrefix(line, "#PL:") {
			continue
		}
		rest := line[len("#PL:"):]
		// Split on first ":" to get name and raw values
		idx := strings.Index(rest, ":")
		if idx < 0 {
			return nil, fmt.Errorf("invalid #PL: line (missing second ':'): %s", line)
		}
		name := strings.TrimSpace(rest[:idx])
		rawValues := StripInlineComment(rest[idx+1:])

		if name == "" {
			return nil, fmt.Errorf("empty placeholder name in: %s", line)
		}
		if seenNames[name] {
			return nil, fmt.Errorf("duplicate placeholder name %q", name)
		}
		seenNames[name] = true

		values, open, err := parsePLValues(rawValues)
		if err != nil {
			return nil, fmt.Errorf("invalid #PL: values for %q: %w", name, err)
		}

		defs = append(defs, PlaceholderDef{Name: name, Values: values, Open: open})
	}

	if err := scanner.Err(); err != nil {
		return nil, fmt.Errorf("error reading script: %w", err)
	}

	return defs, nil
}

// GetTargetFromScript reads a build script and returns the value of the first
// #TARGET: line (trimmed). Returns "" if no #TARGET: line is present.
func GetTargetFromScript(scriptPath string) (string, error) {
	if !FileExists(scriptPath) {
		return "", fmt.Errorf("build script not found at %s", scriptPath)
	}

	file, err := os.Open(scriptPath)
	if err != nil {
		return "", fmt.Errorf("failed to open script: %w", err)
	}
	defer file.Close()

	scanner := bufio.NewScanner(file)
	for scanner.Scan() {
		line := scanner.Text()
		if strings.HasPrefix(line, "#TARGET:") {
			return strings.TrimSpace(line[len("#TARGET:"):]), nil
		}
	}

	if err := scanner.Err(); err != nil {
		return "", fmt.Errorf("error reading script: %w", err)
	}

	return "", nil
}

// InterpolateVars replaces {varname} occurrences in s with the corresponding
// value from vars. Unknown keys are left unchanged.
func InterpolateVars(s string, vars map[string]string) string {
	if len(vars) == 0 {
		return s
	}
	pairs := make([]string, 0, len(vars)*2)
	for k, v := range vars {
		pairs = append(pairs, "{"+k+"}", v)
	}
	return strings.NewReplacer(pairs...).Replace(s)
}

// ExpandPlaceholders returns the Cartesian product of all placeholder value
// combinations. "*" is never used as a concrete iteration value. Returns
// [{}] (a slice with one empty map) when defs is empty.
func ExpandPlaceholders(defs []PlaceholderDef) []map[string]string {
	result := []map[string]string{{}}

	for _, def := range defs {
		// Collect only concrete (non-"*") values for iteration
		var concreteVals []string
		for _, v := range def.Values {
			if v != "*" {
				concreteVals = append(concreteVals, v)
			}
		}
		if len(concreteVals) == 0 {
			// Open-only placeholder with no suggestions — skip expansion
			continue
		}

		var next []map[string]string
		for _, existing := range result {
			for _, val := range concreteVals {
				combo := make(map[string]string, len(existing)+1)
				for k, v := range existing {
					combo[k] = v
				}
				combo[def.Name] = val
				next = append(next, combo)
			}
		}
		result = next
	}

	return result
}

// DefaultVarsFromPlaceholders returns a map of {name: first_concrete_value}
// for each placeholder. If a placeholder is open-only (no concrete values),
// the key is omitted.
func DefaultVarsFromPlaceholders(defs []PlaceholderDef) map[string]string {
	vars := make(map[string]string, len(defs))
	for _, def := range defs {
		for _, v := range def.Values {
			if v != "*" {
				vars[def.Name] = v
				break
			}
		}
	}
	return vars
}
