package utils

import (
	"fmt"
	"os"
	"regexp"
	"sort"
	"strconv"
	"strings"
	"time"

	"github.com/Justype/condatainer/catalog"
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

// FormatSize renders a byte count for display, in binary units to three
// significant figures. Artifacts are gigabytes, so a raw byte count is a number
// nobody reads.
func FormatSize(bytes int64) string {
	const unit = 1024
	if bytes < unit {
		return fmt.Sprintf("%d B", bytes)
	}
	value, exp := float64(bytes), 0
	for value >= unit && exp < 4 {
		value /= unit
		exp++
	}
	suffix := [...]string{"B", "KiB", "MiB", "GiB", "TiB"}[exp]
	if value >= 100 {
		return fmt.Sprintf("%.0f %s", value, suffix)
	}
	if value >= 10 {
		return fmt.Sprintf("%.1f %s", value, suffix)
	}
	return fmt.Sprintf("%.2f %s", value, suffix)
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

// GetDescriptionFromScript returns a script's first #DESC: value, or "".
func GetDescriptionFromScript(scriptPath string) string {
	text, err := os.ReadFile(scriptPath)
	if err != nil {
		return ""
	}
	description, _ := catalog.Find(catalog.ScanAnnotations(text), "#DESC")
	return description
}

// GetDependenciesFromScript parses a script and extracts the dependencies its
// #DEP: annotations declare, normalized and deduplicated.
//
// Position carries no meaning: an annotation counts wherever it is written.
//
// Only #DEP: counts. A `module load` line names the site's own module tree —
// another tool's namespace that merely spells names alike — so it is not a
// CondaTainer dependency and there is no artifact behind it to resolve.
func GetDependenciesFromScript(scriptPath string) ([]string, error) {
	if !FileExists(scriptPath) {
		return nil, fmt.Errorf("build script not found at %s", scriptPath)
	}
	text, err := os.ReadFile(scriptPath)
	if err != nil {
		return nil, fmt.Errorf("failed to open script: %w", err)
	}

	dependencies := []string{}
	seen := make(map[string]bool)
	for _, annotation := range catalog.Select(catalog.ScanAnnotations(text), "#DEP") {
		if annotation.Value == "" {
			continue
		}
		key := annotation.Value
		if !catalog.IsPathDep(key) {
			key = catalog.Normalize(key)
		}
		if !seen[key] {
			dependencies = append(dependencies, key)
			seen[key] = true
		}
	}
	return dependencies, nil
}

// GetTypeFromScript reads an external build script's #TYPE: header and returns
// the payload type it declares. This is the same axis catalog.DeriveType
// resolves — what the payload is — not how it is built.
//
// Only "app" and "data" are accepted, matching catalog.DeriveType exactly, so
// #TYPE: cannot mean one thing to a recipe and another to an external script.
// A script that declares none is an app.
func GetTypeFromScript(scriptPath string) (string, error) {
	if !FileExists(scriptPath) {
		return "", fmt.Errorf("build script not found at %s", scriptPath)
	}
	text, err := os.ReadFile(scriptPath)
	if err != nil {
		return "", fmt.Errorf("failed to open script: %w", err)
	}

	for _, annotation := range catalog.Select(catalog.ScanAnnotations(text), "#TYPE") {
		value := strings.ToLower(annotation.Value)
		if value == "" {
			continue
		}
		switch value {
		case "app", "data":
			return value, nil
		default:
			return "", fmt.Errorf("invalid TYPE value %q: valid values are app or data", value)
		}
	}
	return "app", nil
}

// GetTargetFromScript reads the artifact name an external build script declares
// with #TARGET:, or "" when it declares none.
//
// The name is not cosmetic: it becomes the payload's /cnt/<name> prefix, and
// key.Role reads it to decide which dependencies count toward equivalence — an
// app or OS dependency contributes only when its components occur in the
// artifact's own name. Taking it from the script rather than from the `-p` path
// keeps an identity-relevant input out of the invocation, so the same script
// built to two different paths cannot classify its dependencies two ways.
//
// A {placeholder} is refused. #TARGET: doubles as a template pattern for catalog
// recipes, where #PH: declarations supply the values and a requested name
// selects them; an external build is addressed by path and has neither, so a
// pattern here could never be filled in.
func GetTargetFromScript(scriptPath string) (string, error) {
	if !FileExists(scriptPath) {
		return "", fmt.Errorf("build script not found at %s", scriptPath)
	}
	text, err := os.ReadFile(scriptPath)
	if err != nil {
		return "", fmt.Errorf("failed to open script: %w", err)
	}

	for _, annotation := range catalog.Select(catalog.ScanAnnotations(text), "#TARGET") {
		value := strings.TrimSpace(annotation.Value)
		if value == "" {
			continue
		}
		if strings.ContainsAny(value, "{}") {
			return "", fmt.Errorf("invalid TARGET value %q: an external build takes a plain name, not a {placeholder}", value)
		}
		// Trimmed and checked for empty components, because catalog.Normalize does
		// neither and key.Role splits the name on "/". A stray slash would leave an
		// empty component, which quietly stops a dependency matching and downgrades
		// it to build history — the exact misclassification #TARGET: exists to stop.
		name := strings.Trim(catalog.Normalize(value), "/")
		if name == "" {
			continue
		}
		for _, component := range strings.Split(name, "/") {
			if component == "" {
				return "", fmt.Errorf("invalid TARGET value %q: it has an empty path component", value)
			}
		}
		return name, nil
	}
	return "", nil
}

// SortVersionsDescending sorts version strings in descending natural order.
// Segments are split on ".", "-", or "_" and compared numerically when both
// segments are integers, otherwise lexicographically. The highest version
// comes first.
func SortVersionsDescending(values []string) []string {
	result := make([]string, len(values))
	copy(result, values)
	sort.Slice(result, func(i, j int) bool {
		return catalog.CompareVersions(result[i], result[j]) > 0
	})
	return result
}

// phTokenRe matches {varname} placeholder tokens.
var phTokenRe = regexp.MustCompile(`\{(\w+)\}`)

// imgPackageTokenRe matches {KEY} tokens in an #IMG_PACKAGES: template.
var imgPackageTokenRe = regexp.MustCompile(`\{([A-Z][A-Z0-9_]*)\}`)

// ExtractImgPackageTokens returns the unique {KEY} token names from an #IMG_PACKAGES: template.
func ExtractImgPackageTokens(imgPackages string) []string {
	matches := imgPackageTokenRe.FindAllStringSubmatch(imgPackages, -1)
	seen := map[string]bool{}
	var result []string
	for _, m := range matches {
		if !seen[m[1]] {
			seen[m[1]] = true
			result = append(result, m[1])
		}
	}
	return result
}

// MatchVersion does a partial version match against a list sorted newest-first.
// Input "3.12" matches the first entry beginning with "3.12." (latest patch).
// Exact match always takes priority. Returns input unchanged if nothing matches.
func MatchVersion(input string, versions []string) string {
	for _, v := range versions {
		if v == input {
			return v
		}
	}
	prefix := input + "."
	for _, v := range versions {
		if strings.HasPrefix(v, prefix) {
			return v
		}
	}
	return input
}

// LatestVersion returns the first (newest) entry in a versions list, or "" if empty.
func LatestVersion(versions []string) string {
	if len(versions) == 0 {
		return ""
	}
	return versions[0]
}

// VersionChoicesDisplay formats a newest-first full-patch list into grouped minor
// version summary strings: ["3.13.2","3.12.10","3.11.12"] →
// ["3.13(.2)", "3.12(.10)", "3.11(.12)"].
func VersionChoicesDisplay(versions []string) []string {
	seen := map[string]bool{}
	var result []string
	for _, v := range versions {
		parts := strings.SplitN(v, ".", 3)
		if len(parts) < 2 {
			continue
		}
		minor := parts[0] + "." + parts[1]
		if !seen[minor] {
			seen[minor] = true
			if len(parts) == 3 {
				result = append(result, fmt.Sprintf("%s(.%s)", minor, parts[2]))
			} else {
				result = append(result, minor)
			}
		}
	}
	return result
}

// FormatChoicesInline formats a value list for single-line display.
// Version lists (entries containing a dot) are grouped by minor version with
// the patch shown in dim parentheses: "3.13(.2), 3.12(.10), 3.11(.12)".
// Option lists are joined with " | ": "github | microsoft".
func FormatChoicesInline(vlist []string) string {
	if len(vlist) == 0 {
		return ""
	}
	// Detect version list: first entry must have at least one dot.
	if !strings.Contains(vlist[0], ".") {
		return strings.Join(vlist, " | ")
	}
	seen := map[string]bool{}
	var parts []string
	for _, v := range vlist {
		segs := strings.SplitN(v, ".", 3)
		if len(segs) < 2 {
			continue
		}
		minor := segs[0] + "." + segs[1]
		if seen[minor] {
			continue
		}
		seen[minor] = true
		entry := minor
		if len(segs) == 3 {
			entry += "(." + segs[2] + ")"
		}
		parts = append(parts, entry)
	}
	return strings.Join(parts, ", ")
}
