package utils

import (
	"bufio"
	"fmt"
	"os"
	"regexp"
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
	memStr = strings.ToUpper(strings.TrimSpace(memStr))

	var value int64
	var unit string

	n, err := fmt.Sscanf(memStr, "%d%s", &value, &unit)
	if err != nil && n == 0 {
		return 0, fmt.Errorf("invalid memory format: %s", memStr)
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

// ParseSizeToMB converts strings like "10G", "500M", "1024" into Megabytes (int).
// Default unit is MB if no suffix is provided. Delegates to ParseMemoryMB.
func ParseSizeToMB(sizeStr string) (int, error) {
	mb, err := ParseMemoryMB(sizeStr)
	if err != nil {
		return 0, fmt.Errorf("invalid size format: %s (expected '10G', '500M', etc.)", sizeStr)
	}
	return int(mb), nil
}

// NormalizeNameVersion normalizes package spec formats so that
// "name/version", "name=version", "name@version" are treated the same.
// Converts = and @ to /, and -- to /, then strips whitespace.
func NormalizeNameVersion(nameVersion string) string {
	s := strings.TrimSpace(nameVersion)
	s = strings.ReplaceAll(s, "=", "/")
	s = strings.ReplaceAll(s, "@", "/")
	s = strings.ReplaceAll(s, "--", "/")
	return s
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
					normalized := NormalizeNameVersion(depLine)
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
