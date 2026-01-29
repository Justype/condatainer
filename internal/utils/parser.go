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

// ParseSizeToMB converts strings like "10G", "500M", "1024" into Megabytes (int).
// Default unit is MB if no suffix is provided.
func ParseSizeToMB(sizeStr string) (int, error) {
	s := strings.TrimSpace(strings.ToUpper(sizeStr))

	// Regex to separate number and unit
	re := regexp.MustCompile(`^(\d+)(G|GB|M|MB|T|TB)?$`)
	matches := re.FindStringSubmatch(s)

	if len(matches) < 2 {
		return 0, fmt.Errorf("invalid size format: %s (expected '10G', '500M', etc.)", sizeStr)
	}

	val, err := strconv.Atoi(matches[1])
	if err != nil {
		return 0, fmt.Errorf("invalid number: %s", matches[1])
	}

	unit := matches[2]
	switch unit {
	case "G", "GB":
		return val * 1024, nil
	case "T", "TB":
		return val * 1048576, nil
	case "M", "MB", "":
		return val, nil
	default:
		return 0, fmt.Errorf("unsupported unit: %s", unit)
	}
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

// ParseDuration parses a duration string supporting multiple formats:
//   - Go duration: "2h", "30m", "1h30m", "90s"
//   - HH:MM:SS format: "02:00:00", "2:30:00", "00:30:00"
//   - H:MM format: "2:30" (interpreted as hours:minutes)
//
// Returns the duration in time.Duration format.
func ParseDuration(s string) (time.Duration, error) {
	s = strings.TrimSpace(s)
	if s == "" {
		return 0, fmt.Errorf("empty duration string")
	}

	// Try HH:MM:SS or H:MM:SS or HH:MM format first
	if strings.Contains(s, ":") {
		parts := strings.Split(s, ":")
		switch len(parts) {
		case 2:
			// H:MM or HH:MM format (hours:minutes)
			hours, err := strconv.Atoi(parts[0])
			if err != nil {
				return 0, fmt.Errorf("invalid hours: %s", parts[0])
			}
			minutes, err := strconv.Atoi(parts[1])
			if err != nil {
				return 0, fmt.Errorf("invalid minutes: %s", parts[1])
			}
			return time.Duration(hours)*time.Hour + time.Duration(minutes)*time.Minute, nil
		case 3:
			// HH:MM:SS format
			hours, err := strconv.Atoi(parts[0])
			if err != nil {
				return 0, fmt.Errorf("invalid hours: %s", parts[0])
			}
			minutes, err := strconv.Atoi(parts[1])
			if err != nil {
				return 0, fmt.Errorf("invalid minutes: %s", parts[1])
			}
			seconds, err := strconv.Atoi(parts[2])
			if err != nil {
				return 0, fmt.Errorf("invalid seconds: %s", parts[2])
			}
			return time.Duration(hours)*time.Hour +
				time.Duration(minutes)*time.Minute +
				time.Duration(seconds)*time.Second, nil
		default:
			return 0, fmt.Errorf("invalid time format: %s (use HH:MM:SS or HH:MM)", s)
		}
	}

	// Try Go duration format (2h, 30m, 1h30m, etc.)
	dur, err := time.ParseDuration(s)
	if err != nil {
		return 0, fmt.Errorf("invalid duration: %s (use '2h', '30m', '1h30m', or '02:00:00')", s)
	}
	return dur, nil
}

// GetDependenciesFromScript parses a build script and extracts dependencies
// from #DEP: comments and module load commands.
// Returns a list of normalized dependency names with duplicates removed.
func GetDependenciesFromScript(scriptPath string) ([]string, error) {
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

	moduleLoadRegex := regexp.MustCompile(`^\s*(module\s+load)\s+(.+)$`)
	scanner := bufio.NewScanner(file)

	for scanner.Scan() {
		line := scanner.Text()

		// Check for #DEP: comments
		if strings.HasPrefix(line, "#DEP:") {
			depLine := strings.TrimSpace(line[5:])
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

		// Check for module load commands
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
					// Skip ml subcommands
					if mod == "purge" || mod == "list" || mod == "avail" || mod == "av" || mod == "load" {
						break
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
