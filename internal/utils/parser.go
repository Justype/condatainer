package utils

import (
	"fmt"
	"regexp"
	"strconv"
	"strings"
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
