package conda

import (
	"fmt"
	"os/exec"
	"strings"

	"github.com/Justype/condatainer/internal/utils"
)

// ListCondaPackages reads the conda-meta directory from a writable ext3 image and
// returns a map of package name to installed version for every package found.
// It returns (nil, nil) for non-.img files or when no conda-meta directory exists.
func ListCondaPackages(imgPath string) (map[string]string, error) {
	if !utils.IsImg(imgPath) {
		return nil, nil
	}
	dbg, err := exec.LookPath("debugfs")
	if err != nil {
		return nil, fmt.Errorf("debugfs not found: %w", err)
	}
	cmd := exec.Command(dbg, "-R", "ls -p upper/cnt_env/conda-meta", imgPath)
	out, err := cmd.Output()
	if err != nil {
		return nil, nil
	}
	pkgs := make(map[string]string)
	for _, line := range strings.Split(string(out), "\n") {
		line = strings.TrimSpace(line)
		parts := strings.Split(line, "/")
		if len(parts) < 6 {
			continue
		}
		filename := parts[5]
		if !strings.HasSuffix(filename, ".json") {
			continue
		}
		name, ver := parseCondaMetaFilename(filename)
		if name != "" {
			pkgs[name] = ver
		}
	}
	return pkgs, nil
}

// parseCondaMetaFilename parses a conda-meta filename into its name and version.
func parseCondaMetaFilename(filename string) (name, version string) {
	s := strings.TrimSuffix(filename, ".json")
	for i := 0; i < len(s)-1; i++ {
		if s[i] == '-' && s[i+1] >= '0' && s[i+1] <= '9' {
			name = s[:i]
			rest := s[i+1:]
			if j := strings.Index(rest, "-"); j >= 0 {
				version = rest[:j]
			} else {
				version = rest
			}
			return
		}
	}
	return "", ""
}
