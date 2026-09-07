package conda

import (
	"fmt"
	"os/exec"
	"strings"

	"github.com/Justype/condatainer/internal/toolpath"
	"github.com/Justype/condatainer/internal/utils"
)

// ListCondaPackages reads the conda-meta directory from a writable ext3 image and
// returns a map of package name to installed version for every package found.
// It returns (nil, nil) for non-.img files or when no conda-meta directory exists.
func ListCondaPackages(imgPath string) (map[string]string, error) {
	if !utils.IsImg(imgPath) {
		return nil, nil
	}
	dbg, err := toolpath.Resolve("debugfs")
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

// ListCondaPackagesSqf lists packages recorded in an env-typed .sqf's
// cnt_env/conda-meta directory — the read-only counterpart to
// ListCondaPackages, using `unsquashfs -l` instead of debugfs since a
// SquashFS archive has no ext3 upper/lower distinction to read. Returns (nil,
// nil) for a non-.sqf path or when the directory does not exist, the same
// signal ListCondaPackages gives for "no conda environment here at all."
func ListCondaPackagesSqf(sqfPath string) (map[string]string, error) {
	if !utils.IsSqf(sqfPath) {
		return nil, nil
	}
	bin, err := toolpath.Resolve("unsquashfs")
	if err != nil {
		return nil, nil
	}
	cmd := exec.Command(bin, "-l", "-d", "", "-no-progress", sqfPath, "cnt_env/conda-meta")
	out, err := cmd.Output()
	if err != nil {
		return nil, nil
	}
	var pkgs map[string]string
	for _, line := range strings.Split(string(out), "\n") {
		line = strings.TrimSpace(line)
		if !strings.HasPrefix(line, "/cnt_env/conda-meta/") {
			continue
		}
		filename := line[strings.LastIndex(line, "/")+1:]
		if !strings.HasSuffix(filename, ".json") {
			continue
		}
		name, ver := parseCondaMetaFilename(filename)
		if name == "" {
			continue
		}
		if pkgs == nil {
			pkgs = make(map[string]string)
		}
		pkgs[name] = ver
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
