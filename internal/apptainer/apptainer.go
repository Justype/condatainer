package apptainer

import (
	"bufio"
	"bytes"
	"fmt"
	"os"
	"os/exec"
	"regexp"
	"strconv"
	"strings"

	"github.com/Justype/condatainer/internal/config"
	"github.com/Justype/condatainer/internal/utils"
)

// ExecOptions holds all the configuration for running a container.
type ExecOptions struct {
	BaseImage string
	Command   []string
	Overlays  []string
	Binds     []string
	Env       []string 
	Writable  bool
	Fakeroot  bool
	
	// GPU Flags
	Nv   bool // Enables --nv (NVIDIA)
	Rocm bool // Enables --rocm (AMD)
}

// GetVersion returns the installed Apptainer version string (e.g., "1.5.3")
// Returns an error if not found.
func GetVersion() (string, error) {
	cmd := exec.Command(config.Global.ApptainerBin, "--version")
	output, err := cmd.Output()
	if err != nil {
		return "", fmt.Errorf("apptainer binary not found at %s", config.Global.ApptainerBin)
	}

	outStr := strings.TrimSpace(string(output))
	re := regexp.MustCompile(`(\d+\.\d+(\.\d+)?)`)
	match := re.FindString(outStr)

	if match == "" {
		return "", fmt.Errorf("could not parse version from output: %s", outStr)
	}

	utils.PrintDebug("Detected Apptainer Version: %s", utils.StyleNumber(match))
	return match, nil
}

// CheckZstdSupport verifies if the installed version supports ZSTD (requires >= 1.4)
func CheckZstdSupport(currentVersion string) bool {
	parts := strings.Split(currentVersion, ".")
	if len(parts) < 2 {
		return false 
	}

	major, _ := strconv.Atoi(parts[0])
	minor, _ := strconv.Atoi(parts[1])

	if major > 1 || (major == 1 && minor >= 4) {
		return true
	}
	
	utils.PrintDebug("ZSTD compression disabled (requires Apptainer >= 1.4, found %s)", currentVersion)
	return false
}

// CreateOverlay creates a generic SquashFS overlay image at the specified destination.
func CreateOverlay(path string, sizeMB int, sparse bool, fakeroot bool) error {
	args := []string{"overlay", "create", "--size", fmt.Sprintf("%d", sizeMB)}
	
	if sparse {
		args = append(args, "--sparse")
	}
	if fakeroot {
		args = append(args, "--fakeroot")
	}
	
	args = append(args, path)

	cmd := exec.Command(config.Global.ApptainerBin, args...)
	utils.PrintDebug("Running: %s %s", config.Global.ApptainerBin, strings.Join(args, " "))

	output, err := cmd.CombinedOutput()
	if err != nil {
		return &ApptainerError{
			Op:      fmt.Sprintf("create %dMB overlay", sizeMB),
			Path:    path,
			Cmd:     strings.Join(args, " "),
			Output:  string(output),
			BaseErr: err,
		}
	}
	return nil
}

// Exec runs a command inside the container using the provided options.
func Exec(opts ExecOptions) error {
	args := []string{"exec"}
	
	for _, o := range opts.Overlays {
		if opts.Writable {
			args = append(args, "--overlay", o)
		} else {
			args = append(args, "--overlay", o+":ro")
		}
	}
	for _, b := range opts.Binds {
		args = append(args, "--bind", b)
	}
	for _, e := range opts.Env {
		args = append(args, "--env", e)
	}
	if opts.Fakeroot {
		args = append(args, "--fakeroot")
	}
	if opts.Nv {
		args = append(args, "--nv")
	}
	if opts.Rocm {
		args = append(args, "--rocm")
	}

	args = append(args, opts.BaseImage)
	args = append(args, opts.Command...)

	cmd := exec.Command(config.Global.ApptainerBin, args...)
	cmd.Stdin = os.Stdin
	cmd.Stdout = os.Stdout
	cmd.Stderr = os.Stderr

	utils.PrintDebug("Exec: %s %s", config.Global.ApptainerBin, strings.Join(args, " "))
	
	if err := cmd.Run(); err != nil {
		return &ApptainerError{
			Op:      "execute container command",
			Path:    opts.BaseImage,
			Cmd:     strings.Join(args, " "),
			Output:  "",
			BaseErr: err,
		}
	}
	return nil
}

// BuildSif creates a SIF image from a definition file or source URI.
func BuildSif(destPath string, source string, fakeroot bool) error {
	args := []string{"build"}
	if fakeroot {
		args = append(args, "--fakeroot")
	}
	args = append(args, destPath, source)

	cmd := exec.Command(config.Global.ApptainerBin, args...)
	
	var stderr bytes.Buffer
	if config.Global.Debug {
		cmd.Stdout = os.Stdout
		cmd.Stderr = os.Stderr
	} else {
		cmd.Stderr = &stderr
	}
	
	utils.PrintDebug("Build: %s %s", config.Global.ApptainerBin, strings.Join(args, " "))
	
	if err := cmd.Run(); err != nil {
		return &ApptainerError{
			Op:      "build SIF image",
			Path:    destPath,
			Cmd:     strings.Join(args, " "),
			Output:  stderr.String(),
			BaseErr: err,
		}
	}
	return nil
}

// GetSquashfsSystemID parses 'apptainer sif list' to find the ID of the "System" partition.
func GetSquashfsSystemID(sifPath string) (int, error) {
	cmd := exec.Command(config.Global.ApptainerBin, "sif", "list", sifPath)
	output, err := cmd.Output()
	if err != nil {
		return 0, &ApptainerError{
			Op:      "list SIF content",
			Path:    sifPath,
			Cmd:     "sif list " + sifPath,
			Output:  string(output),
			BaseErr: err,
		}
	}

	scanner := bufio.NewScanner(bytes.NewReader(output))
	for scanner.Scan() {
		line := scanner.Text()
		if strings.Contains(line, "Squashfs") && (strings.Contains(line, "System") || strings.Contains(line, "Container")) {
			parts := strings.Split(line, "|")
			if len(parts) > 0 {
				idStr := strings.TrimSpace(parts[0])
				id, err := strconv.Atoi(idStr)
				if err != nil {
					continue
				}
				utils.PrintDebug("Found System SquashFS partition at ID: %d", id)
				return id, nil
			}
		}
	}

	return 0, &ApptainerError{
		Op:      "find Squashfs system partition",
		Path:    sifPath,
		Cmd:     "sif list " + sifPath,
		Output:  "No partition matching 'Squashfs' and 'System'/'Container' found.",
		BaseErr: fmt.Errorf("invalid or empty SIF image"),
	}
}

// DumpPartition extracts a specific partition ID from a SIF file to the destination.
func DumpPartition(sifPath string, id int, destPath string) error {
	outFile, err := os.Create(destPath)
	if err != nil {
		return fmt.Errorf("failed to create dump file: %w", err)
	}
	defer outFile.Close()

	args := []string{"sif", "dump", strconv.Itoa(id), sifPath}
	cmd := exec.Command(config.Global.ApptainerBin, args...)
	
	cmd.Stdout = outFile
	var stderr bytes.Buffer
	cmd.Stderr = &stderr

	utils.PrintDebug("Dump: %s %s > %s", config.Global.ApptainerBin, strings.Join(args, " "), destPath)

	if err := cmd.Run(); err != nil {
		return &ApptainerError{
			Op:      fmt.Sprintf("dump partition %d", id),
			Path:    sifPath,
			Cmd:     strings.Join(args, " "),
			Output:  stderr.String(),
			BaseErr: err,
		}
	}

	return nil
}

// ExtractSqf extracts the primary SquashFS partition from a SIF image.
// It uses GetSquashfsSystemID and DumpPartition to ensure an exact extraction.
func ExtractSqf(sifPath string, destSqfPath string) error {
	id, err := GetSquashfsSystemID(sifPath)
	if err != nil {
		return err 
	}

	err = DumpPartition(sifPath, id, destSqfPath)
	if err != nil {
		return err 
	}

	return nil
}
