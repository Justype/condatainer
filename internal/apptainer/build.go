package apptainer

import (
	"context"
	"fmt"
	"io"
	"os"
	"os/exec"
	"strconv"
	"strings"

	"github.com/Justype/condatainer/internal/logging"
	"github.com/Justype/condatainer/internal/utils"
)

// BuildOptions contains options for building a container image
type BuildOptions struct {
	Force      bool     // Force overwrite of existing image
	NoCleanup  bool     // Do not clean up bundle after failed build
	Additional []string // Additional flags to pass to apptainer build
}

// Build builds a container image from a definition file
// Always uses --fakeroot (required for non-root users on HPC systems)
func Build(ctx context.Context, imagePath, defFile string, opts *BuildOptions) error {
	if opts == nil {
		opts = &BuildOptions{}
	}

	args := []string{"build"}

	args = append(args, "--fakeroot")

	// Add optional flags
	if opts.Force {
		args = append(args, "--force")
	}
	if opts.NoCleanup {
		args = append(args, "--no-cleanup")
	}

	// Add user-provided extras.
	// args = append(args, DetectGPUFlags()...) // Do not detect the gpu when building
	args = append(args, opts.Additional...)

	// Add image path and definition file
	args = append(args, imagePath, defFile)

	logging.FromContext(ctx).Debug("building container", "image", imagePath, "definition", defFile)

	return runApptainerWithOutput(ctx, "build", imagePath, false, os.Stdin, os.Stdout, os.Stderr, args...)
}

// DumpSifToSquashfs extracts the SquashFS partition from a SIF image to a .sqfs file
// This matches the Python workflow: build SIF -> extract SquashFS -> use as overlay
//
// Workflow:
//  1. apptainer sif list <sifPath> - Find the Squashfs/*System partition ID
//  2. apptainer sif dump <id> <sifPath> > <outputPath> - Dump to .sqfs file
func DumpSifToSquashfs(ctx context.Context, sifPath, outputPath string) error {
	// Step 1: List SIF partitions to find SquashFS ID
	cmd := exec.CommandContext(ctx, apptainerCmd, "sif", "list", sifPath)
	output, err := cmd.Output()
	if err != nil {
		return &ApptainerError{
			Op:      "sif list",
			Cmd:     fmt.Sprintf("%s sif list %s", apptainerCmd, sifPath),
			Path:    sifPath,
			Output:  string(output),
			BaseErr: err,
		}
	}

	// Step 2: Parse output to find Squashfs/*System partition
	var squashfsID int
	found := false
	lines := strings.SplitSeq(string(output), "\n")
	for line := range lines {
		if strings.Contains(line, "Squashfs/*System") {
			// Format: "ID|GROUP|LINK|SIF POSITION|TYPE|..."
			// Example: "1|1|NONE|32768|Squashfs/*System|..."
			parts := strings.Split(line, "|")
			if len(parts) > 0 {
				id, err := strconv.Atoi(strings.TrimSpace(parts[0]))
				if err == nil {
					squashfsID = id
					found = true
					break
				}
			}
		}
	}

	if !found {
		return &ApptainerError{
			Op:      "sif list",
			Cmd:     fmt.Sprintf("%s sif list %s", apptainerCmd, sifPath),
			Path:    sifPath,
			Output:  string(output),
			BaseErr: fmt.Errorf("no Squashfs/*System partition found in SIF"),
		}
	}

	logging.FromContext(ctx).Debug("found SquashFS partition", "id", squashfsID, "sif", sifPath)

	// Step 3: Dump SquashFS partition to output file
	logging.FromContext(ctx).Debug("dumping SquashFS", "output", outputPath)

	// Create output file
	outFile, err := utils.CreateFileWritable(outputPath)
	if err != nil {
		return fmt.Errorf("failed to create output file %s: %w", outputPath, err)
	}
	defer outFile.Close()

	// Run: apptainer sif dump <id> <sifPath>
	dumpCmd := exec.CommandContext(ctx, apptainerCmd, "sif", "dump", strconv.Itoa(squashfsID), sifPath)
	dumpCmd.Stdout = outFile
	dumpCmd.Stderr = io.Discard

	if err := dumpCmd.Run(); err != nil {
		// Clean up partial file on error
		outFile.Close()
		os.Remove(outputPath)
		return &ApptainerError{
			Op:      "sif dump",
			Cmd:     fmt.Sprintf("%s sif dump %d %s", apptainerCmd, squashfsID, sifPath),
			Path:    sifPath,
			Output:  "",
			BaseErr: err,
		}
	}

	utils.ShareWithParentGroup(outputPath)

	logging.FromContext(ctx).Debug("dumped SquashFS", "output", outputPath)
	return nil
}
