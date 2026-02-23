package cmd

import (
	"context"
	"errors"
	"os"
	"path/filepath"
	"strings"

	"github.com/Justype/condatainer/internal/build"
	"github.com/Justype/condatainer/internal/config"
	"github.com/Justype/condatainer/internal/utils"
	"github.com/spf13/cobra"
)

// Variables to hold flag values
var (
	createName      string
	createPrefix    string
	createFile      string
	createBaseImage string
	createSource    string
	createTempSize  string
	createRemote    bool
	createUpdate    bool

	// compression flags are generated dynamically from config.CompressOptions
	compFlags map[string]*bool
)

// compressArgsFromFlags inspects the map of boolean pointers produced by
// flag registration and returns the corresponding mksquashfs arguments.
// If more than one compression flag is set, it returns an error.
func compressArgsFromFlags(flags map[string]*bool) (string, error) {
	selected := ""
	for name, ptr := range flags {
		if ptr != nil && *ptr {
			if selected != "" {
				return "", errors.New("multiple compression options specified")
			}
			selected = name
		}
	}
	if selected == "" {
		return "", nil
	}
	return config.ArgsForCompress(selected), nil
}

var createCmd = &cobra.Command{
	Use:     "create [packages...]",
	Aliases: []string{"install", "i"},
	Short:   "Create a new SquashFS overlay",
	Long: `Create a new SquashFS overlay using available build scripts or Conda packages.

Note: If creation jobs are submitted to a scheduler, exits with code 3.`,
	Example: `  condatainer create samtools/1.22                # Create from build script
  condatainer create python=3.11 numpy -n myenv  # Create conda environment
  condatainer create -f environment.yml -p myenv  # Create from conda file
  condatainer create --source /data -p dataset    # Convert directory to overlay`,
	Run: func(cmd *cobra.Command, args []string) {
		ctx := cmd.Context()
		// 1. Validation Logic
		if len(args) == 0 && createFile == "" && createSource == "" {
			ExitWithError("At least one of [packages], --file, or --source must be provided.")
		}
		if len(args) > 0 && createPrefix != "" {
			ExitWithError("Packages cannot be used with --prefix.")
		}
		if createPrefix != "" && createName != "" {
			ExitWithError("Cannot use both --prefix and --name at the same time.")
		}
		if createPrefix != "" {
			baseName := strings.TrimSuffix(filepath.Base(createPrefix), ".sqf")
			if strings.Contains(baseName, "--") {
				ExitWithError("--prefix name cannot contain '--' (reserved name/version separator)")
			}
		}
		if createSource != "" && createName == "" && createPrefix == "" {
			ExitWithError("When using --source, either --name or --prefix must be provided.")
		}
		if createFile != "" && createPrefix == "" {
			ExitWithError("When using --file, --prefix must be provided.")
		}

		// 2. Ensure base image exists (also checks for apptainer)
		if err := build.EnsureBaseImage(ctx, false); err != nil {
			ExitWithError("Failed to ensure base image: %v", err)
		}

		// 3. Handle Compression Config – consult helper that respects available
		// options and rejects multiple selections.
		if args, err := compressArgsFromFlags(compFlags); err != nil {
			ExitWithError("%v", err)
		} else if args != "" {
			config.Global.Build.CompressArgs = args
		}
		// If no compression flag provided, use config default (auto-detected in root.go)

		// 4. Handle temp size
		if createTempSize != "" {
			sizeMB, err := utils.ParseSizeToMB(createTempSize)
			if err != nil {
				ExitWithError("Invalid temp size: %v", err)
			}
			config.Global.Build.TmpSizeMB = sizeMB
		}

		// 5. Normalize package names (only for build-script mode, not for conda/prefix/source modes)
		normalizedArgs := args
		if createName == "" && createPrefix == "" && createSource == "" {
			normalizedArgs = make([]string, len(args))
			for i, arg := range args {
				normalized := utils.NormalizeNameVersion(arg)
				// Bare name (no slash) → expand to <default_distro>/<name>
				// e.g. "igv" → "ubuntu24/igv"  →  ubuntu24--igv.sqf
				if !strings.Contains(normalized, "/") && config.Global.DefaultDistro != "" {
					expanded := config.Global.DefaultDistro + "/" + normalized
					utils.PrintNote("Expanding '%s' to '%s'", normalized, expanded)
					normalized = expanded
				}
				normalizedArgs[i] = normalized
			}
		}

		// 6. Handle --remote flag (CLI flag or config)
		build.PreferRemote = createRemote || config.Global.PreferRemote

		// 7. Announce update mode
		if createUpdate {
			utils.PrintNote("Update mode: existing overlays will be rebuilt.")
		}

		// 8. Execute create based on mode
		if createSource != "" {
			// Mode: --source (external image like docker://ubuntu)
			runCreateFromSource(ctx)
		} else if createPrefix != "" {
			// Mode: --prefix (with --file for YAML/def/sh)
			runCreateWithPrefix(ctx)
		} else if createName != "" {
			// Mode: --name (multiple packages or YAML into one sqf)
			// Pass original args for conda package names (not normalized)
			runCreateWithName(ctx, args)
		} else {
			// Mode: Default (each package gets its own sqf via BuildObject)
			runCreatePackages(ctx, normalizedArgs)
		}
	},
}

func init() {
	rootCmd.AddCommand(createCmd)

	// Register Flags
	f := createCmd.Flags()
	f.StringVarP(&createName, "name", "n", "", "Custom name for the resulting overlay file")
	f.StringVarP(&createPrefix, "prefix", "p", "", "Custom prefix path for the overlay file")
	f.StringVarP(&createFile, "file", "f", "", "Path to definition file (.yaml, .sh, .def)")
	f.StringVarP(&createBaseImage, "base-image", "b", "", "Base image to use instead of default")
	f.StringVarP(&createSource, "source", "s", "", "Remote source URI (e.g., docker://ubuntu:22.04)")
	f.StringVar(&createTempSize, "temp-size", "20G", "Size of temporary overlay")
	f.BoolVar(&createRemote, "remote", false, "Remote build scripts take precedence over local")
	f.BoolVarP(&createUpdate, "update", "u", false, "Rebuild overlays even if they already exist (atomic .new swap)")

	// Compression flags: create a bool flag for each known option
	compFlags = make(map[string]*bool)
	for _, opt := range config.CompressOptions {
		compFlags[opt.Name] = f.Bool(opt.Name, false, opt.Description)
	}
}

// getWritableImagesDir returns the writable images directory or exits with an error
func getWritableImagesDir() string {
	dir, err := config.GetWritableImagesDir()
	if err != nil {
		ExitWithError("No writable images directory found: %v", err)
	}
	return dir
}

// runCreatePackages creates separate sqf files for each package using BuildGraph
// Example: condatainer create samtools/1.16 bcftools/1.15
func runCreatePackages(ctx context.Context, packages []string) {
	imagesDir := getWritableImagesDir()

	buildObjects := make([]build.BuildObject, 0, len(packages))
	for _, pkg := range packages {
		bo, err := build.NewBuildObject(pkg, false, imagesDir, config.GetWritableTmpDir(), createUpdate)
		if err != nil {
			ExitWithError("Failed to create build object for %s: %v", pkg, err)
		}
		utils.PrintDebug("[CREATE] BuildObject created:\n%s", bo)
		buildObjects = append(buildObjects, bo)
	}

	graph, err := build.NewBuildGraph(buildObjects, imagesDir, config.GetWritableTmpDir(), config.Global.SubmitJob, createUpdate)
	if err != nil {
		ExitWithError("Failed to create build graph: %v", err)
	}

	if err := graph.Run(ctx); err != nil {
		exitOnBuildError(err)
	}
	// If jobs were submitted to the scheduler, exit with a distinct code so downstream tooling
	// can detect that overlays will be created asynchronously by scheduler jobs.
	ExitIfJobsSubmitted(graph)
}

// runCreateWithName creates a single sqf with multiple packages or from YAML
// Example: condatainer create -n myenv nvim nodejs
// Example: condatainer create -n myenv -f environment.yml
func runCreateWithName(ctx context.Context, packages []string) {
	imagesDir := getWritableImagesDir()

	normalizedName := utils.NormalizeNameVersion(createName)
	slashCount := strings.Count(normalizedName, "/")
	if slashCount > 1 {
		ExitWithError("--name cannot contain more than one '/'")
	}

	// Check if already exists (search all paths), skip only when not updating
	if !createUpdate {
		searchName := strings.ReplaceAll(normalizedName, "/", "--") + ".sqf"
		if existingPath, err := config.FindImage(searchName); err == nil {
			utils.PrintMessage("Overlay %s already exists at %s. Skipping creation.",
				utils.StyleName(filepath.Base(existingPath)), utils.StylePath(existingPath))
			return
		}
	}

	utils.PrintDebug("[CREATE] Creating overlay with name: %s", createName)

	// Create a conda BuildObject with buildSource set appropriately
	// The buildSource field will contain either:
	// - Path to YAML file (if -f flag used)
	// - Comma-separated package list (if packages provided)
	var buildSource string
	if createFile != "" {
		// YAML file mode
		if !strings.HasSuffix(createFile, ".yml") && !strings.HasSuffix(createFile, ".yaml") {
			ExitWithError("File must be .yml or .yaml for conda environments")
		}
		buildSource, _ = filepath.Abs(createFile)
	} else if len(packages) > 0 {
		// Multiple packages mode - join with commas
		buildSource = strings.Join(packages, ",")
	}

	// Create CondaBuildObject using the new factory function
	bo, err := build.NewCondaObjectWithSource(normalizedName, buildSource, imagesDir, config.GetWritableTmpDir(), createUpdate)
	if err != nil {
		ExitWithError("Failed to create build object: %v", err)
	}

	if err := bo.Build(ctx, false); err != nil {
		ExitWithError("Build failed: %v", err)
	}
}

// runCreateWithPrefix creates a sqf from external source file (.sh, .def, .yml)
// Example: condatainer create -p myprefix -f environment.yml
// Example: condatainer create -p myprefix -f build.sh
func runCreateWithPrefix(ctx context.Context) {
	absPrefix, _ := filepath.Abs(createPrefix)
	// Use the directory from prefix path as output directory
	outputDir := filepath.Dir(absPrefix)

	if createFile == "" {
		ExitWithError("--prefix requires --file to be specified")
	}

	if !utils.FileExists(createFile) {
		ExitWithError("File %s not found", utils.StylePath(createFile))
	}

	utils.PrintDebug("[CREATE] Creating overlay with prefix: %s", createPrefix)

	// Determine file type and create appropriate BuildObject
	if strings.HasSuffix(createFile, ".yml") || strings.HasSuffix(createFile, ".yaml") {
		// YAML conda environment - use NewCondaObjectWithSource
		absFile, _ := filepath.Abs(createFile)
		bo, err := build.NewCondaObjectWithSource(filepath.Base(absPrefix), absFile, outputDir, config.GetWritableTmpDir(), createUpdate)
		if err != nil {
			ExitWithError("Failed to create build object: %v", err)
		}
		if err := bo.Build(ctx, false); err != nil {
			ExitWithError("Build failed: %v", err)
		}
	} else if strings.HasSuffix(createFile, ".sh") || strings.HasSuffix(createFile, ".bash") || strings.HasSuffix(createFile, ".def") {
		// Shell script or apptainer def file
		isApptainer := strings.HasSuffix(createFile, ".def")
		absFile, _ := filepath.Abs(createFile)
		bo, err := build.FromExternalSource(absPrefix, absFile, isApptainer, outputDir, config.GetWritableTmpDir())
		if err != nil {
			ExitWithError("Failed to create build object from %s: %v", createFile, err)
		}

		buildObjects := []build.BuildObject{bo}
		graph, err := build.NewBuildGraph(buildObjects, outputDir, config.GetWritableTmpDir(), config.Global.SubmitJob, createUpdate)
		if err != nil {
			ExitWithError("Failed to create build graph: %v", err)
		}

		if err := graph.Run(ctx); err != nil {
			exitOnBuildError(err)
		}
		// If jobs were submitted to the scheduler, exit with a distinct code so downstream tooling
		// can detect that overlays will be created asynchronously by scheduler jobs.
		ExitIfJobsSubmitted(graph)
	} else {
		ExitWithError("File must be .yml, .yaml, .sh, .bash, or .def")
	}
}

// runCreateFromSource creates a sqf from an external source (def file or remote URI)
// Example: condatainer create --source docker://ubuntu:22.04 -n myubuntu
func runCreateFromSource(ctx context.Context) {
	imagesDir := getWritableImagesDir()

	if createPrefix == "" && createName == "" {
		ExitWithError("--source requires either --name or --prefix")
	}

	var targetPrefix string
	if createPrefix != "" {
		targetPrefix, _ = filepath.Abs(createPrefix)
	} else {
		normalizedName := utils.NormalizeNameVersion(createName)
		if strings.Count(normalizedName, "/") > 1 {
			ExitWithError("--name cannot contain more than one '/'")
		}
		fileName := strings.ReplaceAll(normalizedName, "/", "--")
		targetPrefix = filepath.Join(imagesDir, fileName)
	}

	source := createSource
	isRemote := strings.Contains(source, "://")
	if !isRemote {
		source, _ = filepath.Abs(source)
		if !utils.FileExists(source) {
			ExitWithError("Source %s not found", utils.StylePath(source))
		}
	}

	isApptainer := strings.HasSuffix(source, ".def") || isRemote
	targetOverlayPath := targetPrefix + ".sqf"

	utils.PrintMessage("Creating overlay %s from %s", filepath.Base(targetOverlayPath), utils.StylePath(source))

	bo, err := build.FromExternalSource(targetPrefix, source, isApptainer, imagesDir, config.GetWritableTmpDir())
	if err != nil {
		ExitWithError("Failed to create build object from %s: %v", source, err)
	}

	buildObjects := []build.BuildObject{bo}
	graph, err := build.NewBuildGraph(buildObjects, imagesDir, config.GetWritableTmpDir(), config.Global.SubmitJob, createUpdate)
	if err != nil {
		ExitWithError("Failed to create build graph: %v", err)
	}

	if err := graph.Run(ctx); err != nil {
		exitOnBuildError(err)
	}
	// If jobs were submitted to the scheduler, exit with a distinct code so downstream tooling
	// can detect that overlays will be created asynchronously by scheduler jobs.
	ExitIfJobsSubmitted(graph)
}

func exitOnBuildError(err error) {
	if errors.Is(err, build.ErrBuildCancelled) ||
		strings.Contains(err.Error(), "signal: killed") ||
		strings.Contains(err.Error(), "signal: interrupt") ||
		strings.Contains(err.Error(), "context canceled") {
		// Suppress output as inner layers handle the "cancelled" messaging
	} else {
		utils.PrintError("Build failed: %v", err)
	}
	os.Exit(ExitCodeError)
}
