package cmd

import (
	"context"
	"errors"
	"fmt"
	"io"
	"os"
	"path/filepath"
	"slices"
	"strings"

	"github.com/Justype/condatainer/internal/build"
	"github.com/Justype/condatainer/internal/config"
	"github.com/Justype/condatainer/internal/container"
	"github.com/Justype/condatainer/internal/utils"
	"github.com/chzyer/readline"
	"github.com/spf13/cobra"
	"github.com/spf13/pflag"
)

// Variables to hold flag values
var (
	createName          string
	createPrefix        string
	createFile          string
	createSource        string
	createTempSize      string
	createBlockSize     string
	createDataBlockSize string
	createChannels      []string
	createRemote        bool
	createUpdate        bool
	createUseTmpOverlay bool

	// compression flags are generated dynamically from config.CompressOptions
	compFlags map[string]*bool

	// buildFlagNames is the set of flags shown under "Build Flags:" in help.
	buildFlagNames = map[string]bool{
		"temp-size": true, "block-size": true, "data-block-size": true, "use-tmp-overlay": true,
	}
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
	Example: `  condatainer create samtools/1.22                       # Create from build script
  condatainer create python=3.11 numpy -n myenv          # Create conda environment
  condatainer create python=3.11 numpy -p /scratch/myenv # Create conda env at custom path
  condatainer create -f environment.yml -p myenv         # Create from conda file with prefix
  condatainer create --source /data -p dataset           # Convert directory to overlay`,
	Run: func(cmd *cobra.Command, args []string) {
		ctx := cmd.Context()
		// 1. Validation Logic
		if len(args) == 0 && createFile == "" && createSource == "" {
			ExitWithUsageError("At least one of [packages], --file, or --source must be provided.")
		}
		if createPrefix != "" && createName != "" {
			ExitWithUsageError("Cannot use both --prefix and --name at the same time.")
		}
		if createPrefix != "" {
			baseName := strings.TrimSuffix(filepath.Base(createPrefix), ".sqf")
			if strings.Contains(baseName, "--") {
				ExitWithUsageError("--prefix name cannot contain '--' (reserved name/version separator)")
			}
		}
		if createSource != "" && createName == "" && createPrefix == "" {
			ExitWithUsageError("When using --source, either --name or --prefix must be provided.")
		}
		if createFile != "" && createPrefix == "" {
			createPrefix = createFile[:len(createFile)-len(filepath.Ext(createFile))]
		}
		if createPrefix != "" && createFile == "" && len(args) == 0 && createSource == "" {
			ExitWithUsageError("--prefix requires either packages, --file, or --source to be specified.")
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

		// 4. Override channels if -c was provided
		if len(createChannels) > 0 {
			config.Global.Build.Channels = createChannels
		}

		// 5. Handle temp size
		if createTempSize != "" {
			sizeMB, err := utils.ParseSizeToMB(createTempSize)
			if err != nil {
				ExitWithError("Invalid temp size: %v", err)
			}
			config.Global.Build.TmpSizeMB = sizeMB
		}

		// 5b. Handle block sizes
		if createBlockSize != "" {
			if !config.IsValidBlockSize(createBlockSize) {
				ExitWithError("Invalid --block-size %q: must be a power of two between 4096 and 1M (e.g. 64k, 128k, 512k, 1m)", createBlockSize)
			}
			config.Global.Build.BlockSize = createBlockSize
		}
		if createDataBlockSize != "" {
			if !config.IsValidBlockSize(createDataBlockSize) {
				ExitWithError("Invalid --data-block-size %q: must be a power of two between 4096 and 1M (e.g. 64k, 128k, 512k, 1m)", createDataBlockSize)
			}
			config.Global.Build.DataBlockSize = createDataBlockSize
		}

		// 5c. Handle use-tmp-overlay
		if createUseTmpOverlay {
			config.Global.Build.UseTmpOverlay = true
		}

		// 6. Normalize package names (only for build-script mode, not for conda/prefix/source modes)
		normalizedArgs := args
		if createName == "" && createPrefix == "" && createSource == "" {
			normalizedArgs = make([]string, len(args))
			for i, arg := range args {
				normalized := utils.NormalizeNameVersion(arg)
				// Bare name (no slash) → expand to <default_distro>/<name> only when
				// a build script exists for that distro/name combination.
				// e.g. "igv" → "ubuntu24/igv"  →  ubuntu24--igv.sqf
				// Without a script, keep the original name so the user gets a clear error.
				if !strings.Contains(normalized, "/") && !strings.Contains(normalized, "::") && config.Global.DefaultDistro != "" {
					candidate := config.Global.DefaultDistro + "/" + normalized
					if _, found := build.FindBuildScript(candidate); found {
						utils.PrintNote("Expanding '%s' to '%s'", normalized, candidate)
						normalized = candidate
					}
				}
				normalizedArgs[i] = normalized
			}
		}

		// 7. Handle --remote flag (CLI flag or config)
		build.PreferRemote = createRemote || config.Global.PreferRemote

		// 8. Announce update mode
		if createUpdate {
			utils.PrintNote("Update mode: existing overlays will be rebuilt.")
		}

		// 9. Execute create based on mode
		if createSource != "" {
			// Mode: --source (external image like docker://ubuntu)
			runCreateFromSource(ctx)
		} else if createPrefix != "" && len(args) > 0 {
			// Mode: --prefix + packages (conda env at custom path, like conda create -p)
			runCreateWithPrefixAndPackages(ctx, args)
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
	f.StringVarP(&createSource, "source", "s", "", "Remote source URI (e.g., docker://ubuntu:22.04)")
	f.StringVar(&createTempSize, "temp-size", "20G", "Size of temporary overlay")
	f.StringVar(&createBlockSize, "block-size", "", "SquashFS block size for app/env/external overlays (e.g. 128k, 512k)")
	f.StringVar(&createDataBlockSize, "data-block-size", "", "SquashFS block size for data overlays (e.g. 512k, 1m)")
	f.StringArrayVarP(&createChannels, "channel", "c", nil, "Conda channel to use (overrides config; repeatable)")
	f.BoolVar(&createRemote, "remote", false, "Remote build scripts take precedence over local")
	f.BoolVarP(&createUpdate, "update", "u", false, "Rebuild overlays even if they already exist (atomic .new swap)")
	f.BoolVar(&createUseTmpOverlay, "use-tmp-overlay", false, "Use a temporary overlay instead of a temp directory")

	// Compression flags: create a bool flag for each known option
	compFlags = make(map[string]*bool)
	for _, opt := range config.CompressOptions {
		compFlags[opt.Name] = f.Bool(opt.Name, false, opt.Description)
	}

	blockSizeCompletion := func(cmd *cobra.Command, args []string, toComplete string) ([]string, cobra.ShellCompDirective) {
		return config.BlockSizeCompletions, cobra.ShellCompDirectiveNoFileComp
	}
	createCmd.RegisterFlagCompletionFunc("block-size", blockSizeCompletion)      //nolint:errcheck
	createCmd.RegisterFlagCompletionFunc("data-block-size", blockSizeCompletion) //nolint:errcheck

	// Mark compression flags as build flags
	for _, opt := range config.CompressOptions {
		buildFlagNames[opt.Name] = true
	}

	// Custom usage: two labeled sections — "Flags:" and "Build Flags:"
	createCmd.SetUsageFunc(func(cmd *cobra.Command) error {
		fmt.Fprintf(cmd.OutOrStderr(), "Usage:\n  %s\n", cmd.UseLine())
		if cmd.HasExample() {
			fmt.Fprintf(cmd.OutOrStderr(), "\nExamples:\n%s\n", cmd.Example)
		}
		general := pflag.NewFlagSet("", pflag.ContinueOnError)
		build := pflag.NewFlagSet("", pflag.ContinueOnError)
		cmd.LocalFlags().VisitAll(func(fl *pflag.Flag) {
			if buildFlagNames[fl.Name] {
				build.AddFlag(fl)
			} else {
				general.AddFlag(fl)
			}
		})
		if general.HasFlags() {
			fmt.Fprintf(cmd.OutOrStderr(), "\nFlags:\n%s", general.FlagUsages())
		}
		if build.HasFlags() {
			fmt.Fprintf(cmd.OutOrStderr(), "\nBuild Flags:\n%s", build.FlagUsages())
		}
		if cmd.HasAvailableInheritedFlags() {
			fmt.Fprintf(cmd.OutOrStderr(), "\nGlobal Flags:\n%s", cmd.InheritedFlags().FlagUsages())
		}
		return nil
	})
}

// getWritableImagesDir returns the writable images directory or exits with an error
func getWritableImagesDir() string {
	dir, err := config.GetWritableImagesDir()
	if err != nil {
		ExitWithError("No writable images directory found: %v", err)
	}
	return dir
}

// readLineWithCompletion reads a line from stdin with tab-completion over completions.
// Ctrl-C / EOF and context cancellation all return context.Canceled.
func readLineWithCompletion(ctx context.Context, prompt string, completions []string) (string, error) {
	items := make([]readline.PrefixCompleterInterface, len(completions))
	for i, v := range completions {
		items[i] = readline.PcItem(v)
	}
	rl, err := readline.NewEx(&readline.Config{
		Prompt:       prompt,
		AutoComplete: readline.NewPrefixCompleter(items...),
	})
	if err != nil {
		return "", err
	}
	defer rl.Close()

	type result struct {
		s   string
		err error
	}
	ch := make(chan result, 1)
	go func() {
		s, err := rl.Readline()
		ch <- result{strings.TrimSpace(s), err}
	}()

	select {
	case <-ctx.Done():
		rl.Close()
		return "", context.Canceled
	case r := <-ch:
		if r.err == readline.ErrInterrupt || r.err == io.EOF {
			return "", context.Canceled
		}
		return r.s, r.err
	}
}

// resolveTemplateInteractively prompts the user to choose a value for each placeholder
// in a PL template script and returns the interpolated concrete name (from #TARGET:).
// When --yes is set, defaults are used without prompting.
func resolveTemplateInteractively(ctx context.Context, info build.ScriptInfo) (string, error) {
	utils.PrintMessage("Placeholder template: %s", info.Name)
	if info.Whatis != "" {
		utils.PrintMessage("%s", utils.StyleHint(info.Whatis))
	}
	if info.TargetTemplate != "" {
		fmt.Fprintf(os.Stdout, "Target: %s\n", utils.HighlightTemplatePlaceholders(info.TargetTemplate))
	}
	chosenVars := make(map[string]string, len(info.PLOrder))

	// Build per-placeholder installed defaults from single-slash tool #DEP: patterns.
	installedDefaults := map[string]string{}
	if overlays, err := container.InstalledOverlays(); err == nil {
		installedVals := map[string][]string{}
		for _, dep := range info.Deps {
			if !strings.Contains(dep, "{") || strings.Count(dep, "/") != 1 {
				continue
			}
			for name := range overlays {
				if vars, ok := utils.MatchTemplateTarget(dep, name); ok {
					for k, v := range vars {
						installedVals[k] = append(installedVals[k], v)
					}
				}
			}
		}
		for k, vals := range installedVals {
			installedDefaults[k] = utils.SortVersionsDescending(vals)[0]
		}
	}

	for _, key := range info.PLOrder {
		vals, ok := info.PL[key]
		if !ok {
			continue
		}

		// Separate concrete values from "*"
		var concrete []string
		hasOpen := false
		for _, v := range vals {
			if v == "*" {
				hasOpen = true
			} else {
				concrete = append(concrete, v)
			}
		}

		// Determine the default: prefer latest installed, fall back to latest available.
		var defaultVal string
		if len(concrete) > 0 {
			defaultVal = concrete[0]
		}
		if iv, ok := installedDefaults[key]; ok {
			if hasOpen || slices.Contains(concrete, iv) {
				defaultVal = iv
			}
		}

		// Build the prompt string
		var prompt string
		n := len(concrete)
		switch {
		case hasOpen && n > 0:
			if n > 8 {
				prompt = fmt.Sprintf("  %s [suggested: %s-%s, or any value] (default: %s): ",
					key, concrete[n-1], concrete[0], defaultVal)
			} else {
				prompt = fmt.Sprintf("  %s [suggested: %s, or any value] (default: %s): ",
					key, strings.Join(concrete, ", "), defaultVal)
			}
		case hasOpen && n == 0:
			prompt = fmt.Sprintf("  %s [any value]: ", key)
		case n > 8:
			prompt = fmt.Sprintf("  %s [%s-%s] (default: %s): ",
				key, concrete[n-1], concrete[0], defaultVal)
		default:
			prompt = fmt.Sprintf("  %s [%s] (default: %s): ",
				key, strings.Join(concrete, ", "), defaultVal)
		}

		if utils.ShouldAnswerYes() {
			if defaultVal != "" {
				fmt.Printf("%s%s\n", prompt, defaultVal)
				chosenVars[key] = defaultVal
			}
			continue
		}

		for {
			input, err := readLineWithCompletion(ctx, prompt, concrete)
			if err != nil {
				return "", err
			}

			// Empty input → use default
			if input == "" {
				if defaultVal == "" {
					utils.PrintWarning("No default available for %q — please enter a value.", key)
					continue
				}
				chosenVars[key] = defaultVal
				break
			}

			// For closed lists, validate against known values
			if !hasOpen {
				valid := false
				for _, v := range concrete {
					if v == input {
						valid = true
						break
					}
				}
				if !valid {
					utils.PrintWarning("Invalid value %q for %s. Valid values: %s",
						input, key, strings.Join(concrete, ", "))
					continue
				}
			}

			chosenVars[key] = input
			break
		}
	}

	concrete := utils.InterpolateVars(info.TargetTemplate, chosenVars)
	fmt.Fprintf(os.Stdout, "  → Creating %s\n", concrete)
	return concrete, nil
}

// runCreatePackages creates separate sqf files for each package using BuildGraph
// Example: condatainer create samtools/1.16 bcftools/1.15
func runCreatePackages(ctx context.Context, packages []string) {
	imagesDir := getWritableImagesDir()

	buildObjects := make([]*build.BuildObject, 0, len(packages))
	for _, pkg := range packages {
		// Resolve template scripts interactively before building
		if info, found := build.FindBuildScript(pkg); found && info.IsTemplate {
			resolved, err := resolveTemplateInteractively(ctx, info)
			if err != nil {
				ExitWithError("Template resolution cancelled for %s: %v", pkg, err)
			}
			pkg = resolved
		}

		bo, err := build.NewBuildObject(ctx, pkg, false, imagesDir, config.GetWritableTmpDir(), createUpdate)
		if err != nil {
			ExitWithError("Failed to create build object for %s: %v", pkg, err)
		}
		utils.PrintDebug("[CREATE] BuildObject created:\n%s", bo)
		buildObjects = append(buildObjects, bo)
	}

	graph, err := build.NewBuildGraph(ctx, buildObjects, imagesDir, config.GetWritableTmpDir(), config.Global.SubmitJob, createUpdate)
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
		if !utils.IsYaml(createFile) {
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

	if !utils.FileExists(createFile) {
		ExitWithError("File %s not found", utils.StylePath(createFile))
	}

	utils.PrintDebug("[CREATE] Creating overlay with prefix: %s", createPrefix)

	// Determine file type and create appropriate BuildObject
	if utils.IsYaml(createFile) {
		// YAML conda environment - use NewCondaObjectWithSource
		absFile, _ := filepath.Abs(createFile)
		bo, err := build.NewCondaObjectWithSource(filepath.Base(absPrefix), absFile, outputDir, outputDir, createUpdate)
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
		bo, err := build.FromExternalSource(ctx, absPrefix, absFile, isApptainer, outputDir)
		if err != nil {
			ExitWithError("Failed to create build object from %s: %v", createFile, err)
		}

		buildObjects := []*build.BuildObject{bo}
		graph, err := build.NewBuildGraph(ctx, buildObjects, outputDir, config.GetWritableTmpDir(), config.Global.SubmitJob, createUpdate)
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

// runCreateWithPrefixAndPackages creates a conda sqf from packages at a custom prefix path.
// Example: condatainer create python=3.11 numpy -p /scratch/myenv
func runCreateWithPrefixAndPackages(ctx context.Context, packages []string) {
	absPrefix, _ := filepath.Abs(createPrefix)
	outputDir := filepath.Dir(absPrefix)
	baseName := filepath.Base(absPrefix)
	buildSource := strings.Join(packages, ",")

	bo, err := build.NewCondaObjectWithSource(baseName, buildSource, outputDir, outputDir, createUpdate)
	if err != nil {
		ExitWithError("Failed to create build object: %v", err)
	}
	if err := bo.Build(ctx, false); err != nil {
		ExitWithError("Build failed: %v", err)
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

	bo, err := build.FromExternalSource(ctx, targetPrefix, source, isApptainer, imagesDir)
	if err != nil {
		ExitWithError("Failed to create build object from %s: %v", source, err)
	}

	buildObjects := []*build.BuildObject{bo}
	graph, err := build.NewBuildGraph(ctx, buildObjects, imagesDir, config.GetWritableTmpDir(), config.Global.SubmitJob, createUpdate)
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
