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
	"sync"

	"github.com/Justype/condatainer/catalog"
	"github.com/Justype/condatainer/internal/build"
	"github.com/Justype/condatainer/internal/config"
	"github.com/Justype/condatainer/internal/runtime/apptainer"
	"github.com/Justype/condatainer/internal/runtime/container"
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
	createFrom          string
	createAppTmpOvlSize string
	createBlockSize     string
	createDataBlockSize string
	createChannels      []string
	createSources       []string
	createUpdate        bool
	createAppTmpOverlay bool
	createAlwaysSubmit  bool

	// compression flags are generated dynamically from config.CompressOptions
	compFlags     map[string]*bool
	compFlagNames map[string]bool

	// buildFlagNames is the set of flags shown under "Build Flags:" in help.
	buildFlagNames = map[string]bool{
		"app-tmp-overlay": true, "app-tmp-overlay-size": true, "block-size": true,
		"data-block-size": true, "always-submit": true, "no-submit": true,
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
	Use:     "create [flags] [packages...]",
	Aliases: []string{"install", "i", "build"},
	Short:   "Create a new SquashFS overlay",
	Long: `Create a new SquashFS overlay using available build scripts or Conda packages.

Submitted build jobs exit with code 3 (useful for scripts).`,
	Example: `  condatainer create orad/2.7.0                   # Create from build script
  condatainer create -n nvim nvim nodejs          # Create conda env
  condatainer create matplotlib pandas  -p /path  # Create conda env at custom path
  condatainer create -f environment.yml -p myenv  # Create from conda file with prefix
  condatainer create --source lab star/2.7.11b     # Resolve recipes only from lab
  condatainer create --from docker://ubuntu:22.04 -p ubuntu  # Build from a container image`,
	PreRunE: func(cmd *cobra.Command, args []string) error {
		return config.SelectSources(createSources)
	},
	Run: func(cmd *cobra.Command, args []string) {
		ctx := cmd.Context()
		// 1. Validation Logic
		if len(args) == 0 && createFile == "" && createFrom == "" {
			ExitWithError("At least one of [packages], --file, or --from must be provided.")
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
		if createFrom != "" && createName == "" && createPrefix == "" {
			ExitWithError("When using --from, either --name or --prefix must be provided.")
		}
		if derived := derivePrefixFromFile(createFile, createPrefix, createName); derived != "" {
			createPrefix = derived
		}
		if createPrefix != "" && createFile == "" && len(args) == 0 && createFrom == "" {
			ExitWithError("--prefix requires either packages, --file, or --from to be specified.")
		}

		// 2. Apptainer runs every build, so fail here rather than after
		// resolution has already fetched recipes. The base image is not checked:
		// it is an implicit prerequisite of the plan, built with everything else.
		if err := apptainer.EnsureApptainer(); err != nil {
			ExitWithError("%v", err)
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

		// 5. Handle the app build overlay size. Only when the flag was given:
		// its default would otherwise outrank build.app_tmp_overlay_size.
		if cmd.Flags().Changed("app-tmp-overlay-size") {
			sizeMB, err := utils.ParseSizeToMB(createAppTmpOvlSize)
			if err != nil {
				ExitWithError("Invalid --app-tmp-overlay-size: %v", err)
			}
			config.Global.Build.AppTmpOverlaySizeMB = sizeMB
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

		// 5c. Handle the app build overlay mode
		if createAppTmpOverlay {
			config.Global.Build.AppTmpOverlay = true
		}

		// 6. Normalize package names (only for build-script mode, not for conda/prefix/source modes)
		normalizedArgs := args
		if createName == "" && createPrefix == "" && createFrom == "" {
			normalizedArgs = make([]string, len(args))
			for i, arg := range args {
				normalized, expanded := expandBareName(cmd.Context(), arg)
				if expanded {
					utils.PrintNote("Expanding '%s' to '%s'", catalog.Normalize(arg), normalized)
				}
				normalizedArgs[i] = normalized
			}
		}

		// 7b. Handle --always-submit flag
		if createAlwaysSubmit {
			config.Global.Build.AlwaysSubmit = true
		}

		// 8. Announce update mode
		if createUpdate {
			utils.PrintNote("Update mode: existing overlays will be rebuilt.")
		}

		// 9. Execute create based on mode
		if createFrom != "" {
			// Mode: --from (external image like docker://ubuntu)
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
	f.StringVarP(&createName, "name", "n", "", "Custom name for the overlay")
	f.StringVarP(&createPrefix, "prefix", "p", "", "Custom prefix path for the overlay")
	f.StringVarP(&createFile, "file", "f", "", "Path to definition file (.yaml, .txt, .sh, .def)")
	f.StringVar(&createFrom, "from", "", "Build from an external image URI (e.g., docker://ubuntu:22.04)")
	f.StringVar(&createBlockSize, "block-size", "", "SquashFS block size of app/external overlays (e.g. 256k)")
	f.StringVar(&createDataBlockSize, "data-block-size", "", "SquashFS block size of data overlays (e.g. 512k, 1m)")
	f.StringArrayVarP(&createChannels, "channel", "c", nil, "Conda channel to use (overrides config; repeatable)")
	f.StringArrayVarP(&createSources, "source", "s", nil,
		"Use only this configured recipe source, in flag order (repeatable)")
	f.BoolVarP(&createUpdate, "update", "u", false, "Rebuild overlays even if they already exist")
	f.BoolVar(&createAppTmpOverlay, "app-tmp-overlay", false, "Assemble an app build in a temporary ext3 overlay instead of host directories")
	f.StringVar(&createAppTmpOvlSize, "app-tmp-overlay-size", "20G", "Size of that temporary overlay")
	f.BoolVar(&createAlwaysSubmit, "always-submit", false, "Submit all builds as scheduler jobs, even no directives")
	f.BoolVar(&noSubmitMode, "no-submit", false, "Disable job submission (build locally)")

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
	createCmd.RegisterFlagCompletionFunc("source", sourceHandleCompletion)       //nolint:errcheck

	// Mark compression flags in their own section
	compFlagNames = make(map[string]bool, len(config.CompressOptions))
	for _, opt := range config.CompressOptions {
		compFlagNames[opt.Name] = true
	}

	// Custom usage: two labeled sections — "Flags:" and "Build Flags:"
	createCmd.SetUsageFunc(func(cmd *cobra.Command) error {
		fmt.Fprintf(cmd.OutOrStderr(), "Usage:\n  %s\n", cmd.UseLine())
		if len(cmd.Aliases) > 0 {
			fmt.Fprintf(cmd.OutOrStderr(), "\nAliases:\n  %s\n", cmd.NameAndAliases())
		}
		if cmd.HasExample() {
			fmt.Fprintf(cmd.OutOrStderr(), "\nExamples:\n%s\n", cmd.Example)
		}
		general := pflag.NewFlagSet("", pflag.ContinueOnError)
		build := pflag.NewFlagSet("", pflag.ContinueOnError)
		compress := pflag.NewFlagSet("", pflag.ContinueOnError)
		cmd.LocalFlags().VisitAll(func(fl *pflag.Flag) {
			if compFlagNames[fl.Name] {
				compress.AddFlag(fl)
			} else if buildFlagNames[fl.Name] {
				build.AddFlag(fl)
			} else {
				general.AddFlag(fl)
			}
		})
		if general.HasFlags() {
			fmt.Fprintf(cmd.OutOrStderr(), "\nFlags:\n%s", flagUsages(general))
		}
		if build.HasFlags() {
			fmt.Fprintf(cmd.OutOrStderr(), "\nBuild Flags:\n%s", flagUsages(build))
		}
		if compress.HasFlags() {
			fmt.Fprintf(cmd.OutOrStderr(), "\nCompression Flags:\n%s", flagUsages(compress))
		}
		if cmd.HasAvailableInheritedFlags() {
			fmt.Fprintf(cmd.OutOrStderr(), "\nGlobal Flags:\n%s", cmd.InheritedFlags().FlagUsages())
		}
		return nil
	})
}

// imagesDirNoteOnce keeps the destination note to one line per command, since the
// three create paths each resolve the directory independently.
var imagesDirNoteOnce sync.Once

// getWritableImagesDir returns the writable images directory or exits with an error.
// The destination and its data layer are reported once: the target depends on which
// directories happen to be writable, so "(app-root)" vs "(user)" is the difference
// between installing for everyone and installing only for yourself.
func getWritableImagesDir() string {
	dir, err := config.GetWritableImagesDir()
	if err != nil {
		ExitWithError("No writable images directory found: %v", err)
	}
	imagesDirNoteOnce.Do(func() {
		utils.PrintNote("Installing to %s (%s)", dir, config.ClassifyDataDir(dir))
	})
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
// in a PH template script and returns the interpolated concrete name (from #TARGET:).
// When --yes is set, defaults are used without prompting.
func resolveTemplateInteractively(ctx context.Context, info *catalog.Entry) (string, error) {
	utils.PrintMessage("Placeholder template: %s", info.Name)
	if info.Description != "" {
		utils.PrintMessage("%s", utils.StyleHint(info.Description))
	}
	if info.TargetTemplate != "" {
		fmt.Fprintf(os.Stdout, "Target: %s\n", utils.HighlightTemplatePlaceholders(info.TargetTemplate))
	}
	tmpl := catalog.NewTemplate(info.TargetTemplate)
	names := tmpl.Names()
	chosenVars := make(map[string]string, len(names))

	// Build per-placeholder installed defaults from single-slash tool #DEP: patterns.
	installedDefaults := map[string]string{}
	if overlays, err := container.InstalledOverlays(); err == nil {
		installedVals := map[string][]string{}
		for _, dep := range info.Deps {
			if !strings.Contains(dep, "{") || strings.Count(dep, "/") != 1 {
				continue
			}
			depTmpl := catalog.NewTemplate(dep)
			for name := range overlays {
				// The dep's tokens are the parent's placeholders, so its declared
				// values are what an installed name has to be one of.
				if vars, ok := depTmpl.Match(name, info.PH); ok {
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

	for _, key := range names {
		vals, ok := info.PH[key]
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

	concrete, err := tmpl.Fill(chosenVars)
	if err != nil {
		return "", err
	}
	fmt.Fprintf(os.Stdout, "  → Creating %s\n", concrete)
	return concrete, nil
}

// expandBareName first preserves an exact catalog name. If no exact entry
// exists, a name with at most one slash is tried below the configured base, so
// both "r" and "r/4.5.3" may find ubuntu24/r[/4.5.3]. If neither resolves, the
// original name is returned for the normal Conda fallback.
func expandBareName(ctx context.Context, nameVersion string) (string, bool) {
	normalized := catalog.Normalize(nameVersion)
	if normalized == "" || strings.Contains(normalized, "::") {
		return normalized, false
	}
	cat, err := config.OpenCatalog(ctx)
	if err != nil {
		return normalized, false
	}
	if _, found, err := cat.Lookup(ctx, normalized); err == nil && found {
		return normalized, false
	}
	if strings.Count(normalized, "/") > 1 {
		return normalized, false
	}
	base := config.ResolvedBase()
	if base == "" {
		return normalized, false
	}
	candidate := base + "/" + normalized
	if _, found, err := cat.Lookup(ctx, candidate); err == nil && found {
		return candidate, true
	}
	return normalized, false
}

// runCreatePackages creates separate sqf files for each package using BuildGraph
// Example: condatainer create samtools/1.16 bcftools/1.15
func runCreatePackages(ctx context.Context, packages []string) {
	imagesDir := getWritableImagesDir()

	buildObjects := make([]*build.BuildObject, 0, len(packages))
	for _, pkg := range packages {
		// A bare template name does not identify an artifact; ask which member.
		if cat, err := config.OpenCatalog(ctx); err == nil {
			if m, found, err := cat.Lookup(ctx, pkg); err == nil && found &&
				m.Entry.IsTemplate && len(m.Vars) == 0 {
				resolved, err := resolveTemplateInteractively(ctx, m.Entry)
				if err != nil {
					ExitWithError("Template resolution cancelled for %s: %v", pkg, err)
				}
				pkg = resolved
			}
		}

		bo, err := build.NewBuildObject(ctx, pkg, false, imagesDir, createUpdate)
		if err != nil {
			ExitWithError("Failed to create build object for %s: %v", pkg, err)
		}
		utils.PrintDebug("[CREATE] BuildObject created:\n%s", bo)
		buildObjects = append(buildObjects, bo)
	}

	graph, err := build.NewBuildGraph(ctx, buildObjects, imagesDir, config.Global.SubmitJob, createUpdate)
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

// derivePrefixFromFile returns the prefix --file implies, or "" when the user
// already named a target. The mode dispatch tests --prefix before --name, so a
// prefix derived while --name is set would silently win over it.
func derivePrefixFromFile(file, prefix, name string) string {
	if file == "" || prefix != "" || name != "" {
		return ""
	}
	return file[:len(file)-len(filepath.Ext(file))]
}

// normalizedTargetName returns --name in catalog form.
//
// Any depth is allowed. Restoring a project from a lockfile has to recreate the
// names the catalog itself uses, and a data image is several levels deep
// (grch38/star/2.7.11b/gencode47-101). The name survives the round trip through
// the filename either way: / becomes -- on the way out, and Normalize turns it
// back on the way in.
func normalizedTargetName() string {
	return catalog.Normalize(createName)
}

// isExternalBuildFile reports whether a file builds through the external-source
// path: a shell recipe or an Apptainer definition.
func isExternalBuildFile(path string) bool {
	return strings.HasSuffix(path, ".sh") ||
		strings.HasSuffix(path, ".bash") ||
		strings.HasSuffix(path, ".def")
}

// buildExternalSource builds one script or definition into targetPrefix, exiting
// on failure. outputDir holds the image and its scratch space.
func buildExternalSource(ctx context.Context, targetPrefix, source string, isApptainer bool, outputDir string) {
	bo, err := build.FromExternalSource(ctx, targetPrefix, source, isApptainer, outputDir, createUpdate)
	if err != nil {
		ExitWithError("Failed to create build object from %s: %v", source, err)
	}

	graph, err := build.NewBuildGraph(ctx, []*build.BuildObject{bo}, outputDir,
		config.Global.SubmitJob, createUpdate)
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

// runCreateWithName creates a single sqf with multiple packages or from a file
// Example: condatainer create -n myenv nvim nodejs
// Example: condatainer create -n myenv -f environment.yml
// Example: condatainer create -n myenv -f build.sh
func runCreateWithName(ctx context.Context, packages []string) {
	imagesDir := getWritableImagesDir()

	normalizedName := normalizedTargetName()

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

	// A script or definition is not a conda input. It builds the same way --prefix
	// builds one, targeting the managed images dir instead of a path the user typed.
	if createFile != "" && !utils.IsCondaFile(createFile) {
		if !isExternalBuildFile(createFile) {
			ExitWithError("File must be .yml, .yaml, .txt, .sh, .bash, or .def")
		}
		absFile, _ := filepath.Abs(createFile)
		targetPrefix := filepath.Join(imagesDir, strings.ReplaceAll(normalizedName, "/", "--"))
		buildExternalSource(ctx, targetPrefix, absFile, strings.HasSuffix(createFile, ".def"), imagesDir)
		return
	}

	// Create a conda BuildObject with buildSource set appropriately
	// The buildSource field will contain either:
	// - Path to YAML file (if -f flag used)
	// - Comma-separated package list (if packages provided)
	var buildSource string
	if createFile != "" {
		buildSource, _ = filepath.Abs(createFile)
	} else if len(packages) > 0 {
		// Multiple packages mode - join with commas
		buildSource = strings.Join(packages, ",")
	}

	// Create CondaBuildObject using the new factory function
	bo, err := build.NewCondaObjectWithSource(normalizedName, buildSource, imagesDir, createUpdate)
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
	if utils.IsCondaFile(createFile) {
		// Conda env/spec file - use NewCondaObjectWithSource
		absFile, _ := filepath.Abs(createFile)
		bo, err := build.NewCondaObjectWithSource(filepath.Base(absPrefix), absFile, outputDir, createUpdate)
		if err != nil {
			ExitWithError("Failed to create build object: %v", err)
		}
		if err := bo.Build(ctx, false); err != nil {
			ExitWithError("Build failed: %v", err)
		}
	} else if isExternalBuildFile(createFile) {
		// Shell script or apptainer def file
		absFile, _ := filepath.Abs(createFile)
		buildExternalSource(ctx, absPrefix, absFile, strings.HasSuffix(createFile, ".def"), outputDir)
	} else {
		ExitWithError("File must be .yml, .yaml, .txt, .sh, .bash, or .def")
	}
}

// runCreateWithPrefixAndPackages creates a conda sqf from packages at a custom prefix path.
// Example: condatainer create python=3.11 numpy -p /scratch/myenv
func runCreateWithPrefixAndPackages(ctx context.Context, packages []string) {
	absPrefix, _ := filepath.Abs(createPrefix)
	outputDir := filepath.Dir(absPrefix)
	baseName := filepath.Base(absPrefix)
	buildSource := strings.Join(packages, ",")

	bo, err := build.NewCondaObjectWithSource(baseName, buildSource, outputDir, createUpdate)
	if err != nil {
		ExitWithError("Failed to create build object: %v", err)
	}
	if err := bo.Build(ctx, false); err != nil {
		ExitWithError("Build failed: %v", err)
	}
}

// runCreateFromSource creates a sqf from an external source (def file or remote URI)
// Example: condatainer create --from docker://ubuntu:22.04 -n myubuntu
func runCreateFromSource(ctx context.Context) {
	imagesDir := getWritableImagesDir()

	if createPrefix == "" && createName == "" {
		ExitWithError("--from requires either --name or --prefix")
	}

	var targetPrefix string
	if createPrefix != "" {
		targetPrefix, _ = filepath.Abs(createPrefix)
	} else {
		fileName := strings.ReplaceAll(normalizedTargetName(), "/", "--")
		targetPrefix = filepath.Join(imagesDir, fileName)
	}

	source := createFrom
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

	buildExternalSource(ctx, targetPrefix, source, isApptainer, imagesDir)
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
