package cmd

import (
	"context"
	"errors"
	"fmt"
	"os"
	"path/filepath"
	"strings"

	"github.com/Justype/condatainer/internal/build"
	"github.com/Justype/condatainer/internal/config"
	"github.com/Justype/condatainer/internal/scheduler"
	"github.com/Justype/condatainer/internal/utils"
	"github.com/spf13/cobra"
)

var (
	checkAutoInstall     bool
	checkParseModuleLoad bool
	checkRemote          bool
)

var scriptCheckCmd = &cobra.Command{
	Use:   "check <script|dir> [script|dir ...]",
	Short: "Check if the dependencies of a script are installed",
	Long: `Check dependencies declared in build scripts using #DEP: tags.

Each argument can be:
  - A local .sh file path
  - A directory (all .sh files under it are checked recursively)
  - A package name (e.g., samtools/1.22) resolved from local or remote build scripts

Dependencies from all scripts are deduplicated and checked together.

Note: If creation jobs are submitted to a scheduler, exits with code 3.`,
	Example: `  condatainer check script.sh             # Check single local script
  condatainer check samtools/1.22         # Check build script by name
  condatainer check script1.sh script2.sh # Check multiple scripts
  condatainer check ./my-scripts/         # Check all .sh files in a directory
  condatainer check -a ./src/*sh          # Check and auto-install missing dependencies from multiple scripts`,
	Args:         cobra.MinimumNArgs(1),
	SilenceUsage: true, // Runtime errors should not show usage
	RunE:         runCheck,
}

func init() {
	rootCmd.AddCommand(scriptCheckCmd)
	scriptCheckCmd.Flags().BoolVarP(&checkAutoInstall, "auto-install", "a", false, "Automatically install missing dependencies")
	scriptCheckCmd.Flags().BoolP("install", "i", false, "Alias for --auto-install")
	scriptCheckCmd.Flags().BoolVar(&checkParseModuleLoad, "module", false, "Also parse 'module load' / 'ml' lines as dependencies")
	scriptCheckCmd.Flags().BoolVar(&checkRemote, "remote", false, "Remote build scripts take precedence over local")
}

func runCheck(cmd *cobra.Command, args []string) error {
	// Check if -i flag was used (alias for -a)
	if installFlag, _ := cmd.Flags().GetBool("install"); installFlag {
		checkAutoInstall = true
	}
	build.PreferRemote = checkRemote || config.Global.PreferRemote

	// Resolve all args to concrete script paths and metadata deps
	scriptPaths, remotePaths, metaDeps, err := resolveAllScriptPaths(args)
	if err != nil {
		return err
	}
	for _, rp := range remotePaths {
		defer os.Remove(rp)
	}

	if len(scriptPaths) == 0 && len(metaDeps) == 0 {
		utils.PrintWarning("No scripts found to check.")
		return nil
	}

	// Collect and deduplicate deps across all scripts
	parseModuleLoad := config.Global.ParseModuleLoad || checkParseModuleLoad
	deps, err := collectDeps(scriptPaths, metaDeps, parseModuleLoad)
	if err != nil {
		return err
	}
	if len(deps) == 0 {
		utils.PrintMessage("No dependencies found in script.")
		return nil
	}

	// Check which deps are installed and print status
	installedOverlays, err := getInstalledOverlaysMap()
	if err != nil {
		return err
	}
	missingDeps := checkDeps(deps, installedOverlays)

	if len(missingDeps) == 0 {
		utils.PrintSuccess("All dependencies are installed.")
		return nil
	}

	if !checkAutoInstall {
		utils.PrintHint("Run the command again with %s or %s to install missing dependencies.",
			utils.StyleAction("-a"), utils.StyleAction("--auto-install"))
		return nil
	}

	// Auto-install missing dependencies
	utils.PrintMessage("Attempting to auto-install missing dependencies...")
	if err := ensureBaseImage(cmd.Context()); err != nil {
		return err
	}

	hasUnresolvable := false
	var packageDeps []string
	for _, dep := range missingDeps {
		if utils.IsOverlay(dep) {
			if !autoCreateExternalOverlay(cmd.Context(), dep) {
				hasUnresolvable = true
			}
		} else {
			packageDeps = append(packageDeps, dep)
		}
	}

	if len(packageDeps) > 0 {
		if err := autoInstallPackages(cmd.Context(), packageDeps); err != nil {
			return err
		}
	}

	if hasUnresolvable {
		return fmt.Errorf("some external overlays could not be auto-created")
	}
	utils.PrintSuccess("All selected overlays installed.")
	return nil
}

// collectDeps gathers deduplicated dependencies from all scripts.
// preSeededDeps are name/version deps from remote metadata (already normalized); they are
// merged first so script-parsed deps deduplicate against them.
// Relative overlay paths are resolved per-script using the scheduler WorkDir (if set) or cwd.
func collectDeps(scriptPaths []string, preSeededDeps []string, parseModuleLoad bool) ([]string, error) {
	multiScript := len(scriptPaths) > 1
	seen := make(map[string]bool)
	var deps []string

	// Pre-seed with metadata deps (name/version strings, never overlay paths)
	for _, dep := range preSeededDeps {
		key := utils.NormalizeNameVersion(dep)
		if !seen[key] {
			seen[key] = true
			deps = append(deps, dep)
		}
	}

	for _, scriptPath := range scriptPaths {
		if multiScript {
			utils.PrintMessage("Checking script: %s", utils.StylePath(scriptPath))
		}
		scriptDeps, err := utils.GetDependenciesFromScript(scriptPath, parseModuleLoad)
		if err != nil {
			return nil, fmt.Errorf("failed to parse dependencies from %s: %w", scriptPath, err)
		}

		// Prefer scheduler specs WorkDir if set, else current working directory.
		workDir, _ := os.Getwd()
		if specs, _ := scheduler.ReadScriptSpecsFromPath(scriptPath); specs != nil && specs.Control.WorkDir != "" {
			workDir = specs.Control.WorkDir
		}

		for _, dep := range scriptDeps {
			key := dep
			if utils.IsOverlay(dep) {
				if !filepath.IsAbs(dep) {
					dep = filepath.Join(workDir, dep)
				}
				key = dep
			} else {
				key = utils.NormalizeNameVersion(dep)
			}
			if !seen[key] {
				seen[key] = true
				deps = append(deps, dep)
			}
		}
	}
	return deps, nil
}

// checkDeps prints the status of each dependency grouped by type and returns the missing ones.
// Uses ✓/✗ symbols with colors. Paths are shown relative to cwd when possible.
func checkDeps(deps []string, installedOverlays map[string]string) []string {
	cwd, _ := os.Getwd()
	depDisplay := func(p string) string {
		if rel, err := filepath.Rel(cwd, p); err == nil && !strings.HasPrefix(rel, "..") {
			return rel
		}
		return p
	}

	check := utils.StyleSuccess("✓")
	cross := utils.StyleError("✗")

	var overlays, packages []string
	for _, dep := range deps {
		if utils.IsOverlay(dep) {
			overlays = append(overlays, dep)
		} else {
			packages = append(packages, dep)
		}
	}

	var missing []string

	if len(packages) > 0 {
		fmt.Fprintf(os.Stdout, "%s\n", utils.StyleTitle("Module Overlays:"))
		for _, dep := range packages {
			normalized := utils.NormalizeNameVersion(dep)
			if _, ok := installedOverlays[normalized]; ok {
				fmt.Fprintf(os.Stdout, "  %s %s\n", check, dep)
			} else {
				fmt.Fprintf(os.Stdout, "  %s %s\n", cross, dep)
				missing = append(missing, dep)
			}
		}
	}

	if len(overlays) > 0 {
		fmt.Fprintf(os.Stdout, "%s\n", utils.StyleTitle("External Overlays:"))
		for _, dep := range overlays {
			if utils.FileExists(dep) {
				fmt.Fprintf(os.Stdout, "  %s %s\n", check, depDisplay(dep))
			} else {
				sibling := findSiblingSourceFile(dep)
				hint := ""
				if sibling != "" {
					ext := filepath.Ext(sibling)
					if strings.HasSuffix(sibling, ".sh") {
						if shDeps, _ := utils.GetDependenciesFromScript(sibling, false); len(shDeps) > 0 {
							hint = "  (" + utils.StyleNote(ext) + " found, but has #DEP - create manually)"
						} else {
							hint = "  (" + utils.StyleNote(ext) + " found)"
						}
					} else {
						hint = "  (" + utils.StyleNote(ext) + " found)"
					}
				}
				fmt.Fprintf(os.Stdout, "  %s %s%s\n", cross, depDisplay(dep), hint)
				missing = append(missing, dep)
			}
		}
	}

	return missing
}

// autoCreateExternalOverlay attempts to create a missing external overlay from a sibling
// source file (.yml, .yaml, .def, .sh). Returns true on success or skip, false if unresolvable.
func autoCreateExternalOverlay(ctx context.Context, dep string) bool {
	sibling := findSiblingSourceFile(dep)
	if sibling == "" {
		utils.PrintError("External overlay %s not found and no source file (.yml/.yaml/.def/.sh) found - cannot auto-create", utils.StyleName(dep))
		return false
	}
	absFile, _ := filepath.Abs(sibling)
	absPrefix, _ := filepath.Abs(dep[:len(dep)-len(filepath.Ext(dep))])
	outputDir := filepath.Dir(absPrefix)
	baseName := filepath.Base(absPrefix)

	if strings.HasSuffix(sibling, ".sh") {
		shDeps, _ := utils.GetDependenciesFromScript(sibling, false)
		if len(shDeps) > 0 {
			utils.PrintError("External overlay %s has .sh with #DEP - create manually:\ncondatainer create -f %s",
				utils.StyleName(filepath.Base(dep)), sibling)
			return false
		}
	}

	utils.PrintMessage("Auto-creating %s from %s...", utils.StyleName(filepath.Base(dep)), utils.StylePath(sibling))

	if utils.IsYaml(sibling) {
		bo, err := build.NewCondaObjectWithSource(baseName, absFile, outputDir, outputDir, false)
		if err != nil {
			utils.PrintError("Failed to create build object for %s: %v", dep, err)
			return false
		}
		if err := bo.Build(ctx, false); err != nil {
			utils.PrintError("Build failed for %s: %v", dep, err)
			return false
		}
		return true
	}

	// .def or .sh
	isApptainer := strings.HasSuffix(sibling, ".def")
	bo, err := build.FromExternalSource(ctx, absPrefix, absFile, isApptainer, outputDir)
	if err != nil {
		utils.PrintError("Failed to create build object for %s: %v", dep, err)
		return false
	}
	graph, err := build.NewBuildGraph(ctx, []*build.BuildObject{bo}, outputDir, config.GetWritableTmpDir(), config.Global.SubmitJob, false)
	if err != nil {
		utils.PrintError("Failed to create build graph for %s: %v", dep, err)
		return false
	}
	if err := graph.Run(ctx); err != nil {
		if !errors.Is(err, build.ErrBuildCancelled) {
			utils.PrintError("Build failed for %s: %v", dep, err)
			return false
		}
	}
	return true
}

// autoInstallPackages installs named package dependencies via BuildGraph.
func autoInstallPackages(ctx context.Context, packages []string) error {
	imagesDir, err := config.GetWritableImagesDir()
	if err != nil {
		return fmt.Errorf("no writable images directory found: %w", err)
	}
	buildObjects := make([]*build.BuildObject, 0, len(packages))
	for _, pkg := range packages {
		bo, err := build.NewBuildObject(ctx, pkg, false, imagesDir, config.GetWritableTmpDir(), false)
		if err != nil {
			return fmt.Errorf("failed to create build object for %s: %w", pkg, err)
		}
		buildObjects = append(buildObjects, bo)
	}

	graph, err := build.NewBuildGraph(ctx, buildObjects, imagesDir, config.GetWritableTmpDir(), config.Global.SubmitJob, false)
	if err != nil {
		return fmt.Errorf("failed to create build graph: %w", err)
	}

	if err := graph.Run(ctx); err != nil {
		if errors.Is(err, build.ErrBuildCancelled) ||
			errors.Is(ctx.Err(), context.Canceled) ||
			strings.Contains(err.Error(), "signal: killed") ||
			strings.Contains(err.Error(), "context canceled") {
			utils.PrintWarning("Installation cancelled.")
			return nil
		}
		return fmt.Errorf("some overlays failed to install: %w", err)
	}

	ExitIfJobsSubmitted(graph)
	return nil
}

// findSiblingSourceFile looks for the first source file (.yml, .yaml, .def, .sh) with the same
// base name as the given overlay path. Returns empty string if none found.
func findSiblingSourceFile(overlayPath string) string {
	base := overlayPath[:len(overlayPath)-len(filepath.Ext(overlayPath))]
	for _, ext := range []string{".yml", ".yaml", ".def", ".sh"} {
		if utils.FileExists(base + ext) {
			return base + ext
		}
	}
	return ""
}

// findScriptsInDir finds all .sh files directly inside dir (non-recursive).
// If no scripts are found, a warning is printed but an empty slice is returned without error.
func findScriptsInDir(dir string) ([]string, error) {
	entries, err := os.ReadDir(dir)
	if err != nil {
		return nil, fmt.Errorf("failed to read directory %s: %w", dir, err)
	}
	var found []string
	for _, e := range entries {
		if !e.IsDir() && strings.HasSuffix(e.Name(), ".sh") {
			found = append(found, filepath.Join(dir, e.Name()))
		}
	}
	if len(found) == 0 {
		utils.PrintWarning("No .sh files found in directory %s", utils.StylePath(dir))
	}
	return found, nil
}

// resolveAllScriptPaths resolves all positional arguments to concrete file paths.
// Directories are expanded to all .sh files found recursively.
// Other args are resolved via resolveScriptPath (file, local build script, or remote download).
// Returns scriptPaths (deduplicated), remotePaths (to defer-remove), metaDeps (deps from
// remote metadata, returned directly without downloading the script), and any error.
func resolveAllScriptPaths(args []string) (scriptPaths []string, remotePaths []string, metaDeps []string, err error) {
	seen := make(map[string]bool)

	for _, arg := range args {
		if utils.DirExists(arg) {
			dirScripts, dirErr := findScriptsInDir(arg)
			if dirErr != nil {
				return nil, remotePaths, metaDeps, dirErr
			}
			for _, p := range dirScripts {
				if !seen[p] {
					seen[p] = true
					scriptPaths = append(scriptPaths, p)
				}
			}
			continue
		}

		resolved, isRemote, argMetaDeps, resolveErr := resolveScriptPath(arg)
		if resolveErr != nil {
			return nil, remotePaths, metaDeps, resolveErr
		}
		if argMetaDeps != nil {
			// Metadata deps available — no script to download
			metaDeps = append(metaDeps, argMetaDeps...)
			continue
		}
		if isRemote {
			remotePaths = append(remotePaths, resolved)
		}
		if !seen[resolved] {
			seen[resolved] = true
			scriptPaths = append(scriptPaths, resolved)
		}
	}

	return scriptPaths, remotePaths, metaDeps, nil
}

// resolveScriptPath resolves a script path or name to an actual file path.
// Returns (path, isRemote, metaDeps, error).
// When metaDeps is non-nil, deps were obtained directly from remote metadata —
// no script was downloaded and path will be empty.
func resolveScriptPath(scriptPathOrName string) (string, bool, []string, error) {
	// If it's a file, use it directly
	if utils.FileExists(scriptPathOrName) {
		return scriptPathOrName, false, nil, nil
	}

	// Try to find as build script (respects build.PreferRemote / --remote flag)
	normalized := utils.NormalizeNameVersion(scriptPathOrName)
	utils.PrintDebug("[CHECK] Checking for build script %s...", normalized)

	info, found := build.FindBuildScript(normalized)
	if !found {
		return "", false, nil, fmt.Errorf("build script for %s not found", normalized)
	}

	if info.IsRemote {
		// If metadata provides deps, use them directly — no download needed.
		if info.Deps != nil {
			utils.PrintDebug("[CHECK] Using metadata deps for %s (no download)", normalized)
			return "", false, info.Deps, nil
		}

		utils.PrintMessage("Downloading build script for %s from remote metadata...", normalized)
		tmpDir := config.GetWritableTmpDir()
		tempPath, err := build.DownloadRemoteScript(info, tmpDir)
		if err != nil {
			return "", false, nil, fmt.Errorf("failed to download build script: %w", err)
		}
		utils.PrintMessage("Downloaded build script to %s", utils.StylePath(tempPath))
		return tempPath, true, nil, nil
	}

	// Local script
	utils.PrintDebug("[CHECK] Found local build script %s", utils.StylePath(info.Path))
	return info.Path, false, nil, nil
}
