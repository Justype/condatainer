package cmd

import (
	"context"
	"errors"
	"fmt"
	"os"
	"path/filepath"
	"sort"
	"strings"

	"github.com/Justype/condatainer/internal/build"
	"github.com/Justype/condatainer/internal/config"
	"github.com/Justype/condatainer/internal/utils"
	"github.com/spf13/cobra"
)

var (
	availRemote  bool
	availInstall bool
	availExpand  bool
	availWhatis  bool
)

var availCmd = &cobra.Command{
	Use:     "avail [terms...]",
	Aliases: []string{"av"},
	Short:   "Check available build scripts (local and remote)",
	Long: `List available build scripts (local and remote).

By default, shows both local and remote build scripts.
When duplicates exist, local scripts take precedence.
Use --remote to make remote scripts take precedence over local.

Search rules (single term):
  Plain string  substring match
  cell*         wildcard (* and ?)
  ^cell.*9\.0   regex (any of ^ $ ( [ + { |)

Search rules (multiple terms, space-separated):
  First term is an exact name  → all terms treated as exact names
  First term not found         → AND substring match

Note: If creation jobs are submitted to a scheduler, exits with code 3.`,
	Example: `  condatainer avail                          # List all
  condatainer avail --remote                 # Remote takes precedence on duplicates
  condatainer avail cellranger               # Substring match
  condatainer avail 'cell*'                  # Wildcard
  condatainer avail 'cellranger|spaceranger' # Regex
  condatainer avail cellranger 9             # AND search (multiple terms)
  condatainer avail cellranger/9.0.1 -i      # Install exact version`,
	SilenceUsage: true, // Runtime errors should not show usage
	RunE:         runAvail,
}

func init() {
	rootCmd.AddCommand(availCmd)
	availCmd.Flags().BoolVar(&availRemote, "remote", false, "Remote build scripts take precedence over local")
	availCmd.Flags().BoolVarP(&availInstall, "install", "i", false, "Install the selected build scripts (used with search terms)")
	availCmd.Flags().BoolP("add", "a", false, "Alias for --install")
	availCmd.Flags().BoolVarP(&availExpand, "expand", "e", false, "Expand template groups to show individual entries")
	availCmd.Flags().BoolVarP(&availWhatis, "whatis", "w", false, "Show description for each build script")
}

// PackageInfo holds information about a build script
type PackageInfo struct {
	Name        string              // name/version format
	Path        string              // full path to build script
	IsContainer bool                // true if .def file
	IsInstalled bool                // true if overlay exists
	IsRemote    bool                // true if from remote repository
	IsTemplate  bool                // true if script has #PL: placeholders
	Whatis      string              // description string
	PL          map[string][]string // placeholder values (all values for templates; single-element for expanded)
	PLOrder     []string            // placeholder key declaration order
}

func runAvail(cmd *cobra.Command, args []string) error {
	// Normalize search terms
	filters := normalizeFilters(args)

	// Get installed overlays for marking
	installedOverlays := getInstalledOverlays()

	// Build package info list
	// When --remote: add remote first, then local (remote wins on duplicates)
	// Default: add local first, then remote (local wins on duplicates)
	packages := make([]PackageInfo, 0)
	seen := make(map[string]bool)

	scriptInfoToPackageInfo := func(name string, info build.ScriptInfo, isRemote bool) PackageInfo {
		return PackageInfo{
			Name:        name,
			Path:        info.Path,
			IsContainer: info.IsContainer,
			IsInstalled: installedOverlays[name],
			IsRemote:    isRemote,
			IsTemplate:  info.IsTemplate,
			Whatis:      info.Whatis,
			PL:          info.PL,
			PLOrder:     info.PLOrder,
		}
	}

	if availRemote {
		// Remote first: remote scripts take precedence
		remoteScripts, err := build.GetRemoteBuildScripts()
		if err != nil {
			utils.PrintDebug("Failed to fetch remote scripts: %v", err)
		} else {
			for name, info := range remoteScripts {
				packages = append(packages, scriptInfoToPackageInfo(name, info, true))
				seen[name] = true
			}
		}

		localScripts, err := build.GetLocalBuildScripts()
		if err != nil {
			return err
		}
		for name, info := range localScripts {
			if seen[name] {
				continue
			}
			packages = append(packages, scriptInfoToPackageInfo(name, info, false))
		}
	} else {
		// Default: local scripts take precedence
		localScripts, err := build.GetLocalBuildScripts()
		if err != nil {
			return err
		}
		for name, info := range localScripts {
			packages = append(packages, scriptInfoToPackageInfo(name, info, false))
			seen[name] = true
		}

		remoteScripts, err := build.GetRemoteBuildScripts()
		if err != nil {
			utils.PrintDebug("Failed to fetch remote scripts: %v", err)
		} else {
			for name, info := range remoteScripts {
				if seen[name] {
					continue
				}
				packages = append(packages, scriptInfoToPackageInfo(name, info, true))
			}
		}
	}

	// Build search query (exact-first: check if first term is a known package name or alias)
	distroLower := strings.ToLower(config.Global.DefaultDistro)
	distroPrefix := distroLower + "/"
	packageNameSet := make(map[string]bool, len(packages))
	for _, pkg := range packages {
		lower := strings.ToLower(pkg.Name)
		packageNameSet[lower] = true
		if distroLower != "" {
			if alias, ok := strings.CutPrefix(lower, distroPrefix); ok {
				packageNameSet[alias] = true
			}
		}
	}
	var availExactLookup func(string) bool
	if len(filters) > 1 {
		availExactLookup = func(term string) bool { return packageNameSet[term] }
	}
	query := NewSearchQuery(filters, availExactLookup)

	// Filter packages; in exact mode also accept distro-prefix alias
	filtered := make([]PackageInfo, 0)
	for _, pkg := range packages {
		var alias string
		if distroLower != "" {
			if a, ok := strings.CutPrefix(strings.ToLower(pkg.Name), distroPrefix); ok {
				alias = a
			}
		}
		if query.MatchesOrAlias(pkg.Name, alias) {
			filtered = append(filtered, pkg)
		}
	}

	if len(filtered) == 0 {
		utils.PrintWarning("No matching build scripts found.")
		return nil
	}

	// Sort by name
	sort.Slice(filtered, func(i, j int) bool {
		return filtered[i].Name < filtered[j].Name
	})

	// Check if -a flag was used (alias for -i)
	if addFlag, _ := cmd.Flags().GetBool("add"); addFlag {
		availInstall = true
	}

	// When --expand is set and the filter matched a template entry, also include all
	// of that template's expanded variants from the full package list.
	if availExpand {
		templatePaths := make(map[string]bool)
		for _, pkg := range filtered {
			if pkg.IsTemplate {
				templatePaths[pkg.Path] = true
			}
		}
		if len(templatePaths) > 0 {
			seen := make(map[string]bool, len(filtered))
			for _, pkg := range filtered {
				seen[pkg.Name] = true
			}
			for _, pkg := range packages {
				if !pkg.IsTemplate && len(pkg.PLOrder) > 0 && templatePaths[pkg.Path] && !seen[pkg.Name] {
					filtered = append(filtered, pkg)
					seen[pkg.Name] = true
				}
			}
			// Re-sort after adding expanded variants
			sort.Slice(filtered, func(i, j int) bool {
				return filtered[i].Name < filtered[j].Name
			})
		}
	}

	// Print results
	// --expand: show individual expanded entries (skip template headers)
	// no filter / no --expand: collapsed view (template headers + non-PL entries)
	// filter without --expand: same collapsed logic, but templates that matched are shown as group headers
	if availExpand {
		for _, pkg := range filtered {
			if !pkg.IsTemplate {
				fmt.Println(formatPackageLine(pkg, availWhatis))
			}
		}
	} else {
		// Collapsed view: show template groups as a single summary line;
		// skip expanded entries whose template is also in the list.
		templateNames := make(map[string]bool)
		for _, pkg := range filtered {
			if pkg.IsTemplate {
				templateNames[pkg.Name] = true
			}
		}
		for _, pkg := range filtered {
			if !pkg.IsTemplate && len(pkg.PLOrder) > 0 {
				// This is an expanded PL entry — suppress it if its template is shown.
				// An expanded entry shares the same script path as its template.
				// Check by seeing whether any template in the list has the same Path.
				suppress := false
				for _, t := range filtered {
					if t.IsTemplate && t.Path == pkg.Path {
						suppress = true
						break
					}
				}
				if suppress {
					continue
				}
			}
			if pkg.IsTemplate {
				fmt.Println(formatTemplateLine(pkg, availWhatis))
			} else {
				fmt.Println(formatPackageLine(pkg, availWhatis))
			}
		}
	}

	// Handle install if requested
	if availInstall && len(filters) > 0 {
		// Resolve templates interactively; collect concrete uninstalled names.
		// Skip expanded PL entries (their template is shown instead).
		uninstalled := []string{}
		for _, pkg := range filtered {
			if pkg.IsInstalled {
				continue
			}
			if pkg.IsTemplate {
				// Re-fetch full ScriptInfo (has TargetTemplate + Deps) then resolve interactively.
				info, found := build.FindBuildScript(pkg.Name)
				if !found {
					utils.PrintWarning("Could not find build script for %s — skipping.", pkg.Name)
					continue
				}
				concrete, err := resolveTemplateInteractively(cmd.Context(), info)
				if err != nil {
					return err
				}
				if !installedOverlays[concrete] {
					uninstalled = append(uninstalled, concrete)
				} else {
					utils.PrintNote("%s is already installed.", concrete)
				}
				continue
			}
			// Skip expanded PL entries (their template handles installation above).
			if len(pkg.PLOrder) > 0 {
				continue
			}
			uninstalled = append(uninstalled, pkg.Name)
		}

		if len(uninstalled) > 0 {
			fmt.Println()
			fmt.Println("==================INSTALL==================")
			fmt.Println("The following overlays will be installed:")
			for _, name := range uninstalled {
				fmt.Printf(" - %s\n", name)
			}
			fmt.Println()

			// Ask for confirmation (skip if --yes flag is set)
			if !utils.ShouldAnswerYes() {
				fmt.Print("Do you want to proceed with the installation? [y/N]: ")
				choice, err := utils.ReadLineContext(cmd.Context())
				if err != nil {
					utils.PrintWarning("Installation cancelled.")
					return nil
				}
				if strings.ToLower(choice) != "y" && strings.ToLower(choice) != "yes" {
					utils.PrintNote("Installation cancelled.")
					return nil
				}
			}

			// Get writable directories
			imagesDir, err := config.GetWritableImagesDir()
			if err != nil {
				return fmt.Errorf("failed to get writable images directory: %w", err)
			}

			// If --remote flag or config is set, ensure build resolution only uses remote scripts
			build.PreferRemote = availRemote || config.Global.PreferRemote

			buildObjects := make([]*build.BuildObject, 0, len(uninstalled))
			for _, pkg := range uninstalled {
				bo, err := build.NewBuildObject(cmd.Context(), pkg, false, imagesDir, config.GetWritableTmpDir(), false)
				if err != nil {
					return fmt.Errorf("failed to create build object for %s: %w", pkg, err)
				}
				buildObjects = append(buildObjects, bo)
			}

			// Build graph and execute
			graph, err := build.NewBuildGraph(cmd.Context(), buildObjects, imagesDir, config.GetWritableTmpDir(), config.Global.SubmitJob, false)
			if err != nil {
				return fmt.Errorf("failed to create build graph: %w", err)
			}

			if err := graph.Run(cmd.Context()); err != nil {
				if errors.Is(err, build.ErrBuildCancelled) ||
					errors.Is(cmd.Context().Err(), context.Canceled) ||
					strings.Contains(err.Error(), "signal: killed") ||
					strings.Contains(err.Error(), "context canceled") {
					utils.PrintWarning("Installation cancelled.")
					return nil
				}
				utils.PrintError("Some overlays failed to install.")
				return err
			}

			// If jobs were submitted to the scheduler, exit with the same distinct code
			// so downstream tooling can detect overlays will be created asynchronously.
			ExitIfJobsSubmitted(graph)
			utils.PrintSuccess("All selected overlays installed.")
		} else {
			utils.PrintNote("All matching packages are already installed.")
		}
		return nil
	}

	// Print summary when showing both
	if !availRemote {
		localCount := 0
		remoteCount := 0
		for _, pkg := range filtered {
			if pkg.IsRemote {
				remoteCount++
			} else {
				localCount++
			}
		}
		if remoteCount > 0 {
			fmt.Println()
			utils.PrintMessage("Found %s local and %s remote build scripts.",
				utils.StyleNumber(localCount), utils.StyleNumber(remoteCount))
		}
	}

	return nil
}

// getInstalledOverlays returns a set of installed overlay names from all search paths
func getInstalledOverlays() map[string]bool {
	installed := make(map[string]bool)

	// Check all image search paths
	for _, imagesDir := range config.GetImageSearchPaths() {
		if !utils.DirExists(imagesDir) {
			continue
		}

		entries, err := os.ReadDir(imagesDir)
		if err != nil {
			utils.PrintWarning("Failed to read directory %s: %v", imagesDir, err)
			continue
		}

		for _, entry := range entries {
			if entry.IsDir() {
				continue
			}
			if !utils.IsOverlay(entry.Name()) {
				continue
			}
			// Convert filename to name/version format
			nameVersion := strings.TrimSuffix(entry.Name(), filepath.Ext(entry.Name()))
			// Convert samtools--1.21 to samtools/1.21
			normalized := strings.ReplaceAll(nameVersion, "--", "/")
			installed[normalized] = true
		}
	}

	return installed
}

// formatTemplateLine formats a template (PL) package as a collapsed group header.
// Shows placeholders and their value lists, plus a variant count.
//
//	grch38/star-gencode  [10 variants]
//	    STAR GRCh38 GENCODE{gencode_version} index for read length {read_length}
//	    star_version:    2.7.11b, 2.7.11a
//	    gencode_version: 22-49  (28 values)
//	    read_length:     151, 101, *
func formatTemplateLine(pkg PackageInfo, showWhatis bool) string {
	// Compute variant count as the Cartesian product of concrete (non-*) values.
	variantCount := 1
	for _, vals := range pkg.PL {
		concrete := 0
		for _, v := range vals {
			if v != "*" {
				concrete++
			}
		}
		if concrete > 0 {
			variantCount *= concrete
		}
	}
	if variantCount == 1 && len(pkg.PL) == 0 {
		variantCount = 0
	}

	var labelParts []string
	if variantCount > 0 {
		labelParts = append(labelParts, utils.StyleDebug(fmt.Sprintf("%d variants", variantCount)))
	} else {
		labelParts = append(labelParts, utils.StyleDebug("template"))
	}
	if pkg.IsRemote {
		labelParts = append(labelParts, utils.StyleHint("remote"))
	}

	line := fmt.Sprintf("%s  %s%s%s", utils.StyleName(pkg.Name),
		utils.StyleDebug("["), strings.Join(labelParts, utils.StyleDebug(", ")), utils.StyleDebug("]"))
	if showWhatis && pkg.Whatis != "" {
		line += "\n  " + utils.StyleHint(pkg.Whatis)
	}

	// Show placeholder value summaries with aligned keys.
	// First pass: compute max key length for alignment.
	maxKeyLen := 0
	for _, key := range pkg.PLOrder {
		if _, ok := pkg.PL[key]; ok && len(key) > maxKeyLen {
			maxKeyLen = len(key)
		}
	}
	for _, key := range pkg.PLOrder {
		vals, ok := pkg.PL[key]
		if !ok {
			continue
		}
		var concrete []string
		hasOpen := false
		for _, v := range vals {
			if v == "*" {
				hasOpen = true
			} else {
				concrete = append(concrete, v)
			}
		}
		n := len(concrete)
		var display string
		if n > 8 {
			display = fmt.Sprintf("%s-%s  %s", concrete[n-1], concrete[0], utils.StyleDebug(fmt.Sprintf("(%d values)", n)))
		} else {
			display = strings.Join(concrete, ", ")
		}
		if hasOpen {
			display += ", *"
		}
		// Pad key to align values column
		padding := strings.Repeat(" ", maxKeyLen-len(key))
		line += fmt.Sprintf("\n  - %s:%s  %s", key, padding, display)
	}

	return line
}

// formatPackageLine formats a package for display.
func formatPackageLine(pkg PackageInfo, showWhatis bool) string {
	line := utils.StyleName(pkg.Name)

	// Compute alias before highlighting (e.g. "ubuntu24/igv" → "[igv]")
	var alias string
	if distro := config.Global.DefaultDistro; distro != "" {
		if a, ok := strings.CutPrefix(pkg.Name, distro+"/"); ok {
			alias = a
		}
	}

	// Build suffix
	var suffixes []string
	if pkg.IsInstalled {
		suffixes = append(suffixes, utils.StyleSuccess("installed"))
	}
	if pkg.IsContainer {
		suffixes = append(suffixes, "container")
	}
	if pkg.IsRemote {
		suffixes = append(suffixes, utils.StyleHint("remote"))
	}

	if len(suffixes) > 0 {
		line = fmt.Sprintf("%s (%s)", line, strings.Join(suffixes, ", "))
	}

	if alias != "" {
		line += "  " + utils.StyleInfo("["+alias+"]")
	}

	if showWhatis && pkg.Whatis != "" {
		line += "  " + utils.StyleHint(pkg.Whatis)
	}

	return line
}
