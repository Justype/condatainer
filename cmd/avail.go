package cmd

import (
	"context"
	"errors"
	"fmt"
	"os"
	"path/filepath"
	"regexp"
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
)

var availCmd = &cobra.Command{
	Use:     "avail [terms...]",
	Aliases: []string{"av"},
	Short:   "Check available build scripts (local and remote)",
	Long: `List available build scripts (local and remote).

By default, shows both local and remote build scripts.
When duplicates exist, local scripts take precedence.
Use --remote to make remote scripts take precedence over local.
Use search terms to filter results (AND logic applied).

Note: If creation jobs are submitted to a scheduler, exits with code 3.`,
	Example: `  condatainer avail              # List all (local takes precedence on duplicates)
  condatainer avail --remote     # List all (remote takes precedence on duplicates)
  condatainer avail sam 1.21     # Search with multiple terms (AND logic)`,
	SilenceUsage: true, // Runtime errors should not show usage
	RunE:         runAvail,
}

func init() {
	rootCmd.AddCommand(availCmd)
	availCmd.Flags().BoolVar(&availRemote, "remote", false, "Remote build scripts take precedence over local")
	availCmd.Flags().BoolVarP(&availInstall, "install", "i", false, "Install the selected build scripts (used with search terms)")
	availCmd.Flags().BoolP("add", "a", false, "Alias for --install")
}

// PackageInfo holds information about a build script
type PackageInfo struct {
	Name        string // name/version format
	Path        string // full path to build script
	IsContainer bool   // true if .def file
	IsInstalled bool   // true if overlay exists
	IsRemote    bool   // true if from remote repository
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

	if availRemote {
		// Remote first: remote scripts take precedence
		remoteScripts, err := build.GetRemoteBuildScripts()
		if err != nil {
			utils.PrintDebug("Failed to fetch remote scripts: %v", err)
		} else {
			for name, info := range remoteScripts {
				packages = append(packages, PackageInfo{
					Name:        name,
					Path:        info.Path,
					IsContainer: info.IsContainer,
					IsInstalled: installedOverlays[name],
					IsRemote:    true,
				})
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
			packages = append(packages, PackageInfo{
				Name:        name,
				Path:        info.Path,
				IsContainer: info.IsContainer,
				IsInstalled: installedOverlays[name],
				IsRemote:    false,
			})
		}
	} else {
		// Default: local scripts take precedence
		localScripts, err := build.GetLocalBuildScripts()
		if err != nil {
			return err
		}
		for name, info := range localScripts {
			packages = append(packages, PackageInfo{
				Name:        name,
				Path:        info.Path,
				IsContainer: info.IsContainer,
				IsInstalled: installedOverlays[name],
				IsRemote:    false,
			})
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
				packages = append(packages, PackageInfo{
					Name:        name,
					Path:        info.Path,
					IsContainer: info.IsContainer,
					IsInstalled: installedOverlays[name],
					IsRemote:    true,
				})
			}
		}
	}

	// Filter packages
	filtered := filterPackages(packages, filters)

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

	// Print results
	for _, pkg := range filtered {
		line := formatPackageLine(pkg, filters)
		fmt.Println(line)
	}

	// Handle install if requested
	if availInstall && len(filters) > 0 {
		// Get uninstalled packages
		uninstalled := []string{}
		for _, pkg := range filtered {
			if !pkg.IsInstalled {
				uninstalled = append(uninstalled, pkg.Name)
			}
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
				var choice string
				fmt.Scanln(&choice)
				choice = strings.ToLower(strings.TrimSpace(choice))

				if choice != "y" && choice != "yes" {
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

			buildObjects := make([]build.BuildObject, 0, len(uninstalled))
			for _, pkg := range uninstalled {
				bo, err := build.NewBuildObject(pkg, false, imagesDir, config.GetWritableTmpDir())
				if err != nil {
					return fmt.Errorf("failed to create build object for %s: %w", pkg, err)
				}
				buildObjects = append(buildObjects, bo)
			}

			// Build graph and execute
			graph, err := build.NewBuildGraph(buildObjects, imagesDir, config.GetWritableTmpDir(), config.Global.SubmitJob)
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

// filterPackages filters packages based on search terms (AND logic)
func filterPackages(packages []PackageInfo, filters []string) []PackageInfo {
	if len(filters) == 0 {
		return packages
	}

	filtered := make([]PackageInfo, 0)
	for _, pkg := range packages {
		nameLower := strings.ToLower(pkg.Name)
		matchAll := true
		for _, filter := range filters {
			if !strings.Contains(nameLower, filter) {
				matchAll = false
				break
			}
		}
		if matchAll {
			filtered = append(filtered, pkg)
		}
	}
	return filtered
}

// formatPackageLine formats a package for display with optional highlighting
func formatPackageLine(pkg PackageInfo, filters []string) string {
	line := pkg.Name

	// Highlight search terms in the name only (before adding suffixes)
	// Use a single-pass approach to avoid ANSI code interference
	if len(filters) > 0 {
		// Build a combined regex pattern for all search terms
		// Sort by length (longest first) to prefer longer matches
		sortedFilters := make([]string, len(filters))
		copy(sortedFilters, filters)
		sort.Slice(sortedFilters, func(i, j int) bool {
			return len(sortedFilters[i]) > len(sortedFilters[j])
		})

		// Combine all terms into one regex with alternation
		patterns := make([]string, len(sortedFilters))
		for i, term := range sortedFilters {
			patterns[i] = regexp.QuoteMeta(term)
		}
		combinedPattern := "(?i)(" + strings.Join(patterns, "|") + ")"
		re := regexp.MustCompile(combinedPattern)

		line = re.ReplaceAllStringFunc(line, func(match string) string {
			return utils.StyleHighlight(match)
		})
	}

	// Build suffix (after highlighting to avoid highlighting text in suffixes)
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

	return line
}
