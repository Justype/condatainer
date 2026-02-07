package cmd

import (
	"fmt"
	"os"
	"path/filepath"
	"regexp"
	"sort"
	"strings"

	"github.com/spf13/cobra"

	"github.com/Justype/condatainer/internal/build"
	"github.com/Justype/condatainer/internal/config"
	"github.com/Justype/condatainer/internal/utils"
)

var (
	listDelete bool
	listExact  bool
)

var listCmd = &cobra.Command{
	Use:     "list [terms...]",
	Aliases: []string{"ls"},
	Short:   "List installed overlays matching search terms",
	Long:    `List installed overlays with optional search filtering.`,
	Example: `  condatainer list              # List all installed overlays
  condatainer list samtools     # Search for samtools
  condatainer list sam 1.22     # Search with multiple terms (AND logic)
  condatainer list -e python    # Exact match only`,
	SilenceUsage: true, // Runtime errors should not show usage
	RunE:         runList,
}

func init() {
	rootCmd.AddCommand(listCmd)
	listCmd.Flags().BoolVarP(&listDelete, "delete", "d", false, "Delete listed overlays after confirmation (used with search terms)")
	listCmd.Flags().BoolP("remove", "r", false, "Alias for --delete")
	listCmd.Flags().BoolVarP(&listExact, "exact", "e", false, "Require exact match instead of substring match")
}

func runList(cmd *cobra.Command, args []string) error {
	filters := normalizeFilters(args)

	appOverlays, err := collectAppOverlays(filters, listExact)
	if err != nil {
		return err
	}
	refOverlays, err := collectReferenceOverlays(filters, listExact)
	if err != nil {
		return err
	}

	// Separate OS overlays (def-built) from app/env overlays
	osOverlays := map[string][]string{}
	moduleOverlays := map[string][]string{}
	for name, versions := range appOverlays {
		var osVers, modVers []string
		for _, v := range versions {
			if v == "(system app)" {
				osVers = append(osVers, v)
			} else {
				modVers = append(modVers, v)
			}
		}
		if len(osVers) > 0 {
			osOverlays[name] = osVers
		}
		if len(modVers) > 0 {
			moduleOverlays[name] = modVers
		}
	}

	printed := false
	if len(osOverlays) > 0 {
		printed = true
		fmt.Println("Available OS overlays:")
		names := make([]string, 0, len(osOverlays))
		for name := range osOverlays {
			names = append(names, name)
		}
		sort.Strings(names)
		for _, name := range names {
			fmt.Printf(" %s\n", utils.StyleName(name))
		}
	}

	if len(moduleOverlays) > 0 {
		if printed {
			fmt.Println()
		}
		printed = true
		fmt.Println("Available app overlays:")
		names := make([]string, 0, len(moduleOverlays))
		for name := range moduleOverlays {
			names = append(names, name)
		}
		sort.Strings(names)
		nameWidth := maxWidth(names)
		for _, name := range names {
			vers := moduleOverlays[name]
			colored := make([]string, 0, len(vers))
			for _, v := range vers {
				if v == "(env)" {
					colored = append(colored, v)
				} else {
					colored = append(colored, utils.StyleInfo(v))
				}
			}
			values := strings.Join(colored, ", ")
			nameField := fmt.Sprintf("%-*s", nameWidth, name)
			fmt.Printf(" %s: %s\n", utils.StyleName(nameField), values)
		}
	}

	if len(refOverlays) > 0 {
		if printed {
			fmt.Println()
		}
		fmt.Println("Available reference overlays:")
		for _, ref := range refOverlays {
			fmt.Printf(" %s\n", utils.StyleName(ref))
		}
		printed = true
	}

	if !printed {
		utils.PrintWarning("No installed overlays match the provided search terms.")
		os.Exit(1)
	}

	// Handle delete if requested
	if removeFlag, _ := cmd.Flags().GetBool("remove"); removeFlag {
		listDelete = true
	}

	if listDelete && len(filters) > 0 && (len(appOverlays) > 0 || len(refOverlays) > 0) {
		fmt.Println()
		fmt.Println("==================REMOVE==================")

		// Collect all matching overlays for deletion
		allMatching := []string{}
		for name, versions := range appOverlays {
			for _, version := range versions {
				if version == "(system app)" || version == "(env)" {
					allMatching = append(allMatching, name)
				} else {
					allMatching = append(allMatching, name+"/"+version)
				}
			}
		}
		for _, ref := range refOverlays {
			allMatching = append(allMatching, ref)
		}

		if len(allMatching) > 0 {
			fmt.Println("Overlays to be removed:")
			for _, name := range allMatching {
				highlighted := name

				// Use single-pass highlighting to avoid ANSI code interference
				if len(filters) > 0 {
					// Sort filters by length (longest first) to prefer longer matches
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

					highlighted = re.ReplaceAllStringFunc(highlighted, func(match string) string {
						return utils.StyleWarning(match)
					})
				}

				fmt.Printf(" - %s\n", highlighted)
			}

			fmt.Print("\n")
			shouldDelete := false
			if utils.ShouldAnswerYes() {
				shouldDelete = true
			} else {
				fmt.Print("Are you sure? Cannot be undone. [y/N]: ")
				var choice string
				fmt.Scanln(&choice)
				choice = strings.ToLower(strings.TrimSpace(choice))
				shouldDelete = (choice == "y" || choice == "yes")
			}

			if shouldDelete {
				// Get installed overlays map for deletion
				installedMap, err := getInstalledOverlaysMap()
				if err != nil {
					return err
				}

				for _, name := range allMatching {
					normalized := utils.NormalizeNameVersion(name)
					if overlayPath, ok := installedMap[normalized]; ok {
						if err := os.Remove(overlayPath); err != nil {
							utils.PrintError("Failed to remove overlay %s: %v", utils.StyleName(name), err)
							continue
						}
						utils.PrintSuccess("Overlay %s removed.", utils.StyleName(name))

						// Also remove .env file if exists
						envPath := overlayPath + ".env"
						if utils.FileExists(envPath) {
							os.Remove(envPath)
						}
					}
				}
			} else {
				utils.PrintNote("Cancelled")
			}
		}
	}

	return nil
}

func collectAppOverlays(filters []string, exactMatch bool) (map[string][]string, error) {
	grouped := map[string]map[string]struct{}{}
	seen := make(map[string]bool) // Track seen files to avoid duplicates

	// Scan all image directories (user → scratch → system)
	for _, imageDir := range config.GetImageSearchPaths() {
		if !utils.DirExists(imageDir) {
			continue
		}

		entries, err := os.ReadDir(imageDir)
		if err != nil {
			utils.PrintDebug("Unable to list overlays in %s: %v", imageDir, err)
			continue
		}

		for _, entry := range entries {
			if entry.IsDir() {
				continue
			}
			if !utils.IsOverlay(entry.Name()) {
				continue
			}

			// Skip if already seen in higher-priority directory
			if seen[entry.Name()] {
				continue
			}
			seen[entry.Name()] = true

			nameVersion := strings.TrimSuffix(entry.Name(), filepath.Ext(entry.Name()))
			normalized := strings.ToLower(utils.NormalizeNameVersion(nameVersion))
			if !matchesFilters(normalized, filters, exactMatch) {
				continue
			}

			var name, version string
			if isDefBuilt(nameVersion) {
				// Def-built overlays are OS overlays regardless of delimiter count
				name = strings.ReplaceAll(nameVersion, "--", "/")
				version = "(system app)"
			} else {
				delimCount := strings.Count(entry.Name(), "--")
				if delimCount > 1 {
					// Not an app overlay
					continue
				}
				if strings.Contains(nameVersion, "--") {
					parts := strings.SplitN(nameVersion, "--", 2)
					name = parts[0]
					version = parts[1]
				} else {
					name = nameVersion
					version = "(env)"
				}
			}
			if name == "" {
				continue
			}

			if grouped[name] == nil {
				grouped[name] = map[string]struct{}{}
			}
			grouped[name][version] = struct{}{}
		}
	}

	result := map[string][]string{}
	for name, versions := range grouped {
		valueList := make([]string, 0, len(versions))
		for version := range versions {
			valueList = append(valueList, version)
		}
		sort.Strings(valueList)
		result[name] = valueList
	}
	return result, nil
}

func collectReferenceOverlays(filters []string, exactMatch bool) ([]string, error) {
	seen := make(map[string]bool) // Track seen files to avoid duplicates
	names := []string{}

	// Scan all image directories (user → scratch → system)
	for _, imageDir := range config.GetImageSearchPaths() {
		if !utils.DirExists(imageDir) {
			continue
		}

		entries, err := os.ReadDir(imageDir)
		if err != nil {
			utils.PrintDebug("Unable to list overlays in %s: %v", imageDir, err)
			continue
		}

		for _, entry := range entries {
			if entry.IsDir() {
				continue
			}
			if !utils.IsOverlay(entry.Name()) {
				continue
			}

			// Skip if already seen in higher-priority directory
			if seen[entry.Name()] {
				continue
			}
			seen[entry.Name()] = true

			delimCount := strings.Count(entry.Name(), "--")
			if delimCount <= 1 {
				// Not a reference overlay
				continue
			}
			nameVersion := strings.TrimSuffix(entry.Name(), filepath.Ext(entry.Name()))
			normalized := strings.ToLower(utils.NormalizeNameVersion(nameVersion))
			if !matchesFilters(normalized, filters, exactMatch) {
				continue
			}
			displayName := strings.ReplaceAll(nameVersion, "--", "/")
			names = append(names, displayName)
		}
	}

	sort.Strings(names)
	return names, nil
}

func normalizeFilters(filters []string) []string {
	normalized := make([]string, 0, len(filters))
	for _, filter := range filters {
		trimmed := strings.TrimSpace(filter)
		if trimmed == "" {
			continue
		}
		normalized = append(normalized, strings.ToLower(utils.NormalizeNameVersion(trimmed)))
	}
	return normalized
}

func matchesFilters(name string, filters []string, exactMatch bool) bool {
	if len(filters) == 0 {
		return true
	}
	if exactMatch {
		// For exact match, the name must equal one of the filters
		for _, filter := range filters {
			if name == filter {
				return true
			}
		}
		return false
	}
	// For substring match, all filters must be present in the name
	for _, filter := range filters {
		if !strings.Contains(name, filter) {
			return false
		}
	}
	return true
}

func maxWidth(names []string) int {
	width := 0
	for _, name := range names {
		if len(name) > width {
			width = len(name)
		}
	}
	return width
}

// isDefBuilt checks if an overlay name is .def-built
func isDefBuilt(nameVersion string) bool {
	normalized := utils.NormalizeNameVersion(nameVersion)
	whitelist := build.GetDefBuiltWhitelist()
	return whitelist[normalized]
}
