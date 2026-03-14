package cmd

import (
	"fmt"
	"os"
	"path/filepath"
	"regexp"
	"sort"
	"strings"

	"github.com/spf13/cobra"

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

	appOverlays, dataOverlays, err := scanOverlays(filters, listExact)
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
		nameWidth := maxWidth(names)
		for _, name := range names {
			nameField := fmt.Sprintf("%-*s", nameWidth, name)
			line := fmt.Sprintf(" %s", utils.StyleName(nameField))
			if distro := config.Global.DefaultDistro; distro != "" {
				if alias, ok := strings.CutPrefix(name, distro+"/"); ok {
					line += "  " + utils.StyleInfo("["+alias+"]")
				}
			}
			fmt.Println(line)
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

	if len(dataOverlays) > 0 {
		if printed {
			fmt.Println()
		}
		fmt.Println("Available data overlays:")
		for _, data := range dataOverlays {
			fmt.Printf(" %s\n", utils.StyleName(data))
		}
		printed = true
	}

	if !printed {
		utils.PrintWarning("No installed overlays match the provided search terms.")
		os.Exit(ExitCodeError)
	}

	// Handle delete if requested
	if removeFlag, _ := cmd.Flags().GetBool("remove"); removeFlag {
		listDelete = true
	}

	if listDelete && len(filters) > 0 && (len(appOverlays) > 0 || len(dataOverlays) > 0) {
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
		for _, data := range dataOverlays {
			allMatching = append(allMatching, data)
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
				choice, choiceErr := utils.ReadLineContext(cmd.Context())
				shouldDelete = choiceErr == nil && (choice == "y" || choice == "yes")
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

// scanOverlays does a single ReadDir pass per image directory and returns both
// app overlays (map[name][]version) and data overlays ([]displayName).
// OS overlays (Apptainer-built SIF/SquashFS) are treated as app overlays
// regardless of their delimiter count.
func scanOverlays(filters []string, exactMatch bool) (map[string][]string, []string, error) {
	appGrouped := map[string]map[string]struct{}{}
	var dataList []string
	seen := make(map[string]bool)
	distroPrefix := strings.ToLower(config.Global.DefaultDistro) + "/"

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
			if entry.IsDir() || !utils.IsOverlay(entry.Name()) || seen[entry.Name()] {
				continue
			}
			seen[entry.Name()] = true

			nameVersion := strings.TrimSuffix(entry.Name(), filepath.Ext(entry.Name()))
			normalized := strings.ToLower(utils.NormalizeNameVersion(nameVersion))
			overlayPath := filepath.Join(imageDir, entry.Name())
			delimCount := strings.Count(entry.Name(), "--")
			osOverlay := isOSOverlay(overlayPath)

			if !osOverlay && delimCount > 1 {
				// Data overlay: filter without alias
				if matchesFilters(normalized, filters, exactMatch) {
					dataList = append(dataList, strings.ReplaceAll(nameVersion, "--", "/"))
				}
				continue
			}

			// App overlay (OS or regular): filter with distro-prefix alias fallback
			if !matchesFilters(normalized, filters, exactMatch) {
				alias := strings.TrimPrefix(normalized, distroPrefix)
				if alias == normalized || !matchesFilters(alias, filters, exactMatch) {
					continue
				}
			}

			var name, version string
			if osOverlay {
				name = strings.ReplaceAll(nameVersion, "--", "/")
				version = "(system app)"
			} else if strings.Contains(nameVersion, "--") {
				parts := strings.SplitN(nameVersion, "--", 2)
				name = parts[0]
				version = parts[1]
			} else {
				name = nameVersion
				version = "(env)"
			}
			if name == "" {
				continue
			}
			if appGrouped[name] == nil {
				appGrouped[name] = map[string]struct{}{}
			}
			appGrouped[name][version] = struct{}{}
		}
	}

	apps := map[string][]string{}
	for name, versions := range appGrouped {
		valueList := make([]string, 0, len(versions))
		for version := range versions {
			valueList = append(valueList, version)
		}
		sort.Strings(valueList)
		apps[name] = valueList
	}
	sort.Strings(dataList)
	return apps, dataList, nil
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
