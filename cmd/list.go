package cmd

import (
	"fmt"
	"os"
	"path/filepath"
	"regexp"
	"sort"
	"strings"

	"github.com/spf13/cobra"
	"golang.org/x/term"

	"github.com/Justype/condatainer/internal/config"
	"github.com/Justype/condatainer/internal/utils"
)

var listDelete bool
var listExact bool

var listCmd = &cobra.Command{
	Use:     "list [terms...]",
	Aliases: []string{"ls"},
	Short:   "List installed overlays matching search terms",
	Long: `List installed overlays grouped by directory, with optional search filtering.

Search rules (single term):
  Plain string  substring match
  cell*         wildcard (* and ?)
  ^cell.*9\.0   regex (any of ^ $ ( [ + { |)
  -e/--exact    force exact full-name match

Search rules (multiple terms, space-separated):
  First term is an exact name  → all terms treated as exact names
  First term not found         → AND substring match`,
	Example: `  condatainer list                        # List all
  condatainer list cellranger             # Substring match
  condatainer list 'cell*'               # Wildcard
  condatainer list cellranger 9          # AND search (multiple terms)
  condatainer list cellranger/9.0.1 -e   # Exact match (single term)`,
	SilenceUsage: true,
	RunE:         runList,
}

func init() {
	rootCmd.AddCommand(listCmd)
	listCmd.Flags().BoolVarP(&listDelete, "delete", "d", false, "Delete listed overlays after confirmation (used with search terms)")
	listCmd.Flags().BoolP("remove", "r", false, "Alias for --delete")
	listCmd.Flags().BoolVarP(&listExact, "exact", "e", false, "Force exact full-name match even for a single term")
	listCmd.Flags().StringP("dir", "D", "", "Limit to a specific image directory (substring match)")
}

// DirOverlays holds the scan results for a single image directory.
type DirOverlays struct {
	Dir       string
	AppGroups map[string][]string // name → sorted []version
	DataList  []string
}

func runList(cmd *cobra.Command, args []string) error {
	filters := normalizeFilters(args)
	filterActive := len(filters) > 0

	// Build exact lookup from installed overlay names (including distro-prefix aliases)
	installedMap, err := getInstalledOverlaysMap()
	if err != nil {
		return err
	}
	distroLower := strings.ToLower(config.Global.DefaultDistro)
	installedLower := make(map[string]bool, len(installedMap))
	for name := range installedMap {
		lower := strings.ToLower(name)
		installedLower[lower] = true
		if distroLower != "" {
			if alias, ok := strings.CutPrefix(lower, distroLower+"/"); ok {
				installedLower[alias] = true
			}
		}
	}
	var listExactLookup func(string) bool
	if listExact {
		listExactLookup = func(string) bool { return true }
	} else if len(filters) > 1 {
		listExactLookup = func(t string) bool { return installedLower[t] }
	}
	query := NewSearchQuery(filters, listExactLookup)

	dirFilter, _ := cmd.Flags().GetString("dir")
	dirs := filterImageDirs(config.GetImageSearchPaths(), dirFilter)
	if dirFilter != "" && len(dirs) == 0 {
		utils.PrintWarning("No image directory matching %q found.", dirFilter)
		os.Exit(ExitCodeError)
	}
	results := scanOverlaysByDir(dirs, query)

	hasAnyMatch := false
	firstSection := true
	for _, d := range results {
		total := len(d.AppGroups) + len(d.DataList)
		// When filtering, skip dirs with no matches
		if filterActive && total == 0 {
			continue
		}
		hasAnyMatch = hasAnyMatch || total > 0

		if !firstSection {
			fmt.Println()
		}
		firstSection = false

		fmt.Println(dirHeader(d.Dir))

		if total == 0 {
			fmt.Println("  (no overlays)")
			continue
		}

		// Separate OS overlays from app overlays
		osOverlays := map[string][]string{}
		moduleOverlays := map[string][]string{}
		for name, versions := range d.AppGroups {
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

		if len(osOverlays) > 0 {
			fmt.Println("Available OS overlays:")
			names := sortedKeys(osOverlays)
			nameWidth := maxWidth(names)
			for _, name := range names {
				nameField := fmt.Sprintf("%-*s", nameWidth, name)
				line := fmt.Sprintf(" %s", utils.StyleName(nameField))
				if distroLower != "" {
					if alias, ok := strings.CutPrefix(name, distroLower+"/"); ok {
						line += "  " + utils.StyleInfo("["+alias+"]")
					}
				}
				fmt.Println(line)
			}
		}

		if len(moduleOverlays) > 0 {
			fmt.Println("Available app overlays:")
			names := sortedKeys(moduleOverlays)
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
				nameField := fmt.Sprintf("%-*s", nameWidth, name)
				fmt.Printf(" %s: %s\n", utils.StyleName(nameField), strings.Join(colored, ", "))
			}
		}

		if len(d.DataList) > 0 {
			fmt.Println("Available data overlays:")
			for _, data := range d.DataList {
				fmt.Printf(" %s\n", utils.StyleName(data))
			}
		}
	}

	if !hasAnyMatch && filterActive {
		utils.PrintWarning("No installed overlays match the provided search terms.")
		os.Exit(ExitCodeError)
	}

	// Handle delete if requested
	if removeFlag, _ := cmd.Flags().GetBool("remove"); removeFlag {
		listDelete = true
	}

	if listDelete && filterActive && hasAnyMatch {
		var allMatching []string
		for _, d := range results {
			for name, versions := range d.AppGroups {
				for _, version := range versions {
					if version == "(system app)" || version == "(env)" {
						allMatching = append(allMatching, name)
					} else {
						allMatching = append(allMatching, name+"/"+version)
					}
				}
			}
			allMatching = append(allMatching, d.DataList...)
		}
		fmt.Print("\n")
		return performDelete(cmd, allMatching)
	}

	return nil
}

// scanOverlaysByDir scans each image directory independently and returns per-dir results.
// Missing directories are skipped; empty directories are included with empty groups.
func scanOverlaysByDir(dirs []string, query *SearchQuery) []DirOverlays {
	distroPrefix := strings.ToLower(config.Global.DefaultDistro) + "/"
	var result []DirOverlays

	for _, imageDir := range dirs {
		if !utils.DirExists(imageDir) {
			continue // skip missing
		}
		d := DirOverlays{Dir: imageDir, AppGroups: map[string][]string{}}

		entries, err := os.ReadDir(imageDir)
		if err != nil {
			result = append(result, d)
			continue
		}

		appGrouped := map[string]map[string]struct{}{}
		for _, entry := range entries {
			if entry.IsDir() || !utils.IsOverlay(entry.Name()) {
				continue
			}
			nameVersion := strings.TrimSuffix(entry.Name(), filepath.Ext(entry.Name()))
			normalized := strings.ToLower(utils.NormalizeNameVersion(nameVersion))
			overlayPath := filepath.Join(imageDir, entry.Name())
			delimCount := strings.Count(entry.Name(), "--")
			osOverlay := isOSOverlay(overlayPath)

			if !osOverlay && delimCount > 1 {
				// Data overlay
				if query.Matches(normalized) {
					d.DataList = append(d.DataList, strings.ReplaceAll(nameVersion, "--", "/"))
				}
				continue
			}

			// App overlay
			if !query.Matches(normalized) {
				alias := strings.TrimPrefix(normalized, distroPrefix)
				if alias == normalized || !query.Matches(alias) {
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

		for name, versions := range appGrouped {
			vers := make([]string, 0, len(versions))
			for v := range versions {
				vers = append(vers, v)
			}
			sort.Strings(vers)
			d.AppGroups[name] = vers
		}
		sort.Strings(d.DataList)
		result = append(result, d)
	}
	return result
}

// dirHeader returns a full-width separator line with the directory path centered.
func dirHeader(path string) string {
	w := terminalWidth()
	inner := " " + path + " "
	pad := w - len(inner)
	if pad < 2 {
		return inner
	}
	left := pad / 2
	right := pad - left
	return strings.Repeat("═", left) + inner + strings.Repeat("═", right)
}

// terminalWidth returns the current terminal width, defaulting to 80.
func terminalWidth() int {
	if w, _, err := term.GetSize(int(os.Stdout.Fd())); err == nil && w > 0 {
		return w
	}
	return 80
}

func sortedKeys(m map[string][]string) []string {
	keys := make([]string, 0, len(m))
	for k := range m {
		keys = append(keys, k)
	}
	sort.Strings(keys)
	return keys
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

// filterImageDirs filters dirs according to dirFilter:
//   - empty          → return all
//   - no leading /   → substring match ("scratch" matches "/scratch/user/images")
//   - starts with /  → exact match, unless it contains * or ? (wildcard via filepath.Match)
func filterImageDirs(dirs []string, dirFilter string) []string {
	if dirFilter == "" {
		return dirs
	}
	var out []string
	for _, d := range dirs {
		if matchDirFilter(d, dirFilter) {
			out = append(out, d)
		}
	}
	return out
}

func matchDirFilter(dir, filter string) bool {
	if !strings.HasPrefix(filter, "/") {
		return strings.Contains(dir, filter)
	}
	if strings.ContainsAny(filter, "*?") {
		// filepath.Match only matches within a single path component (* won't cross /).
		// Convert to regex so * matches across slashes: /scratch/* → ^/scratch/.*$
		pattern := "^" + regexp.QuoteMeta(filter) + "$"
		pattern = strings.ReplaceAll(pattern, `\*`, `.*`)
		pattern = strings.ReplaceAll(pattern, `\?`, `.`)
		matched, _ := regexp.MatchString(pattern, dir)
		return matched
	}
	return dir == filter
}
