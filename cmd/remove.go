package cmd

import (
	"fmt"
	"os"
	"path/filepath"
	"sort"
	"strings"

	"github.com/Justype/condatainer/internal/config"
	"github.com/Justype/condatainer/internal/overlay"
	"github.com/Justype/condatainer/internal/utils"
	"github.com/spf13/cobra"
)

var removeCmd = &cobra.Command{
	Use:     "remove [terms...]",
	Aliases: []string{"rm", "delete", "uninstall"},
	Short:   "Remove installed overlays matching search terms",
	Long: `Remove installed overlays by name, alias, or search terms.

Search rules (single term):
  Plain string  exact match (including distro-prefix alias)
  cell*         wildcard (* and ?)
  ^cell.*9\.0   regex (any of ^ $ ( [ + { |)

Search rules (multiple terms, space-separated):
  First term is an exact name  → all terms treated as exact names
  First term not found         → AND substring match

The command will ask for confirmation before removing overlays.`,
	Example: `  condatainer rm cellranger/9.0.1                   # Remove exact version
  condatainer rm cellranger                         # Exact match: all cellranger
  condatainer rm 'cell*'                            # Wildcard
  condatainer rm cellranger/9.0.1 cellranger/8.0.1  # Remove multiple
  condatainer rm cellranger 9                       # AND search (multiple terms)`,
	Args:         cobra.MinimumNArgs(1),
	SilenceUsage: true,
	RunE:         runRemove,
}

func init() {
	rootCmd.AddCommand(removeCmd)
	removeCmd.Flags().StringP("dir", "D", "", "Limit to a specific image directory (substring match)")
}

func runRemove(cmd *cobra.Command, args []string) error {
	// Split args: file paths (have overlay extension) vs search terms
	var externalPaths []string
	var searchArgs []string
	for _, arg := range args {
		if utils.IsOverlay(arg) {
			externalPaths = append(externalPaths, arg)
		} else {
			searchArgs = append(searchArgs, arg)
		}
	}
	if len(externalPaths) > 0 && len(searchArgs) > 0 {
		return fmt.Errorf("cannot mix file paths and search terms; use separate commands")
	}
	if len(externalPaths) > 0 {
		return removeExternalOverlays(cmd, externalPaths)
	}

	terms := normalizeFilters(searchArgs)

	// Get installed overlays, optionally scoped to a specific dir
	installedOverlays, err := getInstalledOverlaysMap()
	if err != nil {
		return err
	}
	dirFilter, _ := cmd.Flags().GetString("dir")
	if dirFilter != "" {
		known := filterImageDirs(config.GetImageSearchPaths(), dirFilter)
		if len(known) == 0 {
			utils.PrintWarning("No image directory matching %q found.", dirFilter)
			os.Exit(ExitCodeError)
		}
		knownSet := make(map[string]bool, len(known))
		for _, d := range known {
			knownSet[d] = true
		}
		for name, path := range installedOverlays {
			if !knownSet[filepath.Dir(path)] {
				delete(installedOverlays, name)
			}
		}
	}

	if len(installedOverlays) == 0 {
		utils.PrintWarning("No installed overlays found.")
		return nil
	}

	// Build search query (exact match for single term; exact-first for multiple)
	distroLower := strings.ToLower(config.Global.DefaultDistro)
	distroPrefix := distroLower + "/"
	installedLower := make(map[string]bool, len(installedOverlays))
	for name := range installedOverlays {
		lower := strings.ToLower(name)
		installedLower[lower] = true
		if distroLower != "" {
			if alias, ok := strings.CutPrefix(lower, distroPrefix); ok {
				installedLower[alias] = true
			}
		}
	}
	query := NewSearchQuery(terms, func(term string) bool { return installedLower[term] })

	var filtered []string
	for name := range installedOverlays {
		var alias string
		if distroLower != "" {
			if a, ok := strings.CutPrefix(strings.ToLower(name), distroPrefix); ok {
				alias = a
			}
		}
		if query.MatchesOrAlias(name, alias) {
			filtered = append(filtered, name)
		}
	}

	// Warn about terms that matched nothing in exact mode
	if query.Mode == SearchModeExact {
		for _, term := range query.Raw {
			found := false
			for _, name := range filtered {
				lower := strings.ToLower(name)
				if lower == term {
					found = true
					break
				}
				if distroLower != "" {
					if alias, ok := strings.CutPrefix(lower, distroPrefix); ok && alias == term {
						found = true
						break
					}
				}
			}
			if !found {
				utils.PrintWarning("Overlay %s not found among installed overlays.", utils.StyleName(term))
			}
		}
	}

	if len(filtered) == 0 {
		utils.PrintWarning("No matching installed overlays found.")
		return nil
	}

	return performDelete(cmd, filtered)
}

// performDelete handles duplicate detection, confirmation, and deletion of named overlays.
// names should be in "name/version" format. It reads --dir from cmd to restrict to a dir.
// Pass showListing=false to skip reprinting the overlay list (e.g. when list -d already showed it).
func performDelete(cmd *cobra.Command, names []string, showListing ...bool) error {
	printListing := len(showListing) == 0 || showListing[0]
	dirFilter, _ := cmd.Flags().GetString("dir")

	allCopies, err := getAllOverlayCopies()
	if err != nil {
		return err
	}

	installedOverlays, err := getInstalledOverlaysMap()
	if err != nil {
		return err
	}

	// Apply --dir filter
	if dirFilter != "" {
		known := filterImageDirs(config.GetImageSearchPaths(), dirFilter)
		if len(known) == 0 {
			utils.PrintWarning("No image directory matching %q found.", dirFilter)
			os.Exit(ExitCodeError)
		}
		knownSet := make(map[string]bool, len(known))
		for _, d := range known {
			knownSet[d] = true
		}
		for name, path := range installedOverlays {
			if !knownSet[filepath.Dir(path)] {
				delete(installedOverlays, name)
			}
		}
	}

	// Normalize and deduplicate; keep only names present in the (filtered) installed set
	seen := make(map[string]bool)
	var valid []string
	for _, name := range names {
		norm := utils.NormalizeNameVersion(name)
		if seen[norm] {
			continue
		}
		seen[norm] = true
		if _, ok := installedOverlays[norm]; ok {
			valid = append(valid, norm)
		}
	}
	if len(valid) == 0 {
		return nil
	}
	sort.Strings(valid)

	// When not scoped to a specific dir, refuse to remove overlays present in multiple dirs
	if dirFilter == "" {
		var ambiguous []string
		for _, name := range valid {
			if len(allCopies[name]) > 1 {
				ambiguous = append(ambiguous, name)
			}
		}
		if len(ambiguous) > 0 {
			utils.PrintWarning("The following overlays exist in multiple directories. Use --dir to specify which copy to remove:")
			for _, name := range ambiguous {
				fmt.Printf("  %s\n", utils.StyleName(name))
				for _, p := range allCopies[name] {
					fmt.Printf("    %s\n", filepath.Dir(p))
				}
			}
			return nil
		}
	}

	// Display grouped by directory
	byDir := make(map[string][]string)
	for _, name := range valid {
		dir := filepath.Dir(installedOverlays[name])
		byDir[dir] = append(byDir[dir], name)
	}
	if printListing {
		firstDir := true
		for _, dir := range config.GetImageSearchPaths() {
			dirNames := byDir[dir]
			if len(dirNames) == 0 {
				continue
			}
			if !firstDir {
				fmt.Println()
			}
			firstDir = false
			fmt.Println(dirHeader(dir))
			for _, name := range dirNames {
				fmt.Printf("  %s\n", utils.StyleName(name))
			}
		}
	}

	// Confirm
	fmt.Print("\n")
	if !utils.ShouldAnswerYes() {
		fmt.Printf("Remove %d overlay(s)? Cannot be undone. [y/N]: ", len(valid))
		choice, err := utils.ReadLineContext(cmd.Context())
		if err != nil || (choice != "y" && choice != "yes") {
			utils.PrintNote("Cancelled")
			return nil
		}
	}

	for _, name := range valid {
		overlayPath := installedOverlays[name]
		if !utils.CanWriteToDir(filepath.Dir(overlayPath)) {
			utils.PrintWarning("Overlay %s is in a read-only directory and cannot be removed.", utils.StyleName(name))
			continue
		}
		if lock, err := overlay.AcquireLock(overlayPath, true); err != nil {
			utils.PrintError("Cannot remove %s: overlay is currently in use", utils.StyleName(name))
			continue
		} else {
			lock.Close()
		}
		if err := os.Remove(overlayPath); err != nil {
			utils.PrintError("Failed to remove overlay %s: %v", utils.StyleName(name), err)
			continue
		}
		utils.PrintSuccess("Overlay %s removed.", utils.StyleName(name))
		envPath := overlayPath + ".env"
		if utils.FileExists(envPath) {
			if err := os.Remove(envPath); err != nil {
				utils.PrintDebug("Failed to remove env file %s: %v", envPath, err)
			}
		}
	}

	return nil
}

// removeExternalOverlays removes overlay files specified by direct path.
func removeExternalOverlays(cmd *cobra.Command, paths []string) error {
	// Resolve and validate each path
	var resolved []string
	for _, p := range paths {
		abs, err := filepath.Abs(p)
		if err != nil {
			utils.PrintWarning("Cannot resolve path %s: %v", p, err)
			continue
		}
		if !utils.FileExists(abs) {
			utils.PrintWarning("File not found: %s", abs)
			continue
		}
		resolved = append(resolved, abs)
	}
	if len(resolved) == 0 {
		return nil
	}

	// Display files to be removed
	fmt.Println("Overlay files to be removed:")
	for _, p := range resolved {
		fmt.Printf("  %s\n", utils.StyleName(p))
	}

	// Confirm
	fmt.Print("\n")
	if !utils.ShouldAnswerYes() {
		fmt.Printf("Remove %d overlay file(s)? Cannot be undone. [y/N]: ", len(resolved))
		choice, err := utils.ReadLineContext(cmd.Context())
		if err != nil || (choice != "y" && choice != "yes") {
			utils.PrintNote("Cancelled")
			return nil
		}
	}

	for _, p := range resolved {
		if !utils.CanWriteToDir(filepath.Dir(p)) {
			utils.PrintWarning("Overlay %s is in a read-only directory and cannot be removed.", utils.StyleName(filepath.Base(p)))
			continue
		}
		if lock, err := overlay.AcquireLock(p, true); err != nil {
			utils.PrintError("Cannot remove %s: overlay is currently in use", utils.StyleName(filepath.Base(p)))
			continue
		} else {
			lock.Close()
		}
		if err := os.Remove(p); err != nil {
			utils.PrintError("Failed to remove %s: %v", utils.StyleName(filepath.Base(p)), err)
			continue
		}
		utils.PrintSuccess("Overlay %s removed.", utils.StyleName(filepath.Base(p)))
		envPath := p + ".env"
		if utils.FileExists(envPath) {
			if err := os.Remove(envPath); err != nil {
				utils.PrintDebug("Failed to remove env file %s: %v", envPath, err)
			}
		}
	}
	return nil
}

// getAllOverlayCopies returns all paths for each overlay name across all image dirs.
// Unlike getInstalledOverlaysMap, this stores every copy, not just the highest-priority one.
func getAllOverlayCopies() (map[string][]string, error) {
	copies := make(map[string][]string)
	for _, imagesDir := range config.GetImageSearchPaths() {
		if !utils.DirExists(imagesDir) {
			continue
		}
		entries, err := os.ReadDir(imagesDir)
		if err != nil {
			continue
		}
		for _, entry := range entries {
			if entry.IsDir() || !utils.IsOverlay(entry.Name()) {
				continue
			}
			nameVersion := strings.TrimSuffix(entry.Name(), filepath.Ext(entry.Name()))
			normalized := strings.ReplaceAll(nameVersion, "--", "/")
			copies[normalized] = append(copies[normalized], filepath.Join(imagesDir, entry.Name()))
		}
	}
	return copies, nil
}

// getInstalledOverlaysMap returns a map of overlay name to path from all search paths.
// Only the highest-priority (first found) copy per name is kept.
func getInstalledOverlaysMap() (map[string]string, error) {
	overlays := make(map[string]string)
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
			if entry.IsDir() || !utils.IsOverlay(entry.Name()) {
				continue
			}
			nameVersion := strings.TrimSuffix(entry.Name(), filepath.Ext(entry.Name()))
			normalized := strings.ReplaceAll(nameVersion, "--", "/")
			if _, exists := overlays[normalized]; !exists {
				overlays[normalized] = filepath.Join(imagesDir, entry.Name())
			}
		}
	}
	return overlays, nil
}
