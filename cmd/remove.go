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
	Example: `  condatainer rm cellranger/9.0.1    # Remove exact version
  condatainer rm cellranger          # Exact match: all cellranger
  condatainer rm 'cell*'             # Wildcard
  condatainer rm cellranger/9.0.1 cellranger/8.0.1 # Remove multiple
  condatainer rm cellranger 9        # AND search (multiple terms)`,
	Args:         cobra.MinimumNArgs(1),
	SilenceUsage: true, // Runtime errors should not show usage
	RunE:         runRemove,
}

func init() {
	rootCmd.AddCommand(removeCmd)
}

func runRemove(cmd *cobra.Command, args []string) error {
	// Normalize search terms
	terms := normalizeFilters(args)

	// Get installed overlays
	installedOverlays, err := getInstalledOverlaysMap()
	if err != nil {
		return err
	}

	if len(installedOverlays) == 0 {
		utils.PrintWarning("No installed overlays found.")
		return nil
	}

	// Build search query (exact-first detection, including distro-prefix aliases)
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

	// Filter overlays; in exact mode also accept distro-prefix alias
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

	// In exact mode, warn for any term that matched nothing (check full name and alias)
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

	// Sort for consistent display
	sort.Strings(filtered)

	// Display overlays to be removed
	fmt.Println("Overlays to be removed:")
	for _, name := range filtered {
		fmt.Printf(" - %s\n", name)
	}

	// Ask for confirmation (skip if --yes flag is set)
	fmt.Print("\n")
	if !utils.ShouldAnswerYes() {
		fmt.Print("Are you sure? Cannot be undone. [y/N]: ")
		choice, err := utils.ReadLineContext(cmd.Context())
		if err != nil || (choice != "y" && choice != "yes") {
			utils.PrintNote("Cancelled")
			return nil
		}
	}

	// Remove overlays (only if writable)
	writableDir, err := config.GetWritableImagesDir()
	if err != nil {
		return fmt.Errorf("no writable images directory found: %w", err)
	}

	for _, name := range filtered {
		overlayPath := installedOverlays[name]

		// Check if overlay is in writable directory
		if !strings.HasPrefix(overlayPath, writableDir) {
			utils.PrintWarning("Overlay %s is in a read-only directory and cannot be removed.", utils.StyleName(name))
			continue
		}

		// Check if overlay is currently in use (a running exec/run holds a shared lock)
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

		// Also remove .env file if exists
		envPath := overlayPath + ".env"
		if utils.FileExists(envPath) {
			if err := os.Remove(envPath); err != nil {
				utils.PrintDebug("Failed to remove env file %s: %v", envPath, err)
			}
		}
	}

	return nil
}

// getInstalledOverlaysMap returns a map of overlay name to path from all search paths
func getInstalledOverlaysMap() (map[string]string, error) {
	overlays := make(map[string]string)

	// Search all image directories
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
			// Convert samtools--1.21.sqf to samtools/1.21
			normalized := strings.ReplaceAll(nameVersion, "--", "/")

			// Only store first occurrence (higher priority directory wins)
			if _, exists := overlays[normalized]; !exists {
				overlays[normalized] = filepath.Join(imagesDir, entry.Name())
			}
		}
	}

	return overlays, nil
}
