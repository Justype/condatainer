package cmd

import (
	"fmt"
	"os"
	"path/filepath"
	"regexp"
	"sort"
	"strings"

	"github.com/Justype/condatainer/internal/config"
	"github.com/Justype/condatainer/internal/utils"
	"github.com/spf13/cobra"
)

var removeCmd = &cobra.Command{
	Use:     "remove [terms...]",
	Aliases: []string{"rm", "delete", "uninstall"},
	Short:   "Remove installed overlays matching search terms",
	Long: `Remove installed overlays by exact name or search terms.

Multiple terms use AND logic - all terms must match for an overlay to be selected.
The command will ask for confirmation before removing overlays.`,
	Example: `  condatainer remove samtools/1.22      # Remove exact overlay
  condatainer remove samtools           # Remove all samtools versions
  condatainer rm sam 1.22               # Search with multiple terms (AND logic)`,
	Args:         cobra.MinimumNArgs(1),
	SilenceUsage: true, // Runtime errors should not show usage
	RunE:         runRemove,
}

func init() {
	rootCmd.AddCommand(removeCmd)
}

func runRemove(cmd *cobra.Command, args []string) error {
	// Normalize search terms
	terms := make([]string, len(args))
	for i, term := range args {
		terms[i] = utils.NormalizeNameVersion(term)
	}

	// Get installed overlays
	installedOverlays, err := getInstalledOverlaysMap()
	if err != nil {
		return err
	}

	if len(installedOverlays) == 0 {
		utils.PrintWarning("No installed overlays found.")
		return nil
	}

	// Filter overlays
	var filtered []string
	firstTermNormalized := strings.ToLower(terms[0])

	// Check if first term is exact match
	exactMatch := false
	for name := range installedOverlays {
		if strings.ToLower(name) == firstTermNormalized {
			exactMatch = true
			break
		}
	}

	if exactMatch {
		// Exact match mode - treat terms as exact names
		for _, term := range terms {
			termLower := strings.ToLower(term)
			found := false
			for name := range installedOverlays {
				if strings.ToLower(name) == termLower {
					filtered = append(filtered, name)
					found = true
					break
				}
			}
			if !found {
				utils.PrintWarning("Overlay %s not found among installed overlays.", utils.StyleName(term))
			}
		}
	} else {
		// Search mode - match all terms (AND logic)
		for name := range installedOverlays {
			nameLower := strings.ToLower(name)
			matchAll := true
			for _, term := range terms {
				termLower := strings.ToLower(term)
				if !strings.Contains(nameLower, termLower) {
					matchAll = false
					break
				}
			}
			if matchAll {
				filtered = append(filtered, name)
			}
		}
	}

	if len(filtered) == 0 {
		utils.PrintWarning("No matching installed overlays found.")
		return nil
	}

	// Sort for consistent display
	sort.Strings(filtered)

	// Display overlays to be removed with highlighted terms
	fmt.Println("Overlays to be removed:")
	for _, name := range filtered {
		highlighted := name

		// Use single-pass highlighting to avoid ANSI code interference
		if len(terms) > 0 {
			// Sort terms by length (longest first) to prefer longer matches
			sortedTerms := make([]string, len(terms))
			copy(sortedTerms, terms)
			sort.Slice(sortedTerms, func(i, j int) bool {
				return len(sortedTerms[i]) > len(sortedTerms[j])
			})

			// Combine all terms into one regex with alternation
			patterns := make([]string, len(sortedTerms))
			for i, term := range sortedTerms {
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
