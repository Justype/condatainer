package cmd

import (
	"fmt"
	"os"
	"path/filepath"
	"sort"
	"strings"

	"github.com/spf13/cobra"

	"github.com/Justype/condatainer/internal/config"
	"github.com/Justype/condatainer/internal/utils"
)

var listCmd = &cobra.Command{
	Use:     "list [terms...]",
	Aliases: []string{"ls"},
	Short:   "List installed overlays matching search terms",
	RunE:    runList,
}

func init() {
	rootCmd.AddCommand(listCmd)
}

func runList(cmd *cobra.Command, args []string) error {
	filters := normalizeFilters(args)

	appOverlays, err := collectAppOverlays(filters)
	if err != nil {
		return err
	}
	refOverlays, err := collectReferenceOverlays(filters)
	if err != nil {
		return err
	}

	printed := false
	if len(appOverlays) > 0 {
		printed = true
		fmt.Println("Available app overlays:")
		names := make([]string, 0, len(appOverlays))
		for name := range appOverlays {
			names = append(names, name)
		}
		sort.Strings(names)
		nameWidth := maxWidth(names)
		for _, name := range names {
			vers := appOverlays[name]
			colored := make([]string, 0, len(vers))
			for _, v := range vers {
				if v == "(system app overlay)" {
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
	}

	return nil
}

func collectAppOverlays(filters []string) (map[string][]string, error) {
	if !utils.DirExists(config.Global.ImagesDir) {
		return nil, nil
	}

	entries, err := os.ReadDir(config.Global.ImagesDir)
	if err != nil {
		return nil, fmt.Errorf("unable to list app overlays: %w", err)
	}

	grouped := map[string]map[string]struct{}{}
	for _, entry := range entries {
		if entry.IsDir() {
			continue
		}
		if !utils.IsOverlay(entry.Name()) {
			continue
		}
		delimCount := strings.Count(entry.Name(), "--")
		if delimCount > 1 {
			// Not an app overlay
			continue
		}
		nameVersion := strings.TrimSuffix(entry.Name(), filepath.Ext(entry.Name()))
		normalized := strings.ToLower(utils.NormalizeNameVersion(nameVersion))
		if !matchesFilters(normalized, filters) {
			continue
		}

		var name, version string
		if strings.Contains(nameVersion, "--") {
			parts := strings.SplitN(nameVersion, "--", 2)
			name = parts[0]
			version = parts[1]
		} else {
			name = nameVersion
			version = "(system app/env)"
		}
		if name == "" {
			continue
		}

		if grouped[name] == nil {
			grouped[name] = map[string]struct{}{}
		}
		grouped[name][version] = struct{}{}
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

func collectReferenceOverlays(filters []string) ([]string, error) {
	if !utils.DirExists(config.Global.ImagesDir) {
		return nil, nil
	}

	entries, err := os.ReadDir(config.Global.ImagesDir)
	if err != nil {
		return nil, fmt.Errorf("unable to list reference overlays: %w", err)
	}

	names := make([]string, 0, len(entries))
	for _, entry := range entries {
		if entry.IsDir() {
			continue
		}
		if !utils.IsOverlay(entry.Name()) {
			continue
		}
		delimCount := strings.Count(entry.Name(), "--")
		if delimCount <= 1 {
			// Not a reference overlay
			continue
		}
		nameVersion := strings.TrimSuffix(entry.Name(), filepath.Ext(entry.Name()))
		normalized := strings.ToLower(utils.NormalizeNameVersion(nameVersion))
		if !matchesFilters(normalized, filters) {
			continue
		}
		displayName := strings.ReplaceAll(nameVersion, "--", "/")
		names = append(names, displayName)
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

func matchesFilters(name string, filters []string) bool {
	if len(filters) == 0 {
		return true
	}
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
