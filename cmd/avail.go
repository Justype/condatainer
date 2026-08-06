package cmd

import (
	"fmt"
	"os"
	"path/filepath"
	"sort"
	"strings"

	"github.com/Justype/condatainer/catalog"
	"github.com/Justype/condatainer/internal/config"
	"github.com/Justype/condatainer/internal/utils"
	"github.com/spf13/cobra"
)

var (
	availExpand bool
	availWhatis bool
)

var availCmd = &cobra.Command{
	Use:     "avail [flags] [terms...]",
	Aliases: []string{"av"},
	Short:   "List available build scripts",
	Long: `List available build scripts (local and remote).

When duplicates exist, local scripts take precedence.

To install something you find here, use 'condatainer install <name>'.

` + searchSyntaxHint,
	Example: `  condatainer avail                   # List all
  condatainer avail cellranger        # Substring match
  condatainer avail cellranger 9      # AND search (multiple terms)
  condatainer avail 'cell*'           # Wildcard
  condatainer install cellranger/9.0.1 # Install one of the results`,
	SilenceUsage: true, // Runtime errors should not show usage
	RunE:         runAvail,
}

func init() {
	rootCmd.AddCommand(availCmd)
	availCmd.Flags().BoolVarP(&availExpand, "expand", "e", false, "Expand templated script groups into individual entries")
	availCmd.Flags().BoolVarP(&availWhatis, "whatis", "w", false,
		"Show description for each build script (default true, false with --expand); descriptions are searched when shown")
}

// PackageInfo holds information about a build script
type PackageInfo struct {
	Name           string              // name/version format
	Path           string              // full path to build script
	IsContainer    bool                // true if .def file
	IsInstalled    bool                // true if overlay exists
	Source         string              // the handle of the source it came from
	IsTemplate     bool                // true if the recipe has #PH: placeholders
	Whatis         string              // description string
	TargetTemplate string              // raw #TARGET: pattern (e.g. "grch38/star/{star_version}/gencode{gencode_version}-{read_length}")
	PH             map[string][]string // placeholder values (all values for templates; single-element for expanded)
	PHNames        []string            // placeholders in target order
}

func runAvail(cmd *cobra.Command, args []string) error {
	// Normalize search terms
	filters := normalizeFilters(args)

	// Get installed overlays for marking
	installedOverlays := getInstalledOverlays()

	// Sources are already in precedence order; Entries merges them first-wins.
	ctx := cmd.Context()
	cat, err := config.OpenCatalog(ctx)
	if err != nil {
		return err
	}
	entries := cat.Entries(ctx)
	// A missing source changes what the listing means, so say so rather than
	// presenting a short list as if it were complete.
	config.WarnUnreachableSources(ctx, cat)
	sourceOf := map[string]string{}
	for _, src := range cat {
		if found, err := src.Entries(ctx); err == nil {
			for name := range found {
				if _, taken := sourceOf[name]; !taken {
					sourceOf[name] = src.Name
				}
			}
		}
	}

	entryToPackageInfo := func(name string, e *catalog.Entry) PackageInfo {
		var phNames []string
		if e.IsTemplate {
			phNames = catalog.NewTemplate(e.TargetTemplate).Names()
		}
		return PackageInfo{
			Name:           name,
			Path:           e.Path,
			IsContainer:    strings.HasSuffix(e.Path, ".def"),
			IsInstalled:    installedOverlays[name],
			Source:         sourceOf[e.Name],
			IsTemplate:     e.IsTemplate,
			Whatis:         e.Whatis,
			TargetTemplate: e.TargetTemplate,
			PH:             e.PH,
			PHNames:        phNames,
		}
	}

	distroLower := strings.ToLower(config.ResolvedBase())
	distroPrefix := distroLower + "/"
	aliasOf := func(name string) string {
		if distroLower == "" {
			return ""
		}
		a, _ := strings.CutPrefix(strings.ToLower(name), distroPrefix)
		return a
	}

	// Build search query (exact-first: check if the first term names a known script).
	// Templates are unexpanded here, so a concrete variant name misses the map — fall
	// back to Lookup, which matches #TARGET: patterns without expanding.
	var availExactLookup func(string) bool
	if len(filters) > 1 {
		availExactLookup = func(term string) bool {
			for name := range entries {
				lower := strings.ToLower(name)
				if lower == term || aliasOf(name) == term {
					return true
				}
			}
			_, found, err := cat.Lookup(ctx, term)
			return err == nil && found
		}
	}
	query := NewSearchQuery(filters, availExactLookup)

	// Descriptions are shown by default, but not under --expand where they would repeat
	// on every variant. An explicit -w/--whatis=false overrides either default.
	showWhatis := !availExpand
	if cmd.Flags().Changed("whatis") {
		showWhatis = availWhatis
	}

	// Match names always, descriptions only while they are shown, so every hit is
	// visible in the output. Placeholder keys and values are never matched.
	matches := func(name string) bool {
		return query.MatchesOrAlias(name, aliasOf(name))
	}
	matchesEntry := func(name, whatis string) bool {
		if !showWhatis {
			return matches(name)
		}
		return query.MatchesOrAliasWithText(name, aliasOf(name), whatis)
	}

	// Collect matches. Without --expand only templates and plain entries are considered,
	// so templates stay collapsed to a group header and nothing is expanded.
	filtered := make([]PackageInfo, 0)
	seen := make(map[string]bool)
	// Plain entries and template headers.
	for name, e := range entries {
		if seen[name] {
			continue
		}
		seen[name] = true
		if matchesEntry(name, e.Whatis) {
			filtered = append(filtered, entryToPackageInfo(name, e))
		}
	}

	// --expand: list variants individually. A matching template contributes all of
	// its variants; otherwise each variant must match on its own.
	if availExpand {
		for _, e := range entries {
			if !e.IsTemplate {
				continue
			}
			templateMatched := matchesEntry(e.Name, e.Whatis)
			for _, variant := range catalog.NewTemplate(e.TargetTemplate).Enumerate(e.PH) {
				if seen[variant.Name] {
					continue
				}
				seen[variant.Name] = true
				whatis := interpolate(e.Whatis, variant.Vars)
				if templateMatched || matchesEntry(variant.Name, whatis) {
					pi := entryToPackageInfo(variant.Name, e)
					pi.Whatis = whatis
					pi.IsTemplate = false
					pi.PH = singleValues(variant.Vars)
					filtered = append(filtered, pi)
				}
			}
		}
	}

	if len(filtered) == 0 {
		utils.PrintWarning("No matching build scripts found.")
		// A term occurring only inside a variant name (e.g. "gencode47-101") cannot
		// match without --expand, so point at -e instead of leaving a bare miss.
		hasTemplates := false
		for _, e := range entries {
			if e.IsTemplate {
				hasTemplates = true
				break
			}
		}
		if !availExpand && len(filters) > 0 && hasTemplates {
			utils.PrintNote("Templates are listed collapsed — use %s to search individual variants.",
				utils.StyleHint("-e"))
		}
		return nil
	}

	// Sort by name
	sort.Slice(filtered, func(i, j int) bool {
		return filtered[i].Name < filtered[j].Name
	})

	// --expand prints variants only, skipping the template headers they came from;
	// otherwise print template group headers alongside plain entries.
	if availExpand {
		for _, pkg := range filtered {
			if !pkg.IsTemplate {
				fmt.Println(formatPackageLine(pkg, showWhatis))
			}
		}
	} else {
		for _, pkg := range filtered {
			if pkg.IsTemplate {
				fmt.Println(formatTemplateLine(pkg, showWhatis))
			} else {
				fmt.Println(formatPackageLine(pkg, showWhatis))
			}
		}
	}

	// Summarize per source when more than one contributed.
	if len(cat) > 1 {
		counts := map[string]int{}
		for _, pkg := range filtered {
			counts[pkg.Source]++
		}
		if len(counts) > 1 {
			parts := make([]string, 0, len(counts))
			for _, src := range cat {
				if n := counts[src.Name]; n > 0 {
					parts = append(parts, fmt.Sprintf("%s from %s", utils.StyleNumber(n), src.Name))
				}
			}
			fmt.Println()
			utils.PrintMessage("Found %s.", strings.Join(parts, ", "))
		}
	}

	return nil
}

// interpolate substitutes {key} tokens, used to render a variant's description.
func interpolate(text string, vars map[string]string) string {
	for k, v := range vars {
		text = strings.ReplaceAll(text, "{"+k+"}", v)
	}
	return text
}

// singleValues records the one value chosen per placeholder for a variant.
func singleValues(vars map[string]string) map[string][]string {
	out := make(map[string][]string, len(vars))
	for k, v := range vars {
		out[k] = []string{v}
	}
	return out
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

// maxInlinePLValues is the largest number of concrete placeholder values listed
// inline. Longer lists collapse to a "first-last  (n values)" range summary;
// users can <Tab>-complete the individual values when prompted.
const maxInlinePLValues = 5

// formatTemplateLine formats a template (PH) package as a collapsed group header.
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
	for _, vals := range pkg.PH {
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
	if variantCount == 1 && len(pkg.PH) == 0 {
		variantCount = 0
	}

	var labelParts []string
	if variantCount > 0 {
		labelParts = append(labelParts, utils.StyleDebug(fmt.Sprintf("%d variants", variantCount)))
	} else {
		labelParts = append(labelParts, utils.StyleDebug("template"))
	}
	if pkg.Source != "" {
		labelParts = append(labelParts, utils.StyleHint(pkg.Source))
	}

	line := fmt.Sprintf("%s  %s%s%s", utils.StyleName(pkg.Name),
		utils.StyleDebug("["), strings.Join(labelParts, utils.StyleDebug(", ")), utils.StyleDebug("]"))
	if showWhatis && pkg.Whatis != "" {
		line += "\n  " + pkg.Whatis
	}
	if pkg.TargetTemplate != "" {
		line += "\n  " + utils.StyleHint("→ "+pkg.TargetTemplate)
	}

	// Show placeholder value summaries with aligned keys.
	// First pass: compute max key length for alignment.
	maxKeyLen := 0
	for _, key := range pkg.PHNames {
		if _, ok := pkg.PH[key]; ok && len(key) > maxKeyLen {
			maxKeyLen = len(key)
		}
	}
	for _, key := range pkg.PHNames {
		vals, ok := pkg.PH[key]
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
		if n > maxInlinePLValues {
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

	// Compute alias before highlighting (e.g. "ubuntu24/build-essential" → "[build-essential]")
	var alias string
	if distro := config.ResolvedBase(); distro != "" {
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
	if pkg.Source != "" {
		suffixes = append(suffixes, utils.StyleHint(pkg.Source))
	}

	if len(suffixes) > 0 {
		line = fmt.Sprintf("%s (%s)", line, strings.Join(suffixes, ", "))
	}

	if alias != "" {
		line += "  " + utils.StyleInfo("["+alias+"]")
	}

	if showWhatis && pkg.Whatis != "" {
		line += "  " + pkg.Whatis
	}

	return line
}
