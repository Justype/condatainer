package cmd

import (
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
	availRemote bool
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
	availCmd.Flags().BoolVar(&availRemote, "remote", false, "Remote build scripts take precedence over local")
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
	IsRemote       bool                // true if from remote repository
	IsTemplate     bool                // true if script has #PL: placeholders
	Whatis         string              // description string
	TargetTemplate string              // raw #TARGET: pattern (e.g. "grch38/star/{star_version}/gencode{gencode_version}-{read_length}")
	PL             map[string][]string // placeholder values (all values for templates; single-element for expanded)
	PLOrder        []string            // placeholder key declaration order
}

func runAvail(cmd *cobra.Command, args []string) error {
	// Normalize search terms
	filters := normalizeFilters(args)

	// Get installed overlays for marking
	installedOverlays := getInstalledOverlays()

	// Gather script sources in precedence order.
	// When --remote: remote first, then local (remote wins on duplicates)
	// Default: local first, then remote (local wins on duplicates)
	type scriptSource struct {
		scripts  map[string]build.ScriptInfo
		isRemote bool
	}
	loadLocal := func() (scriptSource, error) {
		s, err := build.GetLocalBuildScripts()
		return scriptSource{scripts: s, isRemote: false}, err
	}
	loadRemote := func() scriptSource {
		s, err := build.GetRemoteBuildScripts()
		if err != nil {
			utils.PrintDebug("Failed to fetch remote scripts: %v", err)
			return scriptSource{scripts: map[string]build.ScriptInfo{}, isRemote: true}
		}
		return scriptSource{scripts: s, isRemote: true}
	}

	var sources []scriptSource
	if availRemote {
		local, err := loadLocal()
		if err != nil {
			return err
		}
		sources = []scriptSource{loadRemote(), local}
	} else {
		local, err := loadLocal()
		if err != nil {
			return err
		}
		sources = []scriptSource{local, loadRemote()}
	}

	scriptInfoToPackageInfo := func(name string, info build.ScriptInfo, isRemote bool) PackageInfo {
		return PackageInfo{
			Name:           name,
			Path:           info.Path,
			IsContainer:    info.IsContainer,
			IsInstalled:    installedOverlays[name],
			IsRemote:       isRemote,
			IsTemplate:     info.IsTemplate,
			Whatis:         info.Whatis,
			TargetTemplate: info.TargetTemplate,
			PL:             info.PL,
			PLOrder:        info.PLOrder,
		}
	}

	distroLower := strings.ToLower(config.Global.DefaultDistro)
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
	// back to FindBuildScript, which matches #TARGET: patterns without expanding.
	var availExactLookup func(string) bool
	if len(filters) > 1 {
		availExactLookup = func(term string) bool {
			for _, src := range sources {
				for name := range src.scripts {
					lower := strings.ToLower(name)
					if lower == term || aliasOf(name) == term {
						return true
					}
				}
			}
			_, found := build.FindBuildScript(term)
			return found
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
	for _, src := range sources {
		// Plain entries and template headers.
		for name, info := range src.scripts {
			if seen[name] {
				continue
			}
			seen[name] = true
			if matchesEntry(name, info.Whatis) {
				filtered = append(filtered, scriptInfoToPackageInfo(name, info, src.isRemote))
			}
		}

		if !availExpand {
			continue
		}

		// --expand: list variants individually. A matching template contributes all of
		// its variants; otherwise each variant must match on its own.
		for _, info := range src.scripts {
			if !info.IsTemplate {
				continue
			}
			templateMatched := matchesEntry(info.Name, info.Whatis)
			for name, variant := range build.ExpandTemplate(info) {
				if seen[name] {
					continue
				}
				seen[name] = true
				if templateMatched || matchesEntry(name, variant.Whatis) {
					filtered = append(filtered, scriptInfoToPackageInfo(name, variant, src.isRemote))
				}
			}
		}
	}

	if len(filtered) == 0 {
		utils.PrintWarning("No matching build scripts found.")
		// A term occurring only inside a variant name (e.g. "gencode47-101") cannot
		// match without --expand, so point at -e instead of leaving a bare miss.
		hasTemplates := false
		for _, src := range sources {
			for _, info := range src.scripts {
				if info.IsTemplate {
					hasTemplates = true
					break
				}
			}
			if hasTemplates {
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

// maxInlinePLValues is the largest number of concrete placeholder values listed
// inline. Longer lists collapse to a "first-last  (n values)" range summary;
// users can <Tab>-complete the individual values when prompted.
const maxInlinePLValues = 5

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
		line += "\n  " + pkg.Whatis
	}
	if pkg.TargetTemplate != "" {
		line += "\n  " + utils.StyleHint("→ "+pkg.TargetTemplate)
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
		line += "  " + pkg.Whatis
	}

	return line
}
