package cmd

import (
	"encoding/json"
	"fmt"
	"os"
	"strings"

	"github.com/spf13/cobra"

	"github.com/Justype/condatainer/internal/config"
	"github.com/Justype/condatainer/internal/utils"
)

var searchJSON bool
var searchFuzzy bool
var searchChannels []string
var searchLimit int

var searchCmd = &cobra.Command{
	Use:   "search [flags] <package>",
	Short: "Search conda packages via anaconda.org",
	Long: `Search conda packages using the anaconda.org API.
Results are filtered to the configured channels (or --channel overrides).`,
	Example: `  condatainer search samtools             # Exact match (first channel that has it)
  condatainer search samtools --json      # JSON output
  condatainer search -f samtool           # Fuzzy match (via anaconda search)
  condatainer search -f samtool -l 200    # Fuzzy with higher result limit
  condatainer search samtools -c bioconda # Search specific channel`,
	Args:         cobra.ExactArgs(1),
	SilenceUsage: true,
	RunE:         runSearch,
}

func init() {
	rootCmd.AddCommand(searchCmd)
	searchCmd.Flags().BoolVar(&searchJSON, "json", false, "Output results in JSON format")
	searchCmd.Flags().BoolVarP(&searchFuzzy, "fuzzy", "f", false, "Search package via api.anaconda.org/search")
	searchCmd.Flags().StringArrayVarP(&searchChannels, "channel", "c", nil, "Conda channel to search (overrides config; repeatable)")
	searchCmd.Flags().IntVarP(&searchLimit, "limit", "l", 100, "Number of fuzzy search results before filter")
}

func platformSupported(platform string, platforms []string) bool {
	for _, p := range platforms {
		if p == platform || p == "noarch" {
			return true
		}
	}
	return false
}

func runSearch(cmd *cobra.Command, args []string) error {
	query := args[0]
	channels := config.Global.Build.Channels
	if len(searchChannels) > 0 {
		channels = searchChannels
	}

	results, capped, err := utils.SearchCondaPackages(query, channels, searchFuzzy, searchLimit)
	if err != nil {
		return err
	}

	if len(results) == 0 {
		utils.PrintWarning("No packages found for %q in channels: %s", query, strings.Join(channels, ", "))
		os.Exit(ExitCodeError)
	}

	if searchJSON {
		enc := json.NewEncoder(os.Stdout)
		enc.SetIndent("", "  ")
		return enc.Encode(results)
	}

	currentPlatform := utils.CurrentCondaPlatform()
	for i, r := range results {
		if i > 0 {
			fmt.Println()
		}
		archSuffix := ""
		if currentPlatform != "" && len(r.Platforms) > 0 {
			hasPlatform, hasNoarch := false, false
			for _, p := range r.Platforms {
				switch p {
				case currentPlatform:
					hasPlatform = true
				case "noarch":
					hasNoarch = true
				}
			}
			var parts []string
			if hasPlatform {
				parts = append(parts, currentPlatform)
			}
			if hasNoarch {
				parts = append(parts, "noarch")
			}
			if len(parts) > 0 {
				archSuffix = " " + utils.StyleHint("["+strings.Join(parts, "/")+"]")
			}
		}
		fmt.Printf("%s (%s)%s\n", utils.StyleName(r.Name), utils.StyleHint(r.Channel), archSuffix)
		if r.Summary != "" {
			fmt.Printf("  %s\n", r.Summary)
		}
		if len(r.Versions) > 0 {
			printWrapped("Versions:", strings.Join(r.Versions, " "), 0)
		}
		if currentPlatform != "" && len(r.Platforms) > 0 && !platformSupported(currentPlatform, r.Platforms) {
			utils.PrintWarning("Package %s is not available for %s", r.Name, currentPlatform)
		}
	}

	if capped {
		utils.PrintWarning("Results may be incomplete: reached the limit of %d results. Use -l to increase.", searchLimit)
	}

	return nil
}
