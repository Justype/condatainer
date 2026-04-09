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
	Use:   "search <package>",
	Short: "Search conda packages via anaconda.org",
	Long: `Search conda packages using the anaconda.org API.
Results are filtered to the configured channels (or --channel overrides).`,
	Example: `  condatainer search samtools             # Exact match (first channel that has it)
  condatainer search samtools --json      # JSON output
  condatainer search -f samtool           # Fuzzy/substring match
  condatainer search -f samtool -l 200    # Fuzzy with higher result limit
  condatainer search samtools -c bioconda # Search specific channel`,
	Args:         cobra.MinimumNArgs(1),
	SilenceUsage: true,
	RunE:         runSearch,
}

func init() {
	rootCmd.AddCommand(searchCmd)
	searchCmd.Flags().BoolVar(&searchJSON, "json", false, "Output results in JSON format")
	searchCmd.Flags().BoolVarP(&searchFuzzy, "fuzzy", "f", false, "Substring match instead of exact name match")
	searchCmd.Flags().StringArrayVarP(&searchChannels, "channel", "c", nil, "Conda channel to search (overrides config; repeatable)")
	searchCmd.Flags().IntVarP(&searchLimit, "limit", "l", 100, "Maximum number of fuzzy search results")
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
	query := strings.Join(args, " ")
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
			archSuffix = " " + utils.StyleHint("["+currentPlatform+"]")
		}
		fmt.Printf("%s (%s)%s\n", utils.StyleName(r.Name), utils.StyleHint(r.Channel), archSuffix)
		if r.Summary != "" {
			fmt.Printf("  %s\n", r.Summary)
		}
		if len(r.Versions) > 0 {
			fmt.Printf("  Versions: %s\n", strings.Join(r.Versions, ", "))
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
