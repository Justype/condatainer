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

var searchCmd = &cobra.Command{
	Use:   "search <package>",
	Short: "Search conda packages via anaconda.org",
	Long: `Search conda packages using the anaconda.org API.
Results are filtered to the configured channels (or --channel overrides).`,
	Example: `  condatainer search samtools
  condatainer search samtools --json
  condatainer search samtools -c bioconda`,
	Args:         cobra.MinimumNArgs(1),
	SilenceUsage: true,
	RunE:         runSearch,
}

func init() {
	rootCmd.AddCommand(searchCmd)
	searchCmd.Flags().BoolVar(&searchJSON, "json", false, "Output results in JSON format")
	searchCmd.Flags().BoolVarP(&searchFuzzy, "fuzzy", "f", false, "Substring match instead of exact name match")
	searchCmd.Flags().StringArrayVarP(&searchChannels, "channel", "c", nil, "Conda channel to search (overrides config; repeatable)")
}

func runSearch(cmd *cobra.Command, args []string) error {
	query := strings.Join(args, " ")
	channels := config.Global.Build.Channels
	if len(searchChannels) > 0 {
		channels = searchChannels
	}

	results, capped, err := utils.SearchCondaPackages(query, channels, searchFuzzy)
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

	for i, r := range results {
		if i > 0 {
			fmt.Println()
		}
		fmt.Printf("%s (%s)\n", utils.StyleName(r.Name), utils.StyleHint(r.Channel))
		if r.Summary != "" {
			fmt.Printf("  %s\n", utils.StyleWhatis(r.Summary))
		}
		if len(r.Versions) > 0 {
			fmt.Printf("  Versions: %s\n", strings.Join(r.Versions, ", "))
		}
	}

	if capped {
		utils.PrintWarning("Results may be incomplete: anaconda.org search is capped at 100 entries.")
	}

	return nil
}
