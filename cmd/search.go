package cmd

import (
	"github.com/spf13/cobra"

	"github.com/Justype/condatainer/internal/config"
	"github.com/Justype/condatainer/internal/exec"
)

var searchJSON bool
var searchPretty bool

var searchCmd = &cobra.Command{
	Use:   "search <package>",
	Short: "Search conda packages via micromamba",
	Long: `Search conda packages using micromamba inside the base image.
Runs: micromamba search <package> --no-rc -c <ch1> -c <ch2> ...`,
	Example: `  condatainer search samtools
  condatainer search 'samtools>1.10'
  condatainer search samtools --json
  condatainer search samtools --pretty`,
	Args:         cobra.MinimumNArgs(1),
	SilenceUsage: true,
	RunE:         runSearch,
}

func init() {
	rootCmd.AddCommand(searchCmd)
	searchCmd.Flags().BoolVar(&searchJSON, "json", false, "Output results in JSON format")
	searchCmd.Flags().BoolVar(&searchPretty, "pretty", false, "Pretty-print output")
}

func runSearch(cmd *cobra.Command, args []string) error {
	// Build: micromamba search <args...> --no-rc -c ch1 -c ch2 [--json]
	mmArgs := []string{"micromamba", "search"}
	mmArgs = append(mmArgs, args...)
	mmArgs = append(mmArgs, "--no-rc")
	for _, ch := range config.Global.Build.Channels {
		mmArgs = append(mmArgs, "-c", ch)
	}
	if searchJSON {
		mmArgs = append(mmArgs, "--json")
	}
	if searchPretty {
		mmArgs = append(mmArgs, "--pretty")
	}

	opts := exec.Options{
		BaseImage:    config.GetBaseImage(),
		ApptainerBin: config.Global.ApptainerBin,
		Command:      mmArgs,
		HidePrompt:   true,
	}

	return exec.Run(cmd.Context(), opts)
}
