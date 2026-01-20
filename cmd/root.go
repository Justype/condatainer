package cmd

import (
	"fmt"
	"os"
	"strings"

	"github.com/Justype/condatainer/internal/apptainer"
	"github.com/Justype/condatainer/internal/config"
	"github.com/Justype/condatainer/internal/utils"
	"github.com/spf13/cobra"
)

var (
	debugMode bool
	localMode bool
)

var rootCmd = &cobra.Command{
	Use:          "condatainer",
	Short:        "CondaTainer: Use Apptainer/Conda/Overlays/SquashFS to manage tools/data/env for HPC users.",
	Version:      config.VERSION,
	SilenceErrors: true,

	PersistentPreRun: func(cmd *cobra.Command, args []string) {
		exe, err := os.Executable()
		if err != nil {
			utils.PrintError("Failed to determine executable path: %v", err)
			os.Exit(1)
		}
		config.LoadDefaults(exe)

		if debugMode {
			utils.DebugMode = true
			config.Global.Debug = true
			utils.PrintDebug("Debug mode enabled")
			utils.PrintDebug("CondaTainer Version: %s", utils.StyleInfo(config.VERSION))
		}

		if localMode {
			config.Global.SubmitJob = false
			utils.PrintDebug("Local mode enabled (job submission disabled)")
		}
	},
}

func Execute() {
	if err := rootCmd.Execute(); err != nil {
		// Cobra's automatic error printing is silenced. For Apptainer errors
		// print only the captured output (trimmed) and exit with non-zero
		// status. For other errors, print the default error string.
		if ae, ok := err.(*apptainer.ApptainerError); ok {
			out := strings.TrimSpace(ae.Output)
			if out != "" {
				fmt.Fprintln(os.Stderr, out)
			}
			os.Exit(1)
		}
		fmt.Fprintln(os.Stderr, err)
		os.Exit(1)
	}
}

func init() {
	rootCmd.AddCommand(OverlayCmd)
	rootCmd.PersistentFlags().BoolVar(&debugMode, "debug", false, "Enable debug mode with verbose output")
	rootCmd.PersistentFlags().BoolVar(&localMode, "local", false, "Disable job submission (run locally)")
}
