package cmd

import (
	"os"

	"github.com/Justype/condatainer/internal/config"
	"github.com/Justype/condatainer/internal/utils"
	"github.com/spf13/cobra"
)

var (
	debugMode bool
	localMode bool
)

var rootCmd = &cobra.Command{
	Use:     "condatainer",
	Short:   "CondaTainer: Use Apptainer/Conda/SquashFS to manage tools for HPC users.",
	Version: config.VERSION,

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
			utils.PrintDebug("CondaTainer Version: %s", utils.StyleName(config.VERSION))
		}

		if localMode {
			config.Global.SubmitJob = false
			utils.PrintDebug("Local mode enabled (job submission disabled)")
		}
	},
}

func Execute() {
	if err := rootCmd.Execute(); err != nil {
		os.Exit(1)
	}
}

func init() {
	rootCmd.PersistentFlags().BoolVar(&debugMode, "debug", false, "Enable debug mode with verbose output")
	rootCmd.PersistentFlags().BoolVar(&localMode, "local", false, "Disable job submission (run locally)")
}
