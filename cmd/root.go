package cmd

import (
	"fmt"
	"os"
	"strings"

	"github.com/Justype/condatainer/internal/apptainer"
	"github.com/Justype/condatainer/internal/config"
	"github.com/Justype/condatainer/internal/scheduler"
	"github.com/Justype/condatainer/internal/utils"
	"github.com/spf13/cobra"
)

var (
	debugMode bool
	localMode bool
	quietMode bool
)

var rootCmd = &cobra.Command{
	Use:           "condatainer",
	Short:         "CondaTainer: Use Apptainer/Conda/Overlays/SquashFS to manage tools/data/env for HPC users.",
	Version:       config.VERSION,
	SilenceErrors: true,

	PersistentPreRun: func(cmd *cobra.Command, args []string) {
		exe, err := os.Executable()
		if err != nil {
			utils.PrintError("Failed to determine executable path: %v", err)
			os.Exit(1)
		}

		// Step 1: Load defaults (paths, directories, etc.)
		config.LoadDefaults(exe)

		// Step 2: Initialize Viper (read config file, env vars)
		if err := config.InitViper(); err != nil {
			utils.PrintDebug("Error reading config file: %v", err)
		}

		// Step 3: Auto-detect binaries if needed and save to config
		updated, err := config.AutoDetectAndSave()
		if err != nil {
			utils.PrintDebug("Failed to save config: %v", err)
		} else if updated {
			if configPath, err := config.GetUserConfigPath(); err == nil {
				utils.PrintDebug("Auto-detected binaries saved to: %s", configPath)
			}
		}

		// Step 4: Load detected values from Viper into Global config
		config.LoadFromViper()

		// Warn if apptainer is still not accessible after auto-detection
		if !config.ValidateBinary(config.Global.ApptainerBin) {
			utils.PrintWarning("Apptainer not accessible. The module may have been unloaded or removed.")
			utils.PrintHint("Load the apptainer module: %s", utils.StyleAction("ml apptainer"))
			utils.PrintHint("Then run: %s", utils.StyleAction("condatainer config init"))
		}

		// Step 5: Apply command-line flags (highest priority)
		if debugMode {
			utils.DebugMode = true
			config.Global.Debug = true
			utils.PrintDebug("Debug mode enabled")
			utils.PrintDebug("CondaTainer Version: %s", utils.StyleInfo(config.VERSION))
			utils.PrintDebug("Executable: %s", exe)
			utils.PrintDebug("Base Directory: %s", config.Global.BaseDir)
			utils.PrintDebug("Images Directory: %s", config.GetImagesDir())
			utils.PrintDebug("Base Image: %s", config.Global.BaseImage)
			utils.PrintDebug("Apptainer Binary: %s", config.Global.ApptainerBin)
			if config.Global.SchedulerBin != "" {
				utils.PrintDebug("Scheduler Binary: %s", config.Global.SchedulerBin)
			}
		}

		if localMode {
			config.Global.SubmitJob = false
			utils.PrintDebug("Local mode enabled (job submission disabled)")
		}

		if quietMode {
			utils.QuietMode = true
			utils.PrintDebug("Quiet mode enabled (suppressing verbose messages)")
		}

		// Step 6: Auto-detect compression based on apptainer version
		if err := apptainer.SetBin(config.Global.ApptainerBin); err == nil {
			if version, err := apptainer.GetVersion(); err == nil {
				supportsZstd := apptainer.CheckZstdSupport(version)
				config.AutoDetectCompression(supportsZstd)
			}
		}

		// Step 7: Initialize scheduler if job submission is enabled
		if config.Global.SubmitJob {
			schedType, err := scheduler.Init(config.Global.SchedulerBin)
			if err == nil && schedType != scheduler.SchedulerUnknown {
				utils.PrintDebug("Scheduler initialized: %s", schedType)
			} else if err != nil {
				utils.PrintDebug("Scheduler not available: %v", err)
			}
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
	// Subcommands are attached to rootCmd in their respective init() functions
	rootCmd.PersistentFlags().BoolVar(&debugMode, "debug", false, "Enable debug mode with verbose output")
	rootCmd.PersistentFlags().BoolVar(&localMode, "local", false, "Disable job submission (run locally)")
	rootCmd.PersistentFlags().BoolVarP(&quietMode, "quiet", "q", false, "Suppress verbose overlay creation messages (errors/warnings still shown)")

	// Hide the help command from completions (use -h/--help instead)
	rootCmd.SetHelpCommand(&cobra.Command{Hidden: true})
}
