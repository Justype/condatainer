package cmd

import (
	"context"
	"fmt"
	"os"
	"os/signal"
	"strings"
	"syscall"

	"github.com/Justype/condatainer/internal/apptainer"
	"github.com/Justype/condatainer/internal/config"
	"github.com/Justype/condatainer/internal/scheduler"
	"github.com/Justype/condatainer/internal/utils"
	"github.com/spf13/cobra"
	"github.com/spf13/pflag"
)

var (
	debugMode bool
	localMode bool
	quietMode bool
	yesMode   bool
)

var rootCmd = &cobra.Command{
	Use:           "condatainer",
	Short:         "CondaTainer: Use Apptainer/Conda/Overlays/SquashFS to manage tools/data/env for HPC users.",
	Version:       config.VERSION,
	SilenceErrors: true,

	PersistentPreRun: func(cmd *cobra.Command, args []string) {
		// Skip initialization entirely for completion script generation
		// (no config needed, and any stdout output corrupts the script)
		cmdName := cmd.Name()
		if cmdName == "completion" {
			return
		}

		// For __complete requests (tab-completion), suppress stdout messages
		// so they don't corrupt completion results
		isCompleteRequest := cmdName == "__complete" || cmdName == "__completeNoDesc"
		if isCompleteRequest {
			utils.QuietMode = true
		}

		exe, err := os.Executable()
		if err != nil {
			ExitWithError("Failed to determine executable path: %v", err)
		}

		// Step 1: Load defaults (paths, directories, etc.)
		config.LoadDefaults(exe)

		// Step 2: Initialize Viper (read config file, env vars)
		if err := config.InitViper(); err != nil {
			utils.PrintDebug("Error reading config file: %v", err)
		}

		// Step 3: Load config values into Global (with runtime detection fallback)
		config.LoadFromViper()

		// Warn if apptainer is still not accessible after auto-detection
		// (skip during completion requests to avoid polluting output)
		if !isCompleteRequest && !config.ValidateBinary(config.Global.ApptainerBin) {
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
			if writableDir, _ := config.GetWritableImagesDir(); writableDir != "" {
				utils.PrintDebug("Writable Images Directory: %s", writableDir)
			}
			utils.PrintDebug("Base Image: %s", config.GetBaseImage())
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

		if yesMode {
			utils.YesMode = true
			utils.PrintDebug("Yes mode enabled (automatically answering yes to prompts)")
		}

		// Step 6: Auto-detect compression based on apptainer version
		if err := apptainer.SetBin(config.Global.ApptainerBin); err == nil {
			if version, err := apptainer.GetVersion(); err == nil {
				supportsZstd := apptainer.CheckZstdSupport(version)
				config.AutoDetectCompression(supportsZstd)
			}
		}

		// Step 7: Apply scheduler spec defaults from config
		scheduler.SetSpecDefaults(scheduler.SpecDefaults{
			Ncpus:  config.Global.Scheduler.Ncpus,
			MemMB:  config.Global.Scheduler.MemMB,
			Time:   config.Global.Scheduler.Time,
			Nodes:  config.Global.Scheduler.Nodes,
			Ntasks: config.Global.Scheduler.Ntasks,
		})

		// Step 8: Initialize scheduler if job submission is enabled
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
	ctx, cancel := signal.NotifyContext(context.Background(), os.Interrupt, syscall.SIGTERM)
	defer cancel()

	if err := rootCmd.ExecuteContext(ctx); err != nil {
		// Cobra's automatic error printing is silenced. For Apptainer errors
		// print only the captured output (trimmed) and exit with non-zero
		// status. For other errors, print the default error string.
		if ae, ok := err.(*apptainer.ApptainerError); ok {
			out := strings.TrimSpace(ae.Output)
			if out != "" {
				fmt.Fprintln(os.Stderr, out)
			}
			os.Exit(ExitCodeError)
		}
		ExitWithError("%v", err)
	}
}

func init() {
	// Subcommands are attached to rootCmd in their respective init() functions
	rootCmd.PersistentFlags().BoolVar(&debugMode, "debug", false, "Enable debug mode with verbose output")
	rootCmd.PersistentFlags().BoolVar(&localMode, "local", false, "Disable job submission (run locally)")
	rootCmd.PersistentFlags().BoolVarP(&quietMode, "quiet", "q", false, "Suppress messages (warnings/errors are still shown)")
	rootCmd.PersistentFlags().BoolVarP(&yesMode, "yes", "y", false, "Automatically answer yes to all prompts")

	// Hide the help command from completions (use -h/--help instead)
	rootCmd.SetHelpCommand(&cobra.Command{Hidden: true})

	// Strip short flag shorthands for cleaner completions
	// This is done at init time so __complete also sees the stripped flags
	cobra.OnInitialize(func() {
		stripShortFlags(rootCmd)
	})
}

// stripShortFlags removes short flag shorthands from all commands
// to show only long flags in shell completion
func stripShortFlags(root *cobra.Command) {
	var walk func(c *cobra.Command)
	walk = func(c *cobra.Command) {
		// Strip shorthands from local flags
		c.LocalFlags().VisitAll(func(f *pflag.Flag) {
			f.Shorthand = ""
		})
		// Strip shorthands from persistent flags
		c.PersistentFlags().VisitAll(func(f *pflag.Flag) {
			f.Shorthand = ""
		})
		// Recursively walk child commands
		for _, child := range c.Commands() {
			walk(child)
		}
	}
	walk(root)
}
