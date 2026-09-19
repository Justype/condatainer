package cmd

import (
	"context"
	"errors"
	"fmt"
	"log/slog"
	"os"
	"os/signal"
	"strings"
	"syscall"

	"github.com/Justype/condatainer/cmd/internal/clilog"
	"github.com/Justype/condatainer/internal/build"
	"github.com/Justype/condatainer/internal/config"
	"github.com/Justype/condatainer/internal/libexec"
	"github.com/Justype/condatainer/internal/logging"
	"github.com/Justype/condatainer/internal/runtime/apptainer"
	execpkg "github.com/Justype/condatainer/internal/runtime/exec"
	"github.com/Justype/condatainer/internal/scheduler"
	"github.com/Justype/condatainer/internal/utils"
	"github.com/spf13/cobra"
)

var (
	debugMode    bool
	noSubmitMode bool
	quietMode    bool
	yesMode      bool
)

var rootCmd = &cobra.Command{
	Use:           "condatainer",
	Short:         "Single-file Conda env, tools, and data for HPC — plus app helpers to run RStudio and more.",
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

		// Environment management runs through this same binary inside the
		// container. It needs command IO and cancellation, but none of the host
		// build, scheduler, Apptainer, or persisted-config initialization below.
		if strings.HasPrefix(cmd.CommandPath(), "condatainer env") {
			initializeLightweightCommand(cmd)
			return
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

		// Warn if apptainer is still not accessible after auto-detection, unless
		// the toolchain has its own. Skip for all `config` commands so users can
		// inspect/repair config without seeing contradictory warnings before
		// re-detection runs.
		isConfigCommand := strings.HasPrefix(cmd.CommandPath(), "condatainer config")
		if !isCompleteRequest && !isConfigCommand && !libexec.Installed("apptainer") &&
			!config.ValidateBinary(config.Global.Build.SystemApptainer) {
			utils.PrintWarning("Apptainer not accessible. The module may have been unloaded or removed.")
			utils.PrintHint("Run: %s", utils.StyleAction("condatainer config init"))
		}

		// Step 5: Apply command-line flags (highest priority)
		if debugMode {
			utils.DebugMode = true
			config.Global.Debug = true
			utils.PrintDebug("Debug mode enabled")
			utils.PrintDebug("CondaTainer Version: %s", utils.StyleInfo(config.VERSION))
			utils.PrintDebug("Executable: %s", exe)
			if base, err := config.GetBaseImage(); err == nil {
				utils.PrintDebug("Base Image: %s", base)
			} else {
				utils.PrintDebug("Base Image: %v", err)
			}
			utils.PrintDebug("Apptainer Binary: %s", config.Global.Build.SystemApptainer)
			if config.Global.Scheduler.Bin != "" {
				utils.PrintDebug("Scheduler Binary: %s", config.Global.Scheduler.Bin)
			}
		}

		if noSubmitMode {
			config.Global.SubmitJob = false
			utils.PrintDebug("Job submission disabled (--no-submit)")
		}

		if quietMode {
			utils.QuietMode = true
			utils.PrintDebug("Quiet mode enabled (suppressing verbose messages)")
		}

		if yesMode {
			utils.YesMode = true
			utils.PrintDebug("Yes mode enabled (automatically answering yes to prompts)")
		}

		// Route all slog output (both context-based and slog.Default()) through
		// the CLI's utils.Print* functions so internal packages produce the
		// familiar [CNT] output instead of the raw "2006/01/02 INFO ..." format.
		cliHandler := slog.New(clilog.New())
		slog.SetDefault(cliHandler)
		cmd.SetContext(logging.WithLogger(cmd.Context(), cliHandler))
		cmd.SetContext(execpkg.WithIO(cmd.Context(), execpkg.IO{
			Stdin:  os.Stdin,
			Stdout: os.Stdout,
			Stderr: os.Stderr,
		}))

		// Step 7: Apply debug mode and resource defaults from config
		scheduler.SetDebugMode(config.Global.Debug)
		build.SetBuildDefaults(config.Global.Build.Defaults)
		scheduler.DefaultCommandTimeout = config.Global.Scheduler.Timeout

		// Step 8: Initialize scheduler if job submission is enabled
		if config.Global.SubmitJob {
			schedType, err := scheduler.Init(config.Global.Scheduler.Bin)
			if err == nil && schedType != scheduler.SchedulerUnknown {
				utils.PrintDebug("Scheduler initialized: %s", schedType)
			} else if err != nil {
				utils.PrintDebug("Scheduler not available: %v", err)
			}
		}

	},
}

func sourceHandleCompletion(cmd *cobra.Command, args []string, toComplete string) ([]string, cobra.ShellCompDirective) {
	seen := make(map[string]bool)
	if flag := cmd.Flags().Lookup("source"); flag != nil {
		if values, err := cmd.Flags().GetStringArray("source"); err == nil {
			for _, value := range values {
				seen[value] = true
			}
		}
	}
	var out []string
	for _, source := range config.Global.Sources {
		if !seen[source.Name] && strings.HasPrefix(source.Name, toComplete) {
			out = append(out, source.Name)
		}
	}
	return out, cobra.ShellCompDirectiveNoFileComp
}

func Execute() {
	ctx, cancel := signal.NotifyContext(context.Background(), os.Interrupt, syscall.SIGTERM)
	defer cancel()

	// End an in-place progress line on interrupt, so what follows starts on a
	// fresh line instead of on the same one as the "^C" echo. Routed through the
	// handler rather than writing a newline here: it knows whether a line is
	// actually open, where an unconditional newline leaves a blank one behind
	// whenever a log record follows and closes the line itself.
	interruptCh := make(chan os.Signal, 1)
	signal.Notify(interruptCh, os.Interrupt, syscall.SIGTERM)
	defer signal.Stop(interruptCh)
	go func() {
		if _, ok := <-interruptCh; ok {
			clilog.EndProgressLine()
		}
	}()

	if err := rootCmd.ExecuteContext(ctx); err != nil {
		// Cobra's automatic error printing is silenced. For Apptainer errors
		// print only the captured output (trimmed) and exit with non-zero
		// status. For other errors, print the default error string.
		// Ctrl+C or timeout: ^C on the terminal is feedback enough.
		if errors.Is(err, context.Canceled) || errors.Is(err, context.DeadlineExceeded) {
			os.Exit(ExitCodeError)
		}
		if ae, ok := err.(*apptainer.ApptainerError); ok {
			out := strings.TrimSpace(ae.Output)
			if out != "" {
				fmt.Fprintln(os.Stderr, out)
			}
			os.Exit(ExitCodeError)
		}
		if exitErr, ok := err.(interface{ ExitCode() int }); ok {
			if code := exitErr.ExitCode(); code >= 0 {
				os.Exit(code)
			}
		}
		ExitWithError("%v", err)
	}
}

func initializeLightweightCommand(cmd *cobra.Command) {
	if debugMode {
		utils.DebugMode = true
	}
	if quietMode {
		utils.QuietMode = true
	}
	if yesMode {
		utils.YesMode = true
	}
	cliHandler := slog.New(clilog.New())
	slog.SetDefault(cliHandler)
	cmd.SetContext(logging.WithLogger(cmd.Context(), cliHandler))
	cmd.SetContext(execpkg.WithIO(cmd.Context(), execpkg.IO{
		Stdin:  os.Stdin,
		Stdout: os.Stdout,
		Stderr: os.Stderr,
	}))
}

func init() {
	// Subcommands are attached to rootCmd in their respective init() functions
	rootCmd.PersistentFlags().BoolVar(&debugMode, "debug", false, "Enable debug mode with verbose output")
	rootCmd.PersistentFlags().BoolVarP(&quietMode, "quiet", "q", false, "Suppress messages (warnings/errors are still shown)")
	rootCmd.PersistentFlags().BoolVarP(&yesMode, "yes", "y", false, "Automatically answer yes to all prompts")

	// Hide the help command from completions (use -h/--help instead)
	rootCmd.SetHelpCommand(&cobra.Command{Hidden: true})
}
