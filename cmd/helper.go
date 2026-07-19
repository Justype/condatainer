package cmd

import (
	"context"
	"errors"
	"fmt"
	"os"
	"path/filepath"
	"strconv"
	"strings"
	"time"

	ui "github.com/Justype/condatainer/cmd/internal/ui"
	"github.com/Justype/condatainer/internal/apptainer"
	"github.com/Justype/condatainer/internal/config"
	"github.com/Justype/condatainer/internal/helper"
	"github.com/Justype/condatainer/internal/scheduler"
	"github.com/Justype/condatainer/internal/utils"
	"github.com/spf13/cobra"
)

type postScriptHelperFlags struct {
	cpusSet    bool
	cpus       int
	memSet     bool
	mem        string
	timeSet    bool
	time       string
	gpuSet     bool
	gpu        string
	baseSet    bool
	base       string
	envSet     bool
	env        string
	overlaySet bool
	overlays   []string
	cwdSet     bool
	cwd        string
	newSet     bool
}

var (
	helperPath   bool
	helperList   bool
	helperUpdate bool
	helperStatus bool
)

var helperCmd = &cobra.Command{
	Use:   "helper [flags] [script-name] [script-args...]",
	Short: "Run web apps like RStudio on HPC",
	Long: `Run web apps inside CondaTainer on HPC.

A helper script submits a job, starts the app, and opens an SSH tunnel to your browser.
  - With no arguments, shows the apps you have running.
  - Every helper has its own flags: condatainer helper <name> -h
  - Resources (cpus, mem, time, gpu) can be set per run or saved as defaults.

Note: not available inside a container or a scheduler job.`,
	Example: `  condatainer helper                          # Show running apps
  condatainer helper code-server              # Start code-server
  condatainer helper code-server -h           # Show its flags and defaults
  condatainer helper code-server config       # Show or save its defaults
  condatainer helper stop                     # Stop a running app
  condatainer helper --update                 # Update helper scripts from remote`,
	SilenceUsage:      true,
	RunE:              runHelper,
	ValidArgsFunction: completeHelperScripts,
}

func init() {
	rootCmd.AddCommand(helperCmd)
	helperCmd.Flags().BoolVar(&helperPath, "path", false, "Show path to helper scripts directory")
	helperCmd.Flags().BoolVarP(&helperList, "list", "l", false, "List available helper scripts")
	helperCmd.Flags().BoolVarP(&helperUpdate, "update", "u", false, "Update helper scripts from remote")
	helperCmd.Flags().BoolVar(&helperStatus, "status", false, "Show status of running helpers")
	helperCmd.Flags().MarkHidden("status") //nolint:errcheck

	// Stop flag parsing after the first positional argument so helper-specific flags
	// (e.g. -c/--cpus) are passed through to parsePostScriptHelperFlags rather than cobra.
	helperCmd.Flags().SetInterspersed(false)
}

// completeHelperScripts provides shell completion for helper script names
func completeHelperScripts(cmd *cobra.Command, args []string, toComplete string) ([]string, cobra.ShellCompDirective) {
	// Only complete the first argument (script name), then return to default file completion
	if len(args) > 0 {
		return nil, cobra.ShellCompDirectiveDefault
	}

	// Collect scripts from all search paths (user → system)
	seen := make(map[string]bool)
	var scripts []string

	for _, dir := range config.GetHelperScriptSearchPaths() {
		entries, err := os.ReadDir(dir)
		if err != nil {
			continue
		}

		for _, entry := range entries {
			name := entry.Name()
			// Skip hidden files and directories
			if !entry.IsDir() && !strings.HasPrefix(name, ".") {
				if !seen[name] {
					seen[name] = true
					scripts = append(scripts, name)
				}
			}
		}
	}

	return scripts, cobra.ShellCompDirectiveNoFileComp
}

func parsePostScriptHelperFlags(args []string) (postScriptHelperFlags, []string, error) {
	var flags postScriptHelperFlags
	var passthrough []string

	valueFor := func(name string, i *int, inline string, hasInline bool) (string, error) {
		if hasInline {
			return inline, nil
		}
		*i = *i + 1
		if *i >= len(args) {
			return "", fmt.Errorf("flag %s requires a value", name)
		}
		return args[*i], nil
	}

	for i := 0; i < len(args); i++ {
		arg := args[i]
		if arg == "--" {
			passthrough = append(passthrough, args[i:]...)
			break
		}

		if strings.HasPrefix(arg, "--") {
			name, val, hasVal := strings.Cut(arg, "=")
			switch name {
			case "--cpus":
				v, err := valueFor(name, &i, val, hasVal)
				if err != nil {
					return flags, nil, err
				}
				n, err := parsePositiveIntFlag(name, v)
				if err != nil {
					return flags, nil, err
				}
				flags.cpusSet = true
				flags.cpus = n
			case "--mem":
				v, err := valueFor(name, &i, val, hasVal)
				if err != nil {
					return flags, nil, err
				}
				flags.memSet = true
				flags.mem = v
			case "--time":
				v, err := valueFor(name, &i, val, hasVal)
				if err != nil {
					return flags, nil, err
				}
				flags.timeSet = true
				flags.time = v
			case "--gpu":
				v, err := valueFor(name, &i, val, hasVal)
				if err != nil {
					return flags, nil, err
				}
				flags.gpuSet = true
				flags.gpu = v
			case "--base":
				v, err := valueFor(name, &i, val, hasVal)
				if err != nil {
					return flags, nil, err
				}
				flags.baseSet = true
				flags.base = v
			case "--env":
				v, err := valueFor(name, &i, val, hasVal)
				if err != nil {
					return flags, nil, err
				}
				flags.envSet = true
				flags.env = v
			case "--overlay":
				v, err := valueFor(name, &i, val, hasVal)
				if err != nil {
					return flags, nil, err
				}
				flags.overlaySet = true
				flags.overlays = append(flags.overlays, v)
			case "--cwd":
				v, err := valueFor(name, &i, val, hasVal)
				if err != nil {
					return flags, nil, err
				}
				flags.cwdSet = true
				flags.cwd = v
			case "--new":
				if hasVal {
					return flags, nil, fmt.Errorf("flag %s does not take a value", name)
				}
				flags.newSet = true
			default:
				passthrough = append(passthrough, arg)
			}
			continue
		}

		if strings.HasPrefix(arg, "-") && len(arg) > 1 {
			name := arg[:2]
			attached := arg[2:]
			hasAttached := attached != ""
			switch name {
			case "-c":
				v, err := valueFor(name, &i, attached, hasAttached)
				if err != nil {
					return flags, nil, err
				}
				n, err := parsePositiveIntFlag(name, v)
				if err != nil {
					return flags, nil, err
				}
				flags.cpusSet = true
				flags.cpus = n
			case "-m":
				v, err := valueFor(name, &i, attached, hasAttached)
				if err != nil {
					return flags, nil, err
				}
				flags.memSet = true
				flags.mem = v
			case "-t":
				v, err := valueFor(name, &i, attached, hasAttached)
				if err != nil {
					return flags, nil, err
				}
				flags.timeSet = true
				flags.time = v
			case "-g":
				v, err := valueFor(name, &i, attached, hasAttached)
				if err != nil {
					return flags, nil, err
				}
				flags.gpuSet = true
				flags.gpu = v
			case "-b":
				v, err := valueFor(name, &i, attached, hasAttached)
				if err != nil {
					return flags, nil, err
				}
				flags.baseSet = true
				flags.base = v
			case "-e":
				v, err := valueFor(name, &i, attached, hasAttached)
				if err != nil {
					return flags, nil, err
				}
				flags.envSet = true
				flags.env = v
			case "-o":
				v, err := valueFor(name, &i, attached, hasAttached)
				if err != nil {
					return flags, nil, err
				}
				flags.overlaySet = true
				flags.overlays = append(flags.overlays, v)
			case "-w":
				v, err := valueFor(name, &i, attached, hasAttached)
				if err != nil {
					return flags, nil, err
				}
				flags.cwdSet = true
				flags.cwd = v
			default:
				passthrough = append(passthrough, arg)
			}
			continue
		}

		passthrough = append(passthrough, arg)
	}

	return flags, passthrough, nil
}

func parsePositiveIntFlag(name, value string) (int, error) {
	n, err := strconv.Atoi(value)
	if err != nil || n <= 0 {
		return 0, fmt.Errorf("invalid %s value %q", name, value)
	}
	return n, nil
}

func runHelper(cmd *cobra.Command, args []string) error {
	// Prevent helper commands when inside a container or scheduler job
	// (except --path / --list / --status / config — config writes to NFS and is safe from inside).
	helperConfigMode := len(args) > 1 && args[1] == "config"
	if config.IsInsideContainer() && !helperPath && !helperList && !helperStatus && !helperConfigMode {
		cmd.SilenceUsage = true
		ExitWithError("helper commands are not available inside a container")
	}
	if scheduler.IsInsideJob() && !helperPath && !helperList && !helperStatus && !helperConfigMode {
		cmd.SilenceUsage = true
		ExitWithError("helper commands are not available inside a scheduler job")
	}

	// --- Path Mode ---
	if helperPath {
		if len(args) > 0 {
			// Find specific script in all search paths
			scriptName := args[0]
			scriptPath, err := config.FindHelperScript(scriptName)
			if err != nil {
				cmd.SilenceUsage = true
				ExitWithError("helper script '%s' not found (searched: %s)",
					scriptName, strings.Join(config.GetHelperScriptSearchPaths(), ", "))
			}
			fmt.Println(scriptPath)
			return nil
		}

		// Show all helper script search paths
		fmt.Println("Helper script search paths (priority order):")
		for i, dir := range config.GetHelperScriptSearchPaths() {
			exists := ""
			if !config.DirExists(dir) {
				exists = " (not found)"
			}
			fmt.Printf("  %d. %s%s\n", i+1, dir, exists)
		}

		// Show writable directory
		if writableDir, err := config.GetWritableHelperScriptsDir(); err == nil {
			fmt.Printf("\nWritable directory: %s\n", writableDir)
		}
		return nil
	}

	// --- List Mode ---
	if helperList {
		type scriptInfo struct {
			name   string
			whatis string
		}
		seen := make(map[string]bool)
		var scripts []scriptInfo
		maxNameLen := 0

		for _, dir := range config.GetHelperScriptSearchPaths() {
			entries, err := os.ReadDir(dir)
			if err != nil {
				continue
			}
			for _, entry := range entries {
				name := entry.Name()
				if !entry.IsDir() && !strings.HasPrefix(name, ".") && !seen[name] {
					seen[name] = true
					whatis := utils.GetWhatIsFromScript(filepath.Join(dir, name))
					scripts = append(scripts, scriptInfo{name, whatis})
					if len(name) > maxNameLen {
						maxNameLen = len(name)
					}
				}
			}
		}

		for _, s := range scripts {
			if s.whatis != "" {
				fmt.Printf("  %-*s  %s\n", maxNameLen, s.name, s.whatis)
			} else {
				fmt.Printf("  %s\n", s.name)
			}
		}
		return nil
	}

	// --- Update Mode ---
	if helperUpdate {
		name := ""
		if len(args) > 0 {
			name = args[0]
		}
		if err := helper.UpdateRemoteScripts(cmd.Context(), name, false, os.Stdout); err != nil {
			cmd.SilenceUsage = true
			ExitWithError("%v", err)
		}
		return nil
	}

	// --- Status Mode ---
	if helperStatus {
		name := ""
		if len(args) > 0 {
			name = args[0]
		}
		return showHelperStatus(name)
	}

	// --- Run Mode ---
	if len(args) == 0 {
		helper.EnsureServer()
		return showHelperStatus("")
	}

	// Intercept: condatainer helper stop [--all] — stop across all helper names
	if args[0] == "stop" {
		stopAll := false
		for _, a := range args[1:] {
			if a == "--all" || a == "-a" {
				stopAll = true
			}
		}
		return runHelperStop(cmd.Context(), "", stopAll)
	}

	scriptName := args[0]
	scriptArgs := args[1:]

	// Intercept: condatainer helper <name> config [path|show|get <key>|set <key> <val>]
	if len(scriptArgs) > 0 && scriptArgs[0] == "config" {
		scriptPath, _ := config.FindHelperScript(scriptName)
		return runHelperConfig(scriptName, scriptPath, scriptArgs[1:])
	}

	// Intercept: condatainer helper <name> stop [--all]
	if len(scriptArgs) > 0 && scriptArgs[0] == "stop" {
		stopAll := false
		for _, a := range scriptArgs[1:] {
			if a == "--all" || a == "-a" {
				stopAll = true
			}
		}
		return runHelperStop(cmd.Context(), scriptName, stopAll)
	}

	// Intercept -h/--help in the remaining args to show helper-specific usage.
	for _, a := range scriptArgs {
		if a == "-h" || a == "--help" {
			scriptPath, findErr := config.FindHelperScript(scriptName)
			if findErr != nil {
				return fmt.Errorf("helper script '%s' not found", scriptName)
			}
			return showHelperHelp(scriptName, scriptPath)
		}
	}

	postFlags, scriptArgs, err := parsePostScriptHelperFlags(scriptArgs)
	if err != nil {
		return err
	}

	// startedWithArgs is true when the user supplied explicit resource/path/param flags,
	// meaning we skip the history-reuse menu and go straight to settings.
	startedWithArgs := postFlags.cpusSet || postFlags.memSet ||
		postFlags.timeSet || postFlags.gpuSet || postFlags.baseSet ||
		postFlags.envSet || postFlags.overlaySet || postFlags.cwdSet ||
		postFlags.newSet || len(scriptArgs) > 0

	// Find script in all search paths.
	scriptPath, err := config.FindHelperScript(scriptName)
	if err != nil {
		cmd.SilenceUsage = true
		ExitWithError("helper script '%s' not found\nSearched: %s\nRun '%s' to fetch available helper scripts",
			scriptName, strings.Join(config.GetHelperScriptSearchPaths(), ", "),
			utils.StyleAction("condatainer helper --update"))
	}

	// Ensure executable.
	if err := os.Chmod(scriptPath, utils.PermExec); err != nil {
		utils.PrintDebug("Failed to chmod helper script: %v", err)
	}

	// Check apptainer is available.
	if err := apptainer.EnsureApptainer(); err != nil {
		cmd.SilenceUsage = true
		ExitWithError("apptainer is required to run helper scripts: %v", err)
	}

	// Build resource overrides from explicit flags.
	var overrides *scheduler.ResourceSpec
	if postFlags.cpusSet || postFlags.memSet || postFlags.timeSet || postFlags.gpuSet {
		overrides = &scheduler.ResourceSpec{}
		if postFlags.cpusSet {
			overrides.CpusPerTask = postFlags.cpus
		}
		if postFlags.memSet {
			mb, parseErr := utils.ParseMemoryMB(postFlags.mem)
			if parseErr != nil {
				return fmt.Errorf("invalid --mem value %q: %w", postFlags.mem, parseErr)
			}
			overrides.MemPerNodeMB = mb
		}
		if postFlags.timeSet {
			d, parseErr := utils.ParseWalltime(postFlags.time)
			if parseErr != nil {
				return fmt.Errorf("invalid --time value %q: %w", postFlags.time, parseErr)
			}
			overrides.Time = d
		}
		if postFlags.gpuSet {
			overrides.Gpu = helper.ParseGPUSpec(postFlags.gpu)
		}
	}

	cwd := ""
	if postFlags.cwdSet {
		abs, err := filepath.Abs(postFlags.cwd)
		if err != nil {
			return fmt.Errorf("invalid --cwd path %q: %w", postFlags.cwd, err)
		}
		cwd = abs
	}

	meta, _ := helper.ParseHelperScriptMeta(scriptPath)
	versions := meta.ParamValues

	// -e - means "no env overlay; skip auto-resolve".
	noEnv := postFlags.envSet && postFlags.env == "-"
	envImg := postFlags.env
	if noEnv {
		envImg = ""
	}

	opts := helper.RunOptions{
		ScriptPath: scriptPath,
		ScriptName: scriptName,
		Resources:  overrides,
		BaseImage:  postFlags.base,
		EnvImg:     envImg,
		Overlays:   postFlags.overlays,
		CWD:        cwd,
		ForceNew:   postFlags.newSet,
		FlagArgs:   scriptArgs,
	}

	ctx := cmd.Context()
	cmd.SilenceUsage = true

	// Check for running instances; enforce singleton.
	if !opts.ForceNew {
		running, _ := helper.RunningHelpers(opts.ScriptName)
		if helper.SingletonBlocked(meta, running) {
			return ui.OfferSingletonBlocked(opts.ScriptName, running, config.GetRunningServerPort())
		}
		if !startedWithArgs {
			// Combined running + history session list.
			serverPort := config.GetRunningServerPort()
			chosen, isRunning, startNew, sesErr := ui.OfferRecentSessions(ctx, opts.ScriptName, running, serverPort)
			if sesErr != nil {
				if errors.Is(sesErr, context.Canceled) {
					return nil
				}
				return sesErr
			}
			if isRunning && chosen != nil {
				if chosen.Status == "pending" || chosen.Status == "starting" {
					// Not up yet — re-attach and wait for the service URL.
					return ui.MonitorHelper(ctx, chosen.ID)
				}
				if url := chosen.AccessURL(serverPort); url != "" {
					utils.PrintSuccess("Access at: %s", url)
				}
				return nil
			}
			if !startNew && chosen != nil {
				ui.ApplyHistoryRun(&opts, chosen)
			}
		} else if len(running) > 0 {
			// User supplied explicit flags — just show a brief notice about running jobs.
			serverPort := config.GetRunningServerPort()
			for _, r := range running {
				if url := r.AccessURL(serverPort); url != "" {
					fmt.Printf("%s already running: %s\n", utils.StyleName(opts.ScriptName), url)
				} else if r.Status == "pending" || r.Status == "starting" {
					fmt.Printf("%s already submitted (%s); check with 'condatainer helper --status'\n",
						utils.StyleName(opts.ScriptName), r.Status)
				}
			}
		}
	}

	// Ensure server is running.
	helper.EnsureServer()

	// When -e is unset, search cwd → overlay/ → src/overlay/ (preferring
	// env-$USER.img) and use the first overlay found. An explicit -e value is
	// used literally; -e - opts out.
	if !noEnv && !postFlags.envSet && (opts.EnvImg == "" || opts.EnvImg == "env.img") {
		resolvedCwd := cwd
		if resolvedCwd == "" {
			resolvedCwd, _ = os.Getwd()
		}
		opts.EnvImg = helper.ResolveEnvOverlayInDir(opts.EnvImg, resolvedCwd)
	}

	// Resolve params and resources together in one combined prompt.
	scriptParams, _ := helper.ParseHelperParams(opts.ScriptPath)
	if err := helper.ValidateHelperParamConflicts(scriptParams); err != nil {
		return err
	}
	flagValues, _, err := helper.ApplyParamFlags(scriptParams, opts.FlagArgs)
	if err != nil {
		return fmt.Errorf("parsing param flags: %w", err)
	}
	saved, _ := helper.LoadHelperConfig(opts.ScriptName)
	spec := helper.ResolveHelperSpec(opts.ScriptPath)
	filled, err := ui.PromptSettings(ctx, scriptParams, flagValues, saved, versions, &opts, spec)
	if err != nil {
		if errors.Is(err, context.Canceled) {
			return nil
		}
		return err
	}
	opts.Params = filled

	// Guided overlay creation for helpers that require a writable conda env.
	if meta.ImgPackages != "" && opts.EnvImg == "" {
		created, err := ui.GuidedOverlayCreate(ctx, opts.ScriptName, meta, opts.Params, versions, opts.CWD)
		if err != nil {
			if errors.Is(err, context.Canceled) {
				return nil
			}
			return err
		}
		opts.EnvImg = created
	}

	// Plan: non-interactive validation + resource resolution.
	plan, err := helper.PlanRun(ctx, opts)
	// Offer to install missing conda packages, then retry planning once.
	if err != nil {
		var em *helper.ErrMissingPackages
		if errors.As(err, &em) {
			if installErr := ui.GuidedInstallMissing(ctx, em); installErr != nil {
				if errors.Is(installErr, context.Canceled) {
					return nil
				}
				return installErr
			}
			plan, err = helper.PlanRun(ctx, opts)
		}
	}
	if err != nil {
		var ep *helper.ErrMissingParam
		var es *helper.ErrSingletonRunning
		var ee *helper.ErrEnvInUse
		switch {
		case errors.As(err, &ep):
			return fmt.Errorf("missing required params: %s", ep.Error())
		case errors.As(err, &es):
			return ui.OfferSingletonBlocked(opts.ScriptName, es.Existing, config.GetRunningServerPort())
		case errors.As(err, &ee):
			return fmt.Errorf("env overlay in use: %s", ee.Path)
		}
		if errors.Is(err, context.Canceled) {
			return nil
		}
		return err
	}

	// Print launch spec and wait before submit (when the user supplied explicit flags).
	if startedWithArgs {
		ui.PrintLaunchSpec(plan)
		if err := ui.WaitBeforeStart(ctx, 3*time.Second); err != nil {
			if errors.Is(err, context.Canceled) {
				return nil
			}
			return err
		}
	}

	// Execute: submit to scheduler or run headless.
	helperID, err := helper.ExecutePlan(ctx, plan)
	if err != nil {
		if errors.Is(err, context.Canceled) {
			return nil
		}
		return err
	}

	// Monitor: poll NFS state files until ready or done.
	if err := ui.MonitorHelper(ctx, helperID); err != nil {
		if errors.Is(err, context.Canceled) {
			return nil
		}
		return err
	}
	return nil
}

// showHelperStatus lists running helpers (optionally filtered by name) from JSONL history.
func showHelperStatus(name string) error {
	running, err := helper.RunningHelpers(name)
	if err != nil {
		return err
	}
	if len(running) == 0 {
		if name != "" {
			fmt.Printf("No %s helpers running. Use 'condatainer helper %s' to start.\n", name, name)
		} else {
			fmt.Println("No helpers running. Use 'condatainer helper <name>' to start.")
		}
		return nil
	}
	serverPort := config.GetRunningServerPort()
	for _, r := range running {
		url := r.AccessURL(serverPort)
		age := ui.FormatAge(time.Since(r.StartedAt))
		res := ui.FormatRunResources(r)
		cwd := r.CWD
		if cwd == "" {
			cwd = "(no cwd)"
		}
		tag, verb := ui.RunStatusTag(r)
		fmt.Printf("  %-14s  %-14s  %-40s  %s %s %s ago\n", r.Name, res, cwd, tag, verb, age)
		if url != "" {
			fmt.Printf("    %s\n", utils.StyleDebug("→ "+url))
		}
	}
	return nil
}

// runHelperConfig handles: condatainer helper <name> config [show|get <key>|set <key> <val>|path|-h]
func runHelperConfig(name, scriptPath string, args []string) error {
	sub := ""
	if len(args) > 0 {
		sub = args[0]
	}

	// Load params and spec once (needed for show and -h).
	var params []helper.HelperParam
	var spec *scheduler.ResourceSpec
	if scriptPath != "" {
		params, _ = helper.ParseHelperParams(scriptPath)
		spec = helper.ResolveHelperSpec(scriptPath)
	}

	switch sub {
	case "-h", "--help", "":
		keys := []string{"cpus", "mem", "time", "gpu"}
		for _, p := range params {
			keys = append(keys, p.Key)
		}
		// Laid out to match cobra's help format; this cannot be a real cobra
		// command because the helper name is only known at runtime.
		fmt.Printf("Show or save default settings for the %s helper.\n\n", name)
		fmt.Printf("Usage:\n  condatainer helper %s config <command>\n\n", name)
		fmt.Println("Available Commands:")
		fmt.Println("  get         Print a single saved value")
		fmt.Println("  path        Print the config file path")
		fmt.Println("  set         Save a default value")
		fmt.Println("  show        Show saved defaults")
		fmt.Printf("\nExamples:\n")
		fmt.Printf("  condatainer helper %s config show\n", name)
		fmt.Printf("  condatainer helper %s config set cpus 4\n", name)
		fmt.Printf("  condatainer helper %s config get cpus\n", name)
		fmt.Printf("\nSaveable keys: %s\n", strings.Join(keys, ", "))

	case "path":
		p := helper.HelperConfigPath(name)
		if p == "" {
			return fmt.Errorf("cannot determine config path")
		}
		fmt.Println(p)

	case "get":
		if len(args) < 2 {
			return fmt.Errorf("usage: condatainer helper %s config get <key>", name)
		}
		v, ok := helper.GetHelperConfigKey(name, args[1])
		if !ok {
			return fmt.Errorf("key %q not found in %s config", args[1], name)
		}
		fmt.Println(v)

	case "set":
		if len(args) < 3 {
			return fmt.Errorf("usage: condatainer helper %s config set <key> <value>", name)
		}
		if err := helper.SetHelperConfigKey(name, args[1], args[2]); err != nil {
			return fmt.Errorf("saving config: %w", err)
		}

	case "show":
		cfg, err := helper.LoadHelperConfig(name)
		if err != nil {
			return fmt.Errorf("loading config: %w", err)
		}
		p := helper.HelperConfigPath(name)
		fmt.Printf("# %s\n", p)

		type configRow struct {
			key, def string
		}
		var rows []configRow
		// Resource rows with spec defaults.
		resDefault := func(v string) string {
			if v == "" {
				return "(none)"
			}
			return v
		}
		cpuDef, memDef, timeDef := "", "", ""
		if spec != nil {
			if spec.CpusPerTask > 0 {
				cpuDef = fmt.Sprintf("%d", spec.CpusPerTask)
			}
			if mb := spec.GetMemPerNodeMB(); mb > 0 {
				memDef = utils.FormatMemoryMB(mb)
			}
			if spec.Time > 0 {
				timeDef = utils.FormatDuration(spec.Time)
			}
		}
		rows = append(rows,
			configRow{"cpus", resDefault(cpuDef)},
			configRow{"mem", resDefault(memDef)},
			configRow{"time", resDefault(timeDef)},
			configRow{"gpu", "(none)"},
		)
		for _, pp := range params {
			rows = append(rows, configRow{strings.ToLower(pp.Key), pp.Default})
		}

		maxKey := 4 // len("cpus")
		for _, r := range rows {
			if len(r.key) > maxKey {
				maxKey = len(r.key)
			}
		}
		for _, r := range rows {
			val, ok := cfg[strings.ToLower(r.key)]
			valDisp := "(not set)"
			if ok && val != "" {
				valDisp = val
			}
			defDisp := ""
			if r.def != "" {
				defDisp = utils.StyleDebug(fmt.Sprintf("(default: %s)", r.def))
			}
			fmt.Printf("  %-*s  %-20s  %s\n", maxKey, r.key, valDisp, defDisp)
		}

	default:
		return fmt.Errorf("unknown config command %q — try: condatainer helper %s config -h", sub, name)
	}
	return nil
}

func showHelperHelp(scriptName, scriptPath string) error {
	params, err := helper.ParseHelperParams(scriptPath)
	if err != nil {
		return err
	}
	if err := helper.ValidateHelperParamConflicts(params); err != nil {
		return err
	}
	meta, _ := helper.ParseHelperScriptMeta(scriptPath)
	helper.PrintHelperUsage(scriptName, params, helper.ResolveHelperSpec(scriptPath), meta.ParamValues)
	return nil
}

// runHelperStop stops running instances of a helper.
// name filters to a specific helper script; empty means all helpers.
// If stopAll is true, all matching instances are stopped without prompting.
// When name is set and exactly one instance is running, it is stopped directly.
// When name is empty, the picker is always shown (acts as confirmation).
// When multiple instances are running and stopAll is false, the user is prompted
// to choose one or all via OfferStopPicker.
func runHelperStop(ctx context.Context, name string, stopAll bool) error {
	running, err := helper.RunningHelpers(name)
	if err != nil {
		return err
	}
	if len(running) == 0 {
		if name == "" {
			utils.PrintMessage("No running helpers.")
		} else {
			utils.PrintMessage("No running %s helpers.", name)
		}
		return nil
	}

	var toStop []*helper.HelperRun
	// Skip picker only when a specific name is given and there is exactly one match.
	if stopAll || (name != "" && len(running) == 1) {
		toStop = running
	} else {
		chosen, pickErr := ui.OfferStopPicker(ctx, name, running)
		if pickErr != nil {
			return pickErr
		}
		if len(chosen) == 0 {
			return nil // user cancelled
		}
		toStop = chosen
	}

	for _, r := range toStop {
		if err := stopHelperRun(r); err != nil {
			utils.PrintWarning("Failed to stop %s: %v", r.ID, err)
		} else {
			utils.PrintSuccess("Stopped %s", r.ID)
		}
	}
	return nil
}

// stopHelperRun terminates a single helper run:
//   - Scheduled (JobID != ""): asks the scheduler to cancel
//   - Headless (JobID == ""): sends SIGTERM to the process group via the pid file
//
// Always marks the history entry as "done" so the server watcher and status
// commands see it as finished even if the signal fails (process already dead).
func stopHelperRun(r *helper.HelperRun) error {
	if r.JobID != "" {
		if sched := scheduler.ActiveScheduler(); sched != nil {
			_ = sched.CancelJob(context.Background(), r.JobID)
		}
	} else {
		if err := helper.KillHeadlessProcess(r.ID); err != nil {
			utils.PrintDebug("kill headless process %s: %v", r.ID, err)
		}
	}
	return helper.UpdateHistoryStatus(r.ID, "done", time.Now())
}
