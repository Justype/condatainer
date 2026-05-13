package cmd

import (
	"compress/gzip"
	"context"
	"crypto/sha256"
	"encoding/json"
	"errors"
	"fmt"
	"io"
	"net/http"
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

// ForceRefreshHelpers skips the helper metadata disk cache when set.
var ForceRefreshHelpers bool

// RemoteHelperMetadataCache is the on-disk envelope for cached helper script metadata.
type RemoteHelperMetadataCache struct {
	FetchedAt time.Time                            `json:"fetched_at"`
	SourceURL string                               `json:"source_url"`
	Metadata  map[string]map[string]rawHelperEntry `json:"metadata"`
}

// rawHelperEntry is the wire format for a helper script entry (no SourceURL).
type rawHelperEntry struct {
	Path string `json:"path"`
}

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
	newSet     bool
}

// in-process cache for helper metadata, keyed by base URL
var cachedHelperMetadataByURL = map[string]map[string]map[string]HelperScriptEntry{}

var (
	helperPath   bool
	helperList   bool
	helperUpdate bool
	helperStatus bool
)

// supportedSchedulerTypes lists scheduler types supported for helper scripts.
// Add new supported scheduler types to this slice to enable additional helper categories.
var supportedSchedulerTypes = []scheduler.SchedulerType{
	scheduler.SchedulerSLURM,
}

var helperCmd = &cobra.Command{
	Use:   "helper [flags] [script-name] [script-args...]",
	Short: "Manage and run helper scripts",
	Long: `Manage helper scripts for running services inside CondaTainer on HPC.

Use --list to see available helper scripts.
Use --update to download or refresh helper scripts from remote.

Note: Helper is not available inside a container or a scheduler job.`,
	Example: `  condatainer helper                 # Show running helpers
  condatainer helper code-server     # Run code-server
  condatainer helper code-server -h  # Show code-server flags and options
  condatainer helper --update        # Update all helper scripts
  condatainer helper --list          # List available helper scripts`,
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
				if hasVal {
					return flags, nil, fmt.Errorf("flag %s does not take a value", name)
				}
				flags.cwdSet = true
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
				if hasAttached {
					passthrough = append(passthrough, arg)
					continue
				}
				flags.cwdSet = true
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
	// Prevent helper commands when inside a container or scheduler job (except --path / --list / --status)
	if config.IsInsideContainer() && !helperPath && !helperList && !helperStatus {
		cmd.SilenceUsage = true
		ExitWithError("helper commands are not available inside a container")
	}
	if scheduler.IsInsideJob() && !helperPath && !helperList && !helperStatus {
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
				ExitWithError("helper script '%s' not found (searched: %v)", scriptName, config.GetHelperScriptSearchPaths())
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
		// Get writable directory for updates
		helperScriptsDir, err := config.GetWritableHelperScriptsDir()
		if err != nil {
			cmd.SilenceUsage = true
			ExitWithError("failed to find writable helper scripts directory: %v", err)
		}
		if err := updateHelperScripts(args, helperScriptsDir); err != nil {
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

	scriptName := args[0]
	scriptArgs := args[1:]

	// Intercept: condatainer helper <name> config [path|show|get <key>|set <key> <val>]
	if len(scriptArgs) > 0 && scriptArgs[0] == "config" {
		scriptPath, _ := config.FindHelperScript(scriptName)
		return runHelperConfig(scriptName, scriptPath, scriptArgs[1:])
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
		ExitWithError("helper script '%s' not found\nSearched: %v\nRun '%s' to fetch available helper scripts",
			scriptName, config.GetHelperScriptSearchPaths(), utils.StyleAction("condatainer helper --update"))
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
		cwd, _ = os.Getwd()
	}

	meta, _ := helper.ParseHelperScriptMeta(scriptPath)
	versions := meta.ParamValues

	opts := helper.RunOptions{
		ScriptPath: scriptPath,
		ScriptName: scriptName,
		Resources:  overrides,
		BaseImage:  postFlags.base,
		EnvImg:     postFlags.env,
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
				}
			}
		}
	}

	// Ensure server is running.
	helper.EnsureServer()

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
		created, err := ui.GuidedOverlayCreate(ctx, opts.ScriptName, meta, opts.Params, versions)
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
	helperID, err := helper.ExecutePlan(ctx, plan, helper.IO{Stdout: os.Stdout, Stderr: os.Stderr})
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
		age := formatAge(time.Since(r.StartedAt))
		res := ui.FormatRunResources(r)
		cwd := r.CWD
		if cwd == "" {
			cwd = "(no cwd)"
		}
		fmt.Printf("  %-14s  %-14s  %-40s  %s ago\n", r.Name, res, cwd, age)
		if url != "" {
			fmt.Printf("    %s\n", utils.StyleDebug("→ "+url))
		}
	}
	return nil
}

// formatAge formats an elapsed duration as a compact human-readable string (e.g. "6h11m", "2d3h", "45m").
func formatAge(d time.Duration) string {
	d = d.Round(time.Minute)
	if d < time.Minute {
		return "just now"
	}
	days := int(d.Hours()) / 24
	hours := int(d.Hours()) % 24
	mins := int(d.Minutes()) % 60
	switch {
	case days > 0 && hours > 0:
		return fmt.Sprintf("%dd%dh", days, hours)
	case days > 0:
		return fmt.Sprintf("%dd", days)
	case hours > 0 && mins > 0:
		return fmt.Sprintf("%dh%dm", hours, mins)
	case hours > 0:
		return fmt.Sprintf("%dh", hours)
	default:
		return fmt.Sprintf("%dm", mins)
	}
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
		fmt.Printf("Usage: condatainer helper %s config <command>\n\n", name)
		fmt.Println("  show              Show saved defaults")
		fmt.Println("  get <KEY>         Print a single saved value")
		fmt.Println("  set <KEY> <VALUE> Save a default value")
		fmt.Println("  path              Print config file path")
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

// updateHelperScripts downloads and updates helper scripts from remote metadata
func updateHelperScripts(args []string, helperScriptsDir string) error {
	// Fetch remote metadata
	metadata, err := GetAllRemoteHelperMetadata()
	if err != nil {
		return fmt.Errorf("failed to fetch remote helper metadata: %w", err)
	}

	// Choose category based on scheduler availability
	category := "headless"
	if config.Global.SubmitJob {
		// If user requested submission, ensure a supported scheduler is available
		sched, err := scheduler.DetectSchedulerWithBinary(config.Global.SchedulerBin)
		if err != nil {
			return fmt.Errorf("submit requested but no scheduler found: %v", err)
		}

		if !sched.IsAvailable() || sched.IsInsideJob() {
			return fmt.Errorf("scheduler is not available for submission (available=%v, in_job=%v)", sched.IsAvailable(), sched.IsInsideJob())
		}

		matched := false
		for _, st := range supportedSchedulerTypes {
			if sched.GetType() == st {
				category = strings.ToLower(string(st))
				matched = true
				break
			}
		}

		if !matched {
			return fmt.Errorf("unsupported scheduler type '%s'", sched.GetType())
		}
	}

	category_styled := utils.StyleName(category)
	utils.PrintDebug("Helper category selected: %s", category_styled)

	entries, ok := metadata[category]
	if !ok {
		return fmt.Errorf("no helper scripts found for category '%s'", category_styled)
	}

	// Create helper scripts directory
	if err := os.MkdirAll(helperScriptsDir, utils.PermDir); err != nil {
		return fmt.Errorf("failed to create helper scripts directory: %w", err)
	}

	// Filter to specific script if provided
	if len(args) > 0 {
		scriptName := args[0]
		if entry, ok := entries[scriptName]; ok {
			entries = map[string]HelperScriptEntry{scriptName: entry}
		} else {
			return fmt.Errorf("helper script '%s' not found in remote metadata for category '%s'", scriptName, category_styled)
		}
	} else {
		utils.PrintMessage("Updating all helper scripts for %s...", category_styled)
	}

	// Download each helper script
	for name, entry := range entries {
		if entry.Path == "" {
			continue
		}

		sourceURL := entry.SourceURL
		if sourceURL == "" {
			sourceURL = config.Global.ScriptsLink
		}
		url := fmt.Sprintf("%s/%s", sourceURL, entry.Path)
		destName := filepath.Base(entry.Path)
		dest := filepath.Join(helperScriptsDir, destName)

		utils.PrintMessage("Updating %s/%s", category, name)
		if err := downloadExecutable(url, dest); err != nil {
			utils.PrintWarning("Failed to update %s/%s: %v", category, name, err)
			continue
		}
	}

	utils.PrintSuccess("Helper update finished.")
	return nil
}

// HelperScriptEntry represents a helper script metadata entry.
type HelperScriptEntry struct {
	Path      string `json:"path"`
	SourceURL string `json:"-"` // injected during merge; not stored in remote JSON
}

// helperMetadataCachePathForURL returns the disk cache path for helper metadata from a given base URL.
func helperMetadataCachePathForURL(baseURL string) (string, error) {
	cacheDir, err := config.GetWritableCacheDir()
	if err != nil {
		return "", err
	}
	h := sha256.Sum256([]byte(baseURL))
	name := fmt.Sprintf("helper-scripts-%x.json.gz", h[:6])
	return filepath.Join(cacheDir, name), nil
}

func loadHelperMetadataCacheAny(path, baseURL string, ttl time.Duration, checkTTL bool) (*RemoteHelperMetadataCache, error) {
	var cache RemoteHelperMetadataCache
	if err := utils.ReadGzipJSONFile(path, &cache); err != nil {
		return nil, err
	}
	if cache.SourceURL != baseURL {
		return nil, fmt.Errorf("source URL changed")
	}
	if checkTTL && time.Since(cache.FetchedAt) > ttl {
		return nil, fmt.Errorf("cache expired")
	}
	return &cache, nil
}

// fetchHelperMetadataFromURL fetches (or loads from per-URL cache) helper metadata for one base URL.
// Injects SourceURL into each entry during conversion.
func fetchHelperMetadataFromURL(baseURL string) (map[string]map[string]HelperScriptEntry, error) {
	// In-process cache
	if cached, ok := cachedHelperMetadataByURL[baseURL]; ok {
		return cached, nil
	}

	metaURL := baseURL + "/metadata/helper-scripts.json.gz"
	cachePath, cacheErr := helperMetadataCachePathForURL(baseURL)

	ttl := config.Global.MetadataCacheTTL

	// Try disk cache
	if !ForceRefreshHelpers && cacheErr == nil && ttl > 0 {
		if cache, err := loadHelperMetadataCacheAny(cachePath, baseURL, ttl, true); err == nil {
			utils.PrintDebug("Using cached helper metadata for %s (fetched %s ago)", baseURL, time.Since(cache.FetchedAt).Round(time.Minute))
			result := injectHelperSourceURL(cache.Metadata, baseURL)
			cachedHelperMetadataByURL[baseURL] = result
			return result, nil
		}
	}

	// Network fetch
	utils.PrintDebug("Fetching helper metadata from %s...", metaURL)
	rawMeta, err := fetchRawHelperMetadata(metaURL)
	if err != nil {
		// Stale fallback
		if cacheErr == nil {
			if cache, staleErr := loadHelperMetadataCacheAny(cachePath, baseURL, 0, false); staleErr == nil {
				utils.PrintWarning("Network unavailable for %s; using cached helper metadata (fetched %s ago)",
					metaURL, time.Since(cache.FetchedAt).Round(time.Minute))
				result := injectHelperSourceURL(cache.Metadata, baseURL)
				cachedHelperMetadataByURL[baseURL] = result
				return result, nil
			}
		}
		return nil, err
	}

	// Persist to disk (non-fatal)
	if cacheErr == nil && ttl > 0 {
		envelope := RemoteHelperMetadataCache{
			FetchedAt: time.Now(),
			SourceURL: baseURL,
			Metadata:  rawMeta,
		}
		if err := utils.WriteGzipJSONFileAtomic(cachePath, envelope); err != nil {
			utils.PrintDebug("Failed to write helper metadata cache for %s: %v", baseURL, err)
		}
	}

	result := injectHelperSourceURL(rawMeta, baseURL)
	cachedHelperMetadataByURL[baseURL] = result
	return result, nil
}

// injectHelperSourceURL converts raw entries to HelperScriptEntry, injecting SourceURL.
func injectHelperSourceURL(raw map[string]map[string]rawHelperEntry, baseURL string) map[string]map[string]HelperScriptEntry {
	out := make(map[string]map[string]HelperScriptEntry, len(raw))
	for category, scripts := range raw {
		out[category] = make(map[string]HelperScriptEntry, len(scripts))
		for name, e := range scripts {
			out[category][name] = HelperScriptEntry{Path: e.Path, SourceURL: baseURL}
		}
	}
	return out
}

// fetchRawHelperMetadata performs the HTTP+gzip+JSON fetch for helper metadata.
func fetchRawHelperMetadata(url string) (map[string]map[string]rawHelperEntry, error) {
	resp, err := http.Get(url)
	if err != nil {
		return nil, err
	}
	defer resp.Body.Close()

	if resp.StatusCode != http.StatusOK {
		return nil, fmt.Errorf("failed to fetch metadata: HTTP %d", resp.StatusCode)
	}

	gzReader, err := gzip.NewReader(resp.Body)
	if err != nil {
		return nil, fmt.Errorf("failed to decompress metadata: %w", err)
	}
	defer gzReader.Close()

	data, err := io.ReadAll(gzReader)
	if err != nil {
		return nil, err
	}

	var metadata map[string]map[string]rawHelperEntry
	if err := json.Unmarshal(data, &metadata); err != nil {
		return nil, err
	}
	return metadata, nil
}

// GetAllRemoteHelperMetadata merges helper metadata from all configured remote URLs.
// Earlier URLs in ScriptsLinks take precedence (first match per category/name wins).
func GetAllRemoteHelperMetadata() (map[string]map[string]HelperScriptEntry, error) {
	if ForceRefreshHelpers {
		cachedHelperMetadataByURL = map[string]map[string]map[string]HelperScriptEntry{}
	}

	merged := map[string]map[string]HelperScriptEntry{}
	var firstErr error

	for _, baseURL := range config.Global.ScriptsLinks {
		meta, err := fetchHelperMetadataFromURL(baseURL)
		if err != nil {
			utils.PrintDebug("Failed to fetch helper metadata from %s: %v", baseURL, err)
			if firstErr == nil {
				firstErr = err
			}
			continue
		}
		for category, scripts := range meta {
			if _, ok := merged[category]; !ok {
				merged[category] = map[string]HelperScriptEntry{}
			}
			for name, entry := range scripts {
				if _, exists := merged[category][name]; !exists {
					merged[category][name] = entry // earlier URL wins
				}
			}
		}
	}

	if len(merged) == 0 && firstErr != nil {
		return nil, firstErr
	}
	return merged, nil
}

// downloadExecutable downloads a file and makes it executable
func downloadExecutable(url, destPath string) error {
	resp, err := http.Get(url)
	if err != nil {
		return err
	}
	defer resp.Body.Close()

	if resp.StatusCode != http.StatusOK {
		return fmt.Errorf("HTTP %d", resp.StatusCode)
	}

	// Create parent directory
	if err := os.MkdirAll(filepath.Dir(destPath), utils.PermDir); err != nil {
		return err
	}

	// Write file
	out, err := utils.CreateFileWritable(destPath)
	if err != nil {
		return err
	}
	defer out.Close()

	if _, err := io.Copy(out, resp.Body); err != nil {
		return err
	}

	// Make executable
	return os.Chmod(destPath, utils.PermExec)
}
