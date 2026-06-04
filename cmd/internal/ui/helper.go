package ui

import (
	"context"
	"fmt"
	"os"
	"os/signal"
	"path/filepath"
	"slices"
	"sort"
	"strconv"
	"strings"
	"syscall"
	"time"

	"github.com/Justype/condatainer/internal/config"
	cntexec "github.com/Justype/condatainer/internal/exec"
	"github.com/Justype/condatainer/internal/helper"
	"github.com/Justype/condatainer/internal/overlay"
	"github.com/Justype/condatainer/internal/scheduler"
	"github.com/Justype/condatainer/internal/utils"
)

// OfferRecentSessions shows a combined list of running instances and recent history,
// deduplicated by CWD. Running entries are marked [RUNNING] with their access URL.
//
// Return values:
//   - chosen != nil, isRunning=true  → caller should print URL and exit (no new launch)
//   - chosen != nil, isRunning=false → caller should call ApplyHistoryRun then continue
//   - startNew=true                  → caller continues with defaults
func OfferRecentSessions(ctx context.Context, name string,
	running []*helper.HelperRun, serverPort int,
) (chosen *helper.HelperRun, isRunning bool, startNew bool, err error) {
	// Build a CWD→run map for quick lookup of running jobs.
	runningByCWD := make(map[string]*helper.HelperRun, len(running))
	for _, r := range running {
		runningByCWD[r.CWD] = r
	}

	// Start from recent history (up to 5 unique CWDs).
	history := helper.RecentHistoryByName(name, 5)

	// Build the display list: prefer running entry for any CWD in history.
	type entry struct {
		run       *helper.HelperRun
		isRunning bool
	}
	var entries []entry
	seenCWD := make(map[string]bool, len(history))
	for _, h := range history {
		seenCWD[h.CWD] = true
		if r, ok := runningByCWD[h.CWD]; ok {
			entries = append(entries, entry{r, true})
		} else {
			entries = append(entries, entry{h, false})
		}
	}
	// Append running jobs whose CWD wasn't in the history top-5.
	for _, r := range running {
		if !seenCWD[r.CWD] {
			entries = append(entries, entry{r, true})
		}
	}

	if len(entries) == 0 {
		return nil, false, true, nil
	}

	fmt.Printf("Recent %s sessions:\n", utils.StyleName(name))
	for i, e := range entries {
		r := e.run
		cwd := r.CWD
		if cwd == "" {
			cwd = "(default working directory)"
		}
		res := formatRunResources(r)
		if res == "" {
			res = "default resources"
		}
		age := time.Since(r.StartedAt).Round(time.Minute)
		if e.isRunning {
			url := r.AccessURL(serverPort)
			fmt.Printf("  %d. %-36s  %-14s  [RUNNING] started %s ago\n", i+1, cwd, res, age)
			if url != "" {
				fmt.Printf("     %s\n", utils.StyleDebug("→ "+url))
			}
		} else {
			fmt.Printf("  %d. %-36s  %-14s  last used %s ago\n", i+1, cwd, res, age)
		}
	}
	fmt.Print("  n. Start new\n\n[?] Choose [")
	for i := range entries {
		fmt.Printf("%d/", i+1)
	}
	fmt.Print("n]: ")

	line, readErr := utils.ReadLineContext(ctx)
	if readErr != nil {
		return nil, false, false, readErr
	}
	if line == "" || line == "n" {
		return nil, false, true, nil
	}
	idx, convErr := strconv.Atoi(line)
	if convErr != nil || idx < 1 || idx > len(entries) {
		return nil, false, true, nil
	}
	e := entries[idx-1]
	return e.run, e.isRunning, false, nil
}

// ApplyHistoryRun copies a previous run's settings into opts. All values are
// pre-filled and will appear in the PromptSettings table for the user to review
// and edit before launch.
func ApplyHistoryRun(opts *helper.RunOptions, r *helper.HelperRun) {
	if r.CWD != "" {
		opts.CWD = r.CWD
	}
	if r.CPUs > 0 || r.Mem != "" || r.Walltime > 0 || r.GPU != "" {
		if opts.Resources == nil {
			opts.Resources = &scheduler.ResourceSpec{}
		}
		if r.CPUs > 0 {
			opts.Resources.CpusPerTask = r.CPUs
		}
		if r.Mem != "" {
			if mb, err := utils.ParseMemoryMB(r.Mem); err == nil {
				opts.Resources.MemPerNodeMB = mb
			}
		}
		if r.Walltime > 0 {
			opts.Resources.Time = r.Walltime
		}
		if r.GPU != "" {
			opts.Resources.Gpu = helper.ParseGPUSpec(r.GPU)
		}
	}
	if r.BaseImage != "" {
		opts.BaseImage = r.BaseImage
	}
	if r.EnvOverlay != "" {
		opts.EnvImg = r.EnvOverlay
	}
	if len(r.Overlays) > 0 {
		opts.Overlays = append([]string(nil), r.Overlays...)
	}
	if len(r.Params) > 0 {
		opts.Params = make(map[string]string, len(r.Params))
		for k, v := range r.Params {
			opts.Params[k] = v
		}
	}
}

// PrintDefaultSettings prints the script's default resource settings before prompts.
// OfferSingletonBlocked prints an informative message and returns a non-nil error
// when a singleton helper is already running.
func OfferSingletonBlocked(name string, running []*helper.HelperRun, serverPort int) error {
	utils.PrintMessage("%s is a singleton helper and is already running:", utils.StyleName(name))
	for _, r := range running {
		age := time.Since(r.StartedAt).Round(time.Minute)
		url := r.AccessURL(serverPort)
		utils.PrintMessage("  %s  %s  Started %s ago", r.ID, r.Node, age)
		if url != "" {
			utils.PrintMessage("  Access: %s", url)
		}
	}
	utils.PrintMessage("Stop the existing instance first, or use 'condatainer helper --status'.")
	return fmt.Errorf("singleton helper %s is already running", name)
}

// PromptSettings prints all settings in one combined block, then enters an edit
// loop where the user can type "key value" to update any field. An empty line
// accepts all values and continues. Pass spec=nil to skip resources.
// Empty-default params pass through as empty — the script handles them.
func PromptSettings(ctx context.Context,
	params []helper.HelperParam, values, saved map[string]string, versions map[string][]string,
	opts *helper.RunOptions, spec *scheduler.ResourceSpec,
) (map[string]string, error) {
	// ── resolve params ────────────────────────────────────────────────────────
	result := make(map[string]string, len(params))
	for k, v := range values {
		result[k] = v
	}

	type paramEntry struct {
		p           helper.HelperParam
		displayKey  string
		versionList []string
	}
	paramEntries := make([]paramEntry, 0, len(params))
	for _, p := range params {
		defVal := p.Default
		vlist := versions[p.Key]
		// KEY=? auto-fills with first (latest) value from #VALUE: list.
		if p.Optional && len(vlist) > 0 && defVal == "" {
			defVal = vlist[0]
		}
		// History params seed defaults (below saved config, above script defaults).
		if hp, ok := opts.Params[p.Key]; ok && hp != "" {
			defVal = hp
		}
		if sv, ok := saved[strings.ToLower(p.Key)]; ok && sv != "" {
			defVal = sv
		}
		if flagVal, fromFlag := values[p.Key]; fromFlag && len(vlist) > 0 {
			// Validate flag-provided value against known versions; resolve partials.
			matched := utils.MatchVersion(flagVal, vlist)
			if matched != flagVal {
				result[p.Key] = matched // silently resolve partial (e.g. "4.4" → "4.4.3")
			} else if !slices.Contains(vlist, flagVal) {
				return nil, fmt.Errorf("invalid value %q for --%s: not in known versions (%s)",
					flagVal, strings.TrimPrefix(p.LongFlag, "--"), utils.FormatChoicesInline(vlist))
			}
		} else if _, ok := result[p.Key]; !ok {
			result[p.Key] = defVal
		}
		dk := p.Key
		if p.ShortFlag != "" {
			dk = strings.TrimPrefix(p.ShortFlag, "-")
		}
		paramEntries = append(paramEntries, paramEntry{p, dk, vlist})
	}

	showResources := spec != nil

	// ── resolve resource display values ───────────────────────────────────────
	// Priority: CLI flags (opts.Resources) > saved config > script defaults (spec).
	cpus := 0
	memMB := int64(0)
	var walltime time.Duration
	if opts.Resources != nil {
		cpus = opts.Resources.CpusPerTask
		memMB = opts.Resources.GetMemPerNodeMB()
		walltime = opts.Resources.Time
	}
	// Saved config resource defaults.
	if cpus == 0 {
		if v, ok := saved["cpus"]; ok {
			if n, err := strconv.Atoi(v); err == nil && n > 0 {
				cpus = n
			}
		}
	}
	if memMB == 0 {
		if v, ok := saved["mem"]; ok {
			if mb, err := utils.ParseMemoryMB(v); err == nil && mb > 0 {
				memMB = mb
			}
		}
	}
	if walltime == 0 {
		if v, ok := saved["time"]; ok {
			if d, err := utils.ParseWalltime(v); err == nil && d > 0 {
				walltime = d
			}
		}
	}
	if spec != nil {
		if cpus == 0 && spec.CpusPerTask > 0 {
			cpus = spec.CpusPerTask
		}
		if memMB == 0 {
			if mb := spec.GetMemPerNodeMB(); mb > 0 {
				memMB = mb
			}
		}
		if walltime == 0 && spec.Time > 0 {
			walltime = spec.Time
		}
	}
	cpuDisp := fmt.Sprintf("%d", cpus)
	if cpus == 0 {
		cpuDisp = "(none)"
	}
	memDisp := utils.FormatMemoryMB(memMB)
	if memMB == 0 {
		memDisp = "(none)"
	}
	timeDisp := utils.FormatDuration(walltime)
	if walltime == 0 {
		timeDisp = "(none)"
	}
	gpuDisp := helper.FormatGpuSpec(opts.Resources)
	if gpuDisp == "" {
		if v, ok := saved["gpu"]; ok && v != "" && v != "(none)" {
			if opts.Resources == nil {
				opts.Resources = &scheduler.ResourceSpec{}
			}
			opts.Resources.Gpu = helper.ParseGPUSpec(v)
			gpuDisp = v
		}
	}
	if gpuDisp == "" {
		gpuDisp = "(none)"
	}

	// ── compute column width for path/env/overlay/param rows ─────────────────
	maxDK := 1
	for _, e := range paramEntries {
		if len(e.displayKey) > maxDK {
			maxDK = len(e.displayKey)
		}
	}

	printRow := func(key, val, desc, choices string) {
		valPad := strings.Repeat(" ", max(0, 20-len(val)))
		choicesStr := ""
		if choices != "" {
			choicesStr = "  " + utils.StyleDebug(choices)
		}
		fmt.Printf("  %-*s: %s%s  %s%s\n", maxDK, key, val, valPad, utils.StyleDebug(desc), choicesStr)
	}

	// ── print combined table ──────────────────────────────────────────────────
	utils.PrintMessage("Settings:")

	// Resources: c/m/t on one line, g on its own.
	if showResources {
		fmt.Printf("  c: %s  %s   m: %s  %s   t: %s  %s   g: %s  %s\n",
			cpuDisp, utils.StyleDebug("CPUs"),
			memDisp, utils.StyleDebug("Memory"),
			timeDisp, utils.StyleDebug("Walltime"),
			gpuDisp, utils.StyleDebug("GPU"))
	}

	// Path and env img (always shown).
	effectiveCwd := opts.CWD
	if effectiveCwd == "" {
		effectiveCwd, _ = os.Getwd()
	}
	cwdDisp := opts.CWD
	if cwdDisp == "" {
		cwdDisp = effectiveCwd + " (cwd)"
	}
	envDisp := opts.EnvImg
	if envDisp == "" {
		envDisp = "(none)"
	} else if rel, err := filepath.Rel(effectiveCwd, envDisp); err == nil && !strings.HasPrefix(rel, "..") {
		envDisp = rel
	}
	printRow("w", cwdDisp, "Working directory", "")
	printRow("e", envDisp, "Writable overlay (.img)", "")

	// Read-only overlays (omitted when none).
	if len(opts.Overlays) > 0 {
		printRow("o", strings.Join(opts.Overlays, ", "), "Read-only overlays", "")
	}

	// Script params at the end.
	for _, e := range paramEntries {
		rawVal := result[e.p.Key]
		if rawVal == "" && e.p.Optional {
			rawVal = "(optional)"
		}
		choices := ""
		if len(e.versionList) > 0 {
			choices = utils.FormatChoicesInline(e.versionList)
		}
		printRow(e.displayKey, rawVal, e.p.Desc, choices)
	}

	// ── build param key lookup ────────────────────────────────────────────────
	byParamKey := make(map[string]int, len(paramEntries)*2)
	for i, e := range paramEntries {
		byParamKey[strings.ToLower(e.p.Key)] = i
		byParamKey[strings.ToLower(e.displayKey)] = i
	}

	// ── edit loop ─────────────────────────────────────────────────────────────
	fmt.Print("[Enter to continue, or key value to update]: ")
	for {
		line, err := utils.ReadLineContext(ctx)
		if err != nil {
			return nil, err
		}
		if strings.TrimSpace(line) == "" {
			break
		}
		k, v, ok := strings.Cut(strings.TrimSpace(line), " ")
		if !ok || strings.TrimSpace(v) == "" {
			fmt.Print("[?] Format: key value — try again: ")
			continue
		}
		k, v = strings.ToLower(strings.TrimSpace(k)), strings.TrimSpace(v)

		// Resource keys.
		if showResources {
			if opts.Resources == nil {
				opts.Resources = &scheduler.ResourceSpec{}
			}
			handled := true
			var confirm string
			switch k {
			case "c", "cpus":
				n, err := strconv.Atoi(v)
				if err != nil || n <= 0 {
					fmt.Printf("[?] Invalid CPUs value %q — try again: ", v)
					continue
				}
				opts.Resources.CpusPerTask = n
				confirm = fmt.Sprintf("%d", n)
			case "m", "mem", "memory":
				mb, err := utils.ParseMemoryMB(v)
				if err != nil {
					fmt.Printf("[?] Invalid memory %q — try again: ", v)
					continue
				}
				opts.Resources.MemPerNodeMB = mb
				confirm = utils.FormatMemoryMB(mb)
			case "t", "time", "walltime":
				d, err := utils.ParseWalltime(v)
				if err != nil {
					fmt.Printf("[?] Invalid walltime %q — try again: ", v)
					continue
				}
				opts.Resources.Time = d
				confirm = utils.FormatDuration(d)
			case "g", "gpu":
				opts.Resources.Gpu = helper.ParseGPUSpec(v)
				confirm = v
			default:
				handled = false
			}
			if handled {
				fmt.Printf("  %s updated to %s\n", k, confirm)
				fmt.Print("[Enter to continue, or key value to update]: ")
				continue
			}
		}

		// Path / env / overlay keys.
		switch k {
		case "w":
			if v == "-" {
				opts.CWD = ""
				fmt.Println("  w cleared")
			} else {
				opts.CWD = v
				fmt.Printf("  w updated to %s\n", v)
			}
			fmt.Print("[Enter to continue, or key value to update]: ")
			continue
		case "e":
			if v == "-" {
				opts.EnvImg = ""
				fmt.Println("  e cleared")
			} else {
				opts.EnvImg = v
				fmt.Printf("  e updated to %s\n", v)
			}
			fmt.Print("[Enter to continue, or key value to update]: ")
			continue
		case "o":
			if v == "-" {
				opts.Overlays = nil
				fmt.Println("  o cleared")
			} else {
				opts.Overlays = strings.Fields(v)
				fmt.Printf("  o updated to %s\n", v)
			}
			fmt.Print("[Enter to continue, or key value to update]: ")
			continue
		}

		// Param keys.
		idx, found := byParamKey[k]
		if !found {
			var validKeys []string
			if showResources {
				validKeys = append(validKeys, "c", "m", "t", "g")
			}
			validKeys = append(validKeys, "w", "e", "o")
			for _, e := range paramEntries {
				validKeys = append(validKeys, e.displayKey)
			}
			fmt.Printf("[?] Unknown key %q (valid: %s) — try again: ", k, strings.Join(validKeys, ", "))
			continue
		}
		e := paramEntries[idx]
		if len(e.versionList) > 0 {
			matched := utils.MatchVersion(v, e.versionList)
			if matched != v {
				fmt.Printf("  → resolved to %s\n", matched)
				v = matched
			} else if !slices.Contains(e.versionList, v) {
				fmt.Printf("[?] %q is not a valid value (%s) — try again: ", v, utils.FormatChoicesInline(e.versionList))
				continue
			}
		}
		result[e.p.Key] = v
		fmt.Printf("  %s updated to %s\n", k, v)
		fmt.Print("[Enter to continue, or key value to update]: ")
	}

	return result, nil
}

// PrintLaunchSpec prints the resolved launch settings before a helper starts.
func PrintLaunchSpec(plan *helper.RunPlan) {
	opts := plan.Options
	params := plan.Params
	spec := plan.Spec
	userOverlays := plan.UserOverlays

	utils.PrintMessage("Helper launch settings:")
	utils.PrintMessage("  Helper:   %s", opts.ScriptName)
	if opts.CWD != "" {
		utils.PrintMessage("  CWD:      %s", opts.CWD)
	}
	if spec != nil {
		if spec.CpusPerTask > 0 {
			utils.PrintMessage("  CPUs:     %d", spec.CpusPerTask)
		}
		if mb := spec.GetMemPerNodeMB(); mb > 0 {
			utils.PrintMessage("  Memory:   %s", utils.FormatMemoryMB(mb))
		}
		if spec.Time > 0 {
			utils.PrintMessage("  Walltime: %s", utils.FormatDuration(spec.Time))
		}
	}
	if gpu := helper.FormatGpuSpec(plan.Spec); gpu != "" {
		utils.PrintMessage("  GPU:      %s", gpu)
	}
	base := opts.BaseImage
	if base == "" {
		base = config.GetBaseImage()
	}
	if base != "" {
		utils.PrintMessage("  Base:     %s", base)
	}
	if opts.EnvImg != "" {
		utils.PrintMessage("  Env:      %s", opts.EnvImg)
	}
	if len(userOverlays) > 0 {
		utils.PrintMessage("  Overlay:  %s", strings.Join(userOverlays, ", "))
	}
	if len(params) > 0 {
		keys := make([]string, 0, len(params))
		for k := range params {
			keys = append(keys, k)
		}
		sort.Strings(keys)
		for _, k := range keys {
			utils.PrintMessage("  %s: %s", k, params[k])
		}
	}
}

// WaitBeforeStart prints a countdown and blocks until done or ctx is cancelled.
func WaitBeforeStart(ctx context.Context, d time.Duration) error {
	utils.PrintMessage("Starting in %s. Press Ctrl-C to cancel.", d.Round(time.Second))
	timer := time.NewTimer(d)
	defer timer.Stop()
	select {
	case <-ctx.Done():
		return ctx.Err()
	case <-timer.C:
		return nil
	}
}

// GuidedOverlayCreate runs the interactive wizard for creating a writable conda
// overlay and installing the packages specified by meta.ImgPackages.
// cwd is the working directory used to compute the default overlay path (env.img
// lives next to the project so FindEnvOverlay can auto-detect it on later runs).
// Returns the path of the created overlay.
func GuidedOverlayCreate(ctx context.Context, helperName string, meta helper.HelperScriptMeta, params map[string]string, versions map[string][]string, cwd string) (string, error) {
	if cwd == "" {
		cwd, _ = os.Getwd()
	}
	defaultPath := filepath.Join(cwd, "env.img")

	utils.PrintMessage("%s requires a writable overlay (conda environment).", utils.StyleName(helperName))
	utils.PrintMessage("No overlay specified with -e/--env.")

	imgPath, err := promptDefault(ctx, "[?] Overlay path", defaultPath)
	if err != nil {
		return "", err
	}
	sizeStr, err := promptDefault(ctx, "[?] Overlay size", "20G")
	if err != nil {
		return "", err
	}
	sizeMB, err := utils.ParseSizeToMB(sizeStr)
	if err != nil {
		return "", fmt.Errorf("invalid size %q: %w", sizeStr, err)
	}
	creationVars, err := promptCreationTokens(ctx, meta.ImgPackages, params, versions)
	if err != nil {
		return "", err
	}

	extraPkgs, err := promptDefault(ctx, "[?] Extra packages (space-separated, or Enter to skip)", "")
	if err != nil {
		return "", err
	}

	pkgTemplate := meta.ImgPackages
	for k, v := range params {
		pkgTemplate = strings.ReplaceAll(pkgTemplate, "{"+k+"}", v)
	}
	for k, v := range creationVars {
		pkgTemplate = strings.ReplaceAll(pkgTemplate, "{"+k+"}", v)
	}
	allPkgs := strings.Fields(pkgTemplate)
	if extraPkgs != "" {
		allPkgs = append(allPkgs, strings.Fields(extraPkgs)...)
	}

	if err := os.MkdirAll(filepath.Dir(imgPath), utils.PermDir); err != nil {
		return "", fmt.Errorf("creating overlay directory: %w", err)
	}

	opts := &overlay.CreateOptions{
		Path:    imgPath,
		SizeMB:  sizeMB,
		Profile: overlay.ProfileSmall,
	}
	io := cntexec.IO{Stdin: os.Stdin, Stdout: os.Stdout, Stderr: os.Stderr}
	utils.PrintMessage("Creating overlay and installing packages: %s", cntexec.DescribeInitialCondaPackages(allPkgs))
	if err := cntexec.CreateCondaOverlay(ctx, opts, allPkgs, meta.PostInstallCmd, false, io); err != nil {
		return "", fmt.Errorf("overlay creation failed: %w", err)
	}

	return imgPath, nil
}

// GuidedInstallMissing prompts the user to install packages that failed the
// #IMG_PACKAGES: check. Prints the unsatisfied specs, asks Y/n, then runs
// micromamba install inside the overlay if confirmed.
func GuidedInstallMissing(ctx context.Context, e *helper.ErrMissingPackages) error {
	utils.PrintWarning("conda packages not satisfied in %s:", e.EnvImg)
	for _, m := range e.Messages {
		fmt.Printf("  %s\n", m)
	}

	// For packages with unresolved version tokens, prompt the user to pick a version.
	// Build pkgVersion (package name → chosen version) to inject into the install spec.
	pkgVersion := make(map[string]string)
	for pkg, vlist := range e.VersionChoices {
		choices := utils.VersionChoicesDisplay(vlist)
		utils.PrintMessage("Available %s: %s", pkg, strings.Join(choices, ", "))
		label := strings.ToUpper(pkg[:1]) + pkg[1:] + " version"
		val, err := promptDefault(ctx, "[?] "+label, utils.LatestVersion(vlist))
		if err != nil {
			return err
		}
		if val != "" {
			matched := utils.MatchVersion(val, vlist)
			if matched != val {
				utils.PrintMessage("  → resolved to %s", matched)
				val = matched
			}
			pkgVersion[pkg] = val
		}
	}

	fmt.Print("[?] Install missing packages? [Y/n]: ")
	line, err := utils.ReadLineContext(ctx)
	if err != nil {
		return err
	}
	if line != "" && strings.ToLower(strings.TrimSpace(line)) != "y" {
		return fmt.Errorf("aborted — install missing packages manually, then retry")
	}

	// Build final install specs, injecting chosen versions where applicable.
	specs := make([]string, len(e.Specs))
	for i, spec := range e.Specs {
		if ver, ok := pkgVersion[spec]; ok {
			specs[i] = spec + "=" + ver
		} else {
			specs[i] = spec
		}
	}
	utils.PrintMessage("Installing: %s", strings.Join(specs, " "))
	return cntexec.InstallPackages(ctx, e.EnvImg, specs, false, cntexec.IO{Stdin: os.Stdin, Stdout: os.Stdout, Stderr: os.Stderr})
}

// MonitorHelper polls NFS state files until the job finishes or the user presses
// Ctrl+C. Returns nil once the service is ready or the job completes cleanly.
func MonitorHelper(ctx context.Context, id string) error {
	sig := make(chan os.Signal, 1)
	signal.Notify(sig, syscall.SIGTERM, syscall.SIGINT)
	defer signal.Stop(sig)

	utils.PrintMessage("Waiting for service to start...")
	ticker := time.NewTicker(2 * time.Second)
	defer ticker.Stop()

	var msgOffset int64
	var readySeen bool

	for {
		select {
		case <-sig:
			utils.PrintMessage("Job still running. Use 'condatainer helper --status' to check.")
			return nil

		case <-ctx.Done():
			return nil

		case <-ticker.C:
			r := helper.PollOnce(id, msgOffset)
			msgOffset = r.NewOffset

			for _, m := range r.NewMessages {
				switch m.Level {
				case "warn":
					utils.PrintWarning("%s", m.Text)
				case "error":
					utils.PrintError("%s", m.Text)
				default:
					utils.PrintMessage("%s", m.Text)
				}
			}

			if r.Ready != nil && !readySeen {
				readySeen = true
				rs := r.Ready
				_ = helper.UpdateHistoryRun(id, func(run *helper.HelperRun) {
					run.Status = "running"
					run.Node = rs.Node
					run.Port = rs.Port
					run.URLPath = rs.URLPath
					run.ExternalURL = rs.ExternalURL
					run.StartedAt = rs.Timestamp
					run.Walltime = time.Duration(rs.WalltimeSec) * time.Second
				})
				tmp := &helper.HelperRun{ID: id, Port: rs.Port, URLPath: rs.URLPath, ExternalURL: rs.ExternalURL}
				if u := tmp.AccessURL(config.GetRunningServerPort()); u != "" {
					utils.PrintSuccess("Access at: %s", u)
				}
				if n := config.Global.Notification; n == "terminal" || n == "both" {
					fmt.Print("\a")
					time.Sleep(1100 * time.Millisecond)
					fmt.Print("\a")
				}
				return nil
			}

			if r.Done != nil {
				status := "failed"
				if r.Done.ExitCode == 0 {
					status = "done"
				} else {
					utils.PrintError("Helper %s failed (exit %d). Log: %s",
						id, r.Done.ExitCode, helper.JobLogFilePath(id))
				}
				_ = helper.UpdateHistoryStatus(id, status, r.Done.Ts)
				if r.Done.ExitCode != 0 {
					return fmt.Errorf("helper %s exited with non-zero code", id)
				}
				utils.PrintMessage("Helper %s finished.", id)
				return nil
			}
		}
	}
}

// FormatRunResources returns a compact resource summary like "4c/16GB/12h/a100:1".
func FormatRunResources(r *helper.HelperRun) string { return formatRunResources(r) }

func formatRunResources(r *helper.HelperRun) string {
	var parts []string
	if r.CPUs > 0 {
		parts = append(parts, fmt.Sprintf("%dc", r.CPUs))
	}
	if r.Mem != "" {
		parts = append(parts, r.Mem)
	}
	if r.Walltime > 0 {
		parts = append(parts, utils.FormatDuration(r.Walltime))
	}
	if r.GPU != "" {
		parts = append(parts, r.GPU)
	}
	return strings.Join(parts, "/")
}

func promptDefault(ctx context.Context, prompt, defaultVal string) (string, error) {
	if defaultVal != "" {
		fmt.Printf("%s [%s]: ", prompt, defaultVal)
	} else {
		fmt.Printf("%s: ", prompt)
	}
	line, err := utils.ReadLineContext(ctx)
	if err != nil {
		return "", err
	}
	if line == "" {
		return defaultVal, nil
	}
	return line, nil
}

func promptCreationTokens(ctx context.Context, imgPackages string, params map[string]string, versions map[string][]string) (map[string]string, error) {
	result := make(map[string]string)
	for _, tok := range utils.ExtractImgPackageTokens(imgPackages) {
		if _, ok := params[tok]; ok {
			continue
		}
		vlist := versions[tok]
		defVal := utils.LatestVersion(vlist)
		if len(vlist) > 0 {
			choices := utils.VersionChoicesDisplay(vlist)
			utils.PrintMessage("Available %s: %s", tok, strings.Join(choices, ", "))
		}
		label := formatTokenLabel(tok)
		val, err := promptDefault(ctx, "[?] "+label, defVal)
		if err != nil {
			return nil, err
		}
		if len(vlist) > 0 && val != "" {
			matched := utils.MatchVersion(val, vlist)
			if matched != val {
				utils.PrintMessage("  → resolved to %s", matched)
				val = matched
			}
		}
		result[tok] = val
	}
	return result, nil
}

// OfferStopPicker prompts the user to choose which running instance(s) to stop.
// When name is empty, lists all running helpers with their names shown.
// Returns the runs to stop: one entry when the user picks a number, all when "a" is entered.
// Returns nil slice (no error) when the user cancels with an invalid input.
func OfferStopPicker(ctx context.Context, name string, running []*helper.HelperRun) ([]*helper.HelperRun, error) {
	if name == "" {
		fmt.Println("Running helpers:")
	} else {
		fmt.Printf("Multiple %s instances running:\n", utils.StyleName(name))
	}
	for i, r := range running {
		age := time.Since(r.StartedAt).Round(time.Minute)
		cwd := r.CWD
		if cwd == "" {
			cwd = "(no cwd)"
		}
		res := formatRunResources(r)
		if name == "" {
			fmt.Printf("  %d. %-20s  %-36s  %-14s  started %s ago\n", i+1, utils.StyleName(r.Name), cwd, res, age)
		} else {
			fmt.Printf("  %d. %-36s  %-14s  started %s ago\n", i+1, cwd, res, age)
		}
	}
	fmt.Print("  a. Stop all\n  n. None\n\n[?] Which to stop [")
	for i := range running {
		fmt.Printf("%d/", i+1)
	}
	fmt.Print("a/n]: ")

	line, err := utils.ReadLineContext(ctx)
	if err != nil {
		return nil, err
	}
	line = strings.TrimSpace(line)
	if line == "n" || line == "" {
		return nil, nil // cancelled
	}
	if line == "a" {
		return running, nil
	}
	idx, convErr := strconv.Atoi(line)
	if convErr != nil || idx < 1 || idx > len(running) {
		return nil, nil // cancelled
	}
	return []*helper.HelperRun{running[idx-1]}, nil
}

func formatTokenLabel(tok string) string {
	name := strings.TrimPrefix(tok, "CONDA_")
	name = strings.TrimPrefix(name, "POSIT_")
	words := strings.Fields(strings.ReplaceAll(name, "_", " "))
	for i, w := range words {
		if len(w) > 0 {
			words[i] = strings.ToUpper(w[:1]) + strings.ToLower(w[1:])
		}
	}
	return strings.Join(words, " ") + " version"
}
