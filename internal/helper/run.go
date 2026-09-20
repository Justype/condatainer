package helper

import (
	"context"
	"encoding/json"
	"errors"
	"fmt"
	"log/slog"
	"os"
	"os/exec"
	"path/filepath"
	"strconv"
	"strings"
	"syscall"
	"time"

	"github.com/Justype/condatainer/catalog"
	"github.com/Justype/condatainer/internal/config"
	"github.com/Justype/condatainer/internal/image"
	"github.com/Justype/condatainer/internal/image/ext3"
	"github.com/Justype/condatainer/internal/logging"
	"github.com/Justype/condatainer/internal/project"
	"github.com/Justype/condatainer/internal/runtime/container"
	"github.com/Justype/condatainer/internal/scheduler"
	"github.com/Justype/condatainer/internal/utils"
)

// RunOptions holds all resolved options for running a helper.
type RunOptions struct {
	ScriptPath string
	ScriptName string
	// Resources holds user-specified resource values (nil = use script/config defaults).
	Resources *scheduler.ResourceSpec
	// Container options
	EnvImg   string   // writable overlay (.img)
	Overlays []string // additional read-only overlays
	CWD      string   // working directory inside job
	// NoProject skips resolving CWD's required overlays through a project's
	// lock even when CWD stands inside one, falling back to the ordinary
	// on-disk check and auto-install unconditionally.
	NoProject bool
	// NoSubmit forces the run headless on this node even when a scheduler is
	// available. The CLI also honors config.Global.SubmitJob (set globally by
	// --no-submit); this field is what lets the server apply the override to
	// one request without touching the process-wide config.
	NoSubmit bool
	// Partition to submit under (empty falls back to config.Global.Scheduler.Partition; ignored when headless)
	Partition string
	// Account to submit under (empty falls back to config.Global.Scheduler.Account; ignored when headless)
	Account string
	// Behaviour
	ForceNew bool
	// ExtraBinds is the resolved list of extra bind-mount specs ("src:dest") from
	// #BIND: headers, with {KEY} tokens already substituted from params.
	// $VAR references are left as-is for shell expansion when the wrapper runs.
	ExtraBinds []string
	// #PARAM: values resolved from flags + prompts (raw args passed through)
	FlagArgs []string
	Params   map[string]string
}

func detectScheduler(noSubmit bool) (scheduler.Scheduler, error) {
	if noSubmit || !config.Global.SubmitJob {
		return nil, scheduler.ErrSchedulerNotFound
	}
	if sched := scheduler.ActiveScheduler(); sched != nil {
		if sched.IsInsideJob() {
			return nil, scheduler.ErrSchedulerNotFound
		}
		return sched, nil
	}
	sched, err := scheduler.DetectSchedulerWithBinary(config.Global.Scheduler.Bin)
	if err != nil {
		return nil, scheduler.ErrSchedulerNotFound
	}
	if sched.IsInsideJob() {
		return nil, scheduler.ErrSchedulerNotFound
	}
	return sched, nil
}

// ParseGPUSpec parses a GPU spec string into a GpuSpec.
// Formats:
//
//	"a100:2"  → {Type:"a100", Count:2}
//	"a100"    → {Type:"a100", Count:1}
//	"2"       → {Type:"gpu",  Count:2}
//	""        → nil
func ParseGPUSpec(s string) *scheduler.GpuSpec {
	s = strings.TrimSpace(s)
	if s == "" {
		return nil
	}
	if idx := strings.LastIndex(s, ":"); idx >= 0 {
		typePart := strings.TrimSpace(s[:idx])
		countPart := strings.TrimSpace(s[idx+1:])
		if n, err := strconv.Atoi(countPart); err == nil && n > 0 {
			if typePart == "" {
				typePart = "gpu"
			}
			return &scheduler.GpuSpec{Type: typePart, Count: n}
		}
	}
	if n, err := strconv.Atoi(s); err == nil && n > 0 {
		return &scheduler.GpuSpec{Type: "gpu", Count: n}
	}
	return &scheduler.GpuSpec{Type: s, Count: 1}
}

// FormatGpuSpec formats a ResourceSpec's GPU field back to a display string.
// Returns "" when no GPU is set.
func FormatGpuSpec(spec *scheduler.ResourceSpec) string {
	if spec == nil || spec.Gpu == nil {
		return ""
	}
	if spec.Gpu.Raw != "" {
		return spec.Gpu.Raw
	}
	if spec.Gpu.Count > 1 {
		return fmt.Sprintf("%s:%d", spec.Gpu.Type, spec.Gpu.Count)
	}
	return spec.Gpu.Type
}

// ResolveHelperSpec returns the effective ResourceSpec for a helper script,
// merging script headers with global build defaults. Used for display and prompts.
// Never errors — missing values are left as zero.
func ResolveHelperSpec(scriptPath string) *scheduler.ResourceSpec {
	return resolveSpec(scriptPath, nil)
}

// resolveSpec merges resources with priority: scheduler defaults < script headers < overrides.
// Recognised script headers: #NCPUS:, #MEM:, #TIME:, #GPU:.
func resolveSpec(scriptPath string, overrides *scheduler.ResourceSpec) *scheduler.ResourceSpec {
	base := config.Global.Scheduler.Defaults

	// Script headers (#NCPUS:/#MEM:/#TIME:/#GPU:) override config defaults.
	meta, err := ParseHelperScriptMeta(scriptPath)
	if err != nil {
		slog.Default().Debug("helper: could not parse script meta", "err", err)
	} else {
		scriptDefaults := &scheduler.ResourceSpec{}
		if meta.NCPUs > 0 {
			scriptDefaults.CpusPerTask = meta.NCPUs
		}
		if meta.MemMB > 0 {
			scriptDefaults.MemPerNodeMB = meta.MemMB
		}
		if meta.Walltime > 0 {
			scriptDefaults.Time = meta.Walltime
		}
		if meta.GPU != "" {
			scriptDefaults.Gpu = ParseGPUSpec(meta.GPU)
		}
		base.Override(scriptDefaults)
	}

	// User overrides (CLI flags or server request) take highest priority.
	base.Override(overrides)
	return &base
}

// checkRoot fails a launch whose container root is missing, so the error
// surfaces here rather than in the job log on a compute node. The node's exec
// takes the root from the project cwd stands in, else the default.
func checkRoot(ctx context.Context, cwd string) error {
	standing, err := project.StandingAt(cwd)
	if err != nil {
		return err
	}
	if standing != nil {
		root, err := standing.Base(ctx)
		if err != nil {
			return err
		}
		if root != "" {
			return nil
		}
	}
	_, err = config.GetBaseImage()
	return err
}

// stopGraceSecs is how long the container gets to exit after the wrapper is
// told to stop, so a writable image is unmounted before anything is killed.
// It sits inside the scheduler's own TERM-to-KILL delay.
const stopGraceSecs = 30

// buildCondatainerCmd constructs the condatainer exec command run by the wrapper
// on the compute node. Uses "condatainer exec" (not "e") so that no env.img
// auto-loading occurs — all overlays are explicit. The helper script runs
// directly via its shebang line; no "bash" wrapper is added.
//
// Generated form:
//
//	/abs/condatainer exec [-o overlay…] [-w] [--gpu] \
//	    --bind /abs/condatainer:/usr/local/bin/condatainer \
//	    [--pwd=/cwd] /abs/helper-script
//
// spec is the resolved ResourceSpec (script #GPU: header merged with any
// override); a GPU request there forces --gpu on the exec invocation so a
// helper still gets --nv/--rocm on a node with autoload_gpu disabled.
func buildCondatainerCmd(opts RunOptions, spec *scheduler.ResourceSpec) (string, error) {
	exe, err := os.Executable()
	if err != nil {
		return "", fmt.Errorf("could not locate condatainer binary: %w", err)
	}

	var parts []string
	parts = append(parts, shellQuote(exe), "exec", "--stop-grace", strconv.Itoa(stopGraceSecs))

	// Named read-only overlays (SquashFS) — code-server, igv, etc.
	for _, ol := range opts.Overlays {
		parts = append(parts, "-o", shellQuote(ol))
	}

	// Env overlay: -w only for a writable .img, not a read-only snapshot .sqf.
	if opts.EnvImg != "" {
		parts = append(parts, "-o", shellQuote(opts.EnvImg))
		if utils.IsImg(opts.EnvImg) {
			parts = append(parts, "-w")
		}
	}

	// Script declared a GPU requirement: force detection past autoload_gpu:false.
	if spec != nil && spec.Gpu != nil && spec.Gpu.Count > 0 {
		parts = append(parts, "--gpu")
	}

	// Bind condatainer binary so helpers can call _server_ready/_server_message inside.
	parts = append(parts, "--bind", shellQuote(exe+":/usr/local/bin/condatainer"))

	// Extra binds from #BIND: headers. Double-quoted so $VAR references (e.g.
	// $CNT_HELPER_STATE_DIR, $HOME) are expanded by the shell when the wrapper runs.
	for _, b := range opts.ExtraBinds {
		parts = append(parts, "--bind", `"`+b+`"`)
	}

	// Working directory inside container — passed through to apptainer as --flag=value.
	if opts.CWD != "" {
		parts = append(parts, "--pwd="+shellQuote(opts.CWD))
	}

	// Helper script runs directly via its shebang; no "bash" prefix.
	parts = append(parts, shellQuote(opts.ScriptPath))
	return strings.Join(parts, " "), nil
}

// ResolveEnvOverlayInDir resolves the project's environment overlay in cwd.
// See container.ResolveEnvOverlay, which this delegates to so
// internal/project can share the same resolution without importing this
// package (which itself imports internal/project).
func ResolveEnvOverlayInDir(envImg, cwd string) string {
	return container.ResolveEnvOverlay(envImg, cwd)
}

// FindEnvSnapshot looks in cwd for an env-typed .sqf. See
// container.FindEnvSnapshot.
func FindEnvSnapshot(cwd string) string {
	return container.FindEnvSnapshot(cwd)
}

// checkOverlayIntegrity calls ext3.CheckIntegrity on the given image file.
func checkOverlayIntegrity(ctx context.Context, imgPath string) error {
	return ext3.CheckIntegrity(ctx, imgPath, false)
}

// resolveOverlayTemplate substitutes {KEY} tokens in the #REQUIRED_OVERLAYS: template
// with values from params, then splits on whitespace to return individual overlay names.
func resolveOverlayTemplate(template string, params map[string]string) []string {
	result := template
	for k, v := range params {
		result = strings.ReplaceAll(result, "{"+k+"}", v)
	}
	return strings.Fields(result)
}

// checkAndInstallNamedOverlays ensures each named overlay exists: resolving a
// bare, partial, or distro-unprefixed name (catalog.SolveName) against what's
// already installed first, no network call — then building via "condatainer
// create <name>", passed the name exactly as declared, only once nothing
// local satisfies it. That subprocess runs the same solver plus its
// Conda-search fallback, so a name nothing could ever satisfy fails there
// with a clear message before any build starts.
// Unlike guided overlay creation (#IMG_PACKAGES:), no prompt is shown — named
// overlays are fixed requirements, not user-configurable.
// Returns the resolved absolute paths to be prepended to opts.Overlays.
func checkAndInstallNamedOverlays(ctx context.Context, names []string) ([]string, error) {
	condaBin, err := os.Executable()
	if err != nil {
		condaBin = "condatainer"
	}
	logger := logging.FromContext(ctx)

	distro := config.ResolvedDefaultDistro()

	resolveAll := func(pending []string, installed map[string][]string) (map[string][]string, []string, error) {
		found := make(map[string][]string, len(pending))
		var stillMissing []string
		for _, name := range pending {
			resolvedName, path, ok, err := catalog.SolveInstalled(ctx, installed, distro, name)
			if err != nil {
				return nil, nil, err
			}
			if !ok {
				stillMissing = append(stillMissing, name)
				continue
			}
			if resolvedName != catalog.Normalize(name) {
				logger.Info(fmt.Sprintf("%s -> %s", name, resolvedName), "kind", "note")
			}
			found[name] = []string{path}
		}
		return found, stillMissing, nil
	}

	installed, err := image.ScanOverlays(image.ScanOptions{})
	if err != nil {
		return nil, err
	}
	resolved, missing, err := resolveAll(names, installed)
	if err != nil {
		return nil, err
	}

	// Second pass: build all missing overlays in one condatainer create call.
	if len(missing) > 0 {
		logger.Info("required overlays not found, building now", "overlays", strings.Join(missing, " "))
		args := append([]string{"create"}, missing...)
		cmd := exec.CommandContext(ctx, condaBin, args...)
		cmdOut := logging.WriterFromCtx(ctx)
		if cmdOut == nil {
			cmdOut = os.Stdout
		}
		cmd.Stdout = cmdOut
		cmd.Stderr = cmdOut
		if err := cmd.Run(); err != nil {
			var exitErr *exec.ExitError
			if errors.As(err, &exitErr) && exitErr.ExitCode() == config.ExitCodeJobsSubmitted {
				return nil, fmt.Errorf("overlay build(s) submitted to scheduler: %s — wait for them to finish, then re-run the helper",
					strings.Join(missing, ", "))
			}
			return nil, fmt.Errorf("condatainer create %s failed: %w", strings.Join(missing, " "), err)
		}
		container.InvalidateInstalledOverlaysCache()
		installed, err = image.ScanOverlays(image.ScanOptions{})
		if err != nil {
			return nil, err
		}
		built, stillMissing, err := resolveAll(missing, installed)
		if err != nil {
			return nil, err
		}
		if len(stillMissing) > 0 {
			return nil, fmt.Errorf("overlay(s) not found after build: %s", strings.Join(stillMissing, ", "))
		}
		for name, p := range built {
			resolved[name] = p
		}
	}

	// Reconstruct in original order so overlay precedence matches #REQUIRED_OVERLAYS.
	var out []string
	for _, name := range names {
		out = append(out, resolved[name]...)
	}
	return out, nil
}

// checkPackages verifies that every package declared in meta.ImgPackages is
// installed in the conda environment inside envImg. {KEY} tokens in ImgPackages
// are first substituted from params.
//
// Skipped (returns nil) when ImgPackages is empty or envImg is empty.
// Uses container.PairedPackages (debugfs on envImg, unsquashfs -l on a
// paired snapshot) — no container launch required. Supports conda-style
// version constraints: =, ==, >=, <=, >, <, !=.
func checkPackages(meta HelperScriptMeta, envImg string, params map[string]string) error {
	if meta.ImgPackages == "" || envImg == "" {
		return nil
	}

	// Substitute {KEY} tokens with resolved param values.
	resolved := meta.ImgPackages
	for k, v := range params {
		resolved = strings.ReplaceAll(resolved, "{"+k+"}", v)
	}

	installed, err := container.PairedPackages(envImg)
	if err != nil {
		return fmt.Errorf("reading conda packages from %s: %w", envImg, err)
	}
	if installed == nil {
		return fmt.Errorf("no conda environment found in %s — run guided overlay setup first", envImg)
	}

	var badSpecs, badMsgs []string
	var versionChoices map[string][]string
	for _, spec := range strings.Fields(resolved) {
		name, op, ver := splitPkgConstraint(spec)
		// If the version value still contains an unresolved {TOKEN}, skip the version
		// check and treat it as a name-only requirement. If the token has a known #VALUE:
		// list, record it so the UI can prompt the user to pick a concrete version.
		if strings.Contains(ver, "{") {
			if strings.HasPrefix(ver, "{") && strings.HasSuffix(ver, "}") {
				tokenName := ver[1 : len(ver)-1]
				if vlist := meta.ParamValues[tokenName]; len(vlist) > 0 {
					if versionChoices == nil {
						versionChoices = make(map[string][]string)
					}
					versionChoices[name] = vlist
				}
			}
			op, ver = "", ""
			spec = name
		}
		instVer, ok := installed[name]
		if !ok {
			badSpecs = append(badSpecs, spec)
			badMsgs = append(badMsgs, name+" (not installed)")
			continue
		}
		if op != "" && !checkPkgConstraint(instVer, op, ver) {
			badSpecs = append(badSpecs, spec)
			badMsgs = append(badMsgs, fmt.Sprintf("%s (installed: %s, need: %s%s)", name, instVer, op, ver))
		}
	}
	if len(badSpecs) == 0 {
		return nil
	}
	return &ErrMissingPackages{EnvImg: envImg, Specs: badSpecs, Messages: badMsgs, VersionChoices: versionChoices}
}

// splitPkgConstraint splits a conda package spec (e.g. "python>=3.12", "jupyterlab=4.2.1")
// into (name, op, version). Single "=" is treated as exact match "==".
// Returns (spec, "", "") when no operator is found.
func splitPkgConstraint(spec string) (name, op, ver string) {
	for _, sep := range []string{">=", "<=", "!=", "==", ">", "<"} {
		if idx := strings.Index(spec, sep); idx >= 0 {
			return spec[:idx], sep, spec[idx+len(sep):]
		}
	}
	// conda single "=" means exact version
	if idx := strings.Index(spec, "="); idx >= 0 {
		return spec[:idx], "==", spec[idx+1:]
	}
	return spec, "", ""
}

// checkPkgConstraint reports whether installedVer satisfies op+requiredVer.
func checkPkgConstraint(installedVer, op, requiredVer string) bool {
	cmp := catalog.CompareVersions(installedVer, requiredVer)
	switch op {
	case ">=":
		return cmp >= 0
	case ">":
		return cmp > 0
	case "<=":
		return cmp <= 0
	case "<":
		return cmp < 0
	case "==":
		return cmp == 0
	case "!=":
		return cmp != 0
	}
	return true
}

// buildHelperCommandBody constructs the bash body embedded by CreateScriptWithSpec
// between its writeJobHeader and writeJobFooter calls. The body sets up helper
// identity env vars, redirects output to the per-ID NFS state dir, picks a free
// port, runs the container command, and writes the done file.
//
// The final `(exit $_cnt_exit)` subshell propagates the container exit code to
// the `_EXIT_CODE=$?` line that CreateScriptWithSpec appends after Command.
func buildHelperCommandBody(id, name, cwd, scriptDir, stateDir string, walltime time.Duration,
	params map[string]string, sched scheduler.Scheduler, containerCmd string, bindAll bool) string {

	jobIDExpr := "$$" // headless: shell PID
	if sched != nil {
		jobIDExpr = sched.JobIDEnvExpr()
	}

	var sb strings.Builder

	// Helper identity — ID and state dir are fixed at submission time on the login
	// node, so the compute node never needs to derive them from scheduler env vars.
	fmt.Fprintf(&sb, "export CNT_HELPER_ID=%s\n", shellQuote(id))
	fmt.Fprintf(&sb, "export CNT_HELPER_NAME=%s\n", shellQuote(name))
	fmt.Fprintf(&sb, "export CNT_HELPER_JOB_ID=\"%s\"\n", jobIDExpr) // kept for reference
	fmt.Fprintf(&sb, "export CNT_HELPER_STATE_DIR=%s\n", shellQuote(stateDir))
	fmt.Fprintf(&sb, "export CNT_HELPER_CWD=%s\n", shellQuote(cwd))
	fmt.Fprintf(&sb, "export CNT_HELPER_WALLTIME_SECS=%d\n", int64(walltime.Seconds()))
	fmt.Fprintf(&sb, "export CNT_HELPER_SCRIPT_DIR=%s\n", shellQuote(scriptDir))
	// SCRATCH is not exported: it would leak into nested condatainer calls, where it
	// selects a data layer. Helper scripts use "${SCRATCH:-$HOME}" themselves.

	// CNT_JOB_TMPDIR: per-job node-local scratch dir. Unix sockets (dbus, X11, PulseAudio,
	// rserver RPC) require local storage and do not work over NFS.
	// For schedulers: scheduler-assigned dir (SLURM_TMPDIR etc.) when available, falling
	// back to /tmp/cnt-$USER/$CNT_HELPER_ID. For headless: resolved at generation time.
	// The directory is created here; cleanup is handled by the wrapper trap below.
	if sched != nil {
		if v := sched.TmpDirVar(); v != "" {
			fmt.Fprintf(&sb, "export CNT_JOB_TMPDIR=\"${%s:-/tmp/cnt-${USER:-condatainer}/${CNT_HELPER_ID}}\"\n", v)
		} else {
			fmt.Fprintln(&sb, `export CNT_JOB_TMPDIR="/tmp/cnt-${USER:-condatainer}/${CNT_HELPER_ID}"`)
		}
	} else {
		// Headless: generated and executed on the same login node, so resolve at generation time.
		fmt.Fprintf(&sb, "export CNT_JOB_TMPDIR=%s\n", shellQuote(filepath.Join(utils.GetTmpDir(), id)))
	}
	fmt.Fprintln(&sb, `mkdir -p "$CNT_JOB_TMPDIR"`)

	// #PARAM: values resolved from flags + prompts
	for k, v := range params {
		fmt.Fprintf(&sb, "export %s=%s\n", k, shellQuote(v))
	}

	// State dir is pre-created on the login node; redirect all output there.
	fmt.Fprintln(&sb)
	fmt.Fprintln(&sb, `exec >> "$CNT_HELPER_STATE_DIR/job.log" 2>&1`)
	fmt.Fprintln(&sb)

	// Free port (resolved on the compute node, not the login node)
	fmt.Fprintln(&sb, `export CNT_HELPER_PORT=$(condatainer _pick_port)`)

	// Bind address: 0.0.0.0 when helper_bind_all is set (direct TCP proxy, no SSH tunnel).
	if bindAll {
		fmt.Fprintln(&sb, `export CNT_HELPER_BIND_ADDR=0.0.0.0`)
		fmt.Fprintln(&sb, `export CNT_HELPER_BIND_ALL=1`)
	} else {
		fmt.Fprintln(&sb, `export CNT_HELPER_BIND_ADDR=127.0.0.1`)
	}
	fmt.Fprintln(&sb)

	// Until the container is running there is nothing to forward to: record done
	// and clean up, so a stop during setup still leaves a finished helper.
	// CNT_JOB_TMPDIR is always cleaned up: headless uses rm -rf on the resolved path;
	// scheduler uses the shell variable so it works for both the scheduler-assigned dir
	// and the /tmp fallback (scheduler-assigned dirs being wiped twice is harmless).
	cleanupTrap := `condatainer _server_done --exit-code 130`
	if sched == nil {
		parentDir := utils.GetTmpDir()
		jobTmpDir := filepath.Join(parentDir, id)
		cleanupTrap += `; kill -KILL -$_cnt_watchdog_pid 2>/dev/null`
		cleanupTrap += fmt.Sprintf(`; rm -rf %s; rmdir %s 2>/dev/null || true`,
			shellQuote(jobTmpDir), shellQuote(parentDir))
	} else {
		cleanupTrap += `; rm -rf "$CNT_JOB_TMPDIR"`
	}
	fmt.Fprintf(&sb, "trap '%s; exit 130' TERM INT\n", cleanupTrap)

	// Headless walltime watchdog: background subshell that ignores TERM/INT so it
	// survives the group SIGTERM it sends, then SIGKILLs after a 2-minute grace period.
	// set -m puts it in its OWN process group so it can be torn down with a group kill
	// (kill -KILL -$_cnt_watchdog_pid) that also reaps its sleep child; disown then
	// drops it from the job table so bash prints no "Killed" notice into job.log.
	if sched == nil {
		fmt.Fprintln(&sb, `_cnt_main_pid=$$`)
		fmt.Fprintln(&sb, `set -m`)
		fmt.Fprintf(&sb,
			"( trap '' TERM INT; sleep %d; kill -TERM -$_cnt_main_pid 2>/dev/null; sleep 120; kill -KILL -$_cnt_main_pid 2>/dev/null ) &\n",
			int64(walltime.Seconds()))
		fmt.Fprintln(&sb, `_cnt_watchdog_pid=$!`)
		fmt.Fprintln(&sb, `disown`)
		fmt.Fprintln(&sb, `set +m`)
	}
	fmt.Fprintln(&sb)

	// cd to CWD before invoking condatainer so relative overlay paths resolve correctly.
	if cwd != "" {
		fmt.Fprintf(&sb, "cd %s\n", shellQuote(cwd))
	}

	// Run the helper inside the container in the background and wait on it: bash
	// defers a trap while a foreground command runs, so a scheduler's SIGTERM
	// would never reach it. The trap forwards TERM to the container and the loop
	// keeps waiting, so the container exits (and unmounts its images) before the
	// wrapper does. The trap is armed before the launch, and a TERM that lands
	// before the container has a pid is forwarded right after it starts.
	fmt.Fprintln(&sb, `_cnt_child=`)
	fmt.Fprintln(&sb, `_cnt_stopped=0`)
	fmt.Fprintln(&sb, `trap '_cnt_stopped=1; [ -n "$_cnt_child" ] && kill -TERM $_cnt_child 2>/dev/null' TERM INT`)
	fmt.Fprintf(&sb, "%s &\n", containerCmd)
	fmt.Fprintln(&sb, `_cnt_child=$!`)
	fmt.Fprintln(&sb, `[ $_cnt_stopped = 1 ] && kill -TERM $_cnt_child 2>/dev/null`)
	fmt.Fprintln(&sb, `while :; do`)
	fmt.Fprintln(&sb, `  wait $_cnt_child`)
	fmt.Fprintln(&sb, `  _cnt_exit=$?`)
	fmt.Fprintln(&sb, `  kill -0 $_cnt_child 2>/dev/null || break`)
	fmt.Fprintln(&sb, `done`)
	fmt.Fprintln(&sb, `[ $_cnt_stopped = 1 ] && _cnt_exit=130`)
	fmt.Fprintln(&sb)

	// Disarm the trap, kill the watchdog, then record done with the exit code.
	fmt.Fprintln(&sb, `trap - TERM INT`)
	if sched == nil {
		fmt.Fprintln(&sb, `kill -KILL -$_cnt_watchdog_pid 2>/dev/null`)
	}
	fmt.Fprintln(&sb, `condatainer _server_done --exit-code $_cnt_exit`)
	if sched == nil {
		parentDir := utils.GetTmpDir()
		jobTmpDir := filepath.Join(parentDir, id)
		fmt.Fprintf(&sb, "rm -rf %s\n", shellQuote(jobTmpDir))
		fmt.Fprintf(&sb, "rmdir %s 2>/dev/null || true\n", shellQuote(parentDir))
	} else {
		fmt.Fprintln(&sb, `rm -rf "$CNT_JOB_TMPDIR"`)
	}
	fmt.Fprintln(&sb, `(exit $_cnt_exit)`)

	return sb.String()
}

// pruneStaleTmpdirs removes leftover /tmp/cnt-$USER/{id} dirs for helpers
// that are marked done in history. Handles SIGKILL leaks where the bash
// wrapper's trap never fired. Called in a goroutine — errors are ignored.
func pruneStaleTmpdirs() {
	parentDir := utils.GetTmpDir()
	entries, err := os.ReadDir(parentDir)
	if err != nil {
		return
	}
	for _, e := range entries {
		if !e.IsDir() {
			continue
		}
		id := e.Name()
		run := HistoryEntryForID(id)
		if run != nil && run.Status == "done" {
			os.RemoveAll(filepath.Join(parentDir, id))
		}
	}
	os.Remove(parentDir) // no-op if other dirs remain
}

// newHelperID returns a stable, human-readable helper run ID of the form
// "{name}-YYYYMMDD-HHMMSS". The timestamp is fixed at call time on the login
// node so the ID (and its state directory) are known before job submission.
// Collision within the same second is prevented by probing the state dir:
// if it already exists a "-2", "-3", … suffix is appended.
func newHelperID(name string) string {
	base := name + "-" + time.Now().Format("20060102-150405")
	stateDir := config.GetHelperStateDir(base)
	if err := os.Mkdir(stateDir, utils.PermDir); err == nil {
		return base
	}
	for i := 2; i <= 99; i++ {
		candidate := fmt.Sprintf("%s-%d", base, i)
		if err := os.Mkdir(config.GetHelperStateDir(candidate), utils.PermDir); err == nil {
			return candidate
		}
	}
	return base // fallback: extremely unlikely to reach here
}

// buildHelperScriptSpecs constructs a ScriptSpecs for use with CreateScriptWithSpec.
// Stdout and Stderr are /dev/null because the command body uses `exec >>` to
// redirect output to the per-ID state dir; the scheduler-level log is unused.
// account/partition fall back to config.Global.Scheduler.Account/.Partition when empty —
// there is no script-header tier for either (unlike ncpus/mem/time/gpu), since a script
// author has no way to know a user's cluster account.
func buildHelperScriptSpecs(name, cwd, account, partition string, spec *scheduler.ResourceSpec) *scheduler.ScriptSpecs {
	if account == "" {
		account = config.Global.Scheduler.Account
	}
	if partition == "" {
		partition = config.Global.Scheduler.Partition
	}
	ss := &scheduler.ScriptSpecs{
		Spec: spec,
		Control: scheduler.RuntimeConfig{
			JobName:   "cnt-" + name,
			WorkDir:   cwd,
			Stdout:    "/dev/null",
			Stderr:    "/dev/null",
			Account:   account,
			Partition: partition,
		},
		HasDirectives: spec != nil,
	}
	if config.Global.ProxyPerJob {
		if h, err := os.Hostname(); err == nil && h != "" {
			ss.ProxyVia = h
		}
	}
	return ss
}

// generateWrapper writes the job submission script into the helper's state
// directory (id/stateDir) so the script, logs, and state files all live together.
// Returns the path to pass to scheduler.Submit and any error.
// No cleanup func is needed — the state dir is permanent until the user removes it.
//
// For headless (sched == nil): writes job.sh directly into stateDir.
// For schedulers: delegates to sched.CreateScriptWithSpec using stateDir as outputDir.
func generateWrapper(id, name, cwd, scriptDir, stateDir, account, partition string, walltime time.Duration,
	params map[string]string, spec *scheduler.ResourceSpec, sched scheduler.Scheduler,
	containerCmd string) (string, error) {

	body := buildHelperCommandBody(id, name, cwd, scriptDir, stateDir, walltime, params, sched, containerCmd, config.Global.HelperBindAll)

	if err := utils.MkdirAllShared(stateDir); err != nil {
		return "", fmt.Errorf("creating state dir: %w", err)
	}

	// Headless: plain bash script in the state dir.
	if sched == nil || !sched.IsAvailable() {
		shPath := filepath.Join(stateDir, "job.sh")
		f, err := utils.CreateFileWritable(shPath)
		if err != nil {
			return "", err
		}
		fmt.Fprintln(f, "#!/bin/bash")
		fmt.Fprint(f, body)
		f.Close()
		if err := utils.MakeExecutable(shPath); err != nil {
			return "", err
		}
		return shPath, nil
	}

	jobSpec := &scheduler.JobSpec{
		Name:       "cnt-" + name,
		Command:    body,
		Specs:      buildHelperScriptSpecs(name, cwd, account, partition, spec),
		KeepScript: true,
	}

	scriptPath, err := sched.CreateScriptWithSpec(jobSpec, stateDir)
	if err != nil {
		return "", fmt.Errorf("creating job script: %w", err)
	}
	return scriptPath, nil
}

func shellQuote(s string) string {
	if s == "" {
		return "''"
	}
	if !strings.ContainsAny(s, " \t\n\"'\\$`!") {
		return s
	}
	return "'" + strings.ReplaceAll(s, "'", "'\\''") + "'"
}

// normalizeOverlayForHistory converts a resolved overlay path to a portable form
// for storage in history:
//   - Paths directly inside an image search directory → logical name ("igv/2.19.7")
//   - Other absolute paths → relative to cwd when short enough, else absolute
//   - Non-absolute strings (already logical names or relative paths) → unchanged
//
// A :ro/:rw suffix is preserved through normalization.
func normalizeOverlayForHistory(path, cwd string) string {
	if path == "" {
		return ""
	}
	suffix := ""
	p := path
	for _, s := range []string{":ro", ":rw"} {
		if strings.HasSuffix(p, s) {
			suffix = s
			p = p[:len(p)-len(s)]
			break
		}
	}
	if !filepath.IsAbs(p) {
		return path // already a logical name or relative path
	}
	// Internal: file lives directly inside one of the image search dirs.
	// Derive the logical name from the filename ("igv--2.19.7.sqf" → "igv/2.19.7").
	for _, dir := range config.GetImageSearchPaths() {
		rel, err := filepath.Rel(dir, p)
		if err != nil || strings.Contains(rel, string(filepath.Separator)) || strings.HasPrefix(rel, "..") {
			continue
		}
		name := strings.TrimSuffix(rel, filepath.Ext(rel))
		name = strings.ReplaceAll(name, "--", "/")
		return name + suffix
	}
	// External: store relative to cwd when the path stays nearby (no deep "../../").
	if cwd != "" {
		if rel, err := filepath.Rel(cwd, p); err == nil && !strings.HasPrefix(rel, "..") {
			return rel + suffix
		}
	}
	return path
}

// newHelperRun constructs a HelperRun with all resolved run parameters for history recording.
// userOverlays must be the user-supplied -o overlays only (not required overlays from metadata,
// which are re-derived at rerun time to avoid duplication).
// Overlay paths are normalized to logical names (internal) or CWD-relative paths (external)
// before storage so history entries are portable across installs.
func newHelperRun(id, name, jobID, cwd string, walltime time.Duration,
	opts RunOptions, userOverlays []string, spec *scheduler.ResourceSpec, params map[string]string, status, runner string) *HelperRun {
	normalizedOverlays := make([]string, len(userOverlays))
	for i, ol := range userOverlays {
		normalizedOverlays[i] = normalizeOverlayForHistory(ol, cwd)
	}
	run := &HelperRun{
		ID:         id,
		Name:       name,
		JobID:      jobID,
		Runner:     runner,
		CWD:        cwd,
		Walltime:   walltime,
		GPU:        FormatGpuSpec(spec),
		EnvOverlay: normalizeOverlayForHistory(opts.EnvImg, cwd),
		Overlays:   normalizedOverlays,
		Params:     params,
		StartedAt:  time.Now(),
		Status:     status,
		BindAll:    config.Global.HelperBindAll,
	}
	if spec != nil {
		run.CPUs = spec.CpusPerTask
		if mb := spec.GetMemPerNodeMB(); mb > 0 {
			run.Mem = utils.FormatMemoryMB(mb)
		}
	}
	return run
}

// EnsureServer auto-starts the condatainer server if not running.
// Delegates to "condatainer server start" so the fork logic stays in one place.
// Prints the dashboard URL on first start.
func EnsureServer() {
	pidFile := config.GetServerPidFilePath()
	if pidFile == "" {
		return
	}
	if isServerAlive(pidFile) {
		return
	}

	exe, err := os.Executable()
	if err != nil {
		return
	}

	// "server start" handles port selection, log file, pipe readiness, and printing.
	cmd := exec.Command(exe, "server", "start")
	cmd.Stdout = os.Stdout
	cmd.Stderr = os.Stderr
	if err := cmd.Run(); err != nil {
		slog.Default().Debug("server auto-start failed", "err", err)
	}
}

func isServerAlive(pidFile string) bool {
	ss := readServerState(pidFile)
	if ss == nil || ss.PID == 0 {
		return false
	}
	p, err := os.FindProcess(ss.PID)
	if err != nil {
		return false
	}
	return p.Signal(syscall.Signal(0)) == nil
}

func readServerState(path string) *config.ServerState {
	data, err := os.ReadFile(path)
	if err != nil {
		return nil
	}
	var ss config.ServerState
	if err := json.Unmarshal(data, &ss); err != nil {
		return nil
	}
	return &ss
}

// RecentHistoryByName returns the last N distinct-CWD runs for a helper name.
func RecentHistoryByName(name string, max int) []*HelperRun {
	runs, _ := LoadHistory()
	var result []*HelperRun
	seen := make(map[string]bool)
	for i := len(runs) - 1; i >= 0; i-- {
		r := runs[i]
		if r.Name != name {
			continue
		}
		if seen[r.CWD] {
			continue
		}
		seen[r.CWD] = true
		result = append(result, r)
		if len(result) >= max {
			break
		}
	}
	return result
}

// HistoryEntryForID returns the most recent JSONL record for the given helper ID.
func HistoryEntryForID(id string) *HelperRun {
	runs, _ := LoadHistory()
	for i := len(runs) - 1; i >= 0; i-- {
		if runs[i].ID == id {
			return runs[i]
		}
	}
	return nil
}
