package helper

import (
	"context"
	"errors"
	"fmt"
	"os"
	"os/exec"
	"path/filepath"
	"strings"
	"syscall"

	"github.com/Justype/condatainer/internal/config"
	"github.com/Justype/condatainer/internal/logging"
	"github.com/Justype/condatainer/internal/scheduler"
	"github.com/Justype/condatainer/internal/utils"
)

// RunPlan is the resolved, validated state of a helper invocation.
// PlanRun produces it; ExecutePlan consumes it.
type RunPlan struct {
	Options       RunOptions
	Meta          HelperScriptMeta
	Spec          *scheduler.ResourceSpec
	Scheduler     scheduler.Scheduler
	UserOverlays  []string
	NamedOverlays []string
	Params        map[string]string
	ContainerCmd  string
}

// ErrMissingParam is returned by PlanRun when one or more required #PARAM:
// values are not supplied. The CLI prompts; the server returns HTTP 400.
type ErrMissingParam struct {
	Params []HelperParam
}

func (e *ErrMissingParam) Error() string {
	keys := make([]string, len(e.Params))
	for i, p := range e.Params {
		keys[i] = p.Key
	}
	return "missing required #PARAM: " + strings.Join(keys, ", ")
}

// ErrSingletonRunning indicates the script is marked #SINGLETON: true and at
// least one instance is already running. CLI prints "already running"; server
// returns HTTP 409.
type ErrSingletonRunning struct {
	Existing []*HelperRun
}

func (e *ErrSingletonRunning) Error() string {
	return fmt.Sprintf("singleton helper already running (%d instance(s))", len(e.Existing))
}

// ErrEnvInUse signals that the chosen writable overlay is locked by another
// process. The conflicting path is preserved so callers can render a useful
// message.
type ErrEnvInUse struct {
	Path string
}

// ErrMissingPackages is returned by PlanRun when #IMG_PACKAGES: specs are not
// satisfied in the conda environment inside EnvImg. Specs holds the raw
// unsatisfied specs suitable for micromamba install; Messages holds matching
// human-readable descriptions for display.
// ErrMissingPackages is returned by PlanRun when #IMG_PACKAGES: specs are not
// satisfied in the conda environment inside EnvImg. Specs holds the raw
// unsatisfied specs suitable for micromamba install; Messages holds matching
// human-readable descriptions for display.
// VersionChoices is non-nil when one or more missing packages had an unresolved
// version token ({CONDA_PYTHON} etc.) with a known #VALUE: list — the UI should
// prompt the user to pick a concrete version before installing.
// Keys are package names (e.g. "python"); values are the selectable versions.
type ErrMissingPackages struct {
	EnvImg         string
	Specs          []string            // raw specs for micromamba install, e.g. ["jupyterlab", "python>=3.12"]
	Messages       []string            // human-readable per-spec descriptions
	VersionChoices map[string][]string // package name → selectable versions (from #VALUE:)
}

func (e *ErrMissingPackages) Error() string {
	return fmt.Sprintf("conda packages not satisfied in %s:\n  %s",
		e.EnvImg, strings.Join(e.Messages, "\n  "))
}

func (e *ErrEnvInUse) Error() string {
	return "env overlay in use: " + e.Path
}

// PlanRun validates and resolves opts into a RunPlan without any user I/O.
// All prompts (param menus, resource menus, history menu, guided overlay
// creation) are the caller's responsibility — fill opts first.
//
// Returns typed errors (*ErrMissingParam, *ErrSingletonRunning, *ErrEnvInUse)
// so callers can react contextually. Other errors are plain.
func PlanRun(ctx context.Context, opts RunOptions) (*RunPlan, error) {
	logger := logging.FromContext(ctx)

	meta, err := ParseHelperScriptMeta(opts.ScriptPath)
	if err != nil {
		logger.Debug("planning helper", "warn", "could not parse script meta", "err", err)
	}

	// Singleton scripts forbid a second instance regardless of ForceNew —
	// see README "#SINGLETON: --new/--force prints an error instead".
	if meta.Singleton {
		running, _ := RunningHelpers(opts.ScriptName)
		if len(running) > 0 {
			return nil, &ErrSingletonRunning{Existing: running}
		}
	}

	scriptParams, err := ParseHelperParams(opts.ScriptPath)
	if err != nil {
		return nil, fmt.Errorf("parsing #PARAM: headers: %w", err)
	}
	flagValues, _, err := ApplyParamFlags(scriptParams, opts.FlagArgs)
	if err != nil {
		return nil, fmt.Errorf("parsing param flags: %w", err)
	}
	for k, v := range opts.Params {
		if _, set := flagValues[k]; !set {
			flagValues[k] = v
		}
	}
	for _, p := range scriptParams {
		if _, set := flagValues[p.Key]; set {
			continue
		}
		if p.Optional {
			// KEY=? auto-fills with first #VALUE: entry; falls back to empty.
			if vals := meta.ParamValues[p.Key]; len(vals) > 0 {
				flagValues[p.Key] = vals[0]
			} else {
				flagValues[p.Key] = ""
			}
			continue
		}
		if p.Default != "" {
			flagValues[p.Key] = p.Default
		}
	}
	missing, err := MissingParams(opts.ScriptPath, flagValues)
	if err != nil {
		return nil, err
	}
	if len(missing) > 0 {
		return nil, &ErrMissingParam{Params: missing}
	}
	params := flagValues

	cwd := opts.CWD
	if cwd == "" {
		cwd, _ = os.Getwd()
	}
	if opts.EnvImg != "" {
		logger.Info("Checking env overlay", "path", opts.EnvImg)
		st, err := CheckEnv(ctx, opts.EnvImg)
		if err != nil {
			return nil, fmt.Errorf("checking env overlay: %w", err)
		}
		if st.InUse {
			return nil, &ErrEnvInUse{Path: opts.EnvImg}
		}
		if st.Exists && st.Writable {
			if err := checkOverlayIntegrity(ctx, opts.EnvImg); err != nil {
				return nil, fmt.Errorf("overlay integrity check failed: %w", err)
			}
		}
	}

	// Resolve {KEY} tokens in #BIND: specs from resolved params.
	for _, b := range meta.Binds {
		resolved := b
		for k, v := range params {
			resolved = strings.ReplaceAll(resolved, "{"+k+"}", v)
		}
		opts.ExtraBinds = append(opts.ExtraBinds, resolved)
	}

	userOverlays := opts.Overlays
	if meta.RequiredOverlays != "" {
		logger.Info("Checking required overlays")
	}
	namedOverlays, err := CheckRequiredOverlays(ctx, meta.RequiredOverlays, params)
	if err != nil {
		return nil, err
	}
	if len(namedOverlays) > 0 {
		opts.Overlays = append(append([]string(nil), namedOverlays...), opts.Overlays...)
	}

	if meta.ImgRequired && opts.EnvImg == "" {
		return nil, fmt.Errorf("helper requires a writable overlay (#IMG_PACKAGES set) — create one and pass --env or set EnvImg")
	}

	if err := CheckHelperPackages(meta, opts.EnvImg, params); err != nil {
		return nil, err
	}

	sched, _ := detectScheduler()
	spec := resolveSpec(opts.ScriptPath, opts.Resources)
	if spec.Time == 0 {
		return nil, fmt.Errorf("walltime is required for helpers — set #TIME: in the script or pass walltime explicitly")
	}

	containerCmd, err := buildCondatainerCmd(opts)
	if err != nil {
		return nil, err
	}

	return &RunPlan{
		Options:       opts,
		Meta:          meta,
		Spec:          spec,
		Scheduler:     sched,
		UserOverlays:  userOverlays,
		NamedOverlays: namedOverlays,
		Params:        params,
		ContainerCmd:  containerCmd,
	}, nil
}

// ExecutePlan submits the helper job (or runs it headless). Writes a
// "pending"/"starting" history record and returns the helperID. Status messages
// go through logging.FromContext(ctx); the job's own output goes to the state
// dir's job.log, which the CLI and dashboard tail. ctx bounds only the
// submission/setup — the launched job is detached and outlives the caller.
//
// Does not poll — call helper.Monitor (CLI) or rely on the server watcher to
// observe ready/done events.
func ExecutePlan(ctx context.Context, plan *RunPlan) (string, error) {
	if plan == nil {
		return "", errors.New("nil plan")
	}
	logger := logging.FromContext(ctx)

	cwd := plan.Options.CWD
	if cwd == "" {
		cwd, _ = os.Getwd()
	}

	// Generate the helper ID on the login node before submission so the state
	// directory is known upfront. newHelperID also creates the state dir.
	helperID := newHelperID(plan.Options.ScriptName)
	stateDir := config.GetHelperStateDir(helperID)

	// Pre-create bind source directories that reference the state dir, so
	// Apptainer can mount them before the helper script runs.
	// Only absolute paths (after $CNT_HELPER_STATE_DIR substitution) are created;
	// sources with remaining $VAR references (e.g. $HOME on the compute node) are skipped.
	for _, b := range plan.Options.ExtraBinds {
		src, _, _ := strings.Cut(b, ":")
		src = strings.ReplaceAll(src, "$CNT_HELPER_STATE_DIR", stateDir)
		if filepath.IsAbs(src) {
			_ = utils.MkdirAllShared(src)
		}
	}

	wrapperPath, err := generateWrapper(
		helperID, plan.Options.ScriptName, cwd, filepath.Dir(plan.Options.ScriptPath), stateDir,
		plan.Spec.Time, plan.Params, plan.Spec, plan.Scheduler, plan.ContainerCmd,
	)
	if err != nil {
		os.RemoveAll(stateDir)
		return "", fmt.Errorf("generating wrapper: %w", err)
	}

	if plan.Scheduler != nil && plan.Scheduler.IsAvailable() {
		logger.Info("Submitting helper job", "name", plan.Options.ScriptName)
		jobID, err := plan.Scheduler.Submit(ctx, wrapperPath, nil)
		if err != nil {
			os.RemoveAll(stateDir)
			return "", fmt.Errorf("job submission failed: %w", err)
		}
		_ = AppendHistory(newHelperRun(helperID, plan.Options.ScriptName, jobID, cwd, plan.Spec.Time,
			plan.Options, plan.UserOverlays, plan.Spec, plan.Params, "pending"))
		logger.Info("Helper submitted", "id", helperID)
		return helperID, nil
	}

	logger.Info("Starting helper headless", "name", plan.Options.ScriptName)
	go pruneStaleTmpdirs()
	// Run detached, like a scheduler job: a new session (Setsid) not tied to ctx,
	// so it outlives the server. The job ends only via `helper stop`
	// (KillHeadlessProcess) or the wrapper's walltime watchdog. Output goes to
	// job.log, so stdout/stderr are discarded (nil -> /dev/null).
	cmd := exec.Command("bash", wrapperPath)
	cmd.Stdin = nil
	cmd.Stdout = nil
	cmd.Stderr = nil
	cmd.SysProcAttr = &syscall.SysProcAttr{Setsid: true}
	if err := cmd.Start(); err != nil {
		os.RemoveAll(stateDir)
		return "", fmt.Errorf("starting helper: %w", err)
	}
	// PID == process-group ID (Setsid), so KillHeadlessProcess can signal the
	// whole group. Non-fatal on write error: stop then degrades gracefully
	// (history still marked done by the watcher).
	if err := WriteHelperPid(helperID, cmd.Process.Pid); err != nil {
		logger.Debug("helper: failed to write pid file", "id", helperID, "err", err)
	}
	_ = AppendHistory(newHelperRun(helperID, plan.Options.ScriptName, "", cwd, plan.Spec.Time,
		plan.Options, plan.UserOverlays, plan.Spec, plan.Params, "starting"))
	// Reap the child while the caller lives; init reaps it afterwards.
	go func() { _ = cmd.Wait() }()
	return helperID, nil
}
