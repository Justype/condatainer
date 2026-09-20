package exec

import (
	"context"
	"fmt"
	"strings"

	"github.com/Justype/condatainer/internal/config"
	"github.com/Justype/condatainer/internal/image"
	"github.com/Justype/condatainer/internal/libexec"
	"github.com/Justype/condatainer/internal/logging"
	"github.com/Justype/condatainer/internal/runtime/apptainer"
	"github.com/Justype/condatainer/internal/runtime/container"
	"github.com/Justype/condatainer/internal/runtime/proxy"
	"github.com/Justype/condatainer/internal/utils"
)

// Diagnostic is a non-fatal execution message returned to presentation layers.
type Diagnostic struct {
	Level   string
	Message string
}

// Plan is a prepared container execution. Call Close if RunPrepared is not used.
type Plan struct {
	Options      Options
	Setup        *container.SetupResult
	Fakeroot     bool
	EnvList      []string
	Diagnostics  []Diagnostic
	ExecOptions  *apptainer.ExecOptions
	OverlayLocks []*image.Lock
}

// Close releases resources acquired during Prepare.
func (p *Plan) Close() {
	if p == nil {
		return
	}
	for _, l := range p.OverlayLocks {
		l.Close()
	}
	p.OverlayLocks = nil
}

// Prepare resolves and validates container execution inputs without printing.
func Prepare(ctx context.Context, options Options) (*Plan, error) {
	options, err := options.ensureDefaults()
	if err != nil {
		return nil, err
	}

	// Use shared container setup logic
	setupResult, err := container.Setup(container.SetupConfig{
		Overlays:       options.Overlays,
		WritableImg:    options.WritableImg,
		EnvSettings:    options.EnvSettings,
		BindPaths:      options.BindPaths,
		Fakeroot:       options.Fakeroot,
		ApptainerFlags: options.ApptainerFlags,
		GpuRequested:   options.GpuRequested,
		BindLibexec:    options.BindLibexec,
	})
	if err != nil {
		return nil, err
	}

	options, err = options.resolveBaseImage(setupResult.Root)
	if err != nil {
		return nil, err
	}

	// Auto-enable fakeroot if needed for writable .img
	fakeroot, fakerootDiagnostics := container.AutoEnableFakeroot(setupResult.LastImg, options.WritableImg, setupResult.Fakeroot)

	// Resolved last, once fakeroot is final: fakeroot (explicit or
	// auto-enabled) needs the system apptainer, everything else uses libexec's
	// when installed and the system one otherwise.
	resolve := apptainer.Normal
	if fakeroot {
		resolve = apptainer.Fakeroot
	}
	bin, err := resolve()
	if err != nil {
		if config.IsInsideContainer() {
			return nil, fmt.Errorf("cannot start a nested container: %w", err)
		}
		return nil, err
	}

	logger := logging.FromContext(ctx)
	logger.Debug("prepared exec setup",
		"overlays", setupResult.Overlays,
		"bind_paths", setupResult.BindPaths,
		"overlay_mounts", setupResult.OverlayArgs,
		"env", setupResult.EnvList,
		"command", strings.Join(options.Command, " "),
		"env_overrides", options.EnvSettings,
	)

	diagnostics := make([]Diagnostic, 0, len(setupResult.Diagnostics)+len(fakerootDiagnostics))
	for _, d := range setupResult.Diagnostics {
		diagnostics = append(diagnostics, Diagnostic{Level: d.Level, Message: d.Message})
	}
	for _, d := range fakerootDiagnostics {
		diagnostics = append(diagnostics, Diagnostic{Level: d.Level, Message: d.Message})
	}

	// Acquire read locks on all overlay files for the duration of exec.
	// This prevents concurrent remove/update from deleting files in use.
	var execLocks []*image.Lock
	releaseLocks := func() {
		for _, l := range execLocks {
			l.Close()
		}
	}
	for _, ol := range setupResult.Overlays {
		if !utils.FileExists(ol) || utils.DirExists(ol) {
			continue
		}
		// Skip .img files: Apptainer locks them itself during execution, so
		// acquiring our own lock would conflict. They are pre-checked via
		// container.Setup() → CheckAvailable() before exec starts.
		if utils.IsImg(ol) {
			continue
		}
		lock, err := image.AcquireLock(ol, false) // .sqf: always read-only
		if err != nil {
			releaseLocks()
			return nil, err
		}
		execLocks = append(execLocks, lock)
	}
	if utils.FileExists(options.BaseImage) {
		lock, err := image.AcquireLock(options.BaseImage, false)
		if err != nil {
			releaseLocks()
			return nil, err
		}
		execLocks = append(execLocks, lock)
	}

	// Same reason and duration as the overlay locks above; AcquireUse returns
	// (nil, nil) when nothing is provisioned, so this is a no-op until then.
	if lock, err := libexec.AcquireUse(); err != nil {
		releaseLocks()
		return nil, err
	} else if lock != nil {
		execLocks = append(execLocks, lock)
	}

	// Inject SOCKS5 proxy env vars if a proxy tunnel is active (inside a job only).
	// Checks per-job local proxy first, then shared NFS proxy.
	// Prepend so explicit --env flags from the user take precedence.
	envList := setupResult.EnvList
	if proxyURL, ok := proxy.GetJobProxy(); ok {
		envList = append(proxy.ProxyEnvList(proxyURL), envList...)
	}

	// A mounted overlay's activate.d is shell script, not a static list, so
	// apptainer's own --env can't express it: the command is rewrapped to
	// source it (and define the mm helper) first. Skipped when neither
	// contributes anything for this mount.
	preamble := container.ActivationScript(setupResult.Overlays, setupResult.EnvMounted, options.Activation)
	preamble += container.MMHelperScript(setupResult.EnvMounted)
	if preamble != "" {
		options.Command = wrapWithActivation(preamble, options.Command)
	}

	opts := &apptainer.ExecOptions{
		Bin:        bin,
		Bind:       setupResult.BindPaths,
		Overlay:    setupResult.OverlayArgs,
		Env:        envList,
		Fakeroot:   fakeroot,
		Additional: setupResult.ApptainerFlags,
		StopGrace:  options.StopGrace,
	}

	return &Plan{
		Options:      options,
		Setup:        setupResult,
		Fakeroot:     fakeroot,
		EnvList:      envList,
		Diagnostics:  diagnostics,
		ExecOptions:  opts,
		OverlayLocks: execLocks,
	}, nil
}

// wrapWithActivation rewraps command into one bash -c invocation that runs
// activation first, then execs the original command/argv in its place —
// activation is shell script sourced into the running shell, so it must run
// ahead of the real command inside the same process rather than as a
// separate step.
func wrapWithActivation(activation string, command []string) []string {
	script := activation + "exec \"$@\"\n"
	return append([]string{"bash", "-c", script, "cnt-activate"}, command...)
}

// RunPrepared executes a prepared plan with caller-owned IO streams.
func RunPrepared(ctx context.Context, plan *Plan, ioStreams IO) error {
	defer plan.Close()

	opts := *plan.ExecOptions
	opts.Stdin = ioStreams.Stdin
	opts.Stdout = ioStreams.Stdout
	opts.Stderr = ioStreams.Stderr

	if err := apptainer.Exec(ctx, plan.Options.BaseImage, plan.Options.Command, &opts); err != nil {
		return err
	}
	return nil
}

// Run executes a command inside a configured Apptainer container. It is silent by default;
// callers that want terminal or web output must pass explicit IO writers.
func Run(ctx context.Context, options Options, ioStreams IO) error {
	plan, err := Prepare(ctx, options)
	if err != nil {
		return err
	}
	return RunPrepared(ctx, plan, ioStreams)
}
