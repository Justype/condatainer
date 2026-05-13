package exec

import (
	"context"
	"strings"

	"github.com/Justype/condatainer/internal/apptainer"
	"github.com/Justype/condatainer/internal/container"
	"github.com/Justype/condatainer/internal/logging"
	"github.com/Justype/condatainer/internal/overlay"
	"github.com/Justype/condatainer/internal/proxy"
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
	OverlayLocks []*overlay.Lock
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
	options = options.ensureDefaults()

	if err := apptainer.SetBin(options.ApptainerBin); err != nil {
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
	})
	if err != nil {
		return nil, err
	}

	// Auto-enable fakeroot if needed for writable .img
	fakeroot, fakerootDiagnostics := container.AutoEnableFakeroot(setupResult.LastImg, options.WritableImg, setupResult.Fakeroot)

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
	var execLocks []*overlay.Lock
	releaseLocks := func() {
		for _, l := range execLocks {
			l.Close()
		}
	}
	for _, ol := range setupResult.Overlays {
		if !utils.FileExists(ol) || utils.DirExists(ol) {
			continue
		}
		// Skip .img files: Apptainer flocks them itself during execution, so
		// acquiring our own lock would conflict. They are pre-checked via
		// container.Setup() → CheckAvailable() before exec starts.
		if utils.IsImg(ol) {
			continue
		}
		lock, err := overlay.AcquireLock(ol, false) // .sqf: always read-only
		if err != nil {
			releaseLocks()
			return nil, err
		}
		execLocks = append(execLocks, lock)
	}
	if utils.FileExists(options.BaseImage) {
		lock, err := overlay.AcquireLock(options.BaseImage, false)
		if err != nil {
			releaseLocks()
			return nil, err
		}
		execLocks = append(execLocks, lock)
	}

	// Inject SOCKS5 proxy env vars if a proxy tunnel is active.
	// Checks per-job local proxy first, then shared NFS proxy.
	// Auto-starts a per-job proxy if proxy_perjob=true or CNT_PROXY_PERJOB=1.
	// Prepend so explicit --env flags from the user take precedence.
	envList := setupResult.EnvList
	if proxyURL, ok := proxy.GetJobProxy(); ok {
		envList = append(proxy.ProxyEnvList(proxyURL), envList...)
	}

	opts := &apptainer.ExecOptions{
		Bind:       setupResult.BindPaths,
		Overlay:    setupResult.OverlayArgs,
		Env:        envList,
		Fakeroot:   fakeroot,
		Additional: setupResult.ApptainerFlags,
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
