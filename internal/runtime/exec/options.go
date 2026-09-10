package exec

import (
	"context"
	"fmt"
	"io"

	"github.com/Justype/condatainer/internal/config"
	"github.com/Justype/condatainer/internal/runtime/container"
	"github.com/Justype/condatainer/internal/utils"
)

// Options configures how CondaTainer executes a command inside an Apptainer container.
type Options struct {
	Overlays       []string
	Command        []string
	EnvSettings    []string
	BindPaths      []string
	ApptainerFlags []string // Flags to pass directly to apptainer (e.g., --home=/path, --nv)
	Fakeroot       bool
	WritableImg    bool
	HidePrompt     bool
	// Activation selects which activate.d scripts container.ActivationScript
	// sources before the command — container.ActivationAll (default when
	// left ""), container.ActivationEnv, or container.ActivationNone, e.g. so
	// a hung or misbehaving activation script can be ruled out.
	Activation container.ActivationMode
	// GpuRequested forces GPU flag detection even when autoload_gpu is disabled,
	// for a command that explicitly declared a GPU requirement.
	GpuRequested bool

	BaseImage string

	// PassThruStdin is retained for callers that track whether stdin is expected.
	// Actual stdin is owned by IO.Stdin so internal execution never assumes a terminal.
	PassThruStdin bool
}

// IO contains caller-owned process streams. Nil streams are silent/no input.
type IO struct {
	Stdin  io.Reader
	Stdout io.Writer
	Stderr io.Writer
}

type ioContextKey struct{}

// WithIO returns a context carrying caller-owned process streams.
func WithIO(ctx context.Context, ioStreams IO) context.Context {
	return context.WithValue(ctx, ioContextKey{}, ioStreams)
}

// IOFromContext returns caller-owned process streams from ctx, if present.
func IOFromContext(ctx context.Context) IO {
	if ioStreams, ok := ctx.Value(ioContextKey{}).(IO); ok {
		return ioStreams
	}
	return IO{}
}

// IsZero reports whether no streams are set.
func (ioStreams IO) IsZero() bool {
	return ioStreams.Stdin == nil && ioStreams.Stdout == nil && ioStreams.Stderr == nil
}

// ensureDefaults fills in what the caller left blank: a command to run.
// BaseImage is resolved separately, after overlay setup — see Prepare —
// since the exec root may come from the requested overlays. The apptainer
// binary is not caller-configurable at all: apptainer.ResolveBin decides it
// from whether this invocation ends up needing fakeroot.
func (o Options) ensureDefaults() (Options, error) {
	if len(o.Command) == 0 {
		o.Command = []string{"bash"}
	}
	if o.Activation == "" {
		o.Activation = container.ActivationAll
	}
	return o, nil
}

// resolveBaseImage finalizes BaseImage now that Setup has run: root, an exec
// root pulled out of the requested overlays, wins when present. Otherwise the
// caller's own BaseImage is kept, falling back to the configured default —
// found, never built: this package runs images, it does not make them, so a
// missing default is the caller's job (internal/build.ResolveBase).
//
// There is no overlay-only execution, so a container with no root cannot be
// started at all — this is checked here rather than letting Apptainer report
// a missing file.
func (o Options) resolveBaseImage(root string) (Options, error) {
	switch {
	case root != "":
		o.BaseImage = root
	case o.BaseImage == "":
		base, err := config.GetBaseImage()
		if err != nil {
			return o, err
		}
		o.BaseImage = base
	}
	// An image file, or a sandbox: apptainer runs either as a root, and a
	// definition build packs its own sandbox by running it.
	if !utils.FileExists(o.BaseImage) && !utils.IsSandboxDir(o.BaseImage) {
		return o, fmt.Errorf("base image not found: %s", o.BaseImage)
	}
	return o, nil
}
