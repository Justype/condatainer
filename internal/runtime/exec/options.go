package exec

import (
	"context"
	"fmt"
	"io"

	"github.com/Justype/condatainer/internal/artifact/meta"
	"github.com/Justype/condatainer/internal/config"
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

	BaseImage    string
	ApptainerBin string

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

// ensureDefaults fills in what the caller left blank and requires a base image
// that exists and is one: there is no overlay-only execution, so a container
// with no root cannot be started at all. Building a missing managed base is the
// caller's job — this package runs images, it does not make them.
func (o Options) ensureDefaults() (Options, error) {
	if o.BaseImage == "" {
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
	if err := meta.CheckBase(o.BaseImage); err != nil {
		return o, err
	}
	if o.ApptainerBin == "" {
		o.ApptainerBin = config.Global.ApptainerBin
	}
	if len(o.Command) == 0 {
		o.Command = []string{"bash"}
	}
	return o, nil
}
