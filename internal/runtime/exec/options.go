package exec

import (
	"context"
	"io"

	"github.com/Justype/condatainer/internal/config"
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

func (o Options) ensureDefaults() Options {
	if o.BaseImage == "" {
		o.BaseImage = config.GetBaseImage()
	}
	if o.ApptainerBin == "" {
		o.ApptainerBin = config.Global.ApptainerBin
	}
	if len(o.Command) == 0 {
		o.Command = []string{"bash"}
	}
	return o
}
