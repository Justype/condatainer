package exec

import (
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
	HideOutput     bool
	HidePrompt     bool

	BaseImage    string
	ApptainerBin string

	// PassThruStdin allows the container command to read from stdin
	// This is needed for build scripts that require interactive input (download links, passwords, etc.)
	PassThruStdin bool

	// Stdin specifies a custom input reader (optional, defaults to os.Stdin if PassThruStdin=true)
	Stdin io.Reader
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
