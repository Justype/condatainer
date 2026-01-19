package exec

import "github.com/Justype/condatainer/internal/config"

// Options configures how CondaTainer executes a command inside an Apptainer container.
type Options struct {
	Overlays      []string
	Command       []string
	EnvSettings   []string
	BindPaths     []string
	Fakeroot      bool
	WritableImg   bool
	CaptureOutput bool

	BaseImage    string
	ApptainerBin string
}

func (o Options) ensureDefaults() Options {
	if o.BaseImage == "" {
		o.BaseImage = config.Global.BaseImage
	}
	if o.ApptainerBin == "" {
		o.ApptainerBin = config.Global.ApptainerBin
	}
	if len(o.Command) == 0 {
		o.Command = []string{"bash"}
	}
	return o
}
