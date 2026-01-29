package exec

import "github.com/Justype/condatainer/internal/config"

// Options configures how CondaTainer executes a command inside an Apptainer container.
type Options struct {
	Overlays    []string
	Command     []string
	EnvSettings []string
	BindPaths   []string
	Fakeroot    bool
	WritableImg bool
	HideOutput  bool

	BaseImage    string
	ApptainerBin string
}

func (o Options) ensureDefaults() Options {
	if o.BaseImage == "" {
		if config.Global.BaseImage != "" {
			o.BaseImage = config.Global.BaseImage
		} else {
			o.BaseImage = config.GetBaseImage()
		}
	}
	if o.ApptainerBin == "" {
		o.ApptainerBin = config.Global.ApptainerBin
	}
	if len(o.Command) == 0 {
		o.Command = []string{"bash"}
	}
	return o
}
