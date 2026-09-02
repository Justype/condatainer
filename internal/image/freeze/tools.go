package freeze

import (
	"context"

	"github.com/Justype/condatainer/internal/artifact/meta"
	"github.com/Justype/condatainer/internal/config"
	"github.com/Justype/condatainer/internal/logging"
	"github.com/Justype/condatainer/internal/runtime/apptainer"
)

// buildTools records what produced the artifact. Micromamba is not among them: a
// freeze packs an environment that already exists and never solves one.
func buildTools(ctx context.Context) meta.BuildTools {
	log := logging.FromContext(ctx)
	tools := meta.BuildTools{Condatainer: meta.Tool{Version: config.VERSION}}
	if tools.Condatainer.Version == "" {
		tools.Condatainer.Version = meta.Unrecorded
	}

	implementation := apptainer.Implementation()
	version, err := apptainer.GetVersion()
	if err != nil || version == "" {
		log.Warn("could not record Apptainer version", "err", err)
		version = meta.Unrecorded
	}
	tools.Apptainer = meta.Tool{Name: implementation, Version: version}
	return tools
}
