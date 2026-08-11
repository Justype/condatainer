package build

import (
	"bytes"
	"context"
	"errors"
	"strings"

	"github.com/Justype/condatainer/internal/artifact/meta"
	"github.com/Justype/condatainer/internal/config"
	"github.com/Justype/condatainer/internal/logging"
	"github.com/Justype/condatainer/internal/runtime/apptainer"
	execpkg "github.com/Justype/condatainer/internal/runtime/exec"
)

const maxRecordedToolVersion = 256

// captureCommonBuildTools records the worker-side implementations shared by
// every build. It is called only after the skip check, so an existing image is
// never relabelled with tools from a later invocation.
func (b *BuildObject) captureCommonBuildTools(ctx context.Context) {
	log := logging.FromContext(ctx)

	version, ok := normalizedToolVersion(config.VERSION)
	if !ok {
		version = meta.Unrecorded
		log.Warn("could not record Condatainer version", "name", b.spec.Image.Name)
	}
	b.buildTools.Condatainer = meta.Tool{Version: version}

	if err := apptainer.EnsureApptainer(); err != nil {
		b.buildTools.Apptainer = meta.Tool{Name: "apptainer", Version: meta.Unrecorded}
		log.Warn("could not identify Apptainer version", "name", b.spec.Image.Name, "err", err)
		return
	}

	implementation := apptainer.Implementation()
	version, err := apptainer.GetVersion()
	if err != nil {
		b.buildTools.Apptainer = meta.Tool{Name: implementation, Version: meta.Unrecorded}
		log.Warn("could not record Apptainer version", "name", b.spec.Image.Name, "err", err)
		return
	}
	if version, ok = normalizedToolVersion(version); !ok {
		version = meta.Unrecorded
		log.Warn("could not record Apptainer version", "name", b.spec.Image.Name)
	}
	b.buildTools.Apptainer = meta.Tool{Name: implementation, Version: version}
}

// captureMicromambaVersion asks the Micromamba inside the resolved build base,
// through the same container execution path used for install and export. A host
// Micromamba binary therefore cannot leak into the manifest.
func (b *BuildObject) captureMicromambaVersion(ctx context.Context) {
	log := logging.FromContext(ctx)
	version := meta.Unrecorded

	opts, err := b.condaExecOpts("micromamba --version", nil)
	if err == nil {
		opts.PassThruStdin = false
		var stdout bytes.Buffer
		err = execpkg.Run(ctx, opts, execpkg.IO{Stdout: &stdout})
		if err == nil {
			if captured, ok := normalizedToolVersion(stdout.String()); ok {
				version = captured
			} else {
				err = errEmptyToolVersion
			}
		}
	}

	b.buildTools.Micromamba = meta.Tool{Version: version}
	if err != nil {
		log.Warn("could not record Micromamba version", "name", b.spec.Image.Name, "err", err)
	}
}

// normalizedToolVersion accepts one short, printable line. Tool output is
// diagnostic data, but it is still embedded metadata and must not become an
// unbounded or multiline log fragment.
func normalizedToolVersion(raw string) (string, bool) {
	version := strings.TrimSpace(raw)
	if version == "" || len(version) > maxRecordedToolVersion || strings.ContainsAny(version, "\r\n") {
		return "", false
	}
	return version, true
}

var errEmptyToolVersion = errors.New("version command returned no usable version")
