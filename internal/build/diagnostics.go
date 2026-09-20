package build

import (
	"bytes"
	"context"
	"errors"
	osexec "os/exec"
	"strings"

	"github.com/Justype/condatainer/internal/artifact/meta"
	"github.com/Justype/condatainer/internal/config"
	"github.com/Justype/condatainer/internal/logging"
	"github.com/Justype/condatainer/internal/runtime/apptainer"
)

const maxRecordedToolVersion = 256

// captureCommonBuildTools records the worker-side implementations shared by
// every build. It is called only after the skip check, so an existing image is
// never relabelled with tools from a later invocation.
//
// Apptainer is read back through apptainer.Current(), never resolved here: by
// the time this runs, whichever binary this build's own container step needed
// (ForBuild for a def build, Normal or Fakeroot for a conda/script one) has
// already been decided, and re-resolving independently would risk recording a
// different binary than the one that actually ran.
func (b *BuildObject) captureCommonBuildTools(ctx context.Context) {
	log := logging.FromContext(ctx)

	version, ok := normalizedToolVersion(config.VERSION)
	if !ok {
		version = meta.Unrecorded
		log.Warn("could not record Condatainer version", "name", b.spec.Image.Name)
	}
	b.buildTools.Condatainer = meta.Tool{Version: version}

	implementation, rawVersion, err := apptainer.Current()
	if err != nil {
		b.buildTools.Apptainer = meta.Tool{Name: "apptainer", Version: meta.Unrecorded}
		log.Warn("could not identify Apptainer version", "name", b.spec.Image.Name, "err", err)
		return
	}
	if version, ok = normalizedToolVersion(rawVersion); !ok {
		version = meta.Unrecorded
		log.Warn("could not record Apptainer version", "name", b.spec.Image.Name)
	}
	b.buildTools.Apptainer = meta.Tool{Name: implementation, Version: version}
}

// captureMicromambaVersion runs the self-provisioned Micromamba's own
// --version directly on host.
func (b *BuildObject) captureMicromambaVersion(ctx context.Context) {
	log := logging.FromContext(ctx)
	version := meta.Unrecorded

	mmCmd, err := micromambaCmd()
	if err == nil {
		var out bytes.Buffer
		cmd := osexec.CommandContext(ctx, mmCmd, "--version")
		cmd.Stdout = &out
		if err = cmd.Run(); err == nil {
			if captured, ok := normalizedToolVersion(out.String()); ok {
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

// captureMksquashfsVersion runs mksquashfsBin's own -version and records its
// first line. mksquashfsBin is already resolved by the caller (squashfs.go);
// a failure here is only ever the version parse, never the build.
func (b *BuildObject) captureMksquashfsVersion(ctx context.Context, mksquashfsBin string) {
	log := logging.FromContext(ctx)
	version := meta.Unrecorded

	var out bytes.Buffer
	cmd := osexec.CommandContext(ctx, mksquashfsBin, "-version")
	cmd.Stdout = &out
	err := cmd.Run()
	if err == nil {
		firstLine, _, _ := strings.Cut(out.String(), "\n")
		if captured, ok := normalizedToolVersion(firstLine); ok {
			version = captured
		} else {
			err = errEmptyToolVersion
		}
	}

	b.buildTools.Mksquashfs = meta.Tool{Version: version}
	if err != nil {
		log.Warn("could not record mksquashfs version", "name", b.spec.Image.Name, "err", err)
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
