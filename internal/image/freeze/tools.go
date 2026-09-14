package freeze

import (
	"bytes"
	"context"
	"errors"
	"os/exec"
	"strings"

	"github.com/Justype/condatainer/internal/artifact/meta"
	"github.com/Justype/condatainer/internal/config"
	"github.com/Justype/condatainer/internal/logging"
)

// buildTools records what produced the artifact, apart from what Pack itself
// resolves (Mksquashfs, Fuse2fs) — merged in by the caller once Pack has run.
// Neither Apptainer nor Micromamba is among them: a freeze never runs
// Apptainer (mount.go's "the same apptainer-free way"), and it packs an
// environment that already exists rather than solving one.
func buildTools(ctx context.Context) meta.BuildTools {
	version := config.VERSION
	if version == "" {
		version = meta.Unrecorded
	}
	return meta.BuildTools{Condatainer: meta.Tool{Version: version}}
}

// maxRecordedToolVersion bounds a captured version string the same way
// internal/build's own guard does: diagnostic data, but still embedded
// metadata, and must not become an unbounded or multiline log fragment.
const maxRecordedToolVersion = 256

func normalizedToolVersion(raw string) (string, bool) {
	version := strings.TrimSpace(raw)
	if version == "" || len(version) > maxRecordedToolVersion || strings.ContainsAny(version, "\r\n") {
		return "", false
	}
	return version, true
}

var errEmptyToolVersion = errors.New("version command returned no usable version")

// mksquashfsVersion runs bin's own -version and records its first line. bin is
// already resolved by the caller (Pack); a failure here is only ever the
// version parse, never the pack itself.
func mksquashfsVersion(ctx context.Context, bin string) meta.Tool {
	var out bytes.Buffer
	cmd := exec.CommandContext(ctx, bin, "-version")
	cmd.Stdout = &out
	err := cmd.Run()
	if err == nil {
		firstLine, _, _ := strings.Cut(out.String(), "\n")
		if version, ok := normalizedToolVersion(firstLine); ok {
			return meta.Tool{Version: version}
		}
		err = errEmptyToolVersion
	}
	logging.FromContext(ctx).Warn("could not record mksquashfs version", "err", err)
	return meta.Tool{Version: meta.Unrecorded}
}

// fuse2fsVersion runs bin's own -V and records the first line of its stderr,
// where e2fsprogs tools print their version banner. Same shape as
// mksquashfsVersion.
func fuse2fsVersion(ctx context.Context, bin string) meta.Tool {
	var errOut bytes.Buffer
	cmd := exec.CommandContext(ctx, bin, "-V")
	cmd.Stderr = &errOut
	err := cmd.Run()
	if err == nil {
		firstLine, _, _ := strings.Cut(errOut.String(), "\n")
		if version, ok := normalizedToolVersion(firstLine); ok {
			return meta.Tool{Version: version}
		}
		err = errEmptyToolVersion
	}
	logging.FromContext(ctx).Warn("could not record fuse2fs version", "err", err)
	return meta.Tool{Version: meta.Unrecorded}
}
