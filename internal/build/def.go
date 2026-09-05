package build

import (
	"bytes"
	"context"
	"errors"
	"fmt"
	"os"
	"path/filepath"
	"strings"
	"time"

	"github.com/Justype/condatainer/catalog"
	"github.com/Justype/condatainer/internal/artifact/meta"
	"github.com/Justype/condatainer/internal/config"
	"github.com/Justype/condatainer/internal/logging"
	"github.com/Justype/condatainer/internal/runtime/apptainer"
	execpkg "github.com/Justype/condatainer/internal/runtime/exec"
	"github.com/Justype/condatainer/internal/utils"
)

// buildDef builds an image from an Apptainer definition, real or synthesized
// from a scheme:// URI. Apptainer writes a sandbox, the staged metadata is
// copied into it, and it is packed to .sqf like every other build.
func (b *BuildObject) buildDef(ctx context.Context) error {
	targetPath := b.tgt.Path
	log := logging.FromContext(ctx)

	if skip, err := checkShouldBuild(b); skip || err != nil {
		b.Cleanup(err != nil) //nolint:errcheck
		return err
	}

	if err := b.createBuildLock(); err != nil {
		b.Cleanup(true) //nolint:errcheck
		return err
	}
	defer b.removeBuildLock()
	preparedPath := b.tgt.Prepared
	if pulled, err := b.tryPrebuilt(ctx); err != nil || pulled {
		b.Cleanup(err != nil) //nolint:errcheck
		return err
	}

	if b.skipIfInstalled(ctx) {
		b.Cleanup(false) //nolint:errcheck
		return nil
	}

	log.Info("building image", "kind", "note", "image", filepath.Base(targetPath),
		"mode", buildModeLabel(b), "source", b.buildSource)

	done := watchContext(ctx, "def build")
	defer close(done)

	// Ensure the tmp directory exists before apptainer tries to write the sandbox there.
	if err := ensureWorkspaceRoot(b); err != nil {
		b.Cleanup(true) //nolint:errcheck
		return err
	}

	// Resolve the definition source. A scheme:// source (docker://ubuntu:22.04)
	// has no def file, so synthesize one from the URI; a real .def gets its {key}
	// placeholders substituted (Apptainer reads directives like From: verbatim).
	defSource := b.buildSource
	if strings.Contains(b.buildSource, "://") {
		synthPath, err := synthesizeDefFromURI(b.buildSource, b.ws.Root)
		if err != nil {
			b.Cleanup(true) //nolint:errcheck
			return fmt.Errorf("failed to synthesize def from %s: %w", b.buildSource, err)
		}
		defSource = synthPath
	}

	defData, err := os.ReadFile(defSource)
	if err != nil {
		b.Cleanup(true) //nolint:errcheck
		return fmt.Errorf("failed to read definition %s: %w", defSource, err)
	}

	// A synthesized definition is the only record of what a scheme:// build was.
	b.captureSynthesizedRecipe(defSource, defData)
	b.captureCommonBuildTools(ctx)

	// Before the build: the digest goes into the identity record and into the
	// definition Apptainer is handed, so asking after the pull would be too late.
	b.resolveUpstream(ctx, defData)

	// A def build stages its metadata *into* the definition, so the records are
	// written before Apptainer runs — possible only because the upstream digest
	// was resolved above rather than read off the finished image.
	if err := b.deriveKeys(ctx); err != nil {
		b.Cleanup(true)
		return err
	}

	metaDir, err := stageMetadata(ctx, b)
	if err != nil {
		b.Cleanup(true)
		return err
	}

	// Pin the bootstrap. The metadata is copied into the sandbox after the build.
	buildDefSource, err := writePinnedDef(defSource, b.ws.Root, b.spec.Source.UpstreamDigest())
	if err != nil {
		b.Cleanup(true) //nolint:errcheck
		return fmt.Errorf("failed to prepare pinned def: %w", err)
	}

	log.Info("running apptainer build", "source", b.buildSource)

	buildOpts := &apptainer.BuildOptions{
		Force:     false,
		NoCleanup: false,
		Sandbox:   true,
		TmpDir:    b.ws.BaseRoot,
	}

	if err := apptainer.Build(ctx, b.ws.Sandbox, buildDefSource, buildOpts); err != nil {
		b.Cleanup(true)
		if apptainer.IsBuildCancelled(err) {
			log.Warn("build cancelled, image unchanged", "image", filepath.Base(targetPath))
			return ErrBuildCancelled
		}
		return fmt.Errorf("failed to build sandbox from %s: %w", b.buildSource, err)
	}

	if b.spec.Image.Type == catalog.TypeBase {
		if err := checkBaseTools(ctx, b.ws.Sandbox); err != nil {
			b.Cleanup(true) //nolint:errcheck
			return err
		}
	}

	// The metadata goes in after the build rather than through a %files section:
	// a sandbox is a host directory, so it can simply be written into.
	if err := copyMetaIntoSandbox(metaDir, b.ws.Sandbox); err != nil {
		b.Cleanup(true) //nolint:errcheck
		return err
	}

	if err := createSquashfs(ctx, b, false, sandboxMountPath, "", preparedPath); err != nil {
		os.Remove(preparedPath) //nolint:errcheck
		b.Cleanup(true)
		if errors.Is(err, context.Canceled) || apptainer.IsBuildCancelled(err) {
			log.Warn("build cancelled, image unchanged", "image", filepath.Base(targetPath))
			return ErrBuildCancelled
		}
		return err
	}
	utils.ShareWithParentGroup(preparedPath)

	installed, err := b.publish(preparedPath, targetPath)
	if err != nil {
		b.Cleanup(true) //nolint:errcheck
		return err
	}

	log.Info("image ready", "kind", "success", "path", installed)
	b.Cleanup(false)
	return nil
}

// synthesizeDefFromURI generates a definition from a scheme:// URI such as
// docker://ubuntu:22.04 — the scheme becomes the bootstrap agent, the rest the
// reference, and the URI is written as #DESC:/#URL: headers so it reaches the
// metadata.
func synthesizeDefFromURI(uri, tmpDir string) (string, error) {
	scheme, from, ok := strings.Cut(uri, "://")
	if !ok || scheme == "" || from == "" {
		return "", fmt.Errorf("not a valid source URI: %s", uri)
	}
	// The directive lines are byte for byte what Apptainer synthesizes for the
	// same URI and stores at /.singularity.d/Singularity — lowercase keys and a
	// trailing blank line — so an image built either way derives one recipe
	// digest. The #-headers are stripped before hashing and only feed metadata.
	content := fmt.Sprintf(
		"#DESC:Built from %s\n"+
			"#URL:%s\n"+
			"# Auto-generated by CondaTainer on %s\n"+
			"bootstrap: %s\n"+
			"from: %s\n\n",
		uri, uri, time.Now().Format("2006-01-02"), scheme, from,
	)
	tmpPath := filepath.Join(tmpDir, "cnt-synth.def")
	if err := os.WriteFile(tmpPath, []byte(content), utils.PermFile); err != nil {
		return "", fmt.Errorf("failed to write synthesized def: %w", err)
	}
	return tmpPath, nil
}

// writePinnedDef returns a definition that builds from cleanDefPath with its
// Bootstrap pinned to pinDigest, so an upstream retagged mid-build cannot leave
// the identity record describing bytes the image does not contain. Only this
// transient copy is rewritten; /.cnt/recipe keeps the definition as written. An
// unchanged definition is built in place.
func writePinnedDef(cleanDefPath, tmpDir, pinDigest string) (string, error) {
	absClean, err := filepath.Abs(cleanDefPath)
	if err != nil {
		return "", err
	}
	data, err := os.ReadFile(absClean)
	if err != nil {
		return "", fmt.Errorf("failed to read def file: %w", err)
	}
	pinned, changed := pinDefinition(data, pinDigest)
	if !changed {
		return absClean, nil
	}
	tmpPath := filepath.Join(tmpDir, "cnt-pinned.def")
	if err := os.WriteFile(tmpPath, pinned, utils.PermFile); err != nil {
		return "", fmt.Errorf("failed to write pinned def: %w", err)
	}
	return tmpPath, nil
}

// copyMetaIntoSandbox places the staged metadata at <sandbox>/.cnt, where the
// pack picks it up as part of the tree. An empty metaDir stages nothing.
func copyMetaIntoSandbox(metaDir, sandbox string) error {
	if metaDir == "" {
		return nil
	}
	dest := filepath.Join(sandbox, meta.DirName)
	if err := utils.MkdirAllShared(dest); err != nil {
		return fmt.Errorf("failed to create %s: %w", dest, err)
	}
	entries, err := os.ReadDir(metaDir)
	if err != nil {
		return fmt.Errorf("failed to list staged metadata %s: %w", metaDir, err)
	}
	for _, entry := range entries {
		if entry.IsDir() {
			continue
		}
		data, err := os.ReadFile(filepath.Join(metaDir, entry.Name()))
		if err != nil {
			return fmt.Errorf("failed to read staged metadata %s: %w", entry.Name(), err)
		}
		if err := os.WriteFile(filepath.Join(dest, entry.Name()), data, utils.PermFile); err != nil {
			return fmt.Errorf("failed to stage %s into the sandbox: %w", entry.Name(), err)
		}
	}
	return nil
}

// baseTools are the executables a base must provide on PATH, because these are
// the ones that run inside it: mksquashfs packs every artifact and micromamba
// builds every conda environment. The rest of what CondaTainer shells out to —
// unsquashfs, debugfs, e2fsck, resize2fs, mke2fs — runs on the host and is no
// concern of the base's.
var baseTools = []string{"mksquashfs", "micromamba"}

// baseToolPaths are wanted at an exact path rather than on PATH. Everything
// CondaTainer runs inside a base it launches as /bin/bash — a recipe as
// `bash -euo pipefail`, every pack and unfreeze script the same way — so a base
// carrying bash only at /usr/bin/bash runs none of them.
var baseToolPaths = []string{"/bin/bash"}

// checkBaseTools refuses a base that cannot serve the builds that will run in
// it. Apptainer is only needed for nested mounting, so its absence is reported
// rather than refused.
//
// The question is put to the container, with the PATH the container sets,
// because that is the PATH the scripts using these tools will have. Looking for
// files in the sandbox would answer a different question: Apptainer ships its
// own mksquashfs under /usr/libexec, which is present and not on PATH.
func checkBaseTools(ctx context.Context, sandbox string) error {
	var out bytes.Buffer
	runErr := execpkg.Run(ctx, execpkg.Options{
		BaseImage:    sandbox,
		ApptainerBin: config.Global.ApptainerBin,
		Command:      []string{"/bin/sh", "-c", baseToolsScript},
		HidePrompt:   true,
	}, execpkg.IO{Stdout: &out, Stderr: &out})

	said := out.String()
	report := readBaseTools(said)
	if len(report.missing) > 0 {
		return fmt.Errorf("this base provides no %s; every build runs inside the base and needs them",
			strings.Join(report.missing, ", "))
	}
	if runErr != nil {
		return fmt.Errorf("could not check the base's tools: %w: %s", runErr, strings.TrimSpace(said))
	}
	if report.noApptainer {
		logging.FromContext(ctx).Warn("base provides no apptainer; builds that mount an image inside the container will fail")
	}
	return nil
}

// baseToolsScript reports what the container cannot find, on its own PATH.
//
// It runs under /bin/sh, not bash: a base with no bash is one of the things it
// exists to report, and a checker that cannot start on it reports nothing.
var baseToolsScript = fmt.Sprintf(`
missing=
for t in %s; do
    command -v "$t" >/dev/null 2>&1 || missing="$missing $t"
done
for p in %s; do
    [ -x "$p" ] || missing="$missing $p"
done
command -v apptainer >/dev/null 2>&1 || command -v singularity >/dev/null 2>&1 || echo NOAPPTAINER
[ -z "$missing" ] || { echo "MISSING:$missing"; exit 1; }
`, strings.Join(baseTools, " "), strings.Join(baseToolPaths, " "))

// baseToolsReport is what the container said about itself.
type baseToolsReport struct {
	missing     []string
	noApptainer bool
}

// readBaseTools reads the script's output. Apptainer writes its own greetings
// and warnings to the same streams, so the markers are matched per line rather
// than the output being read whole.
func readBaseTools(said string) baseToolsReport {
	var report baseToolsReport
	for line := range strings.SplitSeq(said, "\n") {
		line = strings.TrimSpace(line)
		if names, ok := strings.CutPrefix(line, "MISSING:"); ok {
			report.missing = strings.Fields(names)
		}
		if line == "NOAPPTAINER" {
			report.noApptainer = true
		}
	}
	return report
}
