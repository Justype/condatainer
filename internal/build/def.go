package build

import (
	"context"
	"fmt"
	"os"
	"path/filepath"
	"strings"
	"time"

	"github.com/Justype/condatainer/internal/apptainer"
	"github.com/Justype/condatainer/internal/config"
	"github.com/Justype/condatainer/internal/logging"
	"github.com/Justype/condatainer/internal/utils"
)

// recordedDefPath is where the build definition is embedded inside every
// def-built overlay, so the overlay carries the recipe that produced it.
var recordedDefPath = "/" + utils.BuildScriptDefName

// buildDef implements the Apptainer .def build workflow on BuildObject.
// Workflow:
//  1. Check if overlay already exists (skip if yes)
//  2. Try downloading prebuilt overlay if available (unless SkipPrebuilt)
//  3. Resolve the definition: a real .def file, or one synthesized from a
//     scheme:// source URI (e.g. docker://ubuntu:22.04)
//  4. Embed a copy of that definition at /.cnt-build-script.def for provenance
//  5. Build SIF from the definition using apptainer build --fakeroot
//  6. Extract SquashFS partition from SIF
//  7. Set permissions
//  8. Cleanup
func (b *BuildObject) buildDef(ctx context.Context) error {
	targetPath, finalPath := buildOverlayPaths(b)
	log := logging.FromContext(ctx)

	if skip, err := checkShouldBuild(b); skip || err != nil {
		return err
	}

	if err := b.createBuildLock(); err != nil {
		return err
	}
	defer b.removeBuildLock()

	log.Info("building overlay", "overlay", filepath.Base(targetPath), "mode", buildModeLabel(b), "source", b.buildSource)

	done := watchContext(ctx, "def build")
	defer close(done)

	// Try to download prebuilt overlay first, fall back to building if not available.
	// Only attempt download if the build script source is remote (skipped with --no-prebuilt).
	if b.isRemote && !SkipPrebuilt {
		if writableDir, err := config.GetWritableImagesDir(); err == nil {
			if filepath.Dir(targetPath) == writableDir {
				downloadPath := buildFinalPath(targetPath, b.update)
				if tryDownloadPrebuiltOverlay(ctx, b.nameVersion, downloadPath, b.prebuiltLink) {
					if err := atomicInstall(downloadPath, targetPath, b.update); err != nil {
						return err
					}
					b.Cleanup(false)
					return nil
				}
				if b.update {
					os.Remove(downloadPath) //nolint:errcheck
				}
			}
		}
	}

	// Ensure the tmp directory exists before apptainer tries to write the SIF there.
	if err := utils.EnsureTmpSubdir(b.tmpDir); err != nil {
		return fmt.Errorf("failed to create tmp dir %s: %w", b.tmpDir, err)
	}

	// Resolve the definition source. A scheme:// source (docker://ubuntu:22.04)
	// has no def file, so synthesize one from the URI; a real .def gets its {key}
	// placeholders substituted (Apptainer reads directives like From: verbatim).
	defSource := b.buildSource
	if strings.Contains(b.buildSource, "://") {
		synthPath, err := synthesizeDefFromURI(b.buildSource, b.tmpDir)
		if err != nil {
			return fmt.Errorf("failed to synthesize def from %s: %w", b.buildSource, err)
		}
		defSource = synthPath
	} else if len(b.vars) > 0 {
		subPath, err := substituteTemplateFile(b.buildSource, b.vars, b.tmpDir)
		if err != nil {
			return fmt.Errorf("failed to substitute placeholders in .def file: %w", err)
		}
		defer os.Remove(subPath) //nolint:errcheck
		defSource = subPath
	}

	// Embed a copy of the definition inside the overlay for provenance.
	buildDefSource, err := writeRecordingDef(defSource, b.tmpDir)
	if err != nil {
		return fmt.Errorf("failed to prepare recording def: %w", err)
	}

	log.Info("running apptainer build", "source", b.buildSource)

	buildOpts := &apptainer.BuildOptions{
		Force:     false,
		NoCleanup: false,
	}

	if err := apptainer.Build(ctx, b.tmpOverlayPath, buildDefSource, buildOpts); err != nil {
		b.Cleanup(true)
		if apptainer.IsBuildCancelled(err) {
			log.Info("build cancelled, overlay unchanged", "overlay", filepath.Base(targetPath))
			return ErrBuildCancelled
		}
		return fmt.Errorf("failed to build SIF from %s: %w", b.buildSource, err)
	}

	log.Info("extracting SquashFS", "path", finalPath)

	if err := apptainer.DumpSifToSquashfs(ctx, b.tmpOverlayPath, finalPath); err != nil {
		os.Remove(finalPath) //nolint:errcheck
		b.Cleanup(true)
		if apptainer.IsBuildCancelled(err) {
			log.Info("build cancelled, overlay unchanged", "overlay", filepath.Base(targetPath))
			return ErrBuildCancelled
		}
		return fmt.Errorf("failed to dump SquashFS from SIF: %w", err)
	}

	utils.ShareWithParentGroup(finalPath)

	if err := atomicInstall(finalPath, targetPath, b.update); err != nil {
		return err
	}

	log.Info("overlay ready", "kind", "success", "path", targetPath)
	b.Cleanup(false)
	return nil
}

// tryDownloadPrebuiltOverlay attempts to download a prebuilt overlay from the given prebuiltLink base URL.
func tryDownloadPrebuiltOverlay(ctx context.Context, nameVersion, destPath, prebuiltLink string) bool {
	return tryDownloadPrebuilt(ctx, nameVersion, destPath, "sqf", prebuiltLink, utils.DownloadFile)
}

// synthesizeDefFromURI generates an Apptainer definition file from a scheme://
// source URI (e.g. docker://ubuntu:22.04). The scheme becomes the Bootstrap
// agent and the remainder becomes From. A header comment records the original
// URI and date so the generated def is self-describing once embedded. Returns
// the path to the generated def file in tmpDir.
func synthesizeDefFromURI(uri, tmpDir string) (string, error) {
	scheme, from, ok := strings.Cut(uri, "://")
	if !ok || scheme == "" || from == "" {
		return "", fmt.Errorf("not a valid source URI: %s", uri)
	}
	content := fmt.Sprintf(
		"# Auto-generated by CondaTainer from %s\n"+
			"# date: %s\n"+
			"Bootstrap: %s\n"+
			"From: %s\n",
		uri, time.Now().Format("2006-01-02"), scheme, from,
	)
	tmpPath := filepath.Join(tmpDir, "cnt-synth.def")
	if err := os.WriteFile(tmpPath, []byte(content), utils.PermFile); err != nil {
		return "", fmt.Errorf("failed to write synthesized def: %w", err)
	}
	return tmpPath, nil
}

// writeRecordingDef returns the path to a definition file that builds from
// cleanDefPath and additionally embeds a copy of that definition inside the
// image at recordedDefPath (via an appended %files section), so every def-built
// overlay carries the recipe that produced it. The embedded copy is the clean
// definition — the appended %files section is not part of what gets recorded.
func writeRecordingDef(cleanDefPath, tmpDir string) (string, error) {
	absClean, err := filepath.Abs(cleanDefPath)
	if err != nil {
		return "", err
	}
	data, err := os.ReadFile(absClean)
	if err != nil {
		return "", fmt.Errorf("failed to read def file: %w", err)
	}

	var sb strings.Builder
	sb.Write(data)
	if len(data) > 0 && data[len(data)-1] != '\n' {
		sb.WriteByte('\n')
	}
	fmt.Fprintf(&sb, "\n%%files\n    %s %s\n", absClean, recordedDefPath)

	tmpPath := filepath.Join(tmpDir, "cnt-record.def")
	if err := os.WriteFile(tmpPath, []byte(sb.String()), utils.PermFile); err != nil {
		return "", fmt.Errorf("failed to write recording def: %w", err)
	}
	return tmpPath, nil
}
