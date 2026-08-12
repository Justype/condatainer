package build

import (
	"context"
	"errors"
	"fmt"
	"os"
	"path/filepath"
	"strings"
	"time"

	"github.com/Justype/condatainer/catalog"
	"github.com/Justype/condatainer/internal/artifact/meta"
	"github.com/Justype/condatainer/internal/image/sif"
	"github.com/Justype/condatainer/internal/logging"
	"github.com/Justype/condatainer/internal/runtime/apptainer"
	"github.com/Justype/condatainer/internal/utils"
)

// buildDef builds an image from an Apptainer definition, real or synthesized
// from a scheme:// URI. A base keeps the .sif Apptainer produced; every other
// type is extracted to .sqf, metadata and all.
func (b *BuildObject) buildDef(ctx context.Context) error {
	targetPath := b.tgt.Path
	log := logging.FromContext(ctx)
	isBase := b.spec.Image.Type == catalog.TypeBase

	if skip, err := checkShouldBuild(b); skip || err != nil {
		return err
	}

	if err := b.createBuildLock(); err != nil {
		return err
	}
	defer b.removeBuildLock()
	preparedPath := b.tgt.Prepared
	if pulled, err := b.tryPrebuilt(ctx); err != nil || pulled {
		return err
	}

	log.Info("building image", "image", filepath.Base(targetPath), "mode", buildModeLabel(b), "source", b.buildSource)

	done := watchContext(ctx, "def build")
	defer close(done)

	// Ensure the tmp directory exists before apptainer tries to write the SIF there.
	if err := utils.EnsureTmpSubdir(b.ws.Root); err != nil {
		return fmt.Errorf("failed to create tmp dir %s: %w", b.ws.Root, err)
	}

	// Resolve the definition source. A scheme:// source (docker://ubuntu:22.04)
	// has no def file, so synthesize one from the URI; a real .def gets its {key}
	// placeholders substituted (Apptainer reads directives like From: verbatim).
	defSource := b.buildSource
	if strings.Contains(b.buildSource, "://") {
		synthPath, err := synthesizeDefFromURI(b.buildSource, b.ws.Root)
		if err != nil {
			return fmt.Errorf("failed to synthesize def from %s: %w", b.buildSource, err)
		}
		defSource = synthPath
	}

	defData, err := os.ReadFile(defSource)
	if err != nil {
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

	// Embed the definition and staged metadata, pinning the bootstrap.
	buildDefSource, err := writeRecordingDef(defSource, metaDir, b.ws.Root, b.spec.Source.UpstreamDigest())
	if err != nil {
		return fmt.Errorf("failed to prepare recording def: %w", err)
	}

	log.Info("running apptainer build", "source", b.buildSource)

	buildOpts := &apptainer.BuildOptions{
		Force:     false,
		NoCleanup: false,
	}

	if err := apptainer.Build(ctx, b.ws.Overlay, buildDefSource, buildOpts); err != nil {
		b.Cleanup(true)
		if apptainer.IsBuildCancelled(err) {
			log.Info("build cancelled, image unchanged", "image", filepath.Base(targetPath))
			return ErrBuildCancelled
		}
		return fmt.Errorf("failed to build SIF from %s: %w", b.buildSource, err)
	}

	if isBase {
		// The container root is executed directly, so the SIF is the product. The
		// build tmp root and the images directory are routinely separate mounts,
		// so this cannot be a bare rename.
		copied, err := utils.MoveFile(ctx, b.ws.Overlay, preparedPath, false)
		if err != nil {
			os.Remove(preparedPath) //nolint:errcheck
			b.Cleanup(true)
			if errors.Is(err, context.Canceled) {
				log.Info("build cancelled, image unchanged", "image", filepath.Base(targetPath))
				return ErrBuildCancelled
			}
			return fmt.Errorf("failed to move SIF to %s: %w", preparedPath, err)
		}
		if copied {
			log.Debug("copied the SIF across filesystems", "from", b.ws.Overlay, "to", preparedPath)
		}
		if err := utils.MakeExecutable(preparedPath); err != nil {
			log.Debug("failed to set permissions", "path", preparedPath, "err", err)
		}
	} else {
		log.Info("extracting SquashFS", "path", preparedPath)
		if err := sif.ExtractPartition(ctx, b.ws.Overlay, preparedPath); err != nil {
			os.Remove(preparedPath) //nolint:errcheck
			b.Cleanup(true)
			if errors.Is(err, context.Canceled) || apptainer.IsBuildCancelled(err) {
				log.Info("build cancelled, image unchanged", "image", filepath.Base(targetPath))
				return ErrBuildCancelled
			}
			return fmt.Errorf("failed to extract SquashFS from SIF: %w", err)
		}
		utils.ShareWithParentGroup(preparedPath)
	}

	if err := atomicInstall(preparedPath, targetPath); err != nil {
		return err
	}

	log.Info("image ready", "kind", "success", "path", targetPath)
	b.Cleanup(false)
	return nil
}

// synthesizeDefFromURI generates a definition from a scheme:// URI such as
// docker://ubuntu:22.04 — the scheme becomes Bootstrap, the rest becomes From,
// and the URI is written as #DESC:/#URL: headers so it reaches the metadata.
func synthesizeDefFromURI(uri, tmpDir string) (string, error) {
	scheme, from, ok := strings.Cut(uri, "://")
	if !ok || scheme == "" || from == "" {
		return "", fmt.Errorf("not a valid source URI: %s", uri)
	}
	content := fmt.Sprintf(
		"#DESC:Built from %s\n"+
			"#URL:%s\n"+
			"# Auto-generated by CondaTainer on %s\n"+
			"Bootstrap: %s\n"+
			"From: %s\n",
		uri, uri, time.Now().Format("2006-01-02"), scheme, from,
	)
	tmpPath := filepath.Join(tmpDir, "cnt-synth.def")
	if err := os.WriteFile(tmpPath, []byte(content), utils.PermFile); err != nil {
		return "", fmt.Errorf("failed to write synthesized def: %w", err)
	}
	return tmpPath, nil
}

// writeRecordingDef returns a definition that builds from cleanDefPath and
// embeds every file staged in metaDir under /.cnt, via an appended %files
// section. An empty metaDir builds the definition as-is.
//
// pinDigest rewrites the From: line to name that digest, so an upstream retagged
// mid-build cannot leave the identity record describing bytes the image does not
// contain. Only this transient copy is rewritten; /.cnt/recipe keeps the
// definition as written.
//
// The files are listed one by one rather than by naming the directory: %files
// copies with cp -a, which would nest the whole directory inside a /.cnt the
// base already has.
func writeRecordingDef(cleanDefPath, metaDir, tmpDir, pinDigest string) (string, error) {
	absClean, err := filepath.Abs(cleanDefPath)
	if err != nil {
		return "", err
	}
	data, err := os.ReadFile(absClean)
	if err != nil {
		return "", fmt.Errorf("failed to read def file: %w", err)
	}
	pinned, changed := pinDefinition(data, pinDigest)
	if metaDir == "" && !changed {
		return absClean, nil
	}
	data = pinned
	if metaDir == "" {
		tmpPath := filepath.Join(tmpDir, "cnt-pinned.def")
		if err := os.WriteFile(tmpPath, data, utils.PermFile); err != nil {
			return "", fmt.Errorf("failed to write pinned def: %w", err)
		}
		return tmpPath, nil
	}
	absMeta, err := filepath.Abs(metaDir)
	if err != nil {
		return "", err
	}
	staged, err := os.ReadDir(absMeta)
	if err != nil {
		return "", fmt.Errorf("failed to list staged metadata %s: %w", absMeta, err)
	}

	var sb strings.Builder
	sb.Write(data)
	if len(data) > 0 && data[len(data)-1] != '\n' {
		sb.WriteByte('\n')
	}
	sb.WriteString("\n%files\n")
	for _, entry := range staged {
		if entry.IsDir() {
			continue
		}
		fmt.Fprintf(&sb, "    %s /%s/%s\n", filepath.Join(absMeta, entry.Name()), meta.DirName, entry.Name())
	}

	tmpPath := filepath.Join(tmpDir, "cnt-metadata.def")
	if err := os.WriteFile(tmpPath, []byte(sb.String()), utils.PermFile); err != nil {
		return "", fmt.Errorf("failed to write recording def: %w", err)
	}
	return tmpPath, nil
}
