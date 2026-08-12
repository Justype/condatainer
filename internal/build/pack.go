package build

import (
	"context"
	"fmt"
	"time"

	"github.com/Justype/condatainer/internal/artifact/meta"
	"github.com/Justype/condatainer/internal/logging"
	"github.com/Justype/condatainer/internal/utils"
)

// stageMetadata validates this build's metadata, writes both documents into the
// workspace, and returns the directory the packer should add as a second archive
// root. Staging is the last point invalid metadata can still stop the build.
func stageMetadata(ctx context.Context, b *BuildObject) (string, error) {
	manifest := b.Manifest()
	// Stamped here, not in Manifest(), which must stay a projection of Spec.
	// Staging is when the payload is final, so this is the build's finish time.
	manifest.Build.Created = time.Now().UTC()
	if err := meta.ValidateManifest(manifest); err != nil {
		return "", fmt.Errorf("refusing to pack %s: %w", b.spec.Image.Name, err)
	}
	rt := b.Runtime()
	if err := meta.ValidateRuntime(rt); err != nil {
		return "", fmt.Errorf("refusing to pack %s: %w", b.spec.Image.Name, err)
	}

	// In ext3 mode the payload is inside the temporary image, so nothing has
	// created the host build directory yet.
	if b.ws.Root != "" {
		if err := utils.EnsureTmpSubdir(b.ws.Root); err != nil {
			return "", fmt.Errorf("failed to create tmp dir %s: %w", b.ws.Root, err)
		}
	}

	dir := b.ws.MetaDir
	if err := meta.StageRuntime(dir, rt); err != nil {
		return "", err
	}
	if err := meta.StageManifest(dir, manifest); err != nil {
		return "", err
	}
	// The sources go in verbatim: a recipe as it was fetched, tokens and all, and
	// a Conda build's exports as captured. manifest.source.files names exactly
	// these, so a reader never has to probe for what is there.
	for _, file := range b.embedded {
		if err := meta.StageBytes(dir, file.Name, file.Data); err != nil {
			return "", err
		}
	}

	logging.FromContext(ctx).Debug("staged metadata",
		"name", b.spec.Image.Name, "dir", dir, "type", manifest.Type, "build_type", manifest.BuildType)
	return dir, nil
}
