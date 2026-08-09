package build

import (
	"context"
	"fmt"

	"github.com/Justype/condatainer/internal/image/meta"
	"github.com/Justype/condatainer/internal/logging"
	"github.com/Justype/condatainer/internal/utils"
)

// metaDirPath is the host .cnt directory a build stages its manifest into.
// It sits beside the payload, not inside it, so the packer can pass both to
// mksquashfs as separate archive roots.
func metaDirPath(b *BuildObject) string {
	return b.ws.MetaDir
}

// stageMetadata validates this build's manifest, writes it into the workspace,
// and returns the directory the packer should add as a second archive root.
// Staging is the last point an invalid manifest can still stop the build.
func stageMetadata(ctx context.Context, b *BuildObject) (string, error) {
	manifest := b.Manifest()
	if err := meta.Validate(manifest); err != nil {
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
	if err := meta.Stage(dir, manifest); err != nil {
		return "", err
	}

	logging.FromContext(ctx).Debug("staged manifest",
		"name", b.spec.Image.Name, "dir", dir, "type", manifest.Type, "build_type", manifest.BuildType)
	return dir, nil
}
