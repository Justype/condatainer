package build

import (
	"context"
	"fmt"
	"path/filepath"

	"github.com/Justype/condatainer/internal/artifact/key"
	"github.com/Justype/condatainer/internal/artifact/meta"
)

// stagePayloadKey hashes the payload the packer is about to read and records the
// key in the staged manifest. sourceDir is the host directory being packed.
//
// Only a recipe's build is keyed: a Conda environment is pinned by its explicit
// export and a definition by its upstream digest, so neither needs it. A build
// that cannot key its payload stops here — an artifact without the key is not
// built.
func (b *BuildObject) stagePayloadKey(ctx context.Context, sourceDir, metaDir string) error {
	if b.buildType != BuildTypeScript {
		return nil
	}
	// mksquashfs names an archive root after its source's basename.
	records, err := key.TreeOf(ctx, sourceDir, "/"+filepath.Base(sourceDir), b.effectiveNcpus())
	if err != nil {
		return fmt.Errorf("cannot key the payload of %s: %w", b.spec.Image.Name, err)
	}
	ref, err := key.PayloadKey(records)
	if err != nil {
		return fmt.Errorf("cannot key the payload of %s: %w", b.spec.Image.Name, err)
	}
	b.keys.Payload = ref
	return meta.StagePayloadKey(metaDir, ref)
}
