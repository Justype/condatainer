package sif

import (
	"context"
	"fmt"
	"io"
	"os"

	"github.com/Justype/condatainer/internal/image/squashfs"
	"github.com/Justype/condatainer/internal/image/tool"
	"github.com/Justype/condatainer/internal/utils"
)

// ReadFile reads one file out of a SIF's SquashFS partition without extracting
// the partition: the partition offset found in the descriptor table is handed
// to unsquashfs -offset, which reads the archive in place.
//
// Errors keep their cause, so a caller can tell a missing file from a missing
// tool from a corrupt image.
func ReadFile(path, innerPath string) ([]byte, error) {
	part, err := PrimarySystemPartition(path)
	if err != nil {
		return nil, err
	}
	return squashfs.CatFile(path, innerPath, part.Offset)
}

// ExtractPartition copies a SIF's primary SquashFS partition out to a standalone
// .sqf. The partition is already a complete SquashFS archive, so this is a byte
// copy of one range — no repacking, no recompression.
//
// Byte-identical to `apptainer sif dump`, but does not depend on apptainer being
// installed or on the wording of its table, and is several times faster.
func ExtractPartition(ctx context.Context, sifPath, outPath string) error {
	part, err := PrimarySystemPartition(sifPath)
	if err != nil {
		return err
	}

	in, err := os.Open(sifPath)
	if err != nil {
		return fmt.Errorf("%w: %s: %w", tool.ErrUnreadable, sifPath, err)
	}
	defer in.Close()

	out, err := utils.CreateFileWritable(outPath)
	if err != nil {
		return fmt.Errorf("failed to create %s: %w", outPath, err)
	}
	// A partial .sqf is worse than none: it looks like an image and fails later.
	success := false
	defer func() {
		out.Close() //nolint:errcheck
		if !success {
			os.Remove(outPath) //nolint:errcheck
		}
	}()

	src := io.NewSectionReader(in, part.Offset, part.Size)
	if err := copyCancellable(ctx, out, src); err != nil {
		return err
	}
	if err := out.Sync(); err != nil {
		return fmt.Errorf("failed to flush %s: %w", outPath, err)
	}

	success = true
	utils.ShareWithParentGroup(outPath)
	return nil
}

// copyCancellable copies src to dst in chunks, checking ctx between them, so a
// cancelled build stops during a multi-gigabyte partition copy instead of after.
func copyCancellable(ctx context.Context, dst io.Writer, src io.Reader) error {
	const chunk = 4 << 20
	for {
		if err := ctx.Err(); err != nil {
			return err
		}
		n, err := io.CopyN(dst, src, chunk)
		if err == io.EOF {
			return nil
		}
		if err != nil {
			return fmt.Errorf("failed to copy SquashFS partition: %w", err)
		}
		if n == 0 {
			return nil
		}
	}
}
