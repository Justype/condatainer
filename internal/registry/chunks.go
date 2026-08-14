package registry

import (
	"context"
	"errors"
	"fmt"
	"io"
	"os"
	"path/filepath"
	"time"

	"github.com/opencontainers/go-digest"
	ocispec "github.com/opencontainers/image-spec/specs-go/v1"
	"oras.land/oras-go/v2/errdef"

	"github.com/Justype/condatainer/internal/logging"
)

// chunkSuffix formats what chunk i appends to the artifact's filename. The
// number is the reassembly order, and the width is fixed so the names sort in
// the order they concatenate.
const chunkSuffix = ".part%06d"

// blobPusher is the part of a blob store a chunked push needs. Narrow on
// purpose: a remote repository's Blobs() satisfies it, and so does a test double
// that never opens a socket.
type blobPusher interface {
	Exists(ctx context.Context, target ocispec.Descriptor) (bool, error)
	Push(ctx context.Context, expected ocispec.Descriptor, content io.Reader) error
}

// pushArtifactLayers uploads path as one blob per chunk of layerSize bytes and
// returns the layer descriptors in offset order, which is the order a pull
// concatenates them in. Each carries its filename as the OCI title annotation,
// naming the range for a reader inspecting the manifest; pull takes the order
// from the manifest itself, not from the names.
//
// layerSize comes from [planLayerSize] and is recorded nowhere: a pull reads the
// boundary off the manifest, so artifacts published at earlier sizes stay
// readable and this is free to change per push.
//
// Nothing is staged — a chunk is pushed straight from an [io.SectionReader].
// Writing every chunk out first would cost a second full copy of the artifact,
// 60 GB of scratch for a 60 GB overlay.
func pushArtifactLayers(ctx context.Context, blobs blobPusher, path, mediaType string, layerSize int64) ([]ocispec.Descriptor, error) {
	f, err := os.Open(path)
	if err != nil {
		return nil, err
	}
	defer f.Close() //nolint:errcheck

	info, err := f.Stat()
	if err != nil {
		return nil, err
	}
	name, size := filepath.Base(path), info.Size()

	if layerSize <= 0 || size <= layerSize {
		desc, _, err := pushRange(ctx, blobs, f, 0, size, name, mediaType, true)
		if err != nil {
			return nil, fmt.Errorf("failed to push %s: %w", name, err)
		}
		return []ocispec.Descriptor{desc}, nil
	}

	count := layerCount(size, layerSize)
	layers := make([]ocispec.Descriptor, 0, count)
	// Ask once. Layers upload in order, so a first layer the registry does not
	// have means nothing is committed and every further probe is a request spent
	// to be told 404 — a third of the request budget on a fresh push. When it
	// *is* present, probing continues for every layer, which is the case that
	// pays for it: a resumed push, or a second architecture with identical
	// content.
	probe := true
	guard := throughputGuardFrom(ctx)
	for i := range count {
		offset := int64(i) * layerSize
		bytes := min(layerSize, size-offset)
		// Asked before the layer is started, because the answer stops being
		// useful once the time has been spent.
		if err := guard.check(bytes); err != nil {
			return nil, fmt.Errorf("cannot push chunk %d of %d for %s: %w", i+1, count, name, err)
		}

		started := time.Now()
		desc, present, err := pushRange(ctx, blobs, f, offset, bytes,
			name+fmt.Sprintf(chunkSuffix, i), mediaType, probe, i+1, count)
		if err != nil {
			// An interrupt otherwise leaves the last progress line, "^C" on the
			// end, as the whole account of what happened. Nothing local needs
			// cleaning up, but the layers already accepted stay in the registry,
			// which is what makes running the command again cheap.
			if errors.Is(err, context.Canceled) {
				logging.FromContext(ctx).Warn(verbUpload+" cancelled",
					"layer", fmt.Sprintf("%d/%d", i+1, count),
					"accepted", fmt.Sprintf("%d/%d", i, count))
			}
			return nil, fmt.Errorf("failed to push chunk %d of %d for %s: %w", i+1, count, name, err)
		}
		if !present {
			// Only what was actually transferred. A layer the registry already
			// had returns at once and would otherwise read as infinite bandwidth.
			guard.observe(bytes, time.Since(started))
		}
		if i == 0 {
			probe = present
		}
		layers = append(layers, desc)
	}
	return layers, nil
}

// pushRange uploads one byte range of f as a blob, reporting whether the
// registry already had it.
//
// The range is read twice, once to digest and once to send: a registry needs the
// digest before the upload begins, so the alternative is a staged copy, and
// reading a file twice is cheap next to writing it once more.
//
// probe asks whether the blob is already present. Skipping it is never a
// correctness risk — a registry deduplicates by digest and ErrAlreadyExists is
// tolerated below — it only costs bandwidth in a case the caller has established
// is unlikely.
//
// The upload is retried under [retryPolicy.run], which is what makes the body
// replayable: each attempt opens a new [io.SectionReader], and ORAS streams from
// a reader with no GetBody, so a reused one would send nothing.
func pushRange(ctx context.Context, blobs blobPusher, f *os.File, offset, size int64, name, mediaType string, probe bool, parts ...int) (desc ocispec.Descriptor, present bool, err error) {
	dgst, err := digestRange(f, offset, size)
	if err != nil {
		return ocispec.Descriptor{}, false, err
	}
	desc = ocispec.Descriptor{
		MediaType:   mediaType,
		Digest:      dgst,
		Size:        size,
		Annotations: map[string]string{ocispec.AnnotationTitle: name},
	}

	// A retried push, or a second architecture with an identical chunk, must not
	// send the bytes again.
	if probe {
		if exists, err := blobs.Exists(ctx, desc); err == nil && exists {
			logging.FromContext(ctx).Debug("chunk already present", "chunk", name, "digest", dgst.String())
			return desc, true, nil
		}
	}

	// The same attribute the progress line carries, so a pause reads as part of
	// the transfer it interrupted rather than as an unrelated event.
	attrs := []any{"blob", name}
	if len(parts) == 2 {
		attrs = []any{"layer", fmt.Sprintf("%d/%d", parts[0], parts[1])}
	}
	err = retryPolicyFrom(ctx).run(ctx, mutation{
		verb:  verbUpload,
		attrs: attrs,
		do: func(ctx context.Context) error {
			// Both rebuilt per attempt: the reader because it is consumed, the
			// progress because a retried layer counts from zero again.
			var reader io.Reader = io.NewSectionReader(f, offset, size)
			var progress *progressReader
			// The last layer always reports, however small, because its final
			// record closes the interactive progress line. A 20 GiB artifact in
			// 2 GiB layers ends with 4 KiB, which would otherwise leave
			// "layer=10/11" on screen with the next message printed onto it.
			last := len(parts) == 2 && parts[0] == parts[1]
			if size >= progressMinSize || last {
				if len(parts) == 2 {
					progress = newLayerProgressReader(ctx, reader, size, verbUpload, parts[0], parts[1])
				} else {
					progress = newProgressReader(ctx, reader, size, verbUpload)
				}
				reader = progress
			}
			if err := blobs.Push(ctx, desc, reader); err != nil && !errors.Is(err, errdef.ErrAlreadyExists) {
				return err
			}
			if progress != nil {
				progress.finish()
			}
			return nil
		},
		committed: func(ctx context.Context) (bool, error) { return blobs.Exists(ctx, desc) },
	})
	if err != nil {
		return ocispec.Descriptor{}, false, err
	}
	return desc, false, nil
}

// digestRange computes the digest of a byte range without loading it, reading
// through the same file handle the push will use.
func digestRange(f *os.File, offset, size int64) (digest.Digest, error) {
	digester := digest.Canonical.Digester()
	if _, err := io.Copy(digester.Hash(), io.NewSectionReader(f, offset, size)); err != nil {
		return "", fmt.Errorf("failed to digest bytes %d-%d: %w", offset, offset+size, err)
	}
	return digester.Digest(), nil
}
