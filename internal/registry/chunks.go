package registry

import (
	"context"
	"errors"
	"fmt"
	"io"
	"os"
	"path/filepath"
	"strings"

	"github.com/opencontainers/go-digest"
	ocispec "github.com/opencontainers/image-spec/specs-go/v1"
	"oras.land/oras-go/v2/errdef"

	"github.com/Justype/condatainer/internal/logging"
)

// artifactChunkSize bounds one blob upload, so an interrupted push loses a chunk
// rather than the whole artifact. An artifact at or under it stays a single
// layer. A variable, not a constant, so tests can chunk a few bytes.
var artifactChunkSize int64 = 512 << 20

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

// pushArtifactLayers uploads path as one blob per chunk and returns the layer
// descriptors in offset order, which is the order [assemblePulledArtifact]
// reassembles by. Each descriptor carries its filename as the OCI title
// annotation, which is what lets a pull's file store write the parts back under
// names it can order.
//
// Nothing is staged. A chunk is a byte range of a file that already exists, so
// it is pushed straight from an [io.SectionReader]; the alternative — writing
// every chunk to a temporary directory first — costs a second full copy of the
// artifact, which for a 60 GB overlay means 60 GB of scratch, usually on a
// node-local disk the user never named.
func pushArtifactLayers(ctx context.Context, blobs blobPusher, path, mediaType string) ([]ocispec.Descriptor, error) {
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

	if artifactChunkSize <= 0 || size <= artifactChunkSize {
		desc, err := pushRange(ctx, blobs, f, 0, size, name, mediaType)
		if err != nil {
			return nil, fmt.Errorf("failed to push %s: %w", name, err)
		}
		return []ocispec.Descriptor{desc}, nil
	}

	count := int((size + artifactChunkSize - 1) / artifactChunkSize)
	layers := make([]ocispec.Descriptor, 0, count)
	for i := range count {
		offset := int64(i) * artifactChunkSize
		desc, err := pushRange(ctx, blobs, f, offset, min(artifactChunkSize, size-offset),
			name+fmt.Sprintf(chunkSuffix, i), mediaType, i+1, count)
		if err != nil {
			return nil, fmt.Errorf("failed to push chunk %d of %d for %s: %w", i+1, count, name, err)
		}
		layers = append(layers, desc)
	}
	return layers, nil
}

// pushRange uploads one byte range of f as a blob.
//
// The range is read twice: once to digest it, once to send it. A registry needs
// the digest before the upload begins, so the alternative is not one pass but a
// staged copy — and reading a file twice is cheap next to writing it once more.
func pushRange(ctx context.Context, blobs blobPusher, f *os.File, offset, size int64, name, mediaType string, parts ...int) (ocispec.Descriptor, error) {
	dgst, err := digestRange(f, offset, size)
	if err != nil {
		return ocispec.Descriptor{}, err
	}
	desc := ocispec.Descriptor{
		MediaType:   mediaType,
		Digest:      dgst,
		Size:        size,
		Annotations: map[string]string{ocispec.AnnotationTitle: name},
	}

	// A retried push, or a second architecture with an identical chunk, must not
	// send the bytes again.
	if exists, err := blobs.Exists(ctx, desc); err == nil && exists {
		logging.FromContext(ctx).Debug("chunk already present", "chunk", name, "digest", dgst.String())
		return desc, nil
	}

	var reader io.Reader = io.NewSectionReader(f, offset, size)
	var progress *progressReader
	if size >= progressMinSize {
		if len(parts) == 2 {
			progress = newLayerProgressReader(ctx, reader, size, verbUpload, parts[0], parts[1])
		} else {
			progress = newProgressReader(ctx, reader, size, verbUpload)
		}
		reader = progress
	}
	if err := blobs.Push(ctx, desc, reader); err != nil && !errors.Is(err, errdef.ErrAlreadyExists) {
		return ocispec.Descriptor{}, err
	}
	if progress != nil {
		progress.finish()
	}
	return desc, nil
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

// assemblePulledArtifact returns the path of the complete artifact in dir: the
// single downloaded layer as it stands, or the ordered chunks concatenated into
// one file.
//
// orderedNames comes from the manifest's layer order, and the names must be the
// unbroken .partNNNNNN sequence that produced them. Checking that is what stops a
// registry from serving a plausible-looking artifact with a chunk missing from
// the middle, which would otherwise assemble into a corrupt file.
func assemblePulledArtifact(dir string, orderedNames []string) (string, error) {
	entries, err := os.ReadDir(dir)
	if err != nil {
		return "", fmt.Errorf("failed to read pulled artifact: %w", err)
	}
	files := make(map[string]string)
	for _, e := range entries {
		if e.IsDir() {
			continue
		}
		files[e.Name()] = filepath.Join(dir, e.Name())
	}
	if len(files) == 0 {
		return "", fmt.Errorf("pulled artifact contained no file")
	}
	if len(files) != len(orderedNames) {
		return "", fmt.Errorf("pulled artifact contained %d files, manifest declares %d layers",
			len(files), len(orderedNames))
	}

	paths := make([]string, 0, len(orderedNames))
	for _, name := range orderedNames {
		path, ok := files[name]
		if !ok {
			return "", fmt.Errorf("manifest layer %q is missing from pulled artifact", name)
		}
		paths = append(paths, path)
	}
	if len(paths) == 1 {
		return paths[0], nil
	}

	prefix := strings.TrimSuffix(orderedNames[0], fmt.Sprintf(chunkSuffix, 0))
	if prefix == orderedNames[0] {
		return "", fmt.Errorf("invalid first chunk name %q", orderedNames[0])
	}
	for i, name := range orderedNames {
		if name != prefix+fmt.Sprintf(chunkSuffix, i) {
			return "", fmt.Errorf("invalid or unordered chunk name %q", name)
		}
	}

	assembled := filepath.Join(dir, "artifact.assembled")
	dst, err := os.OpenFile(assembled, os.O_CREATE|os.O_EXCL|os.O_WRONLY, 0o600)
	if err != nil {
		return "", err
	}
	// A half-written artifact must not survive to be installed.
	ok := false
	defer func() {
		dst.Close() //nolint:errcheck
		if !ok {
			os.Remove(assembled) //nolint:errcheck
		}
	}()
	for _, path := range paths {
		src, err := os.Open(path)
		if err != nil {
			return "", err
		}
		_, copyErr := io.Copy(dst, src)
		closeErr := src.Close()
		if copyErr != nil {
			return "", copyErr
		}
		if closeErr != nil {
			return "", closeErr
		}
	}
	if err := dst.Close(); err != nil {
		return "", err
	}
	ok = true
	return assembled, nil
}
