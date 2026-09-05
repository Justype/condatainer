package registry

import (
	"context"
	"errors"
	"fmt"
	"io"
	"os"
	"path/filepath"
	"strings"

	ocispec "github.com/opencontainers/image-spec/specs-go/v1"
	"golang.org/x/sync/errgroup"
	"golang.org/x/sys/unix"
	"oras.land/oras-go/v2/content"
	"oras.land/oras-go/v2/registry"
	"oras.land/oras-go/v2/registry/remote"

	"github.com/Justype/condatainer/catalog"
	"github.com/Justype/condatainer/internal/artifact/compare"
	"github.com/Justype/condatainer/internal/artifact/meta"
	"github.com/Justype/condatainer/internal/image"
	"github.com/Justype/condatainer/internal/image/producer"
	"github.com/Justype/condatainer/internal/logging"
	"github.com/Justype/condatainer/internal/utils"
)

// Pull downloads the artifact at desc and installs it at destPath.
//
// desc and annotations come from [ResolveArtifact], and whether this is the
// wanted artifact is the caller's question, asked with [Check] before spending
// the bytes. What Pull adds is everything that can only be known once the
// payload is here: that it is the kind of image destPath asked for, that its own
// bytes regenerate the keys the annotations advertised, and that installing it
// displaces nothing.
//
// The install is flat and looks exactly like a local build's: the same
// directory, the same permissions, the same locking and protection rules. It
// creates no store entry — that is acquisition by identity, a different
// operation.
func Pull(ctx context.Context, base, repo string, desc ocispec.Descriptor, annotations map[string]string, destPath string, kind Kind) error {
	guard, err := producer.AcquireLocal(destPath)
	if err != nil {
		return err
	}
	defer guard.Release() //nolint:errcheck
	return pull(ctx, base, repo, desc, annotations, destPath, kind)
}

// PullLocked performs Pull while the caller already holds destPath's producer
// lock. Build uses this after deriving the expected equivalence under its lock;
// callers without a lock use Pull.
func PullLocked(ctx context.Context, base, repo string, desc ocispec.Descriptor, annotations map[string]string, destPath string, kind Kind) error {
	return pull(ctx, base, repo, desc, annotations, destPath, kind)
}

func pull(ctx context.Context, base, repo string, desc ocispec.Descriptor, annotations map[string]string, destPath string, kind Kind) error {
	log := logging.FromContext(ctx)

	if err := checkDistributable(destPath); err != nil {
		return err
	}
	wantArtifactType, wantLayerType, err := kind.types()
	if err != nil {
		return err
	}
	repository, err := newRepository(base, repo)
	if err != nil {
		return err
	}
	manifest, err := fetchManifest(ctx, repository, desc)
	if err != nil {
		return classify(err)
	}
	if err := validatePayloadContract(manifest, wantArtifactType, wantLayerType); err != nil {
		return err
	}
	// Protection is a stable fact — the write bit — so it is worth knowing before
	// a multi-gigabyte download rather than after. A conflicting flock is not:
	// it may well be gone by the time the transfer finishes, so it is left to the
	// probe below and not treated as a reason to refuse now.
	destDir := filepath.Dir(destPath)
	if utils.FileExists(destPath) {
		if err := image.CheckAvailable(destPath, true); errors.Is(err, image.ErrProtected) {
			return err
		}
	}
	if err := utils.MkdirAllShared(destDir); err != nil {
		return fmt.Errorf("cannot prepare %s: %w", destDir, err)
	}
	if err := requireFreeSpace(destDir, downloadSize(desc, manifest)); err != nil {
		return err
	}

	// Staged in the destination directory, not a temp root: rename(2) cannot
	// cross filesystems, and an install that is not one rename is an install
	// that can be seen half-done.
	stageDir, err := os.MkdirTemp(destDir, ".cnt-pull-")
	if err != nil {
		return fmt.Errorf("cannot create staging directory in %s: %w", destDir, err)
	}
	defer os.RemoveAll(stageDir) //nolint:errcheck

	staged := filepath.Join(stageDir, filepath.Base(destPath))
	if err := download(ctx, repository, manifest, staged); err != nil {
		return err
	}
	embedded, err := verifyPayload(staged, annotations)
	if err != nil {
		return err
	}
	if name := annotations[AnnTitle]; name != "" && name != embedded.Name {
		// Reported, never fatal: identity already answered whether these are the
		// wanted bytes, and what an image answers to is decided by where its file
		// sits, not by what it calls itself.
		log.Warn("published name disagrees with the artifact's own",
			"published", name, "embedded", embedded.Name)
	}

	// Immediately before the rename, and not before the download: an exec on
	// another node may hold LOCK_SH and be reading the old file lazily over NFS,
	// where replacing it stales the handle mid-job. A lock taken earlier would
	// block every exec for the length of the transfer and would then be held on
	// an orphaned inode anyway.
	if utils.FileExists(destPath) {
		if err := image.CheckAvailable(destPath, true); err != nil {
			return fmt.Errorf("cannot replace %s: %w", utils.StylePath(destPath), err)
		}
	}
	if err := os.Rename(staged, destPath); err != nil {
		return fmt.Errorf("failed to install pulled artifact: %w", err)
	}
	utils.ShareWithParentGroup(destPath)
	log.Info("installed", "artifact", utils.StylePath(destPath))
	return nil
}

// Fetch downloads the artifact desc addresses to destPath and verifies it,
// installing nothing.
//
// This is Pull with every step that makes an artifact *answer to a name*
// removed: no producer lock, no free-name search, no protection check, no
// rename. The caller already owns destPath and decides what happens to it, which
// is what lets a locked restore hand a store transaction's staging path straight
// to the transport and then publish it under the identity the lock pinned.
//
// destPath must not exist. kind says what is expected of the payload, because a
// staging name cannot: it ends in the producer's suffix rather than in .sqf.
//
// The two checks that survive are the two that are about the bytes: the manifest
// must describe the kind of artifact asked for before anything is transferred,
// and the assembled payload must regenerate the keys its annotations advertised
// before the caller is told it succeeded.
func Fetch(ctx context.Context, base, repo string, desc ocispec.Descriptor, annotations map[string]string, destPath string, kind Kind) error {
	wantArtifactType, wantLayerType, err := kind.types()
	if err != nil {
		return err
	}
	repository, err := newRepository(base, repo)
	if err != nil {
		return err
	}
	manifest, err := fetchManifest(ctx, repository, desc)
	if err != nil {
		return classify(err)
	}
	if err := validatePayloadContract(manifest, wantArtifactType, wantLayerType); err != nil {
		return err
	}
	if err := requireFreeSpace(filepath.Dir(destPath), downloadSize(desc, manifest)); err != nil {
		return err
	}
	if err := download(ctx, repository, manifest, destPath); err != nil {
		return err
	}
	if _, err := verifyPayload(destPath, annotations); err != nil {
		return err
	}
	return nil
}

// SplitCoordinate separates a complete OCI repository coordinate into the
// registry base and the repository path beneath it.
//
// A lock records one string — "ghcr.io/org/cnt/grch38/genome" — while the
// transport wants the two halves apart. The split is at the first slash and
// nowhere cleverer: everything after the host is repository path, whether the
// publisher nested it by name or flattened it into one repository.
func SplitCoordinate(coordinate string) (base, repo string, err error) {
	// Trailing slashes only: TrimBaseScheme drops those. A leading one is not
	// tidied away, because "/lab/p" names no registry and reading it as host
	// "lab" would turn a malformed coordinate into a plausible one.
	host, path, found := strings.Cut(TrimBaseScheme(coordinate), "/")
	if !found || host == "" || path == "" {
		return "", "", fmt.Errorf("%q is not a registry/repository coordinate", coordinate)
	}
	return host, path, nil
}

// Kind is the sort of image an operation expects to move: a container root, or
// something mounted onto one. Every artifact is a `.sqf`, so a filename cannot
// say which, and the manifest type is what does.
//
// Push reads it from the artifact it is publishing. Pull and Fetch are told,
// because neither can read a manifest it has not fetched yet — and being told is
// what gives the check its teeth: a caller that knows an overlay was locked here
// learns at the transport, rather than at mount, that a root was served.
type Kind string

const (
	KindOverlay Kind = "overlay" // mounted onto a root
	KindBase    Kind = "base"    // is the root
)

// KindFor reports the kind an artifact of this type travels as.
func KindFor(typ catalog.Type) Kind {
	if typ == catalog.TypeBase {
		return KindBase
	}
	return KindOverlay
}

// checkDistributable rejects a path that names something no registry moves. It
// says what a manifest cannot: a writable overlay has no identity at all, and a
// file that is not an image was never a candidate.
func checkDistributable(path string) error {
	switch {
	case utils.IsSqf(path):
		return nil
	case utils.IsImg(path):
		return fmt.Errorf("%s is a writable overlay, which has no identity and is never distributed", path)
	}
	return fmt.Errorf("%s is not a distributable image (.sqf)", path)
}

// types reports the artifact and layer media types a kind travels as.
func (k Kind) types() (artifactType, layerType string, err error) {
	switch k {
	case KindOverlay:
		return ArtifactTypeOverlay, MediaTypeOverlayBlob, nil
	case KindBase:
		return ArtifactTypeBase, MediaTypeBaseBlob, nil
	}
	return "", "", fmt.Errorf("unknown image kind %q", k)
}

// imageTypes reports the media types a local artifact travels as, from its own
// embedded manifest.
func imageTypes(path string) (artifactType, layerType string, err error) {
	if err := checkDistributable(path); err != nil {
		return "", "", err
	}
	m, err := meta.ReadManifest(path)
	if err != nil {
		return "", "", fmt.Errorf("%s: %w", path, err)
	}
	return KindFor(m.Type).types()
}

// validatePayloadContract rejects a manifest that does not describe the kind of
// artifact being installed, before any of its bytes are fetched.
func validatePayloadContract(manifest ocispec.Manifest, wantArtifactType, wantLayerType string) error {
	if manifest.ArtifactType != wantArtifactType {
		return fmt.Errorf("%w: it is a %q, not a %q",
			ErrInvalidArtifact, manifest.ArtifactType, wantArtifactType)
	}
	if len(manifest.Layers) == 0 {
		return fmt.Errorf("%w: it carries no payload", ErrInvalidArtifact)
	}
	for i, layer := range manifest.Layers {
		if layer.MediaType != wantLayerType {
			return fmt.Errorf("%w: payload %d is %q, not %q",
				ErrInvalidArtifact, i, layer.MediaType, wantLayerType)
		}
	}
	return nil
}

// download copies the resolved manifest and its blobs into stageDir. It is
// addressed by digest on both sides, so no tag is resolved a second time and the
// artifact verified above is the one that arrives.
func download(ctx context.Context, repository *remote.Repository, manifest ocispec.Manifest, destPath string) error {
	var total int64
	offsets := make([]int64, len(manifest.Layers))
	for i, layer := range manifest.Layers {
		offsets[i] = total
		total += positiveSize(layer.Size)
	}

	f, err := os.OpenFile(destPath, os.O_CREATE|os.O_EXCL|os.O_WRONLY, 0o600)
	if err != nil {
		return fmt.Errorf("cannot create %s: %w", destPath, err)
	}
	// A half-written artifact must not survive to be installed. The staging
	// directory is removed by the caller either way; closing early is what makes
	// the error path deterministic rather than dependent on the defer order.
	closed := false
	defer func() {
		if !closed {
			f.Close() //nolint:errcheck
		}
	}()
	// Sized up front so every layer writes into a range that already exists,
	// which is what lets them be written in any order.
	if err := f.Truncate(total); err != nil {
		return fmt.Errorf("cannot size %s: %w", destPath, err)
	}

	progress := newDownloadProgress(ctx, total)
	blobs := repository.Blobs()
	group, groupCtx := errgroup.WithContext(ctx)
	group.SetLimit(pullConcurrency)
	for i, layer := range manifest.Layers {
		i, layer := i, layer
		group.Go(func() error {
			if err := fetchLayerAt(groupCtx, blobs, layer, f, offsets[i], progress); err != nil {
				return fmt.Errorf("failed to download payload %d of %d: %w",
					i+1, len(manifest.Layers), classify(err))
			}
			return nil
		})
	}
	if err := group.Wait(); err != nil {
		if errors.Is(err, context.Canceled) {
			progress.interrupted()
		}
		return err
	}
	progress.finish()

	closed = true
	if err := f.Close(); err != nil {
		return fmt.Errorf("cannot finish writing %s: %w", destPath, err)
	}
	return nil
}

// fetchLayerAt streams one layer straight into its own range of the artifact.
//
// Writing at an offset removes the second copy: layers are a plain concatenation
// in manifest order, so each one's place is known before anything is fetched.
// Staging them separately and joining them afterwards cost a duplicate of the
// whole artifact — 80 GB of scratch for a 40 GB image — plus a full read-write
// pass. It also makes a retry cheap: the range is fixed, so a second attempt
// overwrites what the first one wrote.
func fetchLayerAt(ctx context.Context, blobs registry.BlobStore, layer ocispec.Descriptor, f *os.File, offset int64, progress *downloadProgress) error {
	return retryPolicyFrom(ctx).run(ctx, mutation{
		verb:  verbDownload,
		attrs: []any{"blob", layer.Digest.String()},
		do: func(ctx context.Context) error {
			reader, err := blobs.Fetch(ctx, layer)
			if err != nil {
				return err
			}
			defer reader.Close() //nolint:errcheck

			// ORAS checks Content-Length and the Docker-Content-Digest header on
			// a fetch but never hashes the body — oras.Copy wraps the verifier,
			// and this path does not use it. Without this a registry could serve
			// anything of the right length.
			verifier := content.NewVerifyReader(reader, layer)
			// Rebuilt per attempt: a retried layer rewrites its range from the
			// start, so an abandoned attempt's bytes must leave the total.
			counted := progress.attempt(verifier)
			if _, err := io.Copy(io.NewOffsetWriter(f, offset), counted); err != nil {
				return err
			}
			return verifier.Verify()
		},
	})
}

// pullConcurrency is how many blobs a pull fetches at once.
//
// Stated rather than inherited from ORAS's identical default: a pull runs once
// per node per artifact, so across a cluster it is the larger of the two
// directions and worth deciding on purpose. Three, because a pull is on the path
// of every job needing the artifact and a rate limit now pauses it rather than
// failing it.
//
// It costs nothing in scratch: every layer writes to its own range, so the peak
// is the artifact's size whatever the order.
const pullConcurrency = 3

// verifyPayload holds the downloaded artifact to what its annotations promised,
// returning what the payload says about itself.
//
// The keys are regenerated from the payload's own files by [compare.Read], never
// read off its manifest, so this catches a publisher whose annotations and bytes
// disagree — which is the whole reason to check after the transfer as well as
// before it. What it does not prove is that the payload matches the keys: two
// images with identical keys and different contents are easy to produce, and
// that needs signatures rather than hashes.
func verifyPayload(path string, annotations map[string]string) (compare.Artifact, error) {
	got, err := compare.Read(path)
	if err != nil {
		return got, fmt.Errorf("%w: %w", ErrInvalidArtifact, err)
	}
	return got, checkRegeneratedKeys(got, Identity(annotations), Equiv(annotations), "published")
}

// checkRegeneratedKeys holds the keys an artifact's own files reproduce against
// the keys something claimed for it. claimant names where that claim came from,
// since both directions use this: push checks the artifact against its own
// recorded manifest, pull checks it against the annotations it was advertised
// under. An absent claim is not checked — only a wrong one is a failure.
func checkRegeneratedKeys(got compare.Artifact, identity, equiv meta.KeyRef, claimant string) error {
	if got.Identity == "" {
		return fmt.Errorf("%w: its payload regenerates no identity, so nothing about it can be verified",
			ErrInvalidArtifact)
	}
	for _, check := range []struct {
		what      string
		want      meta.KeyRef
		gotScheme string
		gotDigest string
	}{
		{"identity", identity, got.IdentityScheme, got.Identity},
		{"equivalence", equiv, got.EquivScheme, got.Equiv},
	} {
		if check.want.Empty() {
			continue
		}
		if check.want.Scheme != check.gotScheme || check.want.Digest() != check.gotDigest {
			return fmt.Errorf("%w: %s %s is %s %s, its payload regenerates %s %s",
				ErrInvalidArtifact, claimant, check.what,
				check.want.Scheme, check.want.Digest(), check.gotScheme, check.gotDigest)
		}
	}
	return nil
}

// downloadSize is what must fit on the destination filesystem: the manifest, the
// config, and every layer once — not twice, because each layer streams into its
// own range of the finished file rather than being staged and joined.
func downloadSize(manifestDesc ocispec.Descriptor, manifest ocispec.Manifest) int64 {
	var payload int64
	for _, layer := range manifest.Layers {
		payload += positiveSize(layer.Size)
	}
	return positiveSize(manifestDesc.Size) + positiveSize(manifest.Config.Size) + payload
}

// positiveSize treats an unset or nonsensical size as zero rather than letting it
// subtract from the total.
func positiveSize(size int64) int64 {
	if size > 0 {
		return size
	}
	return 0
}

// requireFreeSpace refuses before the download when dir cannot hold it. Bavail,
// not Bfree: the blocks reserved for root are not available to this process.
func requireFreeSpace(dir string, required int64) error {
	if required <= 0 {
		return nil
	}
	var stat unix.Statfs_t
	if err := unix.Statfs(dir, &stat); err != nil {
		return fmt.Errorf("cannot check free space in %s: %w", dir, err)
	}
	if available := int64(stat.Bavail) * int64(stat.Bsize); available < required {
		return fmt.Errorf("not enough space in %s: need %s, %s available",
			dir, utils.FormatSize(required), utils.FormatSize(available))
	}
	return nil
}
