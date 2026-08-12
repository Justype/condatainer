package registry

import (
	"context"
	"errors"
	"fmt"
	"os"
	"path/filepath"

	ocispec "github.com/opencontainers/image-spec/specs-go/v1"
	"golang.org/x/sys/unix"
	"oras.land/oras-go/v2"
	"oras.land/oras-go/v2/content/file"

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
func Pull(ctx context.Context, base, repo string, desc ocispec.Descriptor, annotations map[string]string, destPath string) error {
	guard, err := producer.AcquireLocal(destPath)
	if err != nil {
		return err
	}
	defer guard.Release() //nolint:errcheck
	return pull(ctx, base, repo, desc, annotations, destPath)
}

// PullLocked performs Pull while the caller already holds destPath's producer
// lock. Build uses this after deriving the expected equivalence under its lock;
// callers without a lock use Pull.
func PullLocked(ctx context.Context, base, repo string, desc ocispec.Descriptor, annotations map[string]string, destPath string) error {
	return pull(ctx, base, repo, desc, annotations, destPath)
}

func pull(ctx context.Context, base, repo string, desc ocispec.Descriptor, annotations map[string]string, destPath string) error {
	log := logging.FromContext(ctx)

	wantArtifactType, wantLayerType, err := imageTypes(destPath)
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
	layerNames, err := layerFilenames(manifest)
	if err != nil {
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

	if err := download(ctx, repository, desc, stageDir); err != nil {
		return err
	}
	staged, err := stageArtifact(stageDir, layerNames, filepath.Base(destPath))
	if err != nil {
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

// imageTypes reports the artifact and layer media types an image's extension
// implies. Push reads it from the file it is publishing, pull from the
// destination it was asked to write.
//
// On pull that makes the destination the statement of what is expected, which is
// the only place it exists: nothing published says whether an artifact is an
// overlay or a base, and nothing should — the caller already decided what it is
// installing. A `.sif` served where a `.sqf` was wanted then fails at the
// transport instead of at mount, which is where it would otherwise surface.
func imageTypes(path string) (artifactType, layerType string, err error) {
	switch {
	case utils.IsSqf(path):
		return ArtifactTypeOverlay, MediaTypeOverlayBlob, nil
	case utils.IsSif(path):
		return ArtifactTypeBase, MediaTypeBaseBlob, nil
	case utils.IsImg(path):
		return "", "", fmt.Errorf("%s is a writable overlay, which has no identity and is never distributed", path)
	}
	return "", "", fmt.Errorf("%s is not a distributable image (.sqf or .sif)", path)
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

// layerFilenames returns each layer's filename in manifest order, which is the
// order the chunks reassemble in.
func layerFilenames(manifest ocispec.Manifest) ([]string, error) {
	names := make([]string, 0, len(manifest.Layers))
	for i, layer := range manifest.Layers {
		name := filepath.Base(layer.Annotations[ocispec.AnnotationTitle])
		if name == "." || name == "" || name == string(filepath.Separator) {
			return nil, fmt.Errorf("%w: payload %d has no filename", ErrInvalidArtifact, i)
		}
		names = append(names, name)
	}
	return names, nil
}

// download copies the resolved manifest and its blobs into stageDir. It is
// addressed by digest on both sides, so no tag is resolved a second time and the
// artifact verified above is the one that arrives.
func download(ctx context.Context, repository oras.ReadOnlyTarget, desc ocispec.Descriptor, stageDir string) error {
	store, err := file.New(stageDir)
	if err != nil {
		return fmt.Errorf("cannot open staging store: %w", err)
	}
	defer store.Close() //nolint:errcheck

	ref := desc.Digest.String()
	if _, err := oras.Copy(ctx, repository, ref, withTransferProgress(store, verbDownload), ref, oras.DefaultCopyOptions); err != nil {
		return fmt.Errorf("failed to download %s: %w", ref, classify(err))
	}
	return nil
}

// stageArtifact assembles the downloaded layers and returns the complete
// artifact under wantName.
//
// The name matters: an image is read according to its extension, so an artifact
// staged under a publisher's layer title would be unreadable the moment that
// title carried the wrong one. Taking the name from the destination makes the
// staged file the destination in every respect except its directory entry.
func stageArtifact(stageDir string, layerNames []string, wantName string) (string, error) {
	pulled, err := assemblePulledArtifact(stageDir, layerNames)
	if err != nil {
		return "", err
	}
	staged := filepath.Join(stageDir, wantName)
	if pulled == staged {
		return staged, nil
	}
	if err := os.Rename(pulled, staged); err != nil {
		return "", fmt.Errorf("cannot stage pulled artifact: %w", err)
	}
	return staged, nil
}

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
// config, and every layer — plus the payload a second time when it arrives in
// chunks, because reassembly holds the parts and the whole at once.
func downloadSize(manifestDesc ocispec.Descriptor, manifest ocispec.Manifest) int64 {
	var payload int64
	for _, layer := range manifest.Layers {
		payload += positiveSize(layer.Size)
	}
	total := positiveSize(manifestDesc.Size) + positiveSize(manifest.Config.Size) + payload
	if len(manifest.Layers) > 1 {
		total += payload
	}
	return total
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
			dir, utils.FormatBytes(required), utils.FormatBytes(available))
	}
	return nil
}
