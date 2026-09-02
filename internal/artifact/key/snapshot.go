package key

import (
	"crypto/sha256"
	"encoding/hex"
	"errors"
	"fmt"
	"io"
	"os"
	"sort"
	"strings"

	"github.com/Justype/condatainer/internal/artifact/meta"
)

// ErrUnkeyed reports an artifact whose keys cannot be regenerated from vendored
// sources. Only a snapshot is, and callers that can recompute from the payload
// branch on this rather than treating it as a missing key.
var ErrUnkeyed = errors.New("artifact carries no regenerable keys")

// IsSnapshot reports whether a manifest describes a frozen writable overlay.
func IsSnapshot(m meta.Manifest) bool { return m.BuildType == meta.BuildTypeSnapshot }

// SnapshotEnvV1 names how a frozen environment's identity is taken: a sha256
// over one record per archive entry, sorted by path.
//
// It is the artifact's payload that is hashed, never the packed file. A repack
// changes the file — mksquashfs stamps its own creation time into the superblock
// — while an environment that has not changed must keep its identity, which is
// what every other scheme in this package promises and what a lock relies on.
const SnapshotEnvV1 Scheme = "snapshot-env-v1"

// TreeRecord is one archive entry's contribution to a snapshot's preimage.
//
// Only what changes behaviour is carried. Mode is here because the executable
// bit decides whether a binary runs; mtime, uid and gid are not, because
// touching a file does not make a different environment and -all-root flattens
// ownership to a constant. Size is absent as the content hash subsumes it.
type TreeRecord struct {
	// Type is the entry's kind, spelled as find does: f d l c b p s.
	Type byte
	// Mode is the permission bits.
	Mode uint32
	// Path is the entry's full path inside the archive.
	Path string
	// ID identifies the entry's content: the sha256 of the bytes for a regular
	// file, the target for a symlink, "major:minor" for a device, empty for
	// anything with no content of its own.
	ID string
}

// SnapshotIdentity hashes a snapshot's records into its identity.
//
// Records are sorted here rather than trusted in the order they arrive, so the
// ordering is a property of the scheme and not of whatever walked the archive.
// Byte order, because a locale-aware comparison would give one tree different
// identities on two machines.
func SnapshotIdentity(records []TreeRecord) (meta.KeyRef, error) {
	if len(records) == 0 {
		return meta.KeyRef{}, errors.New("snapshot identity: no entries to hash")
	}
	sorted := make([]TreeRecord, len(records))
	copy(sorted, records)
	sort.Slice(sorted, func(i, j int) bool { return sorted[i].Path < sorted[j].Path })

	sum := sha256.New()
	for _, r := range sorted {
		if _, err := io.WriteString(sum, r.preimage()); err != nil {
			return meta.KeyRef{}, err
		}
	}
	return meta.KeyRef{Scheme: string(SnapshotEnvV1), SHA256: hex.EncodeToString(sum.Sum(nil))}, nil
}

// preimage renders one record. The NUL before the content id is what keeps a
// path holding a space or a newline from shifting the framing — paths come from
// whatever the user installed, so they are untrusted input.
func (r TreeRecord) preimage() string {
	return fmt.Sprintf("%c %04o %s\x00%s\n", r.Type, r.Mode&0o7777, r.Path, r.ID)
}

// SnapshotPreimage renders the records the way SnapshotIdentity hashes them.
// Exported so a mismatch can be diffed: two artifacts that disagree are read by
// comparing these, not by staring at two digests.
func SnapshotPreimage(records []TreeRecord) string {
	sorted := make([]TreeRecord, len(records))
	copy(sorted, records)
	sort.Slice(sorted, func(i, j int) bool { return sorted[i].Path < sorted[j].Path })
	var b strings.Builder
	for _, r := range sorted {
		b.WriteString(r.preimage())
	}
	return b.String()
}

// ArtifactDigest is the sha256 of the packed file, in OCI form.
//
// Not an identity: it changes on every repack, so nothing may pin it. It exists
// because a registry descriptor speaks this and nothing else.
func ArtifactDigest(path string) (string, error) {
	sum, err := sumFile(path)
	if err != nil {
		return "", err
	}
	return "sha256:" + sum, nil
}

// sumFile streams a sha256 and returns it as bare hex, which is what meta.KeyRef
// holds and what validSHA256 accepts.
func sumFile(path string) (string, error) {
	f, err := os.Open(path)
	if err != nil {
		return "", fmt.Errorf("digest %s: %w", path, err)
	}
	defer f.Close()

	sum := sha256.New()
	if _, err := io.Copy(sum, f); err != nil {
		return "", fmt.Errorf("digest %s: %w", path, err)
	}
	return hex.EncodeToString(sum.Sum(nil)), nil
}
