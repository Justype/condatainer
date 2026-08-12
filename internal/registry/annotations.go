package registry

import (
	"fmt"
	"path"
	"strconv"
	"strings"
	"time"

	"github.com/Justype/condatainer/internal/artifact/meta"
)

// Annotations builds the manifest-level annotation set published with m.
//
// Everything here is derived from what the artifact itself recorded, never
// accepted from a caller — a mirror must not be able to claim an artifact by
// asserting metadata about it. compression is the one value not in the manifest:
// it is read from the payload's own superblock, and it describes how to mount the
// bytes rather than what they are.
//
// Only facts a consumer must act on *before* fetching the blob belong here.
// Everything else an artifact knows travels inside it, where it is authoritative.
func Annotations(m meta.Manifest, compression string) map[string]string {
	ann := map[string]string{
		AnnTitle:  m.Name,
		AnnSchema: strconv.Itoa(m.SchemaVersion),
	}
	set := func(key, value string) {
		if value = strings.TrimSpace(value); value != "" {
			ann[key] = value
		}
	}

	set(AnnVersion, path.Base(m.Name))
	set(AnnDescription, m.Description)
	set(AnnURL, m.URL)
	set(AnnSource, m.Build.Source)
	set(AnnCompression, compression)

	// Omitted rather than filled with the push time: a build time nobody recorded
	// is not the same fact as the moment somebody uploaded it.
	if !m.Build.Created.IsZero() {
		ann[AnnCreated] = m.Build.Created.UTC().Format(time.RFC3339)
	}

	set(AnnIdentityScheme, m.Keys.Identity.Scheme)
	set(AnnIdentitySHA, m.Keys.Identity.SHA256)
	set(AnnEquivScheme, m.Keys.Equiv.Scheme)
	set(AnnEquivSHA, m.Keys.Equiv.SHA256)

	// Architecture is otherwise structural — the index child's platform
	// descriptor — but "runs anywhere" has no platform spelling.
	if m.Platform.Arch == meta.ArchNone {
		ann[AnnNoarch] = "true"
	}
	return ann
}

// Identity and Equiv read the complete keys back out of a published manifest's
// annotations. A half-present key is reported as absent: a scheme without a
// digest, or the reverse, addresses nothing.
func Identity(ann map[string]string) meta.KeyRef {
	return keyRef(ann[AnnIdentityScheme], ann[AnnIdentitySHA])
}

func Equiv(ann map[string]string) meta.KeyRef {
	return keyRef(ann[AnnEquivScheme], ann[AnnEquivSHA])
}

func keyRef(scheme, sha string) meta.KeyRef {
	if scheme == "" || sha == "" {
		return meta.KeyRef{}
	}
	return meta.KeyRef{Scheme: scheme, SHA256: sha}
}

// IsNoarch reports whether the published artifact declared #ARCH:noarch, and so
// may be pulled by any architecture.
func IsNoarch(ann map[string]string) bool { return ann[AnnNoarch] == "true" }

// Want is what a caller is looking for. Identity is optional: a pull by name has
// no identity to insist on, a restore pinning a lock does.
type Want struct {
	Name     string
	Identity meta.KeyRef
}

// Check reports why the artifact described by ann cannot satisfy want, or nil.
//
// This runs on annotations alone, which [ResolveArtifact] fetches without the
// blob, so a wrong artifact costs one manifest request instead of a download. It
// is not the last word: the payload's own embedded manifest is authoritative and
// is checked again after the transfer. Annotations are a fast index over that,
// never a replacement for it.
func Check(ann map[string]string, want Want) error {
	if len(ann) == 0 || ann[AnnTitle] == "" {
		return ErrNoAnnotations
	}
	if err := checkSchema(ann); err != nil {
		return err
	}
	if want.Name != "" && ann[AnnTitle] != want.Name {
		return fmt.Errorf("%w: published as %q, wanted %q", ErrMismatch, ann[AnnTitle], want.Name)
	}
	if !want.Identity.Empty() {
		if got := Identity(ann); got != want.Identity {
			return fmt.Errorf("%w: identity %s, wanted %s",
				ErrMismatch, describeKey(got), describeKey(want.Identity))
		}
	}
	return nil
}

// checkSchema reports whether this build can read the artifact's metadata.
//
// Both directions are named. A newer artifact needs a newer CondaTainer; an older
// one is a format this build has dropped, and telling someone to upgrade in that
// case sends them the wrong way.
func checkSchema(ann map[string]string) error {
	raw, ok := ann[AnnSchema]
	if !ok {
		return fmt.Errorf("%w: it records no metadata schema", ErrNoAnnotations)
	}
	schema, err := strconv.Atoi(raw)
	if err != nil {
		return fmt.Errorf("%w: metadata schema %q is not a number", ErrIncompatible, raw)
	}
	switch {
	case schema > meta.SchemaVersion:
		return fmt.Errorf("%w: artifact schema is %d, this build reads %d — upgrade CondaTainer",
			ErrIncompatible, schema, meta.SchemaVersion)
	case schema < meta.SchemaVersion:
		return fmt.Errorf("%w: artifact schema is %d, this build reads %d — the artifact predates this format and must be rebuilt",
			ErrIncompatible, schema, meta.SchemaVersion)
	}
	return nil
}

// describeKey renders a key for a message, naming an absent one rather than
// printing an empty string beside a real one.
func describeKey(k meta.KeyRef) string {
	if k.Empty() {
		return "absent"
	}
	return k.Scheme + " " + k.Digest()
}
