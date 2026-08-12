// Package registry talks to OCI registries: reference mapping, an authenticated
// client, and the operations built on them. Resolving a tag to a digest is the
// first; push and pull are the same transport.
//
// Registry is who you talk to; OCI is what you speak. Everything here opens a
// connection or speaks the wire format, and nothing above it needs to know which.
//
// The cargo is a completed immutable artifact — a SquashFS `.sqf` overlay or an
// Apptainer `.sif` base — published as generic OCI blobs through oras-go, so
// CondaTainer stays a static binary with no `oras` CLI on a compute node. A
// writable `.img` has no identity and is never distributed.
package registry

// Artifact and blob media types. The artifact type is manifest-level; the blob
// type is per-layer. Both are checked on pull, because a `.sif` served where a
// `.sqf` is expected must fail at the transport rather than at mount.
const (
	ArtifactTypeOverlay  = "application/vnd.condatainer.overlay.v1"
	MediaTypeOverlayBlob = "application/vnd.condatainer.overlay.squashfs.v1"
	ArtifactTypeBase     = "application/vnd.condatainer.base.v1"
	MediaTypeBaseBlob    = "application/vnd.condatainer.base.sif.v1"
)

// Standard OCI annotation keys. Registries and generic tooling read these, so a
// fact with a predefined slot goes in the slot rather than into a private key.
const (
	// AnnTitle carries the canonical artifact name, e.g. "grch38/genome/gencode".
	// It is the name a pull confirms it received, not a human blurb.
	AnnTitle = "org.opencontainers.image.title"
	// AnnVersion is the name's last segment.
	AnnVersion = "org.opencontainers.image.version"
	// AnnCreated is manifest.build.created in RFC 3339, the same field the
	// version-less date tag derives from, so tag and annotation cannot disagree.
	// Omitted rather than filled with the push time.
	AnnCreated = "org.opencontainers.image.created"
	// AnnSource is the recipe collection's repository. GHCR links a package to a
	// repository with this, so it is functional rather than decorative. Omitted
	// when the collection declares none.
	AnnSource = "org.opencontainers.image.source"
	// AnnDescription carries #DESC:. GHCR renders description but not title.
	AnnDescription = "org.opencontainers.image.description"
	// AnnURL carries #URL: — the upstream project or vendor page.
	AnnURL = "org.opencontainers.image.url"
)

// CondaTainer annotation keys, for facts the image spec has no slot for.
//
// Each one exists because a consumer must act on it *before* fetching the blob.
// Everything else an artifact knows travels inside it, where it is authoritative;
// mirroring the manifest into annotations would be a second copy free to drift.
//
// The org.condatainer.* prefix is fixed forever — changing it fragments metadata
// across old and new artifacts.
const (
	// AnnIdentityScheme and AnnIdentitySHA confirm a fetch is the artifact a lock
	// pinned, before spending the bytes.
	AnnIdentityScheme = "org.condatainer.identity.scheme"
	AnnIdentitySHA    = "org.condatainer.identity.sha256"
	// AnnEquivScheme and AnnEquivSHA gate the prebuilt candidate a name maps to,
	// without downloading it.
	AnnEquivScheme = "org.condatainer.equiv.scheme"
	AnnEquivSHA    = "org.condatainer.equiv.sha256"
	// AnnNoarch marks an artifact that declared #ARCH:noarch. Architecture is
	// otherwise structural — the index child's platform descriptor — but "runs
	// anywhere" has no platform spelling, and inferring it from a missing index
	// is fragile.
	AnnNoarch = "org.condatainer.noarch"
	// AnnSchema is the embedded metadata schema, so a future format is refused
	// before the payload is transferred.
	AnnSchema = "org.condatainer.schema"
	// AnnCompression is the SquashFS compressor. The running kernel must support
	// it to mount the overlay, which is worth learning before the download.
	AnnCompression = "org.condatainer.squashfs.compression"
)
