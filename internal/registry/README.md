# internal/registry

`registry` publishes and fetches immutable CondaTainer artifacts through an OCI
registry. Registry is who the package talks to; OCI is the wire format it speaks.
The implementation uses oras-go directly, so compute nodes do not need an `oras`
binary.

The package has two related entry points:

- `Resolve(ctx, ref)` normalizes an upstream container image reference and
  returns the current machine's platform digest for build provenance.
- `Publish`, `ResolveArtifact`, `Check`, and `Pull` distribute completed `.sqf`
  overlays and `.sif` base images. Writable `.img` overlays have no immutable
  identity and are never distributed.

## Artifact contract

Overlay and base images have distinct manifest and layer media types. Pull checks
both before downloading. Artifacts larger than 512 MiB are pushed as ordered
chunks without writing a second local copy; pull verifies every OCI descriptor,
checks destination free space, reassembles in manifest order, regenerates the
embedded identity/equivalence keys, and atomically renames into place.

Native artifacts are children of a multi-platform OCI image index. A `noarch`
artifact is a plain manifest. `ResolveArtifact` always returns the selected image
manifest descriptor and its annotations, including its artifact type.

## References

- Versioned artifacts map the final name segment to the tag:
  `grch38/genome/gencode49` becomes repository `grch38/genome`, tag
  `gencode49`.
- Version-less base and two-segment OS artifacts use the complete name as the
  repository and publish a `YYYYMMDD` date tag plus `latest`.
- A digest selector (`repository@sha256:...`) addresses exact manifest bytes.

The first tag for a versioned artifact is immutable by default. `Publish` only
allows replacement with `Force`; adding a previously absent architecture is an
addition rather than replacement.

## Verification and errors

Annotations are derived only from the embedded manifest. `Check` rejects an
unsupported metadata schema or a requested name/identity mismatch before payload
transfer. Pull then verifies that the payload regenerates the published keys
before installation.

Callers classify outcomes with `errors.Is`: `ErrNotFound`,
`ErrUnsupportedPlatform`, and `ErrUnavailable` may permit a local build fallback;
`ErrUnauthorized`, `ErrIncompatible`, `ErrInvalidArtifact`, `ErrMismatch`, and
`ErrNoAnnotations` must be surfaced.

## Authentication

Credential precedence is:

1. `CNT_REGISTRY_TOKEN` and optional `CNT_REGISTRY_USER` for any host;
2. `GITHUB_TOKEN` for `ghcr.io` only;
3. Docker/OCI credential store, including configured credential helpers;
4. anonymous access.

`Login` and `Logout` manage the same store. Tokens are never logged.

## CLI

The publisher/admin surface is under one noun:

```text
condatainer registry push|pull|tags|resolve|login|logout
```

The CLI currently takes an explicit `--registry` base such as
`ghcr.io/example/condatainer`. Normal `create` does not need that flag: it keeps
the exact recipe source selected by the catalog and tries that source's ordered
`oci.pull` endpoints before building locally. The candidate must match the
equivalence key derived from the selected recipe. Explicit `registry pull` and
build use the same destination producer lock.
