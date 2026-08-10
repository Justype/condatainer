# internal/registry

Talks to OCI registries: reference normalization, an authenticated client, and
the operations built on them. Resolving a tag to a digest is the first one.

**Transport only.** It knows references, auth, and blobs, not what a recipe or a
record is. Reading `Bootstrap:`/`From:` out of a definition belongs to
`internal/build`; deciding what a CondaTainer artifact *is* as an OCI artifact —
media types, which blob is the payload — belongs above this, the same split
`image` and `artifact` already use.

## Why it exists

A tag is a moving target, so two `os` images built months apart from one
byte-identical definition share a recipe digest while their payloads differ by
every patch upstream shipped. A definition build resolves its bootstrap reference
*before* Apptainer runs: the digest goes into the identity record and into the
definition Apptainer is handed, so the record describes what was actually pulled.

## API

```go
digest, err := registry.Resolve(ctx, "ubuntu:24.04")
// -> "sha256:019e8eb29a85e74d64925745884f2ec79aa27e3feab36353d24656f4d6b89467"
```

`Normalize` follows the container-runtime rule: a bare name is a Docker Hub
official image (`ubuntu` → `docker.io/library/ubuntu`), one dotless component is
a Hub namespace (`myorg/tool` → `docker.io/myorg/tool`), and a component with a
dot or colon, or `localhost`, is a host and is left alone.

## Two decisions worth knowing

**Platform-specific, not the index digest.** A multi-arch tag resolves by default
to an index covering every architecture, which would give an x86_64 and an
aarch64 build one digest for two different payloads. `Resolve` passes the host
platform to `oras.Resolve`, which selects through the index.

**Failure is a normal outcome.** Nothing here fails a build. A login node that
cannot reach Docker Hub still produces a usable image, recording `unrecorded`
instead of a digest.

## Implementation

`oras-go/v2`, for registry auth — credential helpers included — that would
otherwise be hand-rolled per registry. `Normalize` and the authenticated client
are the shared foundation push and pull will reuse; `Resolve` is the first caller.

When they land, the registry's manifest digest supplies the third value in the
identity design — *are these the advertised archive bytes?* — which must come
from outside the archive.
