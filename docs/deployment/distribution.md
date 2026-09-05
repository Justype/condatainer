# Distributing Artifacts

An overlay that took four hours to build is worth uploading once. CondaTainer
publishes finished `.sqf` overlays and base images to any OCI registry —
GHCR, Harbor, Quay, an institutional Artifactory — so the next machine downloads
instead of rebuilding.

There are two distinct reasons to publish, and they get two different layouts:

| You are publishing | Who pulls it | Command | Page |
|---|---|---|---|
| a recipe collection's builds | anyone who configured that collection | `condatainer registry push` | [Sharing Your Scripts](./share_scripts.md) |
| one project's pinned artifacts | whoever checks that project out | `condatainer project push` | [Project](../manuals/condatainer.md#project) |

Both use the same transport, the same credentials, and the same
[publishing rules](../manuals/condatainer.md#publishing-rules). They differ only
in how artifacts are addressed, and that difference is worth understanding before
you look at a package list and wonder what you are seeing.

## Two schemes

```text
catalog scheme                          project scheme
──────────────────────────────          ─────────────────────────────────────────
ghcr.io/my-lab/cnt/                     ghcr.io/my-lab/rnaseq-2026/cnt
├── grch38/genome:gencode49             ├── :grch38--genome--gencode49
├── grch38/genome:gencode48             ├── :grch38--genome--gencode49__a31f902c12ab
├── star:2.7.11b                        ├── :star--2.7.11b
├── samtools:1.23.1                     ├── :star--2.7.11b__6f21c0b81a4d
└── ubuntu24/base:20260721, :latest     └── :samtools--1.23.1__0b4e9c72f118
    one package per artifact name           one package for the whole project
```

The catalog scheme puts the name in the repository path and the version in the
tag, which is what a container registry's UI, `skopeo`, and a human all expect.
Browsing `grch38/genome` shows its versions; each artifact is its own package
with its own visibility and its own retention.

The project scheme puts the whole name in the tag of one repository, with `/`
written as `--`.

## Why a project does not publish per artifact

The obvious design — reuse the catalog scheme in the project's own namespace —
runs into how GitHub Packages counts packages.

**A distinct repository path is a distinct package.** `my-lab/rnaseq-2026/star`
and `my-lab/rnaseq-2026/samtools` are not one package under two names; they are
two packages that happen to share a prefix. A project with fifteen artifacts
published the catalog way is fifteen packages.

**Each package's visibility is its own setting, inherited from nothing.** It does
not follow the code repository the package is linked to, and nothing in
CondaTainer touches it. So fifteen packages means fifteen settings changed by
hand, with nothing that changes them together and nothing that reports the one
you missed.

One repository is one package with fifteen tags: one setting to find, and one to
get right.

It is also the shape the artifacts actually have. A project's set is not a
catalogue anyone browses by name — it is one closure, pinned together, checked
out together, restored together. One package is a truer description of that than
fifteen unrelated ones.

## Why a project does not publish to the collection endpoint

A collection's endpoint carries what that collection's recipes produce, at the
identities the collection vouches for. A project's set is a different claim:
these are the exact builds *this analysis* ran with, including locally built
artifacts, path-addressed one-offs, and versions the collection never shipped.
Writing those into the collection's namespace would mean a package whose contents
nobody maintains and whose names collide with real recipes.

Credentials point the same way. Push access to a lab collection is a small set of
maintainers; push access to a project package is whoever works on that project.

## Why the identity tag

Every artifact a project pushes gets a second tag, `<name>__<12 hex>`, carrying
its identity.

Nothing fetches by tag. A `remotes` entry in the lock records the platform
manifest digest, and restore resolves that. What a tag has to do, then, is keep
the manifest **referenced** — content no tag reaches is unreferenced, and
unreferenced content is the registry's to reclaim, which every housekeeping job
people run against GHCR does on purpose.

A plain name tag can only keep one identity alive:

1. `project push` tags `star--2.7.11b` → manifest A, and the lock records A's digest.
2. Someone re-pins a different STAR build; the next push moves that tag to manifest B.
3. A becomes untagged, therefore unreferenced, therefore reclaimable.
4. A checkout of the older commit asks for A's digest and gets a 404.

The identity tag never moves, because it is derived from the content it names.
The plain tag stays as a human handle that says which identity the project uses
*now* — a moving pointer, like `latest`. Dropping the plain tag would break no
restore.

The cost is one extra tag per artifact: a manifest reference, zero additional
blobs.

The catalog scheme needs no equivalent. A versioned tag there is immutable —
`registry push` refuses to replace one without `--force` — so it never orphans
the manifest it named, and adding an architecture republishes an index that still
references the existing per-platform child. Only a same-day rebuild of a
version-less `base` or `os` replaces a date tag, and those are replaceable by
design.

## When `--all` is worth its bytes

By default `project push` skips any pin a configured collection **already
serves at that exact identity**. Restore tries the upstream location first, so a
second copy buys nothing but the upload — and for a reference genome index that
upload is measured in tens of gigabytes.

`--all` publishes them anyway. Do it when the project has to outlive the
collection:

- **an archived paper**, where the lock must still restore in five years and no
  collection's retention policy promises that;
- **a collection you do not control**, which may retag, prune, or disappear;
- **an air-gapped or restricted cluster**, whose nodes can reach your registry but
  not the collection's;
- **regulatory or institutional archival**, where "it is available upstream" is not
  an acceptable answer.

Do not do it as a matter of routine. For a project rebuilt from a live collection
on a network that can reach it, the upstream copy is the same bytes with someone
else paying for the storage.

`--closure` is a separate axis: it adds the build dependencies an artifact needed
to be produced, which a rebuild needs and a run does not. Push never builds, so
those have to be present first — `condatainer project restore --keep-build-deps`
is what installs them.

## Redistribution

Whether something *may* be published depends on who can pull it, and the answer
is declared, not guessed. Two recipe headers carry it:

```bash
#LICENSE: MIT
#REDISTRIBUTE: yes
```

`#REDISTRIBUTE:` is the answer; the artifact's `#TYPE:` is only the default for
an unanswered question. The full table is in the manual under
[Publishing rules](../manuals/condatainer.md#publishing-rules). The short version:

- `#REDISTRIBUTE: no` is an absolute stop at a public endpoint, and no flag overrides it.
- `#REDISTRIBUTE: yes` publishes anything, whatever the type.
- Undeclared, an `app` is refused at a public endpoint — someone else's software with
  unstated terms, and unknown is not permission.

`#LICENSE:` is an SPDX expression published verbatim as
`org.opencontainers.image.licenses`. It is never parsed and gates nothing:
deriving redistribution permission from a licence expression is a judgement a
tool gets wrong in the permissive direction, and that direction cannot be taken
back.

```{note}
`--audience restricted` on an endpoint lifts every one of these checks, because
it is a claim that only a known set of people can pull from it.

It is deliberately **not** called visibility. A GitHub package's visibility is a
separate per-package setting that CondaTainer neither reads nor changes, and
borrowing GitHub's word — whose values are `public`, `private` and `internal` —
for an unrelated declaration invited exactly that misreading. Nothing verifies
the claim; it is a declaration someone writes down and defends, which is why it
lives in the lock or a collection's `source.json` rather than in a flag typed at
push time.
```

An `app` from a third-party collection is the case with no clean answer today: the
`#REDISTRIBUTE:` header lives in a recipe your project does not own, so you cannot
declare for it. Publish to an endpoint you have declared `internal`.

## Credentials

Both commands read the same sources, in order:

1. `CNT_REGISTRY_TOKEN`, with optional `CNT_REGISTRY_USER`
2. `GITHUB_TOKEN`, for `ghcr.io`
3. the Docker/OCI credential store
4. anonymous

```bash
printf '%s\n' "$TOKEN" | condatainer registry login ghcr.io \
  --username "$USER" --password-stdin
```

For GHCR a classic PAT with `write:packages` is what a push needs; `read:packages`
is enough for a restore from a private package. In CI, the workflow's
`GITHUB_TOKEN` covers packages in the same repository with no extra secret.

See [Registry](../manuals/condatainer.md#registry) for layer planning, rate-limit
behaviour on multi-hour transfers, and the `pull`/`tags`/`resolve` surface.
