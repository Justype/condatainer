# artifact

Everything that understands what CondaTainer embeds in an image: what the image
says it is, what loading it does to the environment, and — as the identity work
lands — how two images are compared.

It reads through [`internal/image`](../image/README.md), which handles archives
and files and knows nothing about metadata. The dependency is one-way: `artifact`
imports `image`, never the reverse.

## Architecture

```
meta/       The two embedded documents: types, staging, archive reads, validation, cache
record/     The identity/equivalence record grammar: build, parse, sort, hash
key/        What reaches a key: source lines, and the dependency projection
capsule/    /.cnt/provenance: compose by union, read, validate
compare/    Gates, verdicts, diffs — how one artifact relates to another
```

## Two documents, two read frequencies

| file | read by | how often |
|---|---|---|
| `/.cnt/runtime.json` | container setup | every mount |
| `/.cnt/manifest.json` | `info`, comparison, restore | on demand |

The split is the point. One document carrying both the mount-time contract and
all provenance means every field added for provenance is paid for by every
`exec`. `runtime.json` is a few hundred bytes and stops growing; the manifest is
free to grow without a runtime cost.

Both are **trusted, not verified**: nothing compares them against the payload
beside them or the build that produced them, and validation checks structure
rather than truth.

### runtime.json

Everything container setup needs and nothing else: name, type, description,
platform, install prefix, and the `#ENV:` contributions.

`Prefix` is not a mount point. An overlay is applied over the container root, so
the payload appears at `/cnt/<name>` because that is where it sits in the
archive, not because anything is mounted there. Readers take the recorded value
as authoritative rather than reconstructing it from the name. `Env` keeps
`{prefix}` intact and substitutes it at load time, since the install prefix is
not known when the image is built — and an `os` artifact's empty prefix is not a
special case, because its files really are at the root, so `{prefix}/bin` is
`/bin`.

`ReadRuntime` distinguishes *no runtime document* from *could not look*: only a
genuinely absent one is `ErrNoRuntime`, so a caller never reports a missing
`unsquashfs` or a corrupt archive as "this image has no metadata" — those keep
`tool.ErrToolMissing`, `tool.ErrUnreadable` or `tool.ErrCorrupt`. An unknown
`SchemaVersion` joins them as `ErrUnsupportedSchema`, handled exactly like a
missing document; a reader ignores unknown fields, so a later schema that only
adds fields stays readable here.

**There is no fallback.** An image built before this layout has no
`runtime.json`, and nothing reads a runtime block out of its manifest to
compensate. It mounts, contributes nothing, and reports the same diagnostic as an
image with no metadata at all. A fallback would mean a second archive read on the
miss path and two code paths that must agree about what a runtime is — reading
exactly one file from exactly one place is worth more than sparing a rebuild.

Runtime reads are cached **across processes**, keyed by absolute path and
validated against size and mtime. Repeated scans — `list`, `avail`, PATH
construction, and shell completion, which runs one process per keystroke — would
otherwise spawn one `unsquashfs` per image every time. Negative verdicts are
cached too, or an image predating the format would be re-probed on every listing,
which is the cost the cache exists to avoid.

Degradation is deliberate and asymmetric, because most images in the wild predate
the format. `CheckBase` accepts an image with no runtime document, warns and
accepts one that is present but unreadable — rejecting it would strand every
build behind a base that is most likely fine — and rejects only metadata that
reads and says it is not a base.

### manifest.json

What the image is and where it came from: identity, build type, description, URL,
and platform. It carries **no runtime block at all** — that lives in
`runtime.json` and nowhere else.

Manifest reads are uncached, because nothing asks for one per `exec`.

### Not handled

A writable `.img` carries no embedded metadata. It is a mutable working overlay
rather than a built image, and its environment comes from its `.env` sidecar.

## Records

`record` is the grammar of `/.cnt/identity.record` and `/.cnt/equiv.record`. It
builds them, parses them back, and hashes them. It knows nothing about images,
recipes, or *what belongs* in a key — that is the key package's job.

A record is UTF-8 text: a format line, then `key=value` lines, LF-terminated, no
trailing whitespace, no blank lines, no comments.

```text
cnt-identity-v1
type=data
env=STAR_INDEX_DIR={prefix}
recipe=sha256:1a4f…
ph=gencode_version=49
ph=star_version=2.7.11b
dep=app samtools/1.23.1 unrecorded
dep=app star/2.7.11b sha256:0c19…
dep=data grch38/gtf-gencode/49 sha256:41ab…
```

Text rather than JSON, deliberately: no canonicalization question, diffable in a
terminal, and `sha256sum identity.record` reproduces the key by hand.

Fields appear in a fixed order — `type`, `env`, source, `ph`, `dep` — and
repeated keys sort within their own group: `env` by key, `ph` by name, `dep` by
the full text following `dep=`, which puts `app` before `data` before `os`.
Repeated `dep` lines are kept rather than collapsed, because a recipe that
mounted two equivalent inputs did mount two.

`Marshal` sorts, so a caller may build a record in whatever order its inputs
arrive. `Parse` does **not**: it accepts only the exact bytes `Marshal` produces,
and rejects an unsorted group, a group out of order, a blank line, trailing
whitespace, or an unknown key. A record that parsed but did not re-marshal
identically would hash to something other than the key it carries, and an unknown
key means a format that should have changed its tag.

`Kind.Tag()` — `cnt-identity-v1`, `cnt-equiv-v1` — versions the whole derivation:
the field set, the ordering, and every rule about what reaches a line. Nothing may
change under a tag without moving every key computed under it.

Three things the grammar enforces, because they are not per-kind: a `sha256:`
token is always a full 64 hex characters; `from=` may appear only under the
identity tag; and `type=base` is accepted as a record's *subject* but rejected on
a `dep=` line. A base is built from a definition like any other image, so it
identifies itself the same way — what it may never be is something another
artifact depends on, because it is the environment a build runs inside rather
than an input the build consumed. What a `dep=` line's remaining fields *mean*
differs between the two kinds, so here they are opaque tokens.

## What reaches a key

`key` is the single definition of every rule about what goes into a record, so a
build and a later comparison cannot disagree about an artifact's keys.

`RecipeDigest` hashes the recipe **with its whole-line comments removed**. The
header is all comments and goes with them, which is the point rather than a side
effect: every header item that belongs in a key already has its own record line,
and hashing the file too would re-admit the ones left out on purpose.

| header | reaches a key as |
|---|---|
| `#TYPE:` | `type=` |
| `#ENV:` | `env=`, via `Env` |
| `#DEP:` | `dep=`, carrying the resolved dependency's own identity |
| selected `#PH:` | `ph=`, via `Placeholders` |
| `#DESC:`, `#URL:`, `#INPUT:`, `#ARCH:`, scheduler directives, `#PH:` menus | nothing |

That last row is not merely noise — `#AUTOUPDATE:` rewrites `#PH:` menus and
repins `#DEP:` with no human involved, so a key that moved for them would mint a
new artifact for a robot's edit and `exact` would be unreachable for any
template. Nothing is lost: a repinned `#DEP:` still moves identity, through the
`dep=` line carrying the new dependency's digest.

Both records take the same `recipe=`, `ph=` and `env=` lines. They diverge at
`dep=` and nowhere else, which is why only data ever has two records that say
different things — and why an artifact with no dependencies has two records
differing only by their format line. That is honest rather than degenerate: a
self-contained artifact built from one recipe has nothing to be loose about.

## The dependency projection

`Records` builds both records; the difference between them is entirely here.

| dependency | identity | equivalence |
|---|---|---|
| any | `dep=<type> <name> <identity>` | — |
| data | ↑ | `dep=data <equiv>` — no name |
| app or os the artifact's name mentions | ↑ | `dep=<type> <name>/<version>` — no digest |
| app or os the name does not mention | ↑ | nothing |

Identity pins everything that was mounted, because anything mounted while the
recipe ran could have shaped the payload. Equivalence keeps only what decides
substitution, and the two forms are deliberate opposites: for a data dependency
the *content* is the contract, so it contributes a digest and its name is
irrelevant; for an app the *name and version* are the contract, so a Conda-built
and a script-built STAR of one version are interchangeable producers of the same
index.

Only direct dependencies are ever processed. Data equivalence composes
transitively without walking, because a data dependency's own `equiv` already
covers its important inputs.

**The name decides which tools count.** `HasComponents` requires the
dependency's `name/version` to appear as a contiguous run of components in the
artifact's own name, so `grch38/star/2.7.11b/gencode49` counts `star/2.7.11b` and
`grch38/star2.7.11b/gencode49` does not. The convention is the contract: a tool
whose version changes the result belongs in the name, and an author who gets that
wrong is already wrong about how the artifact should be named. This cannot drift
for a given artifact, because name is a comparison gate — two differently-named
artifacts are never compared to each other.

## Unrecorded dependencies

Every image built before this format carries no keys, and dependency resolution
satisfies a dependency from an installed image without reading a recipe. Failing
the build is not an option; pretending is worse. So a dependency with no records
gets `dep=<type> <name> unrecorded` in identity, `dep=data unrecorded <name>` in
equivalence where a digest would have gone, and `provenance_complete: false` in
the manifest.

Such a record is still stable and comparable — two builds against the same
unrecorded dependency produce the same records. It simply carries a weaker claim,
stated out loud. When that dependency is rebuilt with records, the artifact above
it gets new keys, which is correct: its inputs became knowable.

## What each build type records

| type | identity | equivalence | `manifest.keys` names |
|---|---|---|---|
| recipe (script app, data, os, base) | `identity.record` | `equiv.record` | those two files |
| Conda app | the explicit export | the environment export | `explicit.txt`, `environment.yml` |

A definition build — an `os` or a `base` — adds one line to its identity record
and nothing to its equivalence record: `from=`, the digest its `From:` reference
resolved to before the build ran. The `Bootstrap:` and `From:` directives sit in
`manifest.build.from` as provenance; only the digest is hashed. The reference as written is already inside the
recipe and therefore already in both. So an OS rebuilt after an upstream security
push is a different recorded build that still substitutes, which is what asking
for `ubuntu:24.04` meant. An image carrying no `keys` block at all was imported
rather than built, or built before this format, and compares `unverifiable`.

`manifest.keys` names the file each key hashes, which is what lets a reader
verify one without knowing which build type produced the image: hash the named
file, compare. A Conda app writes no records because everything one could hold is
ruled out for it — no name or prefix, no `env` since it has no recipe to take one
from, no `ph`, no `dep` — leaving a file whose entire content would be a format
tag, a type line and one digest.

Two consequences, stated rather than left to be rediscovered: a recipe that
writes comment lines into its payload — a launcher heredoc beginning
`#!/bin/bash` — can change that payload without moving a key; and a definition
that bootstraps from a mutable tag hashes text that stays still while what it
pulls moves. Apptainer records what it actually pulled at
`/.singularity.d/labels.json` inside the image, so that question is answerable on
demand; it simply does not enter a key.

## The capsule

`/.cnt/provenance` is the closure an artifact was built from, one flat directory
per artifact in it:

```text
/.cnt/provenance/
  grch38--gtf-gencode--49@41ab1c2d3e4f/
    manifest.json
    identity.record
    equiv.record
    recipe
  star--2.7.11b@0c19a7f34b02/
    manifest.json
    explicit.txt
    environment.yml
```

The directory name joins the artifact's slashes with `--`, exactly as image
filenames and store entries already do, so a dependency's own name never becomes
directory levels. `@` cannot occur in a Conda name or version, so it never
collides with the separator. The identity is truncated to 12 characters for
addressing only — the full digest is inside the entry's `identity.record`, and a
reader verifies against that rather than trusting twelve characters of a
filename.

Name and identity *together* address a record set; identity alone never does,
since one solve published under two names has one identity. Deduplication is by
that pair, so a diamond stores the shared dependency once.

`runtime.json` is never copied. A capsule entry exists to rebuild an artifact,
never to mount one, and a rebuilt dependency derives its runtime from its recipe.
Payloads never come along either: a recipe plus its placeholders is the complete
rebuild input for a recipe dependency, and `explicit.txt` is for a Conda one.

### Composition, not traversal

```text
capsule(root) = ⋃ over direct deps d of  { records(d) }  ∪  capsule(d)
```

`Compose` reads each dependency's own `/.cnt` — one extraction per dependency,
not one archive read per file — copies its records in, and copies its capsule
entries across unchanged. Nothing is re-derived, so there is no recursive
resolution, no catalog access and no network at build time.

Two properties follow, and they are why this shape was chosen. **Completeness is
inductive**: if every dependency's capsule is complete, the union is. And
**cycles cannot occur**, because a dependency image existed before the artifact
that mounts it — `Validate` still rejects a self-reference, since a capsule read
from elsewhere is not trusted to be well-formed.

`provenance_complete` is inherited rather than recomputed: an unrecorded
dependency anywhere below makes everything above it incomplete. That is the
honest answer, and when the dependency is rebuilt with records the artifact above
gets new keys — correct, because its inputs became knowable.

**The closure is shallow.** Only data has dependencies and only data dependencies
extend it, so it is a chain through data with apps and OS artifacts as leaves.
That bound is what keeps a capsule small enough to embed without thinking about
it.

## Comparison

`compare` is where anything decides whether a candidate is the artifact that was
asked for: store selection among variants of one name, `build --update`,
dependency resolution choosing among installed images, an explicit verify.
Everything calls this rather than restating the rules.

`Read` builds one side from an image in **one extraction**, and recomputes every
key from the file `manifest.keys` names rather than trusting the number beside
it. A manifest is what the publisher wrote; a key that does not match its own
file is exactly what verification exists to catch.

### Gates, then keys

1. **architecture** — against the host, not against the other artifact: something
   that cannot run here is not a candidate whatever it contains. `noarch` passes
   anywhere. `platform.os` is never compared.
2. **type** — an app and a dataset are never interchangeable.
3. **name** — a substitution question is always about a requested name.
4. **usability** — an image whose records cannot be read makes no claim.
5. **runtime** — `runtime.json`'s `env` must equal the `env=` lines in the
   identity record. Free integrity: `runtime.json` was never part of the hashed
   preimage, so an edited one cannot agree with the record beside it. `prefix` is
   *not* checked — the image states where its payload sits, and reconstructing
   that from a naming convention would reject any artifact whose prefix was
   chosen differently.

Then identities equal ⇒ `exact`; equivalences equal ⇒ `equivalent`; comparable
and different ⇒ `different` with a diff; anything missing or malformed ⇒
`unverifiable`, which is **not** a synonym for different — it means the question
was not answered.

A gate failure is a verdict with a specific diff, never a fallthrough to
equivalence: artifacts that fail a gate are not comparable, so a key match
between them would be meaningless rather than reassuring.

### Diffs

| artifact | names |
|---|---|
| recipe build | `recipe`, `ph:<name>`, `env:<KEY>`, `dep:<name>` |
| Conda app | `packages` (the explicit export), `environment` |

A dependency's own name is what a human needs, and the equivalence record
deliberately drops it — so diffs read `manifest.dependencies`, which carries
name, identity and equiv side by side. A history-only dependency change appears
in a diff while the verdict stays `equivalent`; that asymmetry is the design.

### The one thing on the execution path

`MountAllowed` compares a runtime document's architecture against the host, and
container setup calls it. It is not really an exception: `runtime.json` is
already in hand and `platform` is in it, so it costs a string comparison — and
without it a wrong-architecture image mounts cleanly and fails somewhere further
downstream where the cause is unrecognizable. Everything after that gate belongs
to comparison, which mounting never runs.

### What comparison does not prove

Records are what the publisher wrote. They do not prove the payload matches them:
two images with byte-identical records and completely different payloads take no
cleverness to produce. A candidate from an untrusted source needs a signature
policy, which this does not provide.
