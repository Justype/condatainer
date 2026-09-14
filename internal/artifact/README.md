# artifact

Everything that understands what CondaTainer embeds in an image: what the image
says it is, what loading it does to the environment, and how two images are
compared.

It reads through [`internal/image`](../image/README.md), which handles archives
and files and knows nothing about metadata. The dependency is one-way: `artifact`
imports `image`, never the reverse.

## Architecture

- meta/: embedded runtime and manifest documents, staging, and validation
- key/: scheme dispatch, one rule file per scheme, and canonical encoding helpers
- capsule/: provenance closure composition and validation
- compare/: gates, verdicts, and human-readable differences

The dispatcher in key/scheme.go only selects a supported
identity/equivalence pair, reconstructs neutral inputs, and invokes the selected
schemes. It does not define which inputs contribute to a key. That policy lives
in the seven versioned scheme files.

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
`toolpath.ErrToolMissing`, `tool.ErrUnreadable` or `tool.ErrCorrupt`. An unknown
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

Degradation is deliberate, because most images in the wild predate the format:
a missing or unreadable runtime document never blocks a mount, it only means
this image contributes no PATH/env entries — a manifest's `type` is never
consulted before accepting something as a container root either.

### manifest.json

What the image is and where it came from: identity, build type, description, URL,
and platform. It carries **no runtime block at all** — that lives in
`runtime.json` and nowhere else.

`type` and `build_type` are both enumerations validated here, and neither derives
from the other: `catalog.Type` says what the payload *is* (`os`, `app`,
`data`, `env`) and `meta.BuildType` says how it was produced (`conda`, `script`,
`def`, `snapshot`). One recipe language produces both an app and a data image.
`internal/build` aliases `meta.BuildType` rather than keeping a second list of
its own — the value the build chooses is the value written here — and has no
member for `snapshot`, which a freeze produces without building anything.

Both enumerations are closed, and a value outside either is refused rather than
mapped to something safe. `Normalize` fills in an absent `type` as `app` and an
absent `os` as `linux`, and does nothing else — a lenient default for a field a
document left out, never a rewrite of one it filled in wrong. The type decides
the install prefix, whether a runtime document must carry one, the tag shape, and
whether a public endpoint accepts the artifact; reading an unknown one as `app`
would answer all four questions confidently and wrongly. `APP` is the case that
settles it: a difference in case is not a spelling to correct, and two names that
collapse into one are how an artifact gets published under a rule that was never
written for it.

The cost is that a type introduced later cannot be read by a build that predates
it. That is the same trade `schema_version` already makes with a strict
inequality two checks earlier, and it is the trade this package makes everywhere:
silence and a wrong answer must not look the same.

The one place the two coincide is a frozen environment: `snapshot` and `env`
imply each other, along with the name `env`, and `ValidateManifest` refuses a
manifest where they disagree. So a rule branches on whichever field it is really
about — the build type for publishing and key derivation, the type for what the
runtime does with the payload.

`dependencies[].identity` and `.equiv` are complete keys — scheme and SHA-256
both, the same `KeyRef` the top-level `keys` uses. An edge held to a bare digest
would be a weaker contract than the artifact it points at, and the schemes hash
both halves, so one digest under two schemes is two different edges.

The `build` block is what the build knew and the recipe does not say: tool
versions, the Conda channel order a solve used, the upstream image a definition
bootstrapped from, `build.source` (the repository of the collection that supplied
the recipe — never the local handle), and `build.created` (when the build
finished, stamped at staging).

`build.source` and the top-level `source` are different keys that share a name:
`source` is the embedded files and how to use them, `build.source` is where the
recipe came from. The name is deliberate — it maps one-to-one onto
`org.opencontainers.image.source`, which is the annotation it exists to derive.

Everything in `build` is diagnostic provenance. **No scheme hashes it**, with one
deliberate exception: `build.from.digest` reaches `definition-identity-v1`, which
is what makes an upstream rebuild a new identity. `build.source` and
`build.created` in particular must never move a key — the same recipe built from
a mirror, or built twice, is the same artifact, and the store deduplicates on
exactly that.

Manifest reads are uncached, because nothing asks for one per `exec`.

### Not handled

A writable `.img` carries no embedded metadata: it is a working overlay rather
than a built image, and its environment comes from its `.env` sidecar.

The boundary is the freeze, not the payload. `overlay freeze` packs one into a
`.sqf` carrying a manifest, a runtime document and a `snapshot-env-v1` key, and
from that point everything here applies to it as to any other artifact.

## Scheme-backed keys

The key package owns seven immutable derivation schemes. Each scheme has one file
that states its question, accepted artifact types, complete preimage contents,
and exclusions:

| scheme | implementation | rule |
|---|---|---|
| script-identity-v1 | [script_identity_v1.go](key/script_identity_v1.go) | recipe inputs plus every exact dependency identity |
| script-equiv-v1 | [script_equiv_v1.go](key/script_equiv_v1.go) | recipe inputs plus only substitution-relevant dependencies |
| definition-identity-v1 | [definition_identity_v1.go](key/definition_identity_v1.go) | definition inputs plus the resolved upstream digest |
| definition-equiv-v1 | [definition_equiv_v1.go](key/definition_equiv_v1.go) | definition inputs without the resolved upstream digest |
| conda-explicit-v1 | [conda_explicit_v1.go](key/conda_explicit_v1.go) | stored explicit.txt bytes exactly |
| conda-environment-v1 | [conda_environment_v1.go](key/conda_environment_v1.go) | stored environment.yml bytes exactly |
| snapshot-env-v1 | [snapshot.go](key/snapshot.go) | one record per packed archive entry, sorted by path |

Valid pairs are fixed by build type:

| build type | identity | equivalence |
|---|---|---|
| script | script-identity-v1 | script-equiv-v1 |
| def | definition-identity-v1 | definition-equiv-v1 |
| conda | conda-explicit-v1 | conda-environment-v1 |
| snapshot | snapshot-env-v1 | snapshot-env-v1 |

A snapshot is the one row whose two keys are a single value: it has no inputs to
abstract away, so "is this the same environment" and "can this substitute for it"
cannot come apart. It is also the one that regenerates from nothing a checkout
holds — the preimage is the packed payload — so verification takes its recorded
keys as they stand rather than re-deriving them.

manifest.keys stores a scheme and SHA-256 for each key. The complete key is the
(scheme, SHA-256) pair: the scheme is not repeated inside the hashed preimage,
and equal SHA values under different schemes remain different keys.

Verification selects the
named pair, loads the required source files, runs those exact implementations,
and compares the regenerated digests. An unknown identity scheme, an unknown
equivalence scheme, a mismatched pair, or an old key without a scheme is an
error.

Recipe-backed schemes construct an in-memory canonical Model. model.go provides
only neutral validation, sorting, line encoding, and hashing; it does not decide
which artifact inputs enter the model. recipe.go provides neutral conversion
helpers such as the comment-stripped recipe digest, environment values, and
sorted placeholders.

For script-identity-v1, every dependency is pinned by name and exact identity.
For script-equiv-v1, a data dependency contributes its equivalence key, an app
or OS named by the artifact contributes name/version, and a history-only app or
OS contributes nothing. The manifest freezes this role, so verification never
reapplies newer policy.

A dependency lacking either key is unrecorded. Script identity includes its name
plus that marker; script equivalence includes an unrecorded data dependency's
marker and name; and the manifest sets provenance_complete to false.

Conda schemes deliberately bypass the canonical model and hash their stored
exports byte for byte. The explicit export pins package URLs and build strings;
the environment export pins channels plus package names and versions.

Changing any rule requires a new scheme constant, a new implementation file, and
a new dispatcher case. Existing scheme files remain immutable.

## The capsule

`/.cnt/provenance` is the source closure an artifact was built from, one flat
directory per artifact:

```text
/.cnt/provenance/
  grch38--gtf-gencode--49@41ab1c2d3e4f/
    manifest.json
    recipe
  star--2.7.11b@0c19a7f34b02/
    manifest.json
    explicit.txt
    environment.yml
```

The directory name joins artifact-name slashes with `--` and appends the first
12 identity characters. The full identity is regenerated from the entry's
manifest and sources; validation requires the resulting name and identity prefix
to reproduce the directory name.

Name and identity together address an entry. Deduplication is by that pair, so a
diamond stores a shared dependency once. `runtime.json` and payload bytes are
never copied: capsule entries exist only to verify and rebuild artifacts.

### Composition, not traversal

```text
capsule(root) = ⋃ over direct deps d of  { sources(d) }  ∪  capsule(d)
```

`Compose` reads each dependency's own `/.cnt` — one extraction per dependency,
not one archive read per file — copies its manifest and rebuild sources in, and copies its capsule
entries across unchanged. Keys are verified from those sources; composition performs there is no recursive
resolution, no catalog access and no network at build time.

Two properties follow, and they are why this shape was chosen. **Completeness is
inductive**: if every dependency's capsule is complete, the union is. And
**cycles cannot occur**, because a dependency image existed before the artifact
that mounts it — `Validate` still rejects a self-reference, since a capsule read
from elsewhere is not trusted to be well-formed.

`provenance_complete` is inherited rather than recomputed: an unrecorded
dependency anywhere below makes everything above it incomplete. That is the
honest answer, and when the dependency is rebuilt with scheme-backed keys the artifact above
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
