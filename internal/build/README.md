# build

Build system for creating overlay images from Conda packages, shell scripts, or Apptainer definition files with dependency resolution and scheduler integration.

## Two axes

`BuildType` is **how** a target is built. `catalog.Type` is **what** the payload is.
They are independent — a script build can be any type. (Go reserves `type`, so
locals and parameters spell it `typ`.)

| BuildType | Source | Builder |
|------|--------|-------------|
| Conda | Package spec, YAML, or `channel::pkg/version` | Micromamba environment |
| Def | Apptainer `.def` file | Apptainer build |
| Script | Recipe (`#DEP`, `#SBATCH`, `#ENV`, `#INPUT`) | Recipe run as a shell script |

| Type | Comes from | Decides |
|---|---|---|
| `os` | a `.def` recipe | — |
| `app` | `#TYPE:app`, else a one-component name | fast scratch tmp root, `BlockSize` |
| `data` | `#TYPE:data`, else 2+ components | stable tmp root, `DataBlockSize` |

Type affects **only** where the build works and how the archive is compressed.
Every build writes to `$CNT_PREFIX` = `/cnt/<name>` — the install prefix recorded
in the image, so anything baked into the payload stays valid at run time.

## Architecture

```
object.go       BuildObject interface, base implementation
spec.go         Build/Spec/Options/Workspace/Target; the build environment
conda.go        Conda package builds via micromamba
script.go       Recipe (script) builds, any type
def.go          Apptainer definition file builds
bootstrap.go    Bootstrap:/From: parsing, upstream resolution, definition pinning
keys.go         identity/equiv schemes, dependency keys, capsule composition
pack.go         Metadata staging, shared by every backend that packs
squashfs.go     The common packer: payload plus staged metadata → one .sqf
base_image.go   Resolving the base image every non-def build runs inside
deps.go         Overlay arguments for a build's dependencies
graph.go        Build planning from a solved catalog plan, and its execution
tmppath.go      Which tmp root a build works under, by type
```

## Key Types

**BuildObject** - Interface for all build types:
- `NameVersion()`, `BuildSource()`, `Dependencies()`, `BuildType()`
- `Spec()`, `Manifest()`, `Runtime()` - What is being built, and the two documents it embeds
- `IsInstalled()`, `RequiresScheduler()`, `Update()`
- `ScriptSpecs()` - Parsed scheduler directives
- `LockPath()` - Path to the build lock file
- `Build(ctx, buildDeps bool)`, `GetMissingDependencies()`
- `CreateTmpOverlay(ctx, force bool)`, `Cleanup(failed bool)`

**Spec** - What the image will be, settled by resolution and never changed after:
`ImageSpec` (`Name`, `Type`, `Description`, `URL`, plus the runtime `Prefix` and
`Env`), the `SourceSpec` whose populated field *is* the `BuildType`, the base, and
dependencies. `Spec.Manifest()` and `Spec.Runtime()` render the two embedded
documents, and both are projections of `Spec` and nothing else — no host, job ID,
local path, or `#INPUT:` answer can reach either.

**BuildLockInfo** - JSON metadata stored in `.lock` files:
- `Runner` - `"local"`, `"slurm"`, `"pbs"`, `"lsf"`, or `"htcondor"`
- `JobID` - Scheduler job ID (empty until submit returns; set by `submitJob`)
- `Node` - Short hostname (domain suffix stripped) where lock was created
- `PID` - OS process ID for local builds; 0 for scheduler builds
- `CreatedAt` - RFC3339 timestamp

**BuildType** - `BuildTypeConda`, `BuildTypeDef`, `BuildTypeScript`

**BuildGraph** - Turns a solved `catalog.Plan` into build work: a BuildObject per missing node, run in dependency order — local builds and prebuilt pulls sequentially, then scheduler submissions with job dependencies.

**LockedSpec** - A rebuild described by a project lock rather than the catalog. See **Locked rebuilds** below.

**ScriptSpecs** - Build resource requirements (alias to `scheduler.ScriptSpecs`)

## Recipes

Recipes support metadata headers:
- `#DEP:name/version` - Exact dependency
- `#DEP:name/version>=min` - Dependency with version constraint (range `[min, version]`)
- `#SBATCH` / `#PBS` / `#BSUB` - Scheduler directives (HTCondor uses native `.sub` files)
- `#ENV:VAR={prefix}/sub` - Environment variables to export; `{prefix}` is filled with the install prefix at load time
- `#SOURCE:name url` / `#SOURCE:name ask:prompt` - A file the build downloads, fetched by the tool before the recipe runs and exposed read-only as `$CNT_SRC_<name>`; its digest enters the identity (see Sources)
- `#INPUT:prompt` - User input, fed to the recipe on stdin in declaration order — read it with `IFS= read -r VAR` (collected locally before scheduler submission; embedded as a heredoc in the job script)

Available variables: `$CNT_NAME` (complete name), `$CNT_TYPE`, `$CNT_PREFIX`, `$CNT_TMP` (also `$TMPDIR`, both `/cnt_tmp`),
plus the scheduler's normalized `$NCPUS`, `$MEM`, `$MEM_GB`.
Run as `bash -euo pipefail <recipe>` top to bottom — no `install()` wrapper.
A script build also gets its `bin/` appended to `PATH` and
`MAMBA_ROOT_PREFIX` set to the scratch path, its package cache pinned beneath it (`CONDA_PKGS_DIRS`,
`MAMBA_PKGS_DIRS`) and `MAMBA_NO_RC=true` — one `micromambaEnv` shared with the Conda build, because
Apptainer passes the host environment through and an inherited cache variable would choose the cache — so `micromamba` is callable and its caches stay out of the user's home. A Conda or script build provisions a micromamba-only toolchain first when the host has none (a
download, once); a failure to do so fails the build. A build submitted to the scheduler provisions it on
the submitting host, before the job is queued, because the compute node may have no outbound access.

### Dependency Resolution

**`#DEP:` is a build dependency and nothing else.** It names what must be mounted
*while this recipe runs* — a tool that unpacks, validates or indexes the payload.
It is not recorded in the image and is not re-expanded when the image is
mounted later: there is no runtime dependency tree. What an overlay needs at run
time is whatever the caller names on the command line.

That is why an `app` is **self-contained** — a conda environment or a prebuilt
package that carries its own libraries — while `data` is the type that normally
has deps, because producing an index needs the tool that produces it.

`#DEP:samtools/1.22.1>=1.10` — accepts any installed version in `[1.10, 1.22.1]`:
- If a satisfying version is installed → skip build, mount the latest satisfying version
- If none installed → build the preferred version (`1.22.1`)
- Versions above preferred (`2.0`) are rejected (implicit upper bound)
- Operators: `>=` (inclusive lower bound) and `>` (exclusive lower bound)
- The dependency edge in the artifact's keys is read from the image that was mounted, resolved
  with the constraint, so it names that version and carries its identity. Reading the preferred
  version instead would record a tool that did not run, or `unrecorded` when it is not installed.

Resolution itself lives in `catalog.Resolve`; see **BuildGraph Execution** below.

## Locked rebuilds

`NewLockedObject` (`locked.go`) builds one artifact from a project lock's vendored
records instead of the catalog. It reuses every build implementation, the
workspace, base resolution and scheduler parsing, and is defined by what it
refuses: no catalog lookup, no fuzzy version choice, no prebuilt pull by name, no
flat install.

- **Dependencies are paths.** `LockedSpec.Deps` supplies one absolute image path
  per manifest dependency edge, in that order, and those paths land in
  `Spec.Dependencies` verbatim. Everything downstream already reads a
  dependency's own manifest for its name and keys, so an edge is recorded the
  same way whether the dependency was located by name or handed over as a path.
- **The output belongs to the caller.** Nothing installs and no flat name is
  claimed. An occupied output is refused: restore builds into a temporary
  sibling and renames.
- **A template runs expanded and embeds the template**, the same split the
  catalog path makes, with the placeholders taken from the manifest — the only
  record of which variant this is.
- **A definition keeps its recorded upstream digest.** `resolveUpstream` returns
  early for a locked build, because re-resolving would follow a tag that has
  since moved.
- **A Conda rebuild replays `explicit.txt`**, whose bytes *are* the recorded
  identity. `environment.yml` would be a fresh solve against whatever the
  channels serve today.
- **The root is supplied, not resolved.** `LockedSpec.Base` sets `Spec.Base`
  directly, so `resolveBase`'s early return (below) skips `ResolveBase`
  entirely — a project restore hands over the artifact its lock's reserved
  base pin names, rather than letting this package fall back to
  `config.GetBaseImage()`. Empty behaves exactly like the catalog path: the
  configured default is resolved and, if missing, built. A `.def` rebuild
  never reads `Spec.Base` either way.

Nothing here checks the result. The rebuild derives its own keys from what it
actually mounted and built, and the caller compares them against the lock — so a
disagreement surfaces as a key mismatch rather than as a plausible artifact under
the right name. Construction therefore refuses only what would make that
comparison uninterpretable: absent keys, a missing vendored source, a dependency
count the recipe does not agree with, or an answer set that does not match the
prompts — `#INPUT:` first, then each `#SOURCE: … ask:`.

## Build Lock

Each image target has a producer lock at the image path plus `.lock`, containing
`BuildLockInfo` JSON. The implementation lives in `internal/image/producer` so
build and registry pull serialize on the same pathname. This is distinct from
the inode lock used while an installed image is mounted: the producer lock
coordinates creation; the inode lock protects current readers at replacement.

**Lifecycle:**
- **Scheduler submit** (`submitJob` in `graph.go`): lock created with `runner=<scheduler>`, `job_id=""`, **no node** before calling the scheduler. Updated with the real job ID after `Submit()` returns. Removed if submit fails.
- **Local / scheduler job start** (`Build()`): `createBuildLock()` writes `runner=local` with current node and PID. If a scheduler lock already exists with a matching job ID (`$SLURM_JOB_ID` / `$PBS_JOBID` / `$LSB_JOBID`), it adopts and updates that lock with the runtime node and PID.
- **Build end**: `defer removeBuildLock()` removes the lock on success or failure.
- **Explicit registry pull**: acquires the same producer lock for resolve-through-install.
- **Automatic prebuilt pull**: runs under the build lock already held, so it does not reacquire it.

**Stale detection** (`clearStaleLock`, run by every constructor):

| Lock state | Action |
|---|---|
| Empty or corrupt file (old format) | Treat as stale — remove and continue |
| `runner=slurm/pbs/lsf/htcondor`, `job_id=""` | Stale (failed before submit returned) — remove |
| `runner=slurm/pbs/lsf/htcondor`, `job_id` set | Call `scheduler.IsJobAlive(job_id)` — remove if not alive |
| `runner=local`, same node, PID dead | Stale — remove and continue |
| `runner=local`, same node, PID alive | Active — error with details |
| `runner=local`, different node | Cannot verify remotely — error, ask user to check |
| Scheduler unavailable for job check | Conservative — error, ask user to check |

## Build Workflow

Conda and script share the same three middle phases — **install → stage → pack** —
because a Conda environment and a recipe payload differ only in how the payload
directory got filled. Both end in `packOutput`.

**Conda:**
1. Check if overlay exists (skip if not updating)
2. If updating existing overlay: probe exclusive lock — fail immediately if in use
3. Create build lock (local); create the host workspace directories
4. `installConda`: `micromamba create` inside container, and nothing else. The
   install run never sees the output path or binds the images directory, so a
   cancelled install cannot leave anything beside the installed images
5. `stageMetadata`, then pack to SquashFS
6. Atomic rename prepared → target; remove lock

**Shell:**
1. Check if overlay exists (skip if not updating)
2. If updating existing overlay: probe exclusive lock — fail immediately if in use
3. Create build lock (local, or adopt scheduler lock)
4. Build missing dependencies (if enabled), derive the recipe equivalence key,
   and try the selected source's ordered pull endpoints. A matching artifact is
   installed and stops here (a planned pull uses the candidate planning found);
   absent, unsupported-platform, or unavailable
   artifacts fall through to the local build, as does one made from a different
   recipe, with a warning. Credential, schema and type failures stop instead.
5. Fetch each `#SOURCE:` into the workspace and record its digest
6. Resolve the build base and create the temporary overlay
7. Run recipe inside container. The payload directory is bound at
   `/cnt/<name>` — the leaf, not `/cnt`, so dependency overlays mounted beside
   it stay visible
8. `stageMetadata`, then pack: verify the payload directory is non-empty, then
   `mksquashfs` the host build dir (basename `cnt`, so the archive is exactly
   `cnt/<name>/…` — the build's `tmp/` is a sibling and never enters it)
9. Atomic rename prepared → target; remove lock

### Sources

A `#SOURCE:` is fetched by the tool, not the recipe, because the tool can only hash what it
downloaded. Each file lands in `<workspace>/sources`, is bound read-only at `/cnt_src`, and is named
to the recipe as `$CNT_SRC_<name>`. Read-only, so a recipe cannot alter a file whose digest is
already recorded. The digests go to `Spec.Source.Fetched`, where the manifest and both key paths
read them.

**Order.** The fetch sits after the prebuilt check and before `skipIfInstalled`. A published
artifact is decided on equivalence, which no source enters, so nothing is downloaded to answer a
question it could have answered. `PredictIdentity` needs the digests, so it cannot run earlier; it is
reached only by `--store`, where fetching first turns a wasted rebuild into a wasted download. The
fetch must stay below `createBuildLock`: on a compute node the lock adopts a scheduler lock and
re-sites the workspace root, and the source directory is resolved from it afterwards.

**Answers.** `ask:` prompts join the `#INPUT:` prompts in one list, `#INPUT:` first, collected
before submission and carried in the job's heredoc like any answer. The recipe's stdin gets only the
first `recipeInputs` of them and the fetch takes the rest, so a submitted build needs no terminal.
An empty or non-http answer fails the fetch, which is what `--yes` produces.

**Failure and secrecy.** A short read fails and deletes the partial file, since the digest of a
truncated file is as good a digest as any and would fix an identity for bytes nobody wants. An
answered link carries an auth token, so a log line names the host alone and the transport error is
stripped of the request URL that `net/http` puts in it. The transport reads the proxy from the
environment, which is how a node with no route out reaches the source. Compression is off, so the file
is the bytes the server sent, and a transfer that receives nothing for five minutes fails rather than
hanging the build; a slow steady one never does.

Sources are transient. They live in the build workspace and go with it on success and failure alike.

### Staging and packing

`stageMetadata` validates `Spec.Runtime()` and `Spec.Manifest()` and writes them,
plus the recipe itself, to `<buildDir>/.cnt/`, *beside* the payload rather than inside it. That placement
is the whole trick: `mksquashfs` makes one archive root per source, so handing it
the payload dir and the `.cnt` dir as two sources yields `/cnt/<name>/…` and
`/.cnt/*.json` at the same level, without the payload ever containing a directory
the recipe did not create. There is no option to rename a source, so the staged
directory must already be called `.cnt` when `mksquashfs` sees it.

The recipe is staged **as it was fetched**. A template keeps its `{placeholder}`
tokens: the expansion goes to the workspace file the build executes and is never
embedded, so every variant of one template shares a recipe digest and is told
apart by `manifest.source.placeholders`, which is their only stored copy.

`deriveKeys` runs immediately before staging in every backend — the last point at
which everything a key depends on is known. A script build calls it after the
recipe has run, so every dependency it needed is installed and its exact identity can enter the
identity scheme. A Conda app derives its two keys directly from the exports it
already embeds.

`createSquashfs` then keys the payload (`stagePayloadKey`) and rewrites the staged manifest with it,
immediately before `mksquashfs` runs — the last moment the payload is final, and the same directory
`mksquashfs` reads, so the key describes what the archive will hold. The payload directory is walked
in Go with one worker per core the build was given. A build that cannot key its payload does
not pack. Only a script build is keyed; a Conda or definition build skips the step. The
payload key is not an input to identity or equivalence; see `internal/artifact`'s README.

Every `.def` build is keyed the same way, `os` included — nothing distinguishes
the configured default root from any other. Its `keys` block identifies the
definition plus the upstream image that definition bootstrapped from — both
known before Apptainer runs, which is what lets a definition build stage its
own keys into the sandbox it packs.

## The upstream digest

`bootstrap.go` reads `Bootstrap:`/`From:` out of a definition's header block — a
`From:` inside `%post` is shell text, not a directive — and `resolveUpstream`
asks `internal/registry` what that reference points at *now*, before the build.
The answer does two jobs: it becomes the `from=` field in the identity preimage,
and `writePinnedDef` rewrites `From:` to name the digest in the definition handed
to Apptainer. Without that pin an upstream retagged mid-build would leave the identity describing bytes the image does not contain. Only the transient copy is
rewritten; `/.cnt/recipe` keeps the definition byte for byte.

Resolution never fails a build. A registry that cannot be reached records
`unrecorded` and warns, because a login node behind a proxy must still be able to
build; a bootstrap with no upstream at all — `scratch`, `localimage`,
`debootstrap` — records nothing, and the two cases stay distinguishable.

A `scheme://` build has no recipe on disk, so the synthesized definition becomes
one: without that a `docker://ubuntu:24.04` image would carry no keys, when it is
exactly the two-line definition it is equivalent to. The generated header carries
a build date and drops out of the key, since whole-line comments never reach a
preimage.

**Its directives are byte for byte what Apptainer writes.** Apptainer synthesizes
the same definition for a bare `scheme://` build and stores it in the root at
`/.singularity.d/Singularity` — lowercase keys, trailing blank line. Matching it
means an image built here and a foreign one built from the same URI strip to one
recipe digest, so they share an equivalence key, and share identity too whenever
both resolved the same upstream digest. Capitalizing the keys or dropping the
blank line would split them for no reason a reader could see; every reader of a
def header is case-insensitive, so nothing else depends on the spelling.

`manifest.source.files` names the stored files from which the selected schemes
can regenerate their keys. `manifest.keys` stores the selected scheme names and
digests; no derived record files are staged. Validation ensures every required
source file is present in `/.cnt`.

A def build has no packing step, so `writeRecordingDef` appends a `%files`
section listing each staged file by name. It lists them individually rather than
naming the directory, since `%files` copies with `cp -a` and would nest the whole
directory inside a `/.cnt` the sandbox already has.

Validation runs at staging because that is the only point where bad metadata can
still stop the build: earlier there is nothing to validate, later the image
already exists.

### Prepared output

Every build writes to `Target.Prepared` — the target path plus the lock owner's
tag and `.part` — and `atomicInstall` renames it over `Target.Path`. Beside the
target, so the rename cannot cross filesystems; every build and not just an
update, because writing straight to the installed path would leave a truncated
image indistinguishable from a finished one.

The installed image is never removed first. A rename over an existing file is
atomic, so a reader sees the old image or the new one and never a gap; removing
first would open a window where the image is missing, and would destroy the
installed copy if the rename then failed.

The name is derived from the lock owner rather than recorded, which is how
`clearStaleLock` finds and removes the orphaned output of a killed build.

**Def:**
1. Check if overlay exists (skip if not updating)
2. If updating existing overlay: probe exclusive lock — fail immediately if in use
3. Create build lock
4. Derive recipe equivalence and try the selected source's ordered pull endpoints
5. If no acceptable prebuilt exists, build a sandbox with Apptainer, copy the staged metadata into it, and pack it to `.sqf`
6. Atomic rename prepared → target; remove lock

**The configured default root:**
Built the same way as any other `.def` — a `.sqf` that packs itself, since a
`.def` build never has `Spec.Base` set (`BuildGraph.resolveBase` skips every
`BuildTypeDef` node). `IsInstalled` searches every image path rather than one
target specifically for the build whose name matches `config.BaseRecipeName`,
so a default root supplied by a shared install is not rebuilt into the user's
own directory. Every other `.def` build's `IsInstalled` checks its own target
path only.

## Importing a foreign root

`FromForeignRoot` (`foreign.go`) is `-f`'s third source shape, beside a recipe
file and a URI: a `.sif` or an Apptainer sandbox directory someone already
built, with no Apptainer invocation in this path at all. It is a separate
constructor and a separate `buildForeign`, not a third case bolted onto
`FromExternalSource`/`buildDef` — the two share every generic stage (the lock,
the prebuilt check, key derivation, metadata staging, publish) but a def build
resolves its upstream and *then* asks Apptainer to fetch and build, where an
import has nothing left to fetch: the root already exists, and the identity
question is what it was built from, read after the fact rather than resolved
before.

**Identity comes from the root's own record, never a fresh resolve.** A
foreign root's `/.singularity.d/Singularity` carries the same `Bootstrap:`/
`From:` lines Apptainer synthesizes for a live `scheme://` build (see *The
upstream digest*), so `readForeignBootstrap` parses it with the same
`parseBootstrap` and the recipe hash matches a fresh build from the same URI.
The digest is `labels.json`'s `org.opencontainers.image.base.digest`, when
present, exactly the reasoning `resolveUpstream` already applies to a locked
rebuild: resolving the reference again would follow a tag that has moved since
this root was built. The label is genuinely optional — an older
Apptainer/Singularity build carries only the legacy `org.label-schema.*` set,
sometimes with no base-digest label at all — so its absence records
`meta.Unrecorded` rather than failing the import.

**Only `docker://`, `oras://`, and `library://` roots import.**
`allowedForeignBootstrap` refuses everything else (`shub`, `yum`, `zypper`,
`debootstrap`, `localimage`, `scratch`, and anything unrecognized) before any
key is derived. This is not a degraded identity tier for the rest — it is a
refusal, because none of them can be rebuilt from a fresh machine: `shub` is
dead, `localimage` names a path on the *original* builder's machine,
`scratch` has no upstream a `.sif` preserves, and `yum`/`zypper`/`debootstrap`
use a mirror-plus-version shape the identity model's single `From` string was
never built to represent. A live `--from shub://...` build fails at
Apptainer's own fetch before anything is published, so it needs no equivalent
guard; an import reads a root that already built successfully, so it has no
such safety net and must refuse explicitly instead.

**Extraction never copies the tree out first.** A sandbox source packs
directly, the same as a `.def` build's own sandbox does, since it is already a
plain host directory. A `.sif` mounts its primary SquashFS partition read-only
with `squashfuse -o offset=N` (`packFromSIF`), inside `freeze.MountedRun`'s
unprivileged namespace — the partition's own byte offset stands in for the
whole-file offset a `.sqf` mount uses, so this is the same primitive
`internal/image/freeze` already uses to read a frozen artifact, not a new one.
Either way, the final `mksquashfs -all-root -no-xattrs` repacks everything: a
`.sif` extracted by dumping the partition and appending the staged metadata
would be cheaper, but would leave every file's original ownership untouched,
and a real service image (not a bare OS base) can carry non-root ownership
that `unfreeze`'s single-identity namespace mapping cannot handle later (see
*Unfreeze can't get correct ownership from the mount — only after it* in
`internal/image/freeze/README.md`). A full repack is what actually
renormalizes it.

**Arch is checked against the SIF's own header, not a label.** `sif.Arch`
reads the binary architecture field every SIF carries, always present unlike
the self-reported `org.label-schema.build-arch` label. A mismatch against
`runtime.GOARCH` refuses the import outright rather than packing an artifact
that would report itself as native and never run.

## The container root

Script and Conda builds run their *install* step inside a container, so each
has an implicit edge to the container root — a build prerequisite, separate
from anything a recipe declares with `#DEP:`. `BuildGraph.Run` resolves it
once before any node runs, and records it on every dependent as `Spec.Base`;
`ResolveBase` builds the configured default root first if none is installed,
resolving it through the exact same catalog/build path — `NewBuildObject` +
`NewBuildGraph` + `Run` — every other name takes. There is no dedicated build
path for it: it is an ordinary `.def` recipe, `catalog.TypeOS` like any other
— `catalog.DeriveType` never derives a distinct type from a `.../base` name.
A project may still name its default root `.../base` for readability; that
name carries no special meaning to the type system, only to
`config.BaseRecipeName`.

A locked rebuild (below) is the second caller of the same convention, and the
one that actually matters inside a project: `internal/project/restore` sets
`LockedSpec.Base` from the lock's reserved base pin, once that pin's own step
has produced a local path, so a conda or script rebuild is never silently
resolved against whichever `default_distro` happens to be configured on the
machine running the restore. See `plan/base-provenance.md`.

A `.def` build has no `Spec.Base` of its own: `BuildGraph.resolveBase` skips
every node whose `BuildType` is `BuildTypeDef`, since a definition bootstraps
its own root and resolves nothing.

**Packing (`createSquashfs`, `squashfs.go`) never runs in a container at
all, for any workspace mode.** A sandbox and a dir-mode payload are already
plain host directories — `mksquashfs`, resolved via `toolpath.Resolve` and
run directly on the host, reads them exactly as it would through a bind
mount, since a bind mount contributes nothing a direct read doesn't already
see.

**The reason a Script/Conda build's *install* step still runs inside a
container is not tool availability — it's that `$CNT_PREFIX` is baked into
what gets installed.** Shebangs, activation scripts, and some RPATHs are not
relocatable, so the install has to happen at the exact absolute path the
artifact will be mounted at later; the container root is what supplies that
path. Conda installs (`micromamba`) run through `internal/libexec`'s
self-provisioned toolchain, not whatever the root happens to carry, so a
`.def` build no longer scans its sandbox for it. Packing needs no such thing
— an archive's content does not depend on what process read it off disk, and
neither do `captureMicromambaVersion`/`condaExport`'s post-install
diagnostics, which read the same host path or `fuse2fs` mount packing itself
uses (`squashfs.go`) rather than running through Apptainer.

Every `.def` build — any of them may end up chosen as someone's root, since
root selection is a per-invocation runtime choice, not a declared type —
still checks its sandbox for `/bin/bash` before it can be used as one: every
`execpkg.Run` a later script or Conda install issues launches that exact
path, so a sandbox carrying bash only elsewhere runs nothing a build inside
it asks of it. This check has nothing to do with packing the `.def` build's
own sandbox, which needs no container and so needs no `/bin/bash` of its
own. `/bin/bash` is the only thing asked of the sandbox: nested mounting,
packing, and Conda installs all run through `internal/libexec`'s
self-provisioned toolchain, never whatever the root happens to carry.

**The question is put to the container, not to the sandbox directory.**
`command -v` inside it answers with the PATH `/bin/bash` will actually run
under. Testing for a file under the sandbox's own directories would ask
something else and could get it wrong: a sandbox is free to put its tools
somewhere a file-existence check never guessed.

## BuildGraph Execution

1. **Solve** - `catalog.Resolve` walks the graph from the index alone: transitive
   deps, cycle detection and dependency-first order, without fetching a recipe.
   Installed versions come from the `Have` callback, so an installed dep is a map
   lookup rather than a build object.
2. **Plan** - Create a BuildObject per missing node, in the solved order
3. **Decide prebuilts** - In solved order, each node derives its equivalence from its recipe,
   placeholders and what its dependencies contribute: an installed dependency's own keys, or, for
   one still to be acquired, the key derived for it earlier in this pass (`plannedDeps`), so no
   dependency has to be installed to decide. A node whose source lists pull endpoints is then looked
   up at them by metadata alone (`planPrebuilt`). An artifact made from the same recipe makes it a
   pull; anything else makes it a build, and one made from a different recipe is warned about. A
   dependency neither installed nor planned — a path, a constraint nothing satisfies — leaves the
   node undecided, and it looks when it builds. A Conda build, a source with no endpoints and
   `--no-prebuilt` are decided as builds without a lookup. It is decided here so the plan can say
   what will happen, and so a pull is never queued behind a scheduler job.

   **A pull opens none of its build dependencies**, so `pruneDependencies` drops the nodes only
   pulled nodes need and the plan lists them as not installed; the user adds them later with a
   `create` of their own. A dependency a root or a node that builds needs stays. A pull that falls
   through to a build acquires them then, which is why `buildScript` tries a planned pull before the
   toolchain and the dependencies, and why the local step passes `buildDeps` for such a node.
4. **Resolve the base** - Every script and Conda build runs inside it, so it is
   built first; a scheduler job cannot build one on the node. A pull, and a node dropped for one,
   does not need it; a pull that falls back to a build resolves it then
5. **Separate local/scheduler** - A pull is always local. A build goes to the scheduler when its recipe carries directives,
   or when it is `data` and `build.always_submit_data` is set. A job re-runs `create`, so only a build
   that command can reproduce is submitted: a catalog name, or an external shell script, re-run as
   `create --name|--prefix … --file …` (`SetJobArgs`). A Conda build from packages or a file, a
   definition and a `.sif` or sandbox import run here.
6. **Execute:**
   - Local builds and pulls: sequentially, in dependency order
   - Scheduler builds: submitted with dependency chains, each waiting on its deps. A job re-runs
     `condatainer create <name>` on the node, so the flags that change what is built or where it
     lands travel with it: `--channel`, `--source`, `--layer`, the block sizes and compression
     flag (`BuildGraph.SetJobFlags`), and `--update`, `--store` (only on the build that was asked
     for, never a dependency) and `--no-prebuilt`, which is also set on a node planning decided to
     build, so the node need not look again. Submission flags are not repeated: the job is
     the submission.

## Environment Variables

Images can export environment variables via `#ENV:` directives:
```bash
#ENV:CELLRANGER_ROOT={prefix}
```

These are captured into `Spec.Image.Env` at resolution and embedded in the
image's `/.cnt/runtime.json`. `{prefix}` survives into the image verbatim and is
substituted with the recorded `prefix` at load time (see
`internal/runtime/container` env handling).

## Workspace Strategy

`workspaceFor` derives every path from `(name, tmp root, scratch extension,
is-definition, producer identity)`. Each producer owns one directory:

```text
<tmp-root>/build_<name>/<local-host-pid|scheduler-jobid>/
  cnt--<name>.sh|.def
  rootfs/             sandbox, a definition only
  work/
    cnt/
    tmp/
    .cnt/
```

**`isDef` is passed, not inferred.** The recipe's extension and the sandbox both
follow it, and nothing in the workspace's own paths distinguishes a definition:
every build produces a `.sqf`, so an inference from the output would name a
definition's recipe `.sh`.

The recipe is materialized in a process-private directory before the target
lock because resolution and scheduler planning need to read it. Acquiring or
adopting the lock re-sites that directory to the lock owner's tag. Generated
definition helpers (`cnt-synth.def`, `cnt-pinned.def`) are written inside the
same private directory, so builds of different targets cannot overwrite each
other.

The final `.part` remains beside the installed target for atomic rename, but it
uses the same owner tag. Cleanup removes only the current owner's directory;
stale-lock cleanup reconstructs the stale owner's directory from the lock.

Every payload is written to host directories: a script or Conda build's under `work/cnt`, bound
into the container, and a definition's in the sandbox apptainer writes. The pack reads either
directly on the host — no container, no bind. A sandbox is unwrapped straight into the archive root
(`packSources`'s `keepAsDirectory=false`) rather than nested under the sandbox directory's own name.
There is no scratch image: a build that must keep its small files off a quota-limited filesystem
points `$CNT_TMPDIR` at local scratch instead.

Each build type also uses a different base directory for build artifacts:

| Build path | root | picked by |
|---|---|---|
| Conda (`name/version`) | fast | `tmpRootForType` |
| Script, type `app` | fast | `tmpRootForType` |
| Script, type `data` | stable | `tmpRootForType` |
| Def (internal), base image | fast | `tmpRootForDef`, applied by `asDefinitionBuild` |
| External `-f`, type `app` | fast | `tmpRootForExternal` |
| External `-f`, type `data` | `filepath.Dir(targetPrefix)` | `tmpRootForExternal` |
| External `-f`, `.def` | fast | `tmpRootForExternal`, via `tmpRootForDef` |

The **fast** root is `utils.GetTmpDir()`: `$CNT_TMPDIR` → scheduler scratch →
`$TMPDIR` → `/tmp`, always plus `cnt-$USER`. The **stable** root is
`config.GetWritableTmpDir()`: the first writable `<data-dir>/tmp`, falling back
to the fast root when no data directory is writable at all.

A definition build takes the fast root because it must. Apptainer builds under
`--fakeroot`, which NFS, Lustre, GPFS and PanFS do not support, and a data
directory on an HPC system is routinely one of those. `tmpRootForDef` warns
through `utils.WarnUnfakerootableScratch` — the build will fail, not merely run
slowly.

**The constraint lands on the workspace, not on `APPTAINER_TMPDIR`.** Apptainer
assembles a sandbox next to where it is going — `filepath.Dir(dest)/build-temp-*`,
renamed into place — and consults `TmpDir` only when the format is *not* sandbox.
So `ws.Root` is what has to hold the tree and support the ownership changes
`--fakeroot` performs. Apptainer has a fallback if it cannot: it warns, builds in
the temporary directory instead, and copies. That path is the degraded one, and
choosing the workspace root well is what avoids it.

`APPTAINER_TMPDIR` is still set, to `ws.BaseRoot` — the fast root itself, which
`$CNT_TMPDIR` selects. Left unset it inherits `TMPDIR`, which on a scheduler is
routinely the network scratch this whole section is avoiding, so one knob has to
answer for both or the two halves of a build disagree. It is the root
`ensureWorkspaceRoot` already creates, so nothing is made for Apptainer's sake and
Apptainer cleans up the children it puts there.

`$CNT_TMPDIR` moves the fast root and nothing else. It used to short-circuit
`GetWritableTmpDir` as well, which switched off this whole table: exporting it to
speed up a conda build silently moved the next data build onto node-local scratch
that the job wipes. An external `data` build keeps its intermediates beside the
target for the same reason — the user picked that location, and a multi-GB
payload is the thing least able to survive a scratch quota. An external `.def`
cannot follow it there: the sandbox has to be somewhere `--fakeroot` works, and
where the user put the target says nothing about that.

`utils.GetTmpDir()` priority: scheduler-assigned scratch (`SLURM_TMPDIR`, `PBS_TMPDIR`, `LSF_TMPDIR`, `_CONDOR_SCRATCH_DIR`) → `TMPDIR`/`TEMP`/`TMP` → `/tmp/cnt-$USER`.

## Resource Allocation for Builds

CPU and memory derived from `scriptSpecs`:

| State | `HasDirectives` | `Spec` | CPU source | Memory source |
|---|---|---|---|---|
| No directives | false | non-nil | `Build.Defaults.CpusPerTask` | `Build.Defaults.MemPerNodeMB` |
| Has directives | true | non-nil | `CpusPerTask × TasksPerNode` | `MemPerNodeMB` |
| Passthrough | true | nil | build fails with error | — |

Build defaults are configured via `build.*` keys in `config.yaml`. See [Config README](../config/README.md).

## Error Handling

- `ErrTmpOverlayExists` - Temporary overlay already exists (stale from previous run)
- `ErrBuildCancelled` - User interrupted build (Ctrl+C)
- `"build lock found for <name> ..."` - Lock exists and is active (job pending/running or local build in progress)
- `"build already queued or running for <name>"` - Duplicate scheduler submission blocked by submit-time lock
- `"cannot update <name>: ..."` - Target overlay is locked (currently used by a running `exec`/`run`)

**A missing self-provisioned tool fails in Go, before any work starts.**
`micromambaCmd` and `condaExecOpts` (`conda.go`) return `libexec.ErrNotProvisioned` before the
install container starts: micromamba has no other legitimate source, so there is nothing to fall
back to and no reason to spend a container launch finding that out. `createSquashfs`
(`squashfs.go`) resolves `mksquashfs` the same way, via `toolpath.Resolve`, before rendering any
script — packing runs on the host now, never in a container, so there is no "the container root
might already carry it" case left to fall back to either.
