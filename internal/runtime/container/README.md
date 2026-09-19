# container

Container setup pipeline providing overlay resolution, bind path deduplication, environment collection, and GPU detection for Apptainer execution.

## Architecture

```
setup.go     Main setup pipeline, SetupConfig/SetupResult
resolve.go   Overlay path resolution (name → absolute path)
bind.go      Bind path deduplication and validation
env.go       Environment variable collection from .env files
gpu.go       GPU detection (NVIDIA, AMD) and flag generation
path.go      Container path utilities
```

## Key Types

**SetupConfig** - Input configuration:
- `Overlays` - Overlay paths or names (resolved automatically)
- `WritableImg` - Whether .img overlays are writable
- `EnvSettings` - User environment variables (`KEY=VALUE`)
- `BindPaths` - User bind paths
- `Fakeroot` - Use fakeroot
- `ApptainerFlags` - Additional flags

**SetupResult** - Processed output ready for execution:
- `Root` - The exec root pulled out of `Overlays`, `""` if none requested (see Root selection)
- `Overlays` - Resolved ordered overlay paths, `Root` excluded
- `OverlayArgs` - Paths with `:ro/:rw` suffixes
- `EnvList` - Complete environment variable list
- `EnvNotes` - Environment descriptions for display
- `BindPaths` - Deduplicated bind paths
- `Fakeroot` - Final fakeroot setting
- `ApptainerFlags` - Flags including GPU detection
- `LastImg` - Path to the writable `.img` overlay (if present)
- `EnvMounted` - A conda environment is mounted at `/cnt_env`, in any form:
  `LastImg`'s `.img` or a read-only env-typed `.sqf` mounted alone

## What Setup refuses

Checks run before anything is locked or mounted, because each describes a
container that cannot be what the caller asked for:

- `ensureSingleImage` — at most one writable `.img`. Only one can be the
  principal image, and only that one takes a write lock.
- `ensureAtMostOneSif` — at most one `.sif`. A `.sif`'s only valid use is the
  exec root (see Root selection), and two of them cannot both be it.
- `ensureDistinctPrefixes` — no two images claiming one `/cnt/<name>` subtree.
  Overlays are disjoint subtrees, not stacked diffs, so two claiming one prefix
  do not combine: the later mount takes the subtree and the earlier contributes
  nothing, while both still reach PATH and the environment. Two builds of one
  name are the case this catches — a project's restored copy and a flat install
  of the same name record the same prefix. An OS image and anything without
  readable metadata record no prefix and are exempt for free; the same file
  named twice is redundant, not a collision.

## Root selection

There is no `-b`/`--base-image` flag: every image a command wants is named the
same way, through `Overlays`. `Setup` pulls the root out of the plain
requested order, before `orderOverlays`' os/app/data layering runs on
whatever is left (see Overlay ordering below). A `.sif` wins unconditionally
when present — `ensureAtMostOneSif` already guarantees there is at most one,
and root is its only valid use, so it wins regardless of where it falls among
the requested overlays. Otherwise the first `os`-typed entry, in the order
requested, wins. `TypeEnv` is excluded even though it merges at root the
same way, since an environment's identity presupposes a chosen root and so
can never supply one.

Found, that entry is pulled out of the list into `SetupResult.Root` and run
as the exec root instead of getting its own `--overlay` mount. Environment
and PATH collection still runs over the *full* requested list (`Root`
included): becoming the exec root changes which Apptainer flag carries an
overlay, not what it contributes — `os` never puts anything on PATH
regardless (see Environment Variables below), so in practice this only
matters for a root that declares its own `#ENV:`.

Not found, `SetupResult.Root` is `""` and the caller supplies the fallback:
`Setup` cannot build a missing configured default itself (`internal/build`
would have to import this package's caller to do that, an import cycle), so a
caller that wants one built calls `HasRequestedRoot(overlays)` *before*
`Setup` runs, and only resolves/builds the configured default
(`internal/build.ResolveBase`) when it answers false — building it
unconditionally would be wasted work, and a needless failure, whenever the
request already names its own root. `exec.Options.resolveBaseImage` is where
`Root` and a caller-supplied `BaseImage` are reconciled: `Root` wins when
present, otherwise the caller's `BaseImage` is kept, falling back to
`config.GetBaseImage()` (found, never built) as a last resort for callers that
never call `HasRequestedRoot` at all (`internal/build`'s own container
invocations, which always set `BaseImage` explicitly and are unaffected by
this scan either way).

## Important Diff from Apptainer Flags

**Writable** in CondaTainer means making the ext3 `.img` overlay writable, not adding `--writable` to Apptainer. The `.img` overlay is writable by default when used as an overlay.

By default, CondaTainer will add `:ro` suffix to all overlays for safety. If `WritableImg` is true, the `.img` overlay will have `:rw` to allow writing, but the `.sqf` overlays will still be `:ro`.

## Usage

### Setup Pipeline

```go
cfg := container.SetupConfig{
    Overlays:    []string{"cellranger/9.0.1", "custom.sqf"},
    WritableImg: true,
    EnvSettings: []string{"DEBUG=1"},
    BindPaths:   []string{"/data"},
    Fakeroot:    false,
}

result, err := container.Setup(cfg)
if err != nil {
    // Handle errors
}

// Use result for execution
opts := &apptainer.ExecOptions{
    Overlay:    result.OverlayArgs,
    Env:        result.EnvList,
    Bind:       result.BindPaths,
    Fakeroot:   result.Fakeroot,
    Additional: result.ApptainerFlags,
}
```

### Overlay Resolution

```go
// Resolve overlay names to absolute paths
// Searches: absolute path → images dirs
overlays, err := container.ResolveOverlayPaths([]string{
    "cellranger/9.0.1",        // → /path/cellranger--9.0.1.sqf
    "/abs/path/custom.sqf",    // → /abs/path/custom.sqf
    "myenv.img",               // → /path/myenv.img
    "./tmp/alpine.sif",        // → /abs/tmp/alpine.sif
})
```

### GPU Detection

```go
// Auto-detect GPUs and generate flags
gpuFlags := container.DetectGPUFlags(requested)
// requested=false: skipped when autoload_gpu is disabled
// requested=true: detected regardless of autoload_gpu (script declared a GPU need)
// Returns: ["--nv"] for NVIDIA, ["--rocm"] for AMD, [] otherwise
```

### Fakeroot Management

```go
// Auto-enable fakeroot for writable .img
fakeroot := container.AutoEnableFakeroot(
    lastImg="/path/env.img",
    writableImg=true,
    currentFakeroot=false,
)
// Returns: true (auto-enabled)
```

## Setup Workflow

1. **Resolve overlays** - Name/path → absolute paths
2. **Validate** - At most one `.img`, at most one `.sif`, no two overlays claiming one prefix
3. **Select root** - Pull it out of the plain requested order into `Root` (see Root selection)
4. **Order** - Layer what's left: os, then app/data, then a paired env snapshot, then `.img` last (see Overlay Ordering)
5. **Lock check** - Verify .img availability (exclusive lock if writable)
6. **Build overlay args** - Add `:ro/:rw` suffixes
7. **Collect environment** - Resolve `#ENV:` for each overlay (embedded build script, sidecar `.env` overrides), `Root` included
8. **Deduplicate binds** - Remove conflicting bind paths
9. **Detect GPU** - Add `--nv` or `--rocm` if available; forced past `autoload_gpu:false` when the caller requested one
10. **Return result** - Ready-to-use configuration

## Overlay Ordering

`orderOverlays` runs on what's left after root selection has already pulled
`Root` out (step 3 above), so it never decides which overlay becomes root —
only how the rest stack once mounted:

1. `os` overlays, in the order requested
2. `app`/`data` overlays, in the order requested
3. the one `env`-typed `.sqf` present (autoloaded or explicit)
4. a writable `.img`, always last — required for it to work as the writable layer

## Environment Variables

Each image's environment comes from its embedded runtime document
(`/.cnt/runtime.json`), with `{prefix}` resolved to the image's install prefix at
load time. Nothing is mounted and no sidecar travels with the image.

That document is the *only* metadata this path reads. The manifest beside it
carries provenance and is never opened here, so provenance can grow without
costing every mount.

An image built for another architecture is refused the same way — mounted,
contributing nothing, with a warning. `compare.MountAllowed` is the one
comparison rule on this path, and only because `runtime.json` is already in hand:
it costs a string comparison, and without it a wrong-architecture image mounts
cleanly and fails somewhere downstream where the cause is unrecognizable. An
artifact is portable only where its recipe said `#ARCH:noarch`.

A second comparison rides along for the same reason, and warns rather than
degrades: inside an image directory the filename *is* the address `-o` resolves,
so an image whose recorded name is not the one its filename encodes is listed and
mounted under a name nothing can look it up by. It still contributes normally —
the prefix and environment come from `runtime.json`, so only the address is
wrong. The check stops at the image directories: a path handed in directly
addresses no name, and a `store/` entry is verified by `store.Scan` instead.

An image with no readable runtime document still mounts and contributes nothing —
no variables, no `PATH` entry, no description. There is no fallback to the
manifest. `resolveImage` reports that once per invocation: informational when the
document is simply absent (an image built before the format, or a plain Apptainer
`.sif`), a warning when one is present but unreadable, or when the host is
missing `unsquashfs`.

A writable `.img` is the exception. It carries no embedded metadata, so it reads a
`<overlay>.env` sidecar instead — one `KEY=value` per line, with an optional
`##` note, and `{prefix}` resolved to `/cnt_env`:

```bash
GOROOT=/cnt_env/go   ## Go Installation path
PATH={prefix}/bin:$PATH
```

A writable `.img`'s `Contribution.Type` (`imgContribution`, `env.go`) is honestly `catalog.TypeEnv`
— the same type an env-typed `.sqf` snapshot reports from its own `runtime.json` — not a `TypeApp`
mislabel. `Contribution.ContributesBin()` is the one predicate `BuildPathEnv` and `ActivationScript`
both call instead of each duplicating a `Type != catalog.TypeApp` check: true for `TypeApp` or
`TypeEnv` with a non-empty `Prefix`, false for everything else (`TypeOS`/`TypeData`, or a degraded
zero-value `Contribution`). Both callers additionally deduplicate by `Prefix` — a writable `.img`
and its paired env-typed `.sqf` snapshot both claim `EnvPrefix`, but Apptainer merges them into the
one physical `/cnt_env` directory at mount time, so each needs exactly one `PATH` entry and one
`activate.d` block, not one per overlay that claims the prefix.

`CollectOverlayEnv` merges every overlay's resolved env into the final list;
`ResolveOverlayEnv` returns a single overlay's description/env/notes for `info`.
When the same variable is set by more than one overlay the later one wins and a
diagnostic is recorded. `CollectOverlayEnv` is the one place degradation is
reported, so an image with no readable metadata is announced once per invocation
rather than once per consumer.

A writable `.img` paired with an env-typed snapshot (`LookupSnapshot`, see
Environment Snapshots below) is not read in isolation by `ResolveOverlayEnv`:
the snapshot's variables are merged in first, the `.img`'s own sidecar on top,
with the same "defined in multiple overlays" diagnostic `CollectOverlayEnv`
already produces for the ordinary multi-overlay case. Without this, `info` on a
`.img` would show none of the variables its paired snapshot silently
contributes at real mount time.

`BuildPathEnv` is silent for the same reason. Anything `ContributesBin()` (an `app`, or a mounted
conda environment — writable `.img` or bare env-typed `.sqf` alike) contributes `<prefix>/bin`
unconditionally — `data` and `os` put nothing on `PATH`, and an image with no readable metadata
contributes nothing at all. There is no check that the directory exists: a nonexistent `PATH` entry
is harmless, while the check would cost one archive read per image per invocation.

**Common Environment:**
- `LC_ALL=C.UTF-8`, `LANG=C.UTF-8`
- `CURL_CA_BUNDLE=`, `SSL_CERT_FILE=` (unset to avoid host interference)

## Activation

`#ENV:` (an `app`'s `runtime.json`) and `.env` sidecars (a writable `.img`)
are recipe-declared and static — they cover the two things `container.Setup`
itself computes ahead of time. A conda-forge package's own
`etc/conda/activate.d/*.sh` is neither: it's shell script the package ships,
and real `conda activate` sources it, not condatainer. `ActivationScript`
replays that one piece — nothing else `conda activate` does, see below.

`ActivationScript(overlays, envMounted, mode)` builds a bash preamble, gated
by `mode` (`ActivationAll`/`ActivationEnv`/`ActivationNone`,
`exec.Options.Activation`, CLI `--activation`): `ActivationAll` sources one
block per overlay that `ContributesBin()` (every `app`, plus the mounted
conda environment's own), deduplicated by `Prefix`; `ActivationEnv` sources
only the `/cnt_env` block, reading `envMounted` directly rather than walking
`overlays` at all, since it deliberately looks at nothing else; `ActivationNone`
sources neither, so a hung or misbehaving activation script can be ruled out.
`envMounted` is true for a conda environment mounted in any form — `Setup`'s
`LastImg` (a writable `.img`) or a read-only env-typed `.sqf` mounted alone
(see Environment Snapshots) — not just `LastImg` specifically;
`AutoEnableFakeroot` is the one place that stays keyed on `LastImg` alone,
since fakeroot is about a writable overlay's UID mismatch, meaningless for a
read-only `.sqf`. Each script sources as
`CONDA_PREFIX=<prefix> . "$script"` — scoped to that one call, not exported,
so a script reading `CONDA_PREFIX` (conda-forge's `libxml2` hook does, for
`XML_CATALOG_FILES`) sees its own overlay's prefix and nothing lingers into
the command this preamble `exec`s. `/cnt_env`'s `CONDA_PREFIX` for that
command comes from `container.Setup`'s env list instead. App blocks run in
overlay order, the same order `BuildPathEnv` iterates before its own prepend
reverses it, so the last-listed overlay's hook still wins a name collision —
one ordering rule, not two.

`exec.Prepare` rewraps `Options.Command` into one `bash -c` invocation running
the activation script then `exec "$@"` into the original command, since this
is shell script sourced into the running shell — apptainer's own `--env` is a
static list and cannot express it. Skipped when `ActivationScript` returns
`""`, which it does for `ActivationNone`, and otherwise whenever nothing
mounted could contribute an `activate.d` directory — the common case for an
`os`-only container under `ActivationAll`.

`MMHelperScript(envMounted)` appends one more bash function to the same preamble: `mm`, wrapping
`condatainer env "$@"` and, for `install`/`update`/`remove`, chaining `eval "$(condatainer env
reactivate --shell bash)"` — `--shell bash` is hardcoded rather than left to `reactivate`'s own
`$SHELL` detection, since `mm` only ever runs inside bash regardless of what `$SHELL` says — so a
package's `activate.d`/`deactivate.d` side effects (an env var like
`JAVA_HOME`) refresh in the same shell instead of staying stale until the user exits and
re-enters. It is gated on `envMounted` directly, not on `Activation` mode — `mm` is condatainer's
own script, not a third-party conda-forge hook, so `ActivationNone` (for ruling out a misbehaving
*conda-forge* activate.d script) has no reason to remove it too. `export -f` is what lets the
function survive the subsequent `exec` into an interactive bash session; it does not survive
`exec` into zsh or fish (confirmed empirically), so `mm` only ever helps a bash session — zsh and
fish both keep working via the plain `condatainer env ...` command, they just don't get the `mm`
name or the automatic reactivate chaining. No file is ever materialized for this: the function
definition is just more text in the same Go string `exec.Prepare` already builds, so there is
nothing to clean up if the process is killed mid-command.

**Only the first-activation branch, once, no restore.** Real `conda
activate`/`deactivate` (`activate.py`'s `build_activate`/`build_deactivate`)
also reads two more env-var sources (`conda-meta/state.json`'s `env_vars`,
`etc/conda/envvars/*.json`) and maintains a `CONDA_SHLVL` stack —
`CONDA_PREFIX_<n>`, `CONDA_STACKED_<n>`, and every clobbered variable backed
up as `__CONDA_SHLVL_<n>_<name>` — so a later `deactivate` can undo exactly
what an `activate` changed. A condatainer run mounts, executes once, and
exits: there is no nested activate/deactivate within one invocation, so none
of that stack has anything to undo. The two extra env-var sources are left
unread deliberately, not by oversight — the mechanism is rare in practice
(this project's own `libexec` toolchain has exactly one `activate.d` script
across four packages and zero `conda-meta/state.json`) and reading them would
need a JSON parser with no guaranteed one in the base image, for a source
`activate.d` already covers for everything that matters today.

## Automatic Bind Paths

`BindPaths()` in `bind.go` collects bind mounts automatically before the user-specified paths are appended:

| Source | Condition |
|--------|-----------|
| Current working directory | Always |
| `$SCRATCH` | When env var is set |
| `$TMPDIR` | When env var is set and the path exists |
| Scheduler node-local tmp | When inside a job exposing one (e.g. `SLURM_TMPDIR`) |
| Base data directories | From config; added `:ro` if not writable |
| `condatainer` executable | For nested calls (always; bound to `/.cnt_bin/condatainer`) |

After collection, `DeduplicateBindPaths()` removes conflicting bind paths:
- Keeps longest/most specific paths
- Removes parent paths when child is bound
- Example: `/home/user/data` removes `/home/user`

**One more bind, added by `Setup` itself rather than `BindPaths()`:** when a conda environment
is actually mounted (`EnvMounted`, `.img` or a bare read-only env-typed `.sqf` alike), `Setup` also
binds `libexec.Dir()` — the resolved self-provisioned toolchain directory, same
host-path-equals-container-path convention as everywhere else `libexec` is bound. This is what lets
`internal/conda`'s in-container `mm`/`env` commands resolve `micromamba` via `toolpath.Resolve` once
running (see `internal/conda/README.md`) instead of assuming the base image carries one. Gated on
`EnvMounted`, not unconditional, because nothing else in an ordinary exec/run needs it. The one
other trigger is `SetupConfig.BindLibexec`, set by the caller when `nested_run` is providing apptainer
from `libexec/` (see `cmd/nested.go`): a container that runs apptainer needs its binary and libraries
at the host path.

## Nested Running

`nested_run` (`auto`, `true`, `false`) decides whether a container can start containers of its own.
The decision is made by `cmd/nested.go`, not `Setup`: the caller adds the provider to the overlay
list and sets `SetupConfig.BindLibexec`, and `Setup` treats both like any other input.

Providers, in order: the apptainer installed in `libexec/` (bound, never an overlay), else the newest
installed `apptainer/<version>` overlay, else — only under `true` — an overlay built on the spot.
`auto` builds nothing and stays silent when there is no provider. It is the same at any depth: each
launch is a fresh container, so one started from inside a container is provided for like any other.

- **The overlay goes first in the list.** `BuildPathEnv` prepends each overlay's `bin/` in list
  order, so the first overlay ends up last on `PATH`. The overlay is a full conda prefix (`jq`,
  `openssl`, `bsdtar`, …); anywhere else in the list it would shadow the same-named tools of the
  overlays the user asked for.
- **The build runs before a job is submitted.** `run` builds a missing overlay ahead of the scheduler
  block, on the host that ran the command, where the network is; the job only looks the overlay up.
  It builds locally whatever `build.always_submit` says, as the default base image does.

## GPU Detection

Flags are added from the one device node unique to each vendor:

- **NVIDIA**: `/dev/nvidiactl` → `--nv`
- **AMD**: `/dev/kfd` → `--rocm`

`/dev/nvidiactl` rather than `/dev/nvidia0`, which is absent on MIG nodes and
whenever the allocated GPU index is not 0. `/dev/kfd` rather than `/dev/dri`,
which is DRM and present for any vendor — an NVIDIA-only node has a populated
`/dev/dri` too, so it would not tell AMD apart from anything else.

A stat says a driver is installed, not that it works. On a node where the
driver is loaded but a GPU is unusable, the node is still there, so detection
fires and `--nv` then fails at container creation — `nvidia-container-cli`
enumerates every visible GPU through NVML before the container exists, and one
bad handle aborts the whole thing. `autoload_gpu: false` skips detection
entirely for that case. `--rocm` has no equivalent helper: it binds libraries
and devices, so a sick AMD GPU surfaces inside the program instead.

`SetupConfig.GpuRequested` overrides that toggle: a script that declared a GPU
requirement (`condatainer run`'s resolved `-g`/scheduler GPU allocation, or a
helper's `#GPU:` header) still gets `--nv`/`--rocm` when the host has the
device node, regardless of `autoload_gpu`. The toggle only ever suppresses
detection nobody asked for — it was never meant to silently drop GPU access
from a workload that explicitly needs one.

## Environment Snapshots

`overlay freeze` produces an `env`-typed `.sqf` beside the writable `.img` it
came from — a "snapshot." `LookupSnapshot` (`snapshot.go`) is the one place
that pairs one back up with an `.img`, reused in both directions: `Setup`
autoloads a snapshot at mount time, and a bare `overlay freeze` asks the same
question in reverse to find what it is replacing.

Naming is derived from the `.img`'s own basename, not hardcoded to `env`:
strip a trailing `-<user>` suffix (if present) to get a stem, then look for
`<stem>-<user>.sqf` before `<stem>.sqf`. The first candidate that exists on
disk decides the outcome — if it is `env`-typed it is the pair, and if it
isn't, the lookup reports `Blocked` rather than falling through to the next
candidate. Falling through would be guessing which of two files, if either,
was meant; `Setup` treats `Blocked` as nothing found (the `.img` mounts
alone), while `overlay freeze` treats it as a hard refusal, because a bare
freeze silently replacing an unrelated file at the derived path is the "wrong
container that looks like a working one" class of mistake, not a warning.

`Setup` wires this in before `ensureDistinctPrefixes` and the final ordering
run, so an autoloaded snapshot participates in both exactly as an
explicitly-listed one would:

- **Collision check.** An `.img` never claims a prefix for collision purposes
  — it has no identity of its own, `ensureSingleImage` already guarantees
  there is at most one, and it is now expected to sit on top of whatever
  `env`-typed `.sqf` is present. `env.sqf` + `env.img` is the expected shape.
  Two `env`-typed `.sqf`s together is still refused: that is still two
  snapshots with no way to tell which one is meant, detected the same way as
  any other prefix collision.
- **Ordering.** `orderOverlays` (see Overlay Ordering) places the writable
  `.img` last and, immediately beneath it, the one `env`-typed `.sqf` present
  — regardless of where either appeared in the original list. The reason
  it matters at all: a directory created with no lower counterpart is
  opaque-marked by `fuse-overlayfs` regardless of intent, so the snapshot has
  to already be mounted underneath by the very first write a fresh `.img`
  makes (e.g. `overlay create`'s conda init) — wiring autoload in only for
  `exec`/`run` would let that first write happen with nothing beneath it, and
  the opaque flag does not get retroactively cleared once the snapshot shows
  up afterward.

**A missing `.img` named directly is refused, never silently substituted.**
`ResolveOverlayPaths` (`resolve.go`) refuses a path that does not exist on
disk the same way it always has — naming a specific `.img` is a request for a
writable mount by that exact name, and guessing that a read-only `.sqf` beside
it was meant instead is exactly the silent substitution this codebase's
refusal-over-guessing stance rules out. `missingOverlayError` only changes
what the refusal says: when the named `.img` would have paired with a real
snapshot (`LookupSnapshot`), the message names it and the two ways to
proceed — `overlay create` to continue from it, or mounting the `.sqf`
directly, read-only.

This is also the third overlay state every naive two-state consumer has to
tell apart from "nothing has ever been created here": no `.img`, but a
fully-populated snapshot right beside where one would go. `LookupSnapshot`
answers it without the `.img` existing at all — it only needs the name to
derive the sibling `.sqf` candidates. `FindEnvSnapshot` builds on this
directly; `ResolveEnvOverlay` treats `.img` and `.sqf` as two forms of one
environment overlay and resolves to whichever exists, `.img` first — never
creating either. Both live here rather than in `internal/helper`, which
otherwise owns every helper-facing wrapper (`helper.FindEnvSnapshot`,
`helper.ResolveEnvOverlayInDir`, kept as thin delegates to these): `internal/helper`
imports `internal/project`, so `internal/project` reusing the same
resolution (checking whether a frozen `env.sqf` is pinned) would cycle back
through `internal/helper` if the logic stayed there. Nothing is ever
fabricated on the strength of
a snapshot alone: a helper or `condatainer e` that resolves to a bare `.sqf`
mounts it read-only, and `helper.CheckEnv`'s `Snapshot` field (from
`PairedSize`, below) is what lets the dashboard tell this state apart from
either of the other two.

**`PairedSize`, `PairedPackages` and `PairedInfo` (`conda_pairing.go`) are the
three places that read *through* a pair rather than just locating it.** Every
caller that needs a `.img`'s size, its installed conda packages, or its
`conda-meta/history` — `helper.CheckEnv`, `#IMG_PACKAGES:` checks,
`condatainer info`, the dashboard's env-info endpoint — goes through one of
these rather than calling `LookupSnapshot` and re-deriving the merge itself.
All three take a bare path and handle every shape internally: a `.img` with a
pair reads both and combines them (size adds, packages union with the
`.img`'s own winning on conflict, history concatenates as one continuous
log); a `.img` with no pair, or a bare `.sqf`, is read alone. `PairedPackages`
and `PairedInfo` live here rather than in `internal/conda` because `conda`
cannot import this package — `container` already depends on
`internal/artifact/compare`, which depends on `internal/artifact/key`, which
depends on `conda`, so the pairing logic has to sit on the `container` side
of that edge even though it is conda-specific.

## Error Handling

- Overlay not found → search all images directories
- Multiple .img overlays → error (only one writable layer allowed)
- Locked .img → error with suggestion to check running instances
- Invalid bind path → validation error
