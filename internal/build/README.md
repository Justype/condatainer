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
| `base` / `os` | a `.def` recipe | — |
| `app` | `#TYPE:app`, else a one-component name | fast scratch tmp root, `BlockSize` |
| `data` | `#TYPE:data`, else 2+ components | stable tmp root, `DataBlockSize` |

Type affects **only** where the build works and how the archive is compressed.
Every build writes to `$CNT_PREFIX` = `/cnt/<name>` — the install prefix recorded
in the manifest, so anything baked into the payload stays valid at run time.

## Architecture

```
object.go       BuildObject interface, base implementation
spec.go         Build/Spec/Options/Workspace/Target; the build environment
conda.go        Conda package builds via micromamba
script.go       Recipe (script) builds, any type
def.go          Apptainer definition file builds
pack.go         Manifest staging, shared by every backend that packs
squashfs.go     The common packer: payload plus staged metadata → one .sqf
base_image.go   Resolving the base image every non-def build runs inside
deps.go         Overlay arguments for a build's dependencies
graph.go        Build planning from a solved catalog plan, and its execution
tmppath.go      Which tmp root a build works under, by type
```

## Key Types

**BuildObject** - Interface for all build types:
- `NameVersion()`, `BuildSource()`, `Dependencies()`, `BuildType()`
- `Spec()`, `Manifest()` - What is being built, and the metadata it embeds
- `IsInstalled()`, `RequiresScheduler()`, `Update()`
- `ScriptSpecs()` - Parsed scheduler directives
- `LockPath()` - Path to the build lock file
- `Build(ctx, buildDeps bool)`, `GetMissingDependencies()`
- `CreateTmpOverlay(ctx, force bool)`, `Cleanup(failed bool)`

**Spec** - What the image will be, settled by resolution and never changed after:
identity (`Name`, `Type`, `Description`, `URL`), the `SourceSpec` whose populated
field *is* the `BuildType`, the base, dependencies, and the `meta.Runtime` block.
`Spec.Manifest()` is the embedded metadata, and it is a projection of `Spec` and
nothing else — no host, job ID, local path, or `#INPUT:` answer can reach it.

**BuildLockInfo** - JSON metadata stored in `.lock` files:
- `Runner` - `"local"`, `"slurm"`, `"pbs"`, `"lsf"`, or `"htcondor"`
- `JobID` - Scheduler job ID (empty until submit returns; set by `submitJob`)
- `Node` - Short hostname (domain suffix stripped) where lock was created
- `PID` - OS process ID for local builds; 0 for scheduler builds
- `CreatedAt` - RFC3339 timestamp

**BuildType** - `BuildTypeConda`, `BuildTypeDef`, `BuildTypeScript`

**BuildGraph** - Turns a solved `catalog.Plan` into build work: a BuildObject per missing node, run in dependency order — local builds sequentially, then scheduler submissions with job dependencies.

**ScriptSpecs** - Build resource requirements (alias to `scheduler.ScriptSpecs`)

## Usage

```go
// Auto-resolve from name/version
obj, err := build.NewBuildObject(ctx, "cellranger/9.0.1", false, imagesDir, tmpDir, false)

// Channel-annotated package: skips build script lookup, version required
obj, err := build.NewBuildObject(ctx, "bioconda::star/2.7.11b", false, imagesDir, tmpDir, false)
// → sqf: star--2.7.11b.sqf, micromamba spec: bioconda::star=2.7.11b

// Conda with custom source (YAML or package list)
obj, err := build.NewCondaObjectWithSource("myenv/1.0", "/path/env.yml", imagesDir, tmpDir, false)

// From external file (script or def)
obj, err := build.FromExternalSource(ctx, "myapp/1.0", "/path/script.sh", false, imagesDir, tmpDir)

// Single build
obj.Build(ctx, true)

// Build with dependency graph
graph, err := build.NewBuildGraph(ctx, objects, imagesDir, tmpDir, true, false)
graph.Run(ctx)

// Base image: an existing usable one, or built from its definition first
base, err := build.ResolveBase(ctx)

// Recipe lookup goes through the catalog
cat, err := config.OpenCatalog(ctx)
match, found, err := cat.Lookup(ctx, "cellranger/9.0.1")
```

## Recipes

Recipes support metadata headers:
- `#DEP:name/version` - Exact dependency
- `#DEP:name/version>=min` - Dependency with version constraint (range `[min, version]`)
- `#SBATCH` / `#PBS` / `#BSUB` - Scheduler directives (HTCondor uses native `.sub` files)
- `#ENV:VAR={prefix}/sub` - Environment variables to export; `{prefix}` is filled with the install prefix at load time
- `#INPUT:prompt` - User input, fed to the recipe on stdin in declaration order — read it with `IFS= read -r VAR` (collected locally before scheduler submission; embedded as a heredoc in the job script)

Available variables: `$CNT_NAME` (complete name), `$CNT_TYPE`, `$CNT_PREFIX`, `$CNT_TMP` (also `$TMPDIR`, both `/cnt_tmp`),
plus the scheduler's normalized `$NCPUS`, `$MEM`, `$MEM_GB`.
Run as `bash -euo pipefail <recipe>` top to bottom — no `install()` wrapper.

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

Resolution itself lives in `catalog.Resolve`; see **BuildGraph Execution** below.

## Build Lock

Each build target has a lock file at `Target.Lock` (the image path plus `.lock`) containing `BuildLockInfo` JSON.

**Lifecycle:**
- **Scheduler submit** (`submitJob` in `graph.go`): lock created with `runner=<scheduler>`, `job_id=""`, **no node** before calling the scheduler. Updated with the real job ID after `Submit()` returns. Removed if submit fails.
- **Local / scheduler job start** (`Build()`): `createBuildLock()` writes `runner=local` with current node and PID. If a scheduler lock already exists with a matching job ID (`$SLURM_JOB_ID` / `$PBS_JOBID` / `$LSB_JOBID`), it adopts and updates that lock with the runtime node and PID.
- **Build end**: `defer removeBuildLock()` removes the lock on success or failure.

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
3. Create build lock (local); create the workspace (ext3 image, or host dirs)
4. `installConda`: `micromamba create` inside container, and nothing else. The
   install run never sees the output path or binds the images directory, so a
   cancelled install cannot leave anything beside the installed images
5. `stageMetadata`, then pack to SquashFS
6. Atomic rename prepared → target; remove lock

**Shell:**
1. Check if overlay exists (skip if not updating)
2. If updating existing overlay: probe exclusive lock — fail immediately if in use
3. Create build lock (local, or adopt scheduler lock); create temporary overlay
4. Build missing dependencies (if enabled)
5. Run recipe inside container. The payload directory is bound at
   `/cnt/<name>` — the leaf, not `/cnt`, so dependency overlays mounted beside
   it stay visible
6. `stageMetadata`, then pack: verify the payload directory is non-empty, then
   `mksquashfs` the host build dir (basename `cnt`, so the archive is exactly
   `cnt/<name>/…` — the build's `tmp/` is a sibling and never enters it)
7. Atomic rename prepared → target; remove lock

### Staging and packing

`stageMetadata` validates `Spec.Manifest()` and writes it to `<buildDir>/.cnt/`,
*beside* the payload rather than inside it. That placement is the whole trick:
`mksquashfs` makes one archive root per source, so handing it the payload dir and
the `.cnt` dir as two sources yields `/cnt/<name>/…` and `/.cnt/manifest.json` at
the same level, without the payload ever containing a directory the recipe did
not create. There is no option to rename a source, so the staged directory must
already be called `.cnt` when `mksquashfs` sees it — which is why the ext3 mode,
where the payload is inside the temporary image, binds the host dir to `/.cnt`.

Validation runs at staging because that is the only point where a bad manifest
can still stop the build: earlier there is nothing to validate, later the image
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
4. Build SIF with Apptainer; extract SquashFS partition
5. Atomic rename prepared → target; remove lock

**Base image (`.sif`):**
The definition backend with one branch changed: the SIF is the product, so step 4
keeps it instead of extracting a partition. `IsInstalled` searches every image
path rather than one target, so a base supplied by a shared install is not
rebuilt into the user's own directory.

## The base image

Script and Conda builds run their install *and* their packing step inside a
container, so each has an implicit edge to the base image — a build prerequisite,
separate from anything a recipe declares with `#DEP:`. `BuildGraph.Run` resolves
it once before any node runs, and records it on every dependent as `Spec.Base`;
`ResolveBase` builds the configured base first if none is installed.

A definition bootstraps its own root, so it resolves nothing — which is what lets
the base's own build run when no base exists yet.

An installed base passes `meta.CheckBase`: no manifest is fine (bases predate the
format), an unreadable one warns, and one that reads has to say `type: base`.

## BuildGraph Execution

1. **Solve** - `catalog.Resolve` walks the graph from the index alone: transitive
   deps, cycle detection and dependency-first order, without fetching a recipe.
   Installed versions come from the `Have` callback, so an installed dep is a map
   lookup rather than a build object.
2. **Plan** - Create a BuildObject per missing node, in the solved order
3. **Resolve the base** - Every script and Conda build runs inside it, so it is
   built first; a scheduler job cannot build one on the node
4. **Separate local/scheduler** - Based on resource requirements
5. **Execute:**
   - Local builds: sequentially, in dependency order
   - Scheduler builds: submitted with dependency chains, each waiting on its deps

## Environment Variables

Images can export environment variables via `#ENV:` directives:
```bash
#ENV:CELLRANGER_ROOT={prefix}
#ENV:PATH={prefix}/bin:$PATH
```

These are captured into `Spec.Runtime` at resolution and embedded in the image's
manifest. `{prefix}` survives into the manifest verbatim and is substituted with
`runtime.prefix` at load time (see `internal/runtime/container` env handling).

## Workspace Strategy

`workspaceFor` derives every path from `(name, tmp root, scratch extension,
payload location)`. The tmp root and the mode are both functions of the type:

| type | scratch image | payload |
|---|---|---|
| `app` | `.img` under `use_tmp_overlay`, else none | in the image, or host in directory mode |
| `data` | **never** — see below | always host |
| `os`, `base` | `.sif`, always | apptainer's rootfs |

The ext3 image is an app-build optimisation: it keeps a conda environment's
thousands of small files off the host's inode budget. A data payload is a few
large files staged on the host, so an image would be created, mounted and
discarded holding nothing; a definition build has apptainer's own rootfs. So
`appExt3ScratchExt` returns `""` for every type but `app`, whatever the config
says, and `Workspace.UsesImage()` is the mode everything downstream reads —
never `config.Global` a second time.

Each build type also uses a different base directory for build artifacts:

| Build path | `tmpDir` source | Rationale |
|---|---|---|
| Conda (`name/version`) | `utils.GetTmpDir()` | Fast local node storage (scheduler TMPDIR → TMPDIR → `/tmp/cnt-$USER`) |
| Script, type `app` | `utils.GetTmpDir()` | Fast local node storage |
| Script, type `data` | `config.GetWritableTmpDir()` | Stable condatainer data path (large datasets) |
| Def (internal) | `config.GetWritableTmpDir()` | Stable path; set in `createConcreteType` after type is resolved |
| External sh | `filepath.Dir(targetPrefix)` | Next to output target (user controls location) |
| External def | `filepath.Dir(targetPrefix)` | Next to output target (user controls location) |
| Base image | `config.GetWritableTmpDir()` | Stable condatainer data path |

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
