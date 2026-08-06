# build

Build system for creating overlay images from Conda packages, shell scripts, or Apptainer definition files with dependency resolution and scheduler integration.

## Two axes

`BuildType` is **how** a target is built. `catalog.Kind` is **what** the payload is.
They are independent — a script build can be any kind.

| BuildType | Source | Builder |
|------|--------|-------------|
| Conda | Package spec, YAML, or `channel::pkg/version` | Micromamba environment |
| Def | Apptainer `.def` file | Apptainer build |
| Script | Recipe (`#DEP`, `#SBATCH`, `#ENV`, `#INPUT`) | Recipe run as a shell script |

| Kind | Comes from | Decides |
|---|---|---|
| `base` / `os` | a `.def` recipe | — |
| `app` | `#TYPE:app`, else a one-component name | fast scratch tmp root, `BlockSize` |
| `data` | `#TYPE:data`, else 2+ components | stable tmp root, `DataBlockSize` |

Kind affects **only** where the build works and how the archive is compressed.
Every build writes to `$CNT_PREFIX` = `/cnt/<name>/<version>` — the path the
artifact is mounted at, so anything baked into the payload stays valid at run time.

## Architecture

```
object.go       BuildObject interface, base implementation
conda.go        Conda package builds via micromamba
script.go       Recipe (script) builds, any kind
def.go          Apptainer definition file builds
base_image.go   Base image provisioning (download/build)
graph.go        Build planning from a solved catalog plan, parallel execution
fetch.go        Remote build script downloading
env.go          Environment variable extraction from overlays
```

## Key Types

**BuildObject** - Interface for all build types:
- `NameVersion()`, `BuildSource()`, `Dependencies()`, `Type()`
- `IsInstalled()`, `RequiresScheduler()`, `Update()`
- `ScriptSpecs()` - Parsed scheduler directives
- `LockPath()` - Path to the build lock file
- `Build(ctx, buildDeps bool)`, `GetMissingDependencies()`
- `CreateTmpOverlay(ctx, force bool)`, `Cleanup(failed bool)`

**BuildLockInfo** - JSON metadata stored in `.lock` files:
- `Type` - `"local"`, `"slurm"`, `"pbs"`, `"lsf"`, or `"htcondor"`
- `JobID` - Scheduler job ID (empty until submit returns; set by `submitJob`)
- `Node` - Short hostname (domain suffix stripped) where lock was created
- `PID` - OS process ID for local builds; 0 for scheduler builds
- `CreatedAt` - RFC3339 timestamp

**BuildType** - `BuildTypeConda`, `BuildTypeDef`, `BuildTypeScript`

**BuildGraph** - Turns a solved `catalog.Plan` into build work: a BuildObject per missing node, then parallel execution (local worker pool + scheduler job submission).

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

// Base image
build.EnsureBaseImage(ctx, false) // download prebuilt or build from def

// Recipe lookup goes through the catalog
cat, err := config.OpenCatalog(ctx)
match, found, err := cat.Lookup(ctx, "cellranger/9.0.1")
```

## Recipes

Recipes support metadata headers:
- `#DEP:name/version` - Exact dependency
- `#DEP:name/version>=min` - Dependency with version constraint (range `[min, version]`)
- `#SBATCH` / `#PBS` / `#BSUB` - Scheduler directives (HTCondor uses native `.sub` files)
- `#ENV:VAR={prefix}/sub` - Environment variables to export; `{prefix}` is filled with the mount root at load time
- `#INPUT:prompt` - User input, fed to the recipe on stdin in declaration order — read it with `IFS= read -r VAR` (collected locally before scheduler submission; embedded as a heredoc in the job script)

Available variables: `$CNT_NAME`, `$CNT_VERSION`, `$CNT_KIND`, `$CNT_PREFIX`, `$CNT_TMP` (also `$TMPDIR`),
plus the scheduler's normalized `$NCPUS`, `$MEM`, `$MEM_GB`.
Run as `bash -euo pipefail <recipe>` top to bottom — no `install()` wrapper.

### Dependency Resolution

**`#DEP:` is a build dependency and nothing else.** It names what must be mounted
*while this recipe runs* — a tool that unpacks, validates or indexes the payload.
It is not recorded in the artifact and is not re-expanded when the artifact is
mounted later: there is no runtime dependency tree. What an overlay needs at run
time is whatever the caller names on the command line.

That is why an `app` is **self-contained** — a conda environment or a prebuilt
package that carries its own libraries — while `data` is the kind that normally
has deps, because producing an index needs the tool that produces it.

`#DEP:samtools/1.22.1>=1.10` — accepts any installed version in `[1.10, 1.22.1]`:
- If a satisfying version is installed → skip build, mount the latest satisfying version
- If none installed → build the preferred version (`1.22.1`)
- Versions above preferred (`2.0`) are rejected (implicit upper bound)
- Operators: `>=` (inclusive lower bound) and `>` (exclusive lower bound)

Resolution itself lives in `catalog.Resolve`; see **BuildGraph Execution** below.

## Build Lock

Each build target has a lock file at `<targetOverlayPath>.lock` containing `BuildLockInfo` JSON.

**Lifecycle:**
- **Scheduler submit** (`submitJob` in `graph.go`): lock created with `type=scheduler`, `job_id=""`, **no node** before calling the scheduler. Updated with the real job ID after `Submit()` returns. Removed if submit fails.
- **Local / scheduler job start** (`Build()`): `createBuildLock()` writes `type=local` with current node and PID. If a scheduler lock already exists with a matching job ID (`$SLURM_JOB_ID` / `$PBS_JOBID` / `$LSB_JOBID`), it adopts and updates that lock with the runtime node and PID.
- **Build end**: `defer removeBuildLock()` removes the lock on success or failure.

**Stale detection** (in `NewBuildObject`):

| Lock state | Action |
|---|---|
| Empty or corrupt file (old format) | Treat as stale — remove and continue |
| `type=slurm/pbs/lsf/htcondor`, `job_id=""` | Stale (failed before submit returned) — remove |
| `type=slurm/pbs/lsf/htcondor`, `job_id` set | Call `scheduler.IsJobAlive(job_id)` — remove if not alive |
| `type=local`, same node, PID dead | Stale — remove and continue |
| `type=local`, same node, PID alive | Active — error with details |
| `type=local`, different node | Cannot verify remotely — error, ask user to check |
| Scheduler unavailable for job check | Conservative — error, ask user to check |

## Build Workflow

**Conda:**
1. Check if overlay exists (skip if not updating)
2. If updating existing overlay: probe exclusive lock — fail immediately if in use
3. Create build lock (local); create temporary ext3 overlay
4. `micromamba create` inside container
5. Set permissions, pack to SquashFS
6. Atomic rename `.new` → target; remove lock

**Shell:**
1. Check if overlay exists (skip if not updating)
2. If updating existing overlay: probe exclusive lock — fail immediately if in use
3. Create build lock (local, or adopt scheduler lock); create temporary overlay
4. Build missing dependencies (if enabled)
5. Run recipe inside container. The payload directory is bound at
   `/cnt/<name>/<version>` — the leaf, not `/cnt`, so dependency overlays
   mounted beside it stay visible
6. Pack, as a separate container run: verify the payload directory is non-empty,
   then `mksquashfs` the host build dir (basename `cnt`, so the archive is
   exactly `cnt/<name>/<version>/…` — the build's `tmp/` is a sibling and never
   enters it)
7. Atomic rename `.new` → target; remove lock

The build script (with its `#ENV:` directives) is embedded inside the overlay at `/cnt/<name>/<version>/.cnt-build-script`, so the overlay is self-contained; env is resolved from it at load time rather than written to a sidecar `.env` at build time.

**Def:**
1. Check if overlay exists (skip if not updating)
2. If updating existing overlay: probe exclusive lock — fail immediately if in use
3. Create build lock; try prebuilt download (if remote source)
4. Build SIF with Apptainer; extract SquashFS partition
5. Atomic rename `.new` → target; remove lock

**Base image (`.sif`):**
1. Check if already installed (skip if not updating)
2. If updating existing image: probe exclusive lock — fail immediately if in use
3. Create build lock; try prebuilt `.sif` download (if remote source)
4. Build SIF with Apptainer
5. Atomic rename `.new` → target; remove lock

## BuildGraph Execution

1. **Solve** - `catalog.Resolve` walks the graph from the index alone: transitive
   deps, cycle detection and dependency-first order, without fetching a recipe.
   Installed versions come from the `Have` callback, so an installed dep is a map
   lookup rather than a build object.
2. **Plan** - Create a BuildObject per missing node, in the solved order
3. **Separate local/scheduler** - Based on resource requirements
4. **Parallel execution:**
   - Local builds: Run concurrently with worker pool
   - Scheduler builds: Submit with dependency chains

## Environment Variables

Overlays can export environment variables via `#ENV:` directives:
```bash
#ENV:CELLRANGER_ROOT=$app_root
#ENV:PATH=$app_root/bin:$PATH
```

These live in the build script embedded inside the overlay. `$app_root` is resolved to the overlay's mount root at load time (see `internal/runtime/container` env handling).

## Tmp Directory Strategy

Each build type uses a different base directory for build artifacts (`$TMPDIR`, build dir, tmp overlay):

| Build path | `tmpDir` source | Rationale |
|---|---|---|
| Conda (`name/version`) | `utils.GetTmpDir()` | Fast local node storage (scheduler TMPDIR → TMPDIR → `/tmp/cnt-$USER`) |
| Script, kind `app` | `utils.GetTmpDir()` | Fast local node storage |
| Script, kind `data` | `config.GetWritableTmpDir()` | Stable condatainer data path (large datasets) |
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
