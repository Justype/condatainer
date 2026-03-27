# build

Build system for creating overlay images from Conda packages, shell scripts, or Apptainer definition files with dependency resolution and scheduler integration.

## Build Types

| Type | Source | Description |
|------|--------|-------------|
| Conda | Package spec, YAML, or `channel::pkg/version` | Micromamba environment |
| Shell | Build script (`#DEP`, `#SBATCH`, `#ENV`, `#INTERACTIVE`) | Custom installation |
| Def | Apptainer `.def` file | Apptainer build |
| Ref | Shell script (large datasets) | Reference data (genomes, indices) |

## Architecture

```
object.go       BuildObject interface, base implementation
conda.go        Conda package builds via micromamba
script.go       Shell script builds (apps and ref data)
def.go          Apptainer definition file builds
base_image.go   Base image provisioning (download/build)
graph.go        Dependency graph, topological sort, parallel execution
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

**BuildType** - `IsConda`, `IsDef`, `IsShell`, `IsRef`

**BuildGraph** - Dependency graph with topological sort, cycle detection, and parallel execution (local worker pool + scheduler job submission).

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

// Script lookup
info, found := build.FindBuildScript("cellranger/9.0.1")
```

## Build Scripts

Shell scripts support metadata headers:
- `#DEP:name/version` - Dependencies
- `#SBATCH` / `#PBS` / `#BSUB` - Scheduler directives (HTCondor uses native `.sub` files)
- `#ENV:VAR=$app_root` - Environment variables to export
- `#INTERACTIVE:prompt` - User input prompts

Available variables: `$NCPUS`, `$target_dir`, `$tmp_dir`, `$app_name`, `$version`
Required function: `install()`

### Dependency Resolution

```go
// Check missing dependencies
missing, err := obj.GetMissingDependencies()

// Find and download remote build script
info, found := build.FindBuildScript("cellranger/9.0.1")
if found && info.IsRemote {
    localPath, err := build.DownloadRemoteScript(info, tmpDir)
}
```

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
5. Run script inside container with overlays
6. For ref: verify files, create SquashFS from `$cnt_dir`
7. For apps: create SquashFS from `/cnt`
8. Extract and save ENV variables
9. Atomic rename `.new` → target; remove lock

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

1. **Expand dependencies** - Resolve all transitive deps
2. **Detect cycles** - Fail on circular dependencies
3. **Topological sort** - Order by dependencies
4. **Separate local/scheduler** - Based on resource requirements
5. **Parallel execution:**
   - Local builds: Run concurrently with worker pool
   - Scheduler builds: Submit with dependency chains

## Environment Variables

Overlays can export environment variables via `#ENV:` directives:
```bash
#ENV:CELLRANGER_ROOT=$app_root
#ENV:PATH=$app_root/bin:$PATH
```

Extracted to `.env` file alongside overlay, loaded at runtime.

## Tmp Directory Strategy

Each build type uses a different base directory for build artifacts (`$TMPDIR`, build dir, tmp overlay):

| Build path | `tmpDir` source | Rationale |
|---|---|---|
| Conda (`name/version`) | `utils.GetTmpDir()` | Fast local node storage (scheduler TMPDIR → TMPDIR → `/tmp/cnt-$USER`) |
| App script (`name/version`) | `utils.GetTmpDir()` | Fast local node storage |
| Ref script (`a/b/c`, 2+ slashes) | `config.GetWritableTmpDir()` | Stable condatainer data path (large datasets) |
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
