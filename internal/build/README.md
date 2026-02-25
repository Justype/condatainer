# build

Build system for creating overlay images from Conda packages, shell scripts, or Apptainer definition files with dependency resolution and scheduler integration.

## Build Types

| Type | Source | Description |
|------|--------|-------------|
| Conda | Package spec or YAML | Micromamba environment |
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
- `Build(ctx, buildDeps bool)`, `GetMissingDependencies()`
- `CreateTmpOverlay(ctx, force bool)`, `Cleanup(failed bool)`

**BuildType** - `IsConda`, `IsDef`, `IsShell`, `IsRef`

**BuildGraph** - Dependency graph with topological sort, cycle detection, and parallel execution (local worker pool + scheduler job submission).

**ScriptSpecs** - Build resource requirements (alias to `scheduler.ScriptSpecs`)

## Usage

```go
// Auto-resolve from name/version
obj, err := build.NewBuildObject(ctx, "cellranger/9.0.1", false, imagesDir, tmpDir, false)

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

## Build Workflow

**Conda:**
1. Check if overlay exists (skip if not updating)
2. If updating existing overlay: probe exclusive lock — fail immediately if in use
3. Create temporary ext3 overlay
4. `micromamba create` inside container
5. Set permissions, pack to SquashFS
6. Atomic rename `.new` → target; cleanup

**Shell:**
1. Check if overlay exists (skip if not updating)
2. If updating existing overlay: probe exclusive lock — fail immediately if in use
3. Create temporary overlay
4. Build missing dependencies (if enabled)
5. Run script inside container with overlays
6. For ref: verify files, create SquashFS from `$cnt_dir`
7. For apps: create SquashFS from `/cnt`
8. Extract and save ENV variables
9. Atomic rename `.new` → target; cleanup

**Def:**
1. Check if overlay exists (skip if not updating)
2. If updating existing overlay: probe exclusive lock — fail immediately if in use
3. Try prebuilt download (if remote source)
4. Build SIF with Apptainer; extract SquashFS partition
5. Atomic rename `.new` → target; cleanup

**Base image (`.sif`):**
1. Check if already installed (skip if not updating)
2. If updating existing image: probe exclusive lock — fail immediately if in use
3. Try prebuilt `.sif` download (if remote source)
4. Build SIF with Apptainer
5. Atomic rename `.new` → target; cleanup

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

## Resource Allocation for Builds

CPU and memory derived from `scriptSpecs`:

| State | `HasDirectives` | `Spec` | CPU source | Memory source |
|---|---|---|---|---|
| No directives | false | non-nil | `Build.Defaults.CpusPerTask` | `Build.Defaults.MemPerNodeMB` |
| Has directives | true | non-nil | `CpusPerTask × TasksPerNode` | `MemPerNodeMB` |
| Passthrough | true | nil | build fails with error | — |

Build defaults are configured via `build.*` keys in `config.yaml`. See [Config README](../config/README.md).

## Error Handling

- `ErrTmpOverlayExists` - Temporary overlay already exists (build in progress)
- `ErrBuildCancelled` - User interrupted build (Ctrl+C)
- `"cannot update <name>: ..."` - Target overlay is locked (currently used by a running `exec`/`run`)
