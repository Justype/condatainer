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
graph.go        Dependency graph, topological sort, parallel execution
fetch.go        Remote build script downloading
env.go          Environment variable extraction from overlays
whitelist.go    Package name validation
```

## Key Types

**BuildObject** - Interface for all build types:
- `Type()`, `NameVersion()`, `Build(ctx, buildDeps)`, `CreateTmpOverlay()`
- `GetMissingDependencies()`, `RequiresScheduler()`, `IsInstalled()`

**BuildType** - `IsConda`, `IsDef`, `IsShell`, `IsRef`

**BuildGraph** - Manages dependency resolution:
- Expands transitive dependencies
- Detects cycles
- Topological sort (dependencies → dependents)
- Separates local vs scheduler builds
- Parallel execution with job dependencies

**ScriptSpecs** - Build resource requirements (alias to `scheduler.ScriptSpecs`)

## Usage

### Creating BuildObjects

```go
// Auto-resolve from name/version
obj, err := build.NewBuildObject("cellranger/9.0.1", external=false, imagesDir, tmpDir)

// Conda with custom source (YAML or package list)
obj, err := build.NewCondaObjectWithSource("myenv/1.0", "/path/env.yml", imagesDir, tmpDir)

// From external file (script or def)
obj, err := build.FromExternalSource("myapp/1.0", "/path/script.sh", isApptainer=false, imagesDir, tmpDir)
```

### Building

```go
// Single build
if err := obj.Build(ctx, buildDeps=true); err != nil {
    // Handle build errors
}

// Build with dependency graph
graph, err := build.NewBuildGraph(objects, imagesDir, tmpDir, submitJobs=true)
if err := graph.Run(ctx); err != nil {
    // Handle errors
}
```

### Build Scripts

Shell scripts support metadata headers:
- `#DEP:name/version` - Dependencies
- `#SBATCH` / `#PBS` / `#BSUB` - Scheduler directives (HTCondor uses native `.sub` files)
- `#ENV:VAR=$app_root` - Environment variables to export
- `#INTERACTIVE:prompt` - User input prompts

Available variables in scripts: `$NCPUS`, `$target_dir`, `$tmp_dir`, `$app_name`, `$version`

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
1. Check if overlay exists
2. Create temporary ext3 overlay
3. `micromamba create` inside container
4. Set permissions, pack to SquashFS
5. Cleanup

**Shell:**
1. Check if overlay exists
2. Create temporary overlay
3. Build missing dependencies (if enabled)
4. Run script inside container with overlays
5. For ref: verify files, create SquashFS from `$cnt_dir`
6. For apps: create SquashFS from `/cnt`
7. Extract and save ENV variables
8. Cleanup

**Def:**
1. Check if overlay exists
2. Build SIF with Apptainer
3. Extract SquashFS partition from SIF
4. Cleanup

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

CPU and memory are derived from `scriptSpecs` via `effectiveNcpus()` / `effectiveMemMB()`:

| State | `HasDirectives` | `Spec` | CPU source | Memory source |
|---|---|---|---|---|
| No directives | false | non-nil | `Build.DefaultCPUs` | `Build.DefaultMemMB` |
| Has directives | true | non-nil | `CpusPerTask × TasksPerNode` | `MemPerNodeMB` |
| Passthrough | true | nil | build fails with error | — |

**Passthrough mode** (unsupported scheduler flags).

Build defaults are configured via `build.*` keys in `config.yaml`. See [Config README](../config/README.md) for details.

## Build Script Search

Searches in order:
1. Local absolute path
2. `CNT_EXTRA_BASE_DIRS`
3. Portable dir (`<install>/build-scripts/`)
4. Scratch dir (`$SCRATCH/condatainer/build-scripts/`)
5. User dir (`~/.local/share/condatainer/build-scripts/`)
6. Remote fetch from GitHub (if `PreferRemote` or not found locally)

## Error Handling

- `ErrTmpOverlayExists` - Temporary overlay already exists (build in progress)
- `ErrBuildCancelled` - User interrupted build (Ctrl+C)
- Uses `isCancelledByUser()` to detect interrupts gracefully
