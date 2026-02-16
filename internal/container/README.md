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
- `Overlays` - Resolved ordered overlay paths
- `OverlayArgs` - Paths with `:ro/:rw` suffixes
- `EnvList` - Complete environment variable list
- `EnvNotes` - Environment descriptions for display
- `BindPaths` - Deduplicated bind paths
- `Fakeroot` - Final fakeroot setting
- `ApptainerFlags` - Flags including GPU detection
- `LastImg` - Path to .img overlay (if present)

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
})
```

### GPU Detection

```go
// Auto-detect GPUs and generate flags
gpuFlags := container.DetectGPUFlags()
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
2. **Validate** - Ensure at most one .img overlay
3. **Order** - Place .img last (required for writable overlay layer)
4. **Lock check** - Verify .img availability (exclusive lock if writable)
5. **Build overlay args** - Add `:ro/:rw` suffixes
6. **Collect environment** - Load `.env` files for each overlay
7. **Deduplicate binds** - Remove conflicting bind paths
8. **Detect GPU** - Add `--nv` or `--rocm` if available
9. **Return result** - Ready-to-use configuration

## Overlay Ordering

- SquashFS (`.sqf`) overlays are read-only, order matters for file precedence
- ext3 (`.img`) overlay must be last (writable layer)
- Automatic reordering ensures correct layering

## Environment Variables

Each overlay can have an associated `.env` file (e.g., `cellranger--9.0.1.sqf.env`):
```bash
CELLRANGER_ROOT=/ext3/cnt/cellranger/9.0.1
PATH=/ext3/cnt/cellranger/9.0.1/bin:$PATH
```

Loaded automatically and merged into final environment list.

**Common Environment:**
- `LC_ALL=C.UTF-8`, `LANG=C.UTF-8`
- `CURL_CA_BUNDLE=`, `SSL_CERT_FILE=` (unset to avoid host interference)

## Bind Path Deduplication

Removes conflicting bind paths:
- Keeps longest/most specific paths
- Removes parent paths when child is bound
- Example: `/home/user/data` removes `/home/user`

## GPU Detection

Auto add flags based on available devices:

- **NVIDIA**: Checks for `/dev/nvidia*`, adds `--nv`
- **AMD**: Checks for `/dev/kfd`, adds `--rocm`

## Error Handling

- Overlay not found → search all images directories
- Multiple .img overlays → error (only one writable layer allowed)
- Locked .img → error with suggestion to check running instances
- Invalid bind path → validation error
