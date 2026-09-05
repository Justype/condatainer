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

## What Setup refuses

Two checks run before anything is locked or mounted, because both describe a
container that cannot be what the caller asked for:

- `ensureSingleImage` — at most one writable `.img`. Only one can be the
  principal image, and only that one takes a write lock.
- `ensureDistinctPrefixes` — no two images claiming one `/cnt/<name>` subtree.
  Overlays are disjoint subtrees, not stacked diffs, so two claiming one prefix
  do not combine: the later mount takes the subtree and the earlier contributes
  nothing, while both still reach PATH and the environment. Two builds of one
  name are the case this catches — a project's restored copy and a flat install
  of the same name record the same prefix. A base, an OS image and anything
  without readable metadata record no prefix and are exempt for free; the same
  file named twice is redundant, not a collision.

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
6. **Collect environment** - Resolve `#ENV:` for each overlay (embedded build script, sidecar `.env` overrides)
7. **Deduplicate binds** - Remove conflicting bind paths
8. **Detect GPU** - Add `--nv` or `--rocm` if available
9. **Return result** - Ready-to-use configuration

## Overlay Ordering

- SquashFS (`.sqf`) overlays are read-only, order matters for file precedence
- ext3 (`.img`) overlay must be last (writable layer)
- Automatic reordering ensures correct layering

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

`CollectOverlayEnv` merges every overlay's resolved env into the final list;
`ResolveOverlayEnv` returns a single overlay's description/env/notes for `info`.
When the same variable is set by more than one overlay the later one wins and a
diagnostic is recorded. `CollectOverlayEnv` is the one place degradation is
reported, so an image with no readable metadata is announced once per invocation
rather than once per consumer.

`BuildPathEnv` is silent for the same reason. Only an `app` contributes, and it
contributes `<prefix>/bin` unconditionally — `data`, `os` and `base` put nothing
on `PATH`, and an image with no readable metadata contributes nothing at all.
There is no check that the directory exists: a nonexistent `PATH` entry is
harmless, while the check would cost one archive read per image per invocation.

**Common Environment:**
- `LC_ALL=C.UTF-8`, `LANG=C.UTF-8`
- `CURL_CA_BUNDLE=`, `SSL_CERT_FILE=` (unset to avoid host interference)

## Automatic Bind Paths

`BindPaths()` in `bind.go` collects bind mounts automatically before the user-specified paths are appended:

| Source | Condition |
|--------|-----------|
| Current working directory | Always |
| `$SCRATCH` | When env var is set |
| `$TMPDIR` | When env var is set and the path exists |
| Scheduler node-local tmp | When inside a job exposing one (e.g. `SLURM_TMPDIR`) |
| Base data directories | From config; added `:ro` if not writable |
| `condatainer` executable | For nested calls (always; bound to `/usr/bin/condatainer`) |

After collection, `DeduplicateBindPaths()` removes conflicting bind paths:
- Keeps longest/most specific paths
- Removes parent paths when child is bound
- Example: `/home/user/data` removes `/home/user`

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

## Error Handling

- Overlay not found → search all images directories
- Multiple .img overlays → error (only one writable layer allowed)
- Locked .img → error with suggestion to check running instances
- Invalid bind path → validation error
