# overlay

Overlay image file operations providing creation, resizing, ownership changes, integrity checking, and locking for ext3/SquashFS images.

## Architecture

```
create.go   Image creation with filesystem profiles
resize.go   ext3 image resizing (grow/shrink)
chown.go    Ownership changes (requires fakeroot)
check.go    Integrity checking and repair
lock.go     File locking (shared/exclusive)
info.go     Image information queries
squashfs.go SquashFS packing operations
error.go    Structured error types
```

## Key Types

**Predefined Profiles:**
- `ProfileSmall` - 4KB per inode (Conda, Python packages)
- `ProfileDefault` - 16KB per inode (general purpose)
- `ProfileLarge` - 1MB per inode (genomes, databases)

**Stats** - `Usage() (usedBytes, percent)`, `InodeUsage() percent`

**Lock** - `Close() error`

## Creation API

There are three creation paths depending on the caller's needs:

| Function | When to use |
|---|---|
| `CreateWithOptions(ctx, opts)` | Simple create + move in one call (e.g. `internal/build`) |
| `CreateInTmp(ctx, opts)` → `MoveOverlayCopied` | Caller needs to run work on the tmp image first (e.g. conda init in `cmd/overlay`) |
| `CreateDirectly(ctx, opts)` | Target is already on fast local storage (`--no-tmp`) |

**Default tmp-staging flow** (used by `condatainer o` / `overlay create`):

1. `CreateInTmp` — builds the overlay **sparse** at `utils.GetTmpDir()` (local SSD, fast random I/O for `dd`/`mke2fs`/`debugfs`/conda)
2. Caller runs additional work (e.g. conda environment init) on the tmp path
3. `MoveOverlayCopied(tmpPath, finalPath, sparse)` — moves to destination:
   - Same filesystem → `os.Rename` (instant); may need `AllocateOverlay` afterward
   - Cross-filesystem (e.g. local `/tmp` → LustreFS) → copy then remove src:
     - `sparse=false`: `io.Copy` writes all zeros, destination is fully allocated
     - `sparse=true`: `sparseAwareCopy` via `SEEK_DATA`/`SEEK_HOLE`, holes preserved

Override the tmp location with `CNT_TMPDIR` (takes priority over `SLURM_TMPDIR`, `TMPDIR`, `/tmp`).

## Usage

```go
// --- Simple create (internal/build) ---
overlay.CreateWithOptions(ctx, &overlay.CreateOptions{
    Path: "env.img", SizeMB: 3000,
    UID: os.Getuid(), GID: os.Getgid(),
    Profile: overlay.ProfileDefault, Sparse: true, FilesystemType: "ext3", Quiet: true,
})

// --- Create at tmp, do work, then move (cmd/overlay) ---
opts := &overlay.CreateOptions{Path: "env.img", SizeMB: 10240, ...}
tmpPath, err := overlay.CreateInTmp(ctx, opts)
// ... run conda init on tmpPath ...
copied, err := overlay.MoveOverlayCopied(tmpPath, opts.Path, opts.Sparse)
if !opts.Sparse && !copied {
    overlay.AllocateOverlay(ctx, opts.Path, opts.SizeMB)
}

// Resize (grow or shrink)
overlay.Resize(ctx, "env.img", 5000)

// Ownership (requires fakeroot; internalPath = path inside overlay)
overlay.ChownRecursively(ctx, "env.img", 0, 0, "/ext3")       // → root
overlay.ChownRecursively(ctx, "env.img", 1000, 1000, "/ext3")  // → user

// Integrity
overlay.CheckIntegrity(ctx, "env.img", false) // check
overlay.CheckIntegrity(ctx, "env.img", true)  // check and repair
overlay.CheckAvailable("env.img", false)      // shared lock check
overlay.CheckAvailable("env.img", true)       // exclusive lock check

// Stats
stats, err := overlay.GetStats("env.img")

// Locking
lock, err := overlay.AcquireLock("env.img", false) // shared
lock, err := overlay.AcquireLock("env.img", true)  // exclusive
defer lock.Close()
```

## Error Types

All overlay errors use `overlay.Error` with an `Op` field (e.g. `"lock"`, `"check"`, `"resize"`, `"chown"`):

```go
var overlayErr *overlay.Error
if errors.As(err, &overlayErr) {
    fmt.Println(overlayErr.Op)
}
```

## Locking Strategy

Locks use `syscall.Flock` on the image file itself (no separate `.lock` file), non-blocking (`LOCK_NB`):

- **Shared** (`LOCK_SH`): multiple readers can hold concurrently; acquired read-only (`O_RDONLY`)
- **Exclusive** (`LOCK_EX`): single writer; blocks all other locks; requires write permission (`O_RDWR`)
- Released automatically when the file descriptor is closed

**Who holds locks:**
- `exec`/`run`: acquire and hold shared read locks on all `.sqf` overlays and the base `.sif` for the entire duration of `apptainer exec`. `.img` overlays are skipped — Apptainer flocks them itself; acquiring our own lock conflicts with Apptainer's locking.
- `overlay resize/check/chown`: acquire and hold an exclusive lock for the duration of the operation.
- `remove`: probe-and-release exclusive lock before `os.Remove()` — fails if shared lock is held.
- `build --update`: probe-and-release exclusive lock before starting any build work — fails if shared lock is held.
