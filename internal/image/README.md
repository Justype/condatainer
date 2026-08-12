# image

Image file operations: creation, resizing, ownership changes, integrity
checking, locking, reading bytes back out, and finding what is installed.

This package is **file operations and nothing else**. It opens an archive, reads
a path out of one, locks, and resolves paths — it does not know what a manifest
is, what a key is, or that `/.cnt` means anything. Everything that understands
CondaTainer metadata lives under [`internal/artifact`](../artifact/README.md) and
reads through here. The dependency is one-way: `artifact` imports `image`, never
the reverse.

## Architecture

```
lock.go     File locking (shared/exclusive)
path.go     PathExists inside an image, .sqf or .img
read.go     ReadFile from inside an image; ExtractDir, one directory in one call
scan.go     ScanOverlays: what is installed, across the image search paths
ext3/       Writable .img: create, resize, chown, check, info
squashfs/   Read-only .sqf: stat, cat, pack
sif/        .sif headers and partition offsets
tool/       Shared external-tool invocation and structured errors
```

## The overlay index

`ScanOverlays` is the single walk of the image search paths. It keys every
overlay by normalized name (`samtools--1.21.sqf` → `samtools/1.21`) and keeps
*every* copy in search order, so the first is the highest-priority one and the
shadowed ones stay reachable — `remove --layer` needs a copy the priority map
would have hidden. `FirstPaths` and `Names` reduce a scan for callers that want
one path or just the set.

Two things vary by caller, so they are options rather than separate scanners:

- **`Aliases`** adds a bare-name key for each `<base>/<name>` image, which is
  what lets `run -o build-essential` find `ubuntu24/build-essential`. Runtime
  resolution wants it; a `#DEP:` or a `remove` names an image exactly and must
  not hit it.
- **An unreadable directory** comes back as an error *alongside* a usable map.
  Container launch treats it as fatal, since resolving to a shorter list would
  silently drop an overlay; the CLI listings warn and carry on.

A `.sif` is not an overlay and never appears: it is the container root.

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

## External tools

Every archive read shells out, and two of those calls have non-obvious shapes.

`ExtractDir` decides by its **output, not the exit code**: `unsquashfs` exits 0
when nothing matched, and exits non-zero over warnings about files it extracted
perfectly well — an unprivileged user cannot restore `security.*` xattrs, which
is the common case on a shared filesystem. It passes `-no-xattrs` to avoid that
and then stats the extracted path to decide what happened. The single-file
extraction fallback passes it for the same reason.

`squashfs.PathExists` lists the archive with `unsquashfs -lc -d ""` and requires
an *exact* line match on `/entry` or a prefix match on `/entry/`. `-d ""` lists
matches as `/entry` rather than `squashfs-root/entry`; `-lc` lists only files and
empty directories, so a populated directory shows up through its children; and
unsquashfs 4.4 (Ubuntu 20.04) prints banner lines even when nothing matches, so
"there was output" is not a usable signal.

`ext3.crossFsCopy` shells out to `cp` so a long copy stays cancellable, and uses
`cmd.Start` with a Wait goroutine rather than `CombinedOutput`. On Lustre/NFS a
`cp` in uninterruptible I/O sleep does not honour SIGKILL until the I/O resolves,
and `CombinedOutput` would block in `cmd.Wait` behind it; the goroutine cleans up
whenever the process finally exits.

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
- **Exclusive** (`LOCK_EX`): single writer; blocks all other locks; opened `O_RDWR`
- Released automatically when the file descriptor is closed

**Protection.** The `O_RDWR` above is not what `flock` needs — a lock can be placed on any descriptor.
It is the rule: **an image with its write bit clear is protected and is never modified or removed**,
including for its owner, who can unlink it through the directory anyway and can restore the bit.
`chmod a-w <image>` is how an artifact is pinned, and it stays readable while pinned.

A failed attempt reports which of three things happened, and callers may branch on the first two:
`ErrProtected` (write bit clear), `ErrInUse` (a conflicting flock), or a missing file.

**Who holds locks:**
- `exec`/`run`: acquire and hold shared read locks on all `.sqf` overlays and the base `.sif` for the entire duration of `apptainer exec`. `.img` overlays are skipped — Apptainer flocks them itself; acquiring our own lock conflicts with Apptainer's locking.
- `overlay chown`: acquires and holds an exclusive lock for the duration of the operation.
- `overlay resize/check`: probe-and-release exclusive lock (via `CheckIntegrity` → `CheckAvailable`) — no lock is held across the resize2fs/e2fsck run, so the caller must not pre-acquire one (a held lock collides with the probe).
- `remove`: probe-and-release exclusive lock before `os.Remove()` — fails if shared lock is held.
- `build --update`: probe-and-release exclusive lock before starting any build work — fails if shared lock is held.
