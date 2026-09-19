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

`ext3.Resize` takes the same `sparse` decision and defaults to allocating: growing
extends the file with `os.Truncate`, so the added space is a hole whatever the image
was built as.

## External tools

Finding a runnable path for a bare tool name is `internal/toolpath`'s job,
not this package's — see that package's README for why it's split out
(short version: `internal/conda` needs the same lookup and isn't under
`internal/image`, and `internal/image/squashfs` resolving back through
anything that reaches `internal/libexec`'s consumers would cycle). Every
call site here and in `ext3/`, `freeze/`, and `internal/conda` that shells
out to `debugfs`, `e2fsck`, `tune2fs`, `mke2fs`, `resize2fs`, `fuse2fs`,
`mksquashfs`, `unsquashfs`, or `squashfuse`/`squashfuse_ll` goes through
`toolpath.Resolve`/`toolpath.Command` (directly, or via `tool.
CheckDependencies`/`tool.RunCommand`, which call them internally) rather
than handing a bare name to `exec.Command`/`exec.CommandContext`, which is
`PATH`-only with no fallback. `dd`, `fallocate`, and `unshare` (coreutils/
util-linux, not e2fsprogs) are deliberately excluded — core-OS tools even
more universal than e2fsprogs, never expected to need the fallback.

`toolpath.Resolve` checks `internal/libexec`'s self-provisioned copy first —
required over whatever the host has, when it has one at all — then `PATH`,
then the FHS fallback directories. `mksquashfs`/`unsquashfs`/`squashfuse` are
provisioned into `libexec` alongside `apptainer`; `debugfs`/`e2fsck`/
`tune2fs`/`mke2fs`/`resize2fs`/`fuse2fs` (e2fsprogs) are not, so a name
`libexec` never provisions simply isn't found there and falls through to
`PATH`/FHS unaffected.

e2fsprogs (everything but `fuse2fs`) is host-only, never self-provisioned
into `internal/libexec`: it is practically guaranteed present on any Linux
host already, so bundling it would only duplicate what is already there.
`fuse2fs` itself has no independent source condatainer can bundle either — it
stays a documented prerequisite.

What stays in *this* package's own `tool` subpackage: `Error`/`analyze()`
(ext3/squashfs-specific failure→hint text) and the archive sentinels
`ErrFileNotFound`/`ErrUnreadable`/`ErrCorrupt` — genuinely specific to
interpreting an image-archive read, unlike lookup itself.
`toolpath.ErrToolMissing` is what a caller checks for "couldn't find the
binary at all," not an `internal/image/tool` sentinel — finding a tool in
the first place isn't an image-archive concern.

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

Locks are `fcntl` open-file-description locks over the whole image file (no separate `.lock` file),
non-blocking. They are the kind Apptainer takes on a mounted ext3 image, so a running container
holds a conflicting lock on every filesystem — a `flock` would not see it where `flock` and `fcntl`
locks are separate (NFSv4). A lock belongs to the descriptor, not the process, so two acquisitions in
one process conflict and closing the file releases it.
The bare open+lock mechanics live in `internal/utils.AcquireFileLock` — shared with
`internal/libexec`'s own toolchain-generation lock, which has no write-protection concept and so
calls it directly rather than through this package (see `internal/utils/README.md`, *Locking*).
`Lock` here is a type alias for `*utils.FileLock`; what this package adds on top is the
protection semantics below:

- **Shared**: multiple readers can hold concurrently; acquired read-only (`O_RDONLY`)
- **Exclusive**: single writer; blocks all other locks; opened `O_RDWR`
- Released automatically when the file descriptor is closed

**Protection.** An exclusive `fcntl` lock needs a descriptor open for writing, so a write lock opens
`O_RDWR` — and that is also the rule: **an image with its write bit clear is protected and is never modified or removed**,
including for its owner, who can unlink it through the directory anyway and can restore the bit.
`chmod a-w <image>` is how an artifact is pinned, and it stays readable while pinned.

A failed attempt reports which of three things happened, and callers may branch on the first two:
`ErrProtected` (write bit clear), `ErrInUse` (a conflicting lock), or a missing file.

**Who holds locks:**
- `exec`/`run`: acquire and hold shared read locks on all `.sqf` overlays and the base image for the entire duration of `apptainer exec`. `.img` overlays are skipped — Apptainer locks them itself; acquiring our own lock conflicts with Apptainer's locking.
- `overlay chown`: acquires and holds an exclusive lock for the duration of the operation.
- `overlay resize/check`: probe-and-release exclusive lock (via `CheckIntegrity` → `CheckAvailable`) — no lock is held across the resize2fs/e2fsck run, so the caller must not pre-acquire one (a held lock collides with the probe).
- `remove`: probe-and-release exclusive lock before `os.Remove()` — fails if shared lock is held.
- `build --update`: probe-and-release exclusive lock before starting any build work — fails if shared lock is held.
- `overlay freeze`: acquires and holds a **shared** lock on the source `.img` for the whole pack. Held rather than probed because the pack takes minutes: a probe would only prove the overlay was idle when the command started. Shared is the entire requirement — it conflicts with the exclusive lock a writer takes — and `O_RDONLY` needs no write bit, so a pinned overlay is still freezable.

**An exclusive lock on an image a container will mount is probed, never held.** Apptainer acquires that
lock itself, so holding our own would collide with it — which is why `exec`/`run` skip `.img` overlays
entirely and why `resize`/`check`/`remove`/`build --update` probe and release. The probe exists to put
the error before the container starts instead of inside it. This is settled: none of these becomes a
held lock.

**An exec-path probe asks for the lock the mount will take.** `container.Setup` passes its `writeLock`
and `run` passes `--writable`, so a writable overlay is probed exclusive and a read-only one shared.
The probe answers the question the mount is about to ask, which it can only do by asking the same one.
