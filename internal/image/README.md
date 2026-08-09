# image

Image file operations: creation, resizing, ownership changes, integrity
checking, locking, embedded metadata, and finding what is installed.

## Architecture

```
lock.go     File locking (shared/exclusive)
path.go     PathExists inside an image, .sqf or .img
read.go     ReadFile from inside an image, .sqf or .img
scan.go     ScanOverlays: what is installed, across the image search paths
ext3/       Writable .img: create, resize, chown, check, info
squashfs/   Read-only .sqf: stat, cat, pack
sif/        .sif headers and partition offsets
meta/       The embedded manifest: stage, read, validate
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

## Manifests

`meta` handles `/.cnt/manifest.json`, written at build time and read back from
the packed image. It is **trusted, not verified**: nothing compares it against
the payload beside it or the build that produced it, and `Validate` checks
structure rather than truth — a known schema, a name and install prefix where
they are load-bearing, and an environment that applies without collisions.

`Runtime.Prefix` is not a mount point. An overlay is applied over the container
root, so the payload appears at `/cnt/<name>` because that is where it sits in
the archive, not because anything is mounted there. `Runtime.Env` keeps
`{prefix}` intact and substitutes it at load time, since the install prefix is
not known when the image is built.

`Read` distinguishes *no manifest* from *could not look*: only a genuinely
absent one is `ErrNoManifest`, so a caller never reports a missing `unsquashfs`
or a corrupt archive as "this image has no metadata" — those keep
`tool.ErrToolMissing`, `tool.ErrUnreadable` or `tool.ErrCorrupt`. An unknown
`SchemaVersion` joins them as `ErrUnsupportedSchema`, handled exactly like a
missing manifest; a reader ignores unknown fields, so a later schema that only
adds fields stays readable here.

Reads are cached **across processes**, keyed by absolute path and validated
against size and mtime. Repeated scans — `list`, `avail`, PATH construction, and
shell completion, which runs one process per keystroke — would otherwise spawn
one `unsquashfs` per image every time. Negative verdicts are cached too, or an
image predating the format would be re-probed on every listing, which is the
cost the cache exists to avoid.

Degradation is deliberate and asymmetric, because most images in the wild predate
the format. `CheckBase` accepts an image with no manifest, warns and accepts one
that is present but unreadable — rejecting it would strand every build behind a
base that is most likely fine — and rejects only a manifest that reads and says
it is not a base.

A writable `.img` is not handled at all. It is a mutable working overlay rather
than a built image, and its environment comes from its `.env` sidecar.

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
- **Exclusive** (`LOCK_EX`): single writer; blocks all other locks; requires write permission (`O_RDWR`)
- Released automatically when the file descriptor is closed

**Who holds locks:**
- `exec`/`run`: acquire and hold shared read locks on all `.sqf` overlays and the base `.sif` for the entire duration of `apptainer exec`. `.img` overlays are skipped — Apptainer flocks them itself; acquiring our own lock conflicts with Apptainer's locking.
- `overlay chown`: acquires and holds an exclusive lock for the duration of the operation.
- `overlay resize/check`: probe-and-release exclusive lock (via `CheckIntegrity` → `CheckAvailable`) — no lock is held across the resize2fs/e2fsck run, so the caller must not pre-acquire one (a held lock collides with the probe).
- `remove`: probe-and-release exclusive lock before `os.Remove()` — fails if shared lock is held.
- `build --update`: probe-and-release exclusive lock before starting any build work — fails if shared lock is held.
