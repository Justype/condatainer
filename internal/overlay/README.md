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

## Usage

```go
// Create
overlay.CreateForCurrentUser(ctx, "env.img", 3000, "small", true, "ext3", false)
overlay.CreateForRoot(ctx, "env.img", 3000, "default", false, "ext3", false)

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
