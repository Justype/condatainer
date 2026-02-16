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
utils.go    Helpers (is.sqf check, etc.)
error.go    Structured error types
```

## Key Types

**CreateOptions:**
- `Path`, `SizeMB`, `UID`, `GID`, `Profile`, `Sparse`, `FilesystemType`, `Quiet`

**Profile** - Filesystem tuning:
- `InodeRatio` - Bytes per inode
- `ReservedPerc` - Reserved space percentage

**Predefined Profiles:**
- `ProfileSmall` - 4KB per inode (Conda, Python packages)
- `ProfileDefault` - 16KB per inode (general purpose)
- `ProfileLarge` - 1MB per inode (genomes, databases)

## Usage

### Create Overlay

```go
import "github.com/Justype/condatainer/internal/overlay"

ctx := context.Background()

// For current user
err := overlay.CreateForCurrentUser(
    ctx, "env.img", 3000, "small", sparse=true, "ext3", quiet=false)

// For root (requires fakeroot)
err := overlay.CreateForRoot(
    ctx, "env.img", 3000, "default", sparse=false, "ext3", quiet=false)

// Custom options
opts := &overlay.CreateOptions{
    Path:           "custom.img",
    SizeMB:         5000,
    UID:            1000,
    GID:            1000,
    Profile:        overlay.ProfileLarge,
    Sparse:         true,
    FilesystemType: "ext3",
    Quiet:          false,
}
err = overlay.CreateWithOptions(ctx, opts)
```

### Resize

```go
// Grow to 5000 MB
err := overlay.Resize(ctx, "env.img", 5000)

// Shrink also supported (checks free space first)
err = overlay.Resize(ctx, "env.img", 2000)
```

### Change Ownership

```go
// Change to root (requires fakeroot)
err := overlay.ChownToRoot(ctx, "env.img")

// Change to specific user
err = overlay.ChownToUser(ctx, "env.img", 1000, 1000)
```

### Integrity Check

```go
// Check integrity
err := overlay.Check(ctx, "env.img", repair=false)

// Check and repair
err = overlay.Check(ctx, "env.img", repair=true)

// Verify availability (lock check)
err = overlay.CheckAvailable("env.img", writeLock=false)
err = overlay.CheckAvailable("env.img", writeLock=true)
```

### Information

```go
// Get overlay stats
stats, err := overlay.GetStats("env.img")
// stats.Usage() returns (usedBytes, percent)
// stats.InodeUsage() returns percent

// Check if SquashFS (from utils package)
isSqf := utils.IsSqf("file.sqf")
```

### Locking

```go
// Acquire shared lock (multiple readers)
lock, err := overlay.AcquireLock("env.img", exclusive=false)
defer lock.Release()

// Acquire exclusive lock (single writer)
lock, err := overlay.AcquireLock("env.img", exclusive=true)
defer lock.Release()
```

## Filesystem Profiles

Choose profile based on use case:

| Profile | Inode Ratio | Use Case |
|---------|-------------|----------|
| Small | 4KB | Conda environments, Python site-packages |
| Default | 16KB | General applications |
| Large | 1MB | Reference genomes, large databases |

## Resize Algorithm

1. Check current size via `stat`
2. If growing: `truncate` + `e2fsck` + `resize2fs`
3. If shrinking: Check free space first, then `e2fsck` + `resize2fs` + `truncate`
4. Verify final size

## Ownership Changes

Uses fakeroot + Apptainer to:
1. Mount overlay as writable
2. Run `chown -R UID:GID /ext3` inside container
3. Unmount

## Locking Strategy

- `.img` files are locked during use to prevent concurrent writes
- Shared locks (read-only): Multiple instances can mount same overlay
- Exclusive locks (writable): Only one process can write
- Lock files: `.img.lock` in same directory

## Error Types

- `LockError` - Cannot acquire lock (overlay in use)
- `IntegrityError` - Filesystem corruption detected
- `ResizeError` - Resize operation failed
- `OwnershipError` - Chown operation failed

## Integration

- `internal/container` - Checks overlay availability before mounting
- `internal/build` - Creates temporary overlays for builds
- `cmd/overlay.go` - CLI operations (create, resize, chown, check)
