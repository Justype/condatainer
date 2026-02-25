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

- Shared locks: multiple readers can hold concurrently
- Exclusive locks: single writer; blocks all other locks
- Lock files stored as `<image>.lock` alongside the image
