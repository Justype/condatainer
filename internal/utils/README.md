# utils

Shared utilities for console output, file operations, downloads, and script parsing.

## Architecture

```
console.go   Styled console output, logging, user prompts
download.go  HTTP downloads with progress bars
files.go     File and directory utilities
filelock.go  Bare non-blocking file locking primitive
parser.go    Build script metadata parsing
```

## Locking

`AcquireFileLock(path, write)` is the one non-blocking file-lock primitive in the
codebase — `open` with a mode picked by `write` (`O_RDWR`/exclusive vs.
`O_RDONLY`/shared, since some callers use the open mode itself as a
permission check ahead of and independent from the lock), then a non-blocking
whole-file `fcntl` open-file-description lock (`F_OFD_SETLK`, Linux 3.15+). It is that kind, not
`flock(2)`, because its job is to see an ext3 image a running container has mounted: Apptainer takes a
plain `fcntl` lock there, which a `flock` conflicts with only on filesystems that fold the two together
(NFSv3 does, NFSv4 does not). The descriptor owns the lock, so two acquisitions in one process conflict
and closing the file releases it, as with `flock`. It carries no domain meaning of its own: `internal/image.Lock`/
`AcquireLock` wraps it to add `ErrProtected`/`ErrInUse` classification for
overlay images (a lock sentinel that opens `O_RDWR` and fails as
"write-protected" specifically means an artifact was `chmod a-w`'d to pin
it); `internal/libexec` calls it directly for its own toolchain-generation
sentinel, which has no such "protected" concept. It lives here, not in
`internal/image`, because `internal/libexec` needs the identical mechanism
and must not depend on `internal/image` (or anything that depends on it) —
see `internal/libexec/README.md` and `internal/image/lock.go`'s own comment
for why. `internal/image.Lock` is a type alias for `*utils.FileLock`
rather than a wrapping struct, so both packages' locks are the exact same
type wherever a caller holds them together (`internal/runtime/
exec.Prepare`'s combined overlay + toolchain lock slice).

## Console Output

**Print:** `PrintMessage`, `PrintSuccess`, `PrintWarning`, `PrintError`, `PrintHint`, `PrintNote`, `PrintDebug`

**Style:** `StyleName` (yellow), `StylePath`, `StyleAction`, `StyleCommand`, `StyleHint` (cyan), `StyleError` (red), `StyleSuccess` (green), `StyleWarning` (yellow)

**Modes:** `DebugMode`, `QuietMode`, `YesMode`

```go
utils.PrintMessage("Installing %s", utils.StyleName("package"))
utils.PrintSuccess("Build completed")
utils.ShouldAnswerYes()                      // true when YesMode or non-interactive
utils.ReadLineContext(ctx context.Context)    // read line with cancellation
```

## File Operations

**Checks:** `FileExists`, `DirExists`, `IsImg`, `IsSqf`, `IsSif`, `IsOverlay`
**Operations:** `EnsureDir`
**Permissions:** `PermFile` (0664), `PermDir` (0775), `PermExec` (0775)

`ShareWithParentGroup` only ever fixes the one path it's given — every wired creation helper
(`CreateFileWritable`, `MkdirAllShared`, `MakeExecutable`) calls it on exactly what that helper
itself just created, nothing more. A tree written by something else entirely — an external tool's
own installer, not this package's own creation helpers — needs `ShareTreeWithParentGroup(root)`
instead: it walks top-down and applies the single-path fix at every level, since
`ShareWithParentGroup` only acts once a path's own parent is already shared. `internal/libexec`'s
`provisionAt` is the one caller — `micromamba create` writes every file under its own prefix
directly, bypassing every helper in this package, so nothing below the top level comes out
group-writable on its own.

## Downloads

```go
utils.DownloadFile(url, destPath)       // with progress bar
utils.DownloadExecutable(url, destPath) // sets exec permissions
```

## Script Parsing

These read *user* scripts and helper scripts, not recipes — a recipe is parsed by
`catalog.ParseRecipe`.

```go
// #DEP: only — a "module load" line names the site's module tree, not an artifact
deps, err := utils.GetDependenciesFromScript(scriptPath)

description := utils.GetDescriptionFromScript(scriptPath)

// #TYPE: from an external build script, defaulting to "app"
typ, err := utils.GetTypeFromScript(scriptPath)

// extract scheduler directives and apply defaults (scheduler package)
specs, err := scheduler.ReadScriptSpecsFromPath(scriptPath)
```

### Version Helpers

```go
// Sort version strings descending (newest first). Returns a new slice.
sorted := utils.SortVersionsDescending(versions)

// Render a #TARGET: pattern with {placeholder} tokens in bold-yellow.
styled := utils.HighlightTemplatePlaceholders(pattern)
```

Names, versions and dependency constraints live in `catalog`, so one string
resolves one way everywhere: `catalog.Normalize`, `catalog.CompareVersions`,
`catalog.ParseDep` and `Dep.Satisfies`. Template matching and `#PH:`/`#VALUE:`
value lists likewise: `catalog.NewTemplate` and `catalog.ParseValues`.
