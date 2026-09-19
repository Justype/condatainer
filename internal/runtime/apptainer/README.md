# apptainer

Go wrapper for Apptainer (formerly Singularity) container runtime providing binary management, ephemeral execution, and image building.

## Architecture

```
apptainer.go    Binary setup, version detection, command execution
error.go        Structured error types with exit codes and hints
exec.go         Ephemeral container execution (exec)
build.go        Sandbox and image building
```

## Key Types

**ApptainerError** - Structured execution error with `Op`, `Cmd`, `Path`, `Output`, `BaseErr`
- `ExitCode()` - Returns shell exit code or -1
- `Error()` - Auto-analyzes output and provides hints (corrupted image, no space, permissions, etc.)

**BuildOptions** - `Force`, `NoCleanup`, `Sandbox`, `TmpDir`, `Additional`
**ExecOptions** - `Bind`, `Overlay`, `Fakeroot`, `Env`, `HideOutput`, `Additional`, `Stdin`

## Usage

```go
// Binary setup
apptainer.SetBin("")        // detect from PATH (tries apptainer, then singularity)
apptainer.EnsureApptainer() // error if not found
apptainer.IsSingularity()   // true when the configured binary is singularity

// Ephemeral execution
apptainer.Exec(ctx, "base.sqf", []string{"cmd"}, &apptainer.ExecOptions{
    Overlay:  []string{"overlay.sqf"},
    Fakeroot: true,
})

// Building. A definition writes a sandbox directory, which the caller then packs.
apptainer.Build(ctx, "work/rootfs", "def.def", &apptainer.BuildOptions{
    Sandbox: true,
    TmpDir:  "/tmp/cnt-user/build",
})
```

## Which binary: `ResolveBin`

There is no single "the" apptainer binary — which one runs is decided per
invocation, by `internal/runtime/exec.Prepare`, once fakeroot is final (an
explicit `--fakeroot`, or `container.AutoEnableFakeroot`'s auto-enable for a
writable `.img` with root-owned files):

- **Fakeroot** — `ResolveBin(true)` — always the system/module binary
  (`config.Global.Build.SystemApptainer`). Only its setuid starter (or a module's) can
  escalate privilege; libexec's own apptainer is deliberately non-setuid. It
  is also version-checked: condatainer packs every artifact as zstd
  unconditionally (`config.LoadDefaults`'s `CompressArgs`), so a system
  apptainer below the zstd floor (`>= 1.4`, `CheckZstdSupport`) or Singularity
  (assumed zstd-incapable regardless of version) is refused outright, with a
  message naming the requirement — rather than failing later, unrecognizably,
  inside Apptainer's own mount step.
- **Everything else** — `ResolveBin(false)` — `internal/libexec.ApptainerPath()`
  when an apptainer is installed there, so installing one is the user's way to
  choose it over the host's. Otherwise the system/module binary, under the same
  zstd-floor and Singularity checks as fakeroot. When neither works,
  `ResolveBin` refuses and names both ways out: `condatainer update --libexec
  apptainer`, or loading an apptainer module.

A `.def` (`os`) build's own `apptainer build --fakeroot` does not go
through `ResolveBin` at all — `internal/build/def.go` resolves the system
binary directly (`EnsureApptainer`) with no zstd check, because its output is
a sandbox: an `os` build never mounts a zstd-compressed artifact
during the build itself (`#DEP:` is data-only). Packing that sandbox
afterward never runs in a container at all any more —
`internal/build/squashfs.go`'s `createSquashfs` runs `mksquashfs` directly on
the host (`toolpath.Resolve`), so it does not go through `ResolveBin` either.

### `Current`: reading back what already ran

`internal/build`'s `captureCommonBuildTools` records which binary produced a
build, but it must never decide that itself — resolving independently could
name a different binary than the one the build's own container step actually
used. `Current` reports the implementation and version of whichever binary is
already configured (`Implementation`/`GetVersion` under the hood) and errors
if nothing has been resolved yet, rather than falling back to PATH the way
`SetBin("")` does.

### Apptainer needs squashfs-tools in PATH

Apptainer needs `unsquashfs`/`mksquashfs` in `PATH` for some of its own operations (e.g. extracting
a `.sqf` to build a sandbox) — this is Apptainer's own subprocess lookup, unrelated to the `--env`/
`EnvSettings` a launched container's PATH is set from (`internal/libexec/README.md`'s "not a `PATH`
prepend" rule is about that, different environment). `runApptainerWithOutput` adds
`internal/libexec`'s `bin/` to `PATH` whenever the resolved binary is libexec's own copy, so
Apptainer finds the squashfs-tools libexec provisioned right next to it. Without it: `exec:
"unsquashfs": executable file not found in $PATH`.

## Error Handling

```go
var apErr *apptainer.ApptainerError
if errors.As(err, &apErr) {
    apErr.ExitCode()
}

if apptainer.IsBuildCancelled(err) {
    // User Ctrl+C or declined overwrite
}
```

## Implementation Notes

- Build operations unset `SINGULARITY_BIND`/`APPTAINER_BIND` to avoid `%post` mount conflicts
- Definition builds pass `--fix-perms`, so an OCI base's owner-unreadable files and no-write directories
  are packed with their content and removable at cleanup
- Context cancellation: SIGTERM (5s wait) → SIGKILL
- Version cached after first `GetVersion()`, invalidated on `SetBin()` change
