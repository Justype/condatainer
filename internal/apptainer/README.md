# apptainer

Go wrapper for Apptainer (formerly Singularity) container runtime providing binary management, ephemeral execution, persistent instances, and image building.

## Architecture

```
apptainer.go    Binary setup, version detection, command execution
error.go        Structured error types with exit codes and hints
exec.go         Ephemeral container execution (exec)
instance.go     Persistent instance management (start/stop/exec)
build.go        Image building and SIF-to-SquashFS conversion
```

## Key Types

**ApptainerError** - Structured execution error with `Op`, `Cmd`, `Path`, `Output`, `BaseErr`
- `ExitCode()` - Returns shell exit code or -1
- `Error()` - Auto-analyzes output and provides hints (corrupted image, no space, permissions, etc.)

**BuildOptions** - `Force`, `NoCleanup`, `Additional`
**ExecOptions** - `Bind`, `Overlay`, `Fakeroot`, `Env`, `HideOutput`, `Additional`, `Stdin`
**InstanceStartOptions** - `Bind`, `Overlay`, `Fakeroot`, `Additional`
**InstanceExecOptions** - `Env`, `Additional`, `Stdin`

## Usage

```go
// Binary setup
apptainer.SetBin("")        // detect from PATH (tries apptainer, then singularity)
apptainer.EnsureApptainer() // error if not found
apptainer.IsSingularity()   // true when binary is singularity (use gzip compression)

// Ephemeral execution
apptainer.Exec(ctx, "image.sif", []string{"cmd"}, &apptainer.ExecOptions{
    Overlay:  []string{"overlay.sqf"},
    Fakeroot: true,
})

// Persistent instance
apptainer.InstanceStart(ctx, "base.sif", "name", startOpts)
apptainer.InstanceExec(ctx, "name", []string{"cmd"}, execOpts)
apptainer.InstanceStop(ctx, "name", nil)
apptainer.InstanceList(ctx)

// Building
apptainer.Build(ctx, "out.sif", "def.def", opts)
apptainer.DumpSifToSquashfs(ctx, "in.sif", "out.sqf")
```

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
- Context cancellation: SIGTERM (5s wait) → SIGKILL
- Version cached after first `GetVersion()`, invalidated on `SetBin()` change
