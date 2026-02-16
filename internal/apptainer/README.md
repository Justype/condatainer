# apptainer

Go wrapper for Apptainer (formerly Singularity) container runtime providing binary management, ephemeral execution, persistent instances, and image building with automatic base image provisioning.

## Architecture

```
apptainer.go    Binary setup, version detection, command execution
error.go        Structured error types with exit codes and hints
exec.go         Ephemeral container execution (exec)
instance.go     Persistent instance management (start/stop/exec)
build.go        Image building and SIF-to-SquashFS conversion
ensure.go       Base image provisioning (download/build)
```

## Key Types

**ApptainerNotFoundError** - Binary cannot be located  
**ApptainerError** - Structured execution error with `Op`, `Cmd`, `Path`, `Output`, `BaseErr`
- `ExitCode()` - Returns shell exit code or -1
- `Error()` - Auto-analyzes output and provides hints (corrupted image, no space, permissions, etc.)

**BuildOptions** - `Force`, `NoCleanup`, `Additional`  
**ExecOptions** - `Bind`, `Overlay`, `Fakeroot`, `Env`, `HideOutput`, `Additional`, `Stdin`  
**InstanceStartOptions** - `Bind`, `Overlay`, `Fakeroot`, `Additional`  
**InstanceExecOptions** - `Env`, `Additional`, `Stdin`

## Usage

### Binary and Version

```go
apptainer.SetBin("")              // Detect from PATH
apptainer.Which()                 // Current binary path
apptainer.GetVersion()            // Cached version
apptainer.CheckZstdSupport(ver)   // >= 1.4 check
```

### Execution

```go
// Ephemeral
opts := &apptainer.ExecOptions{
    Bind: []string{"/host:/container"},
    Overlay: []string{"overlay.sqf"},
    Fakeroot: true,
}
apptainer.Exec(ctx, "image.sif", []string{"cmd"}, opts)

// Persistent instance
apptainer.InstanceStart(ctx, "base.sif", "name", startOpts)
apptainer.InstanceExec(ctx, "name", []string{"cmd"}, execOpts)
apptainer.InstanceStop(ctx, "name", nil)
apptainer.InstanceList(ctx)
apptainer.InstanceStats(ctx, "name")
```

### Building

```go
apptainer.Build(ctx, "out.sif", "def.def", opts)
apptainer.IsBuildCancelled(err)           // Ctrl+C detection
apptainer.DumpSifToSquashfs(ctx, "in.sif", "out.sqf")
```

### Base Image Provisioning

```go
// Ensure exists (download or build)
apptainer.EnsureBaseImage(ctx, update=false, downloadOnly=false)
```

**Logic:**  
1. Check Apptainer binary availability
2. If update mode: download to `.new`, replace old
3. Otherwise: search all image paths
4. If not found: download prebuilt (unless debug mode)
5. If download fails and !downloadOnly: build locally from def

Prebuilt platforms: `x86_64`  
Download URL: `{config.PrebuiltBaseURL}/base_image_{arch}.sif`

## Error Handling

```go
var apErr *apptainer.ApptainerError
if errors.As(err, &apErr) {
    apErr.ExitCode()  // Shell exit code
}

if apptainer.IsBuildCancelled(err) {
    // User Ctrl+C or declined overwrite
}
```

## Implementation Notes

**Command Execution:**
- `runApptainer()` - Standard execution
- `runApptainerWithOutput()` - Custom output control (stream/capture/hide)
- Build operations unset `SINGULARITY_BIND`/`APPTAINER_BIND` to avoid `%post` mount conflicts
- Context cancellation: SIGTERM (5s wait) → SIGKILL

**Version Caching:** Cached after first `GetVersion()`, invalidated on `SetBin()` change

**HPC Considerations:** Always use `--fakeroot` for non-root users
