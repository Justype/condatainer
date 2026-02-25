# exec

Ephemeral container execution orchestrating Apptainer with the container setup pipeline.

## Architecture

```
options.go   Options struct with defaults
run.go       Main Run() execution function
```

## Key Types

**Options:**
- `BaseImage` - Base container image path
- `Overlays` - Overlay paths/names
- `WritableImg` - Allow writing to .img overlay
- `EnvSettings` - Environment variables
- `BindPaths` - Bind mounts
- `Fakeroot` - Use fakeroot
- `ApptainerFlags` - Additional flags
- `HideOutput` - Suppress output
- `HidePrompt` - Hide environment notes
- `PassThruStdin` - Forward stdin to container
- `Command` - Command to execute
- `ApptainerBin` - Path to apptainer binary

## Usage

### Basic Execution

```go
import "github.com/Justype/condatainer/internal/exec"

ctx := context.Background()
opts := exec.Options{
    BaseImage: "/path/to/base.sif",
    Overlays:  []string{"cellranger/9.0.1"},
    Command:   []string{"cellranger", "--version"},
}

if err := exec.Run(ctx, opts); err != nil {
    // Handle error
}
```

### Interactive Shell

```go
opts := exec.Options{
    Overlays:    []string{"python/3.11", "env.img"},
    WritableImg: true,
    Command:     []string{"bash"},
}

exec.Run(ctx, opts)
```

### With Environment and Binds

```go
opts := exec.Options{
    Overlays:    []string{"app/1.0"},
    EnvSettings: []string{"DEBUG=1", "WORKERS=8"},
    BindPaths:   []string{"/data:/mnt/data", "/scratch"},
    Command:     []string{"python", "train.py"},
}

exec.Run(ctx, opts)
```

## Execution Flow

1. **Apply defaults** - Fill in missing configuration
2. **Initialize Apptainer** - Set binary path
3. **Container setup** - Call `container.Setup()` for:
   - Overlay resolution and ordering
   - Environment collection
   - Bind path deduplication
   - GPU detection
4. **Auto-enable fakeroot** - If writable .img overlay
5. **Debug output** - Print configuration if debug mode
6. **Print environment** - Show overlay environments (if interactive and not hidden)
7. **Acquire file locks** - Hold shared locks on all `.sqf` overlays and the base `.sif` for the duration of execution. `.img` overlays are skipped — Apptainer flocks them itself and our lock would conflict. Prevents concurrent `remove` or `build --update` from deleting `.sqf`/`.sif` files in use.
8. **Execute** - Call `apptainer.Exec()` with processed configuration
9. **Release locks** - All file locks released after `apptainer.Exec()` returns

## Environment Display

When running interactively, displays environment variables from overlays:
```
Overlay envs:
  CELLRANGER_ROOT: /ext3/cnt/cellranger/9.0.1
  PATH: /ext3/cnt/cellranger/9.0.1/bin:$PATH
```

Can be hidden with `HidePrompt: true`.

## Defaults

Missing fields are filled from config:
- `BaseImage` → `config.GetBaseImage()`
- `ApptainerBin` → `config.Global.ApptainerBin`
- `Fakeroot` → `false` (auto-enabled if needed)

## Stdin Forwarding

When `PassThruStdin: true`, forwards stdin to the container for interactive scripts:
```go
opts := exec.Options{
    Command:       []string{"bash", "script.sh"},
    PassThruStdin: true,
}
```

Used by build system for interactive build scripts.

## Integration

Used by:
- `cmd/run.go` - Direct command execution
- `cmd/exec.go` - Alias to run
- `internal/build` - Build script execution
