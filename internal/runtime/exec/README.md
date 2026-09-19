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

## Usage

### Basic Execution

```go
import "github.com/Justype/condatainer/internal/runtime/exec"

ctx := context.Background()
opts := exec.Options{
    BaseImage: "/path/to/base.sqf",
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

1. **Apply defaults** - Fill in missing configuration (command; not the base
   image or the apptainer binary — both need overlay setup first, see below)
2. **Container setup** - Call `container.Setup()` for:
   - Overlay resolution and ordering, and pulling out an exec root if requested
   - Environment collection
   - Bind path deduplication
   - GPU detection
3. **Resolve base image** - `Root` from step 2 wins if present, else the
   caller's `BaseImage`, else `config.GetBaseImage()`
4. **Auto-enable fakeroot** - If writable .img overlay
5. **Resolve the apptainer binary** - `apptainer.Fakeroot()` or
   `apptainer.Normal()`, now that fakeroot is final — see that package's README for which binary and why
6. **Debug output** - Print configuration if debug mode
7. **Print environment** - Show overlay environments (if interactive and not hidden)
8. **Acquire file locks** - Hold shared locks on all `.sqf` overlays and the base image for the duration of execution. `.img` overlays are skipped — Apptainer flocks them itself and our lock would conflict. Prevents concurrent `remove` or `build --update` from deleting files in use.
9. **Inject proxy env** - If an active SOCKS5 proxy is found via `proxy.FindActiveProxy()`, prepend `http_proxy`/`https_proxy`/`all_proxy` (and uppercase variants) to the container environment so tools inside the container use the tunnel.
10. **Wrap for activation** - If `container.ActivationScript()` is non-empty (an overlay's own `etc/conda/activate.d` needs sourcing — see that package's README, *Activation*), `wrapWithActivation` replaces `Command` with `bash -c <activation + exec "$@"> cnt-activate <original command...>`. A static `--env` list can't express this: it's shell script, not values.
11. **Execute** - Call `apptainer.Exec()` with processed configuration
12. **Release locks** - All file locks released after `apptainer.Exec()` returns

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
- `BaseImage` → `Setup`'s `Root` if the requested overlays name one, else
  `config.GetBaseImage()`
- `Fakeroot` → `false` (auto-enabled if needed)

There is no `ApptainerBin` field at all: which binary runs is never
caller-configurable, only decided from the final `Fakeroot` value, by
`apptainer.Fakeroot` / `apptainer.Normal` — see that package's README.

The base image is required: there is no overlay-only execution, so a container
with no root cannot start, and `Prepare` fails rather than letting Apptainer
report a missing file. It must exist — that is the only check; nothing here or
in `internal/build` reads a manifest's type before accepting something as
root. Building a missing default root is the caller's job — this package runs
images, it does not make them, and `internal/build` is what knows how.

## Stdin Forwarding

When `PassThruStdin: true`, forwards stdin to the container for interactive scripts:
```go
opts := exec.Options{
    Command:       []string{"bash", "script.sh"},
    PassThruStdin: true,
}
```

Used by build system for interactive build scripts.

## Conda Overlay Creation (conda.go)

`CreateCondaOverlay` builds a new writable overlay by installing packages into
a scratch ext3 image (`ext3.CreateInTmp`) and only moving it to its final
destination once installation succeeds. That scratch path is disconnected
from wherever the final destination's paired snapshot lives, so
`container.Setup`'s own autoload (which looks beside the path it is given)
can never find it there. `CreateCondaOverlay` looks it up itself
(`container.LookupSnapshot(opts.Path)` — against the *final* path, not the
scratch one) and passes it through to `InitCondaEnv`/`RunPostInstall`, which
mount it alongside the scratch image. Without this, installing into a `.img`
that will end up paired with an existing snapshot would reinstall everything
the snapshot already has instead of writing only the incremental diff.

`InstallPackages`/`RemovePackages` (adding or removing packages in an
already-placed `.img`) need no such lookup: they run against the real final
path directly, so `container.Setup`'s ordinary autoload already finds the
pairing on its own.

## Integration

Used by:
- `cmd/run.go` - Direct command execution
- `cmd/exec.go` - Alias to run
- `internal/build` - Build script execution
