# instance

Persistent Apptainer instance management with state file persistence for configuration tracking.

## Architecture

```
start.go    Instance startup with container setup
exec.go     Command execution in running instances
stop.go     Instance termination
state.go    State file persistence (JSON)
```

## Key Types

**Options** - `Name`, `Overlays`, `WritableImg`, `EnvSettings`, `BindPaths`, `Fakeroot`, `BaseImage`  
**ExecOptions** - `InstanceName`, `Command`, `ApptainerBin`  
**StateFile** - Persisted config: `Name`, `BaseImage`, `Overlays`, `StartTime`, `PID`

## Usage

```go
// Start instance
opts := instance.Options{
    Name:        "myinstance",
    Overlays:    []string{"python/3.11", "env.img"},
    WritableImg: true,
}
instance.Start(ctx, opts)

// Execute in instance
instance.Exec(ctx, instance.ExecOptions{
    InstanceName: "myinstance",
    Command:      []string{"python", "script.py"},
})

// Stop instance
instance.Stop(ctx, "myinstance")
```

## State Persistence

State files: `$XDG_STATE_HOME/condatainer/instances/<name>.json`

Stores instance configuration for recovery and tracking. Auto-cleaned on normal stop.

## Lifecycle

**Start:** Container setup → `apptainer.InstanceStart()` → save state  
**Exec:** Load state → `apptainer.InstanceExec()`  
**Stop:** `apptainer.InstanceStop()` → delete state

## Instance List/Stats

For listing and stats operations, commands call `internal/apptainer` directly:
- `apptainer.InstanceList(ctx)` - List all running instances
- `apptainer.InstanceStats(ctx, name)` - Show instance statistics

These operations don't require state file access.
