# Custom Helper Scripts

Helper scripts are plain bash files with metadata headers that CondaTainer reads before running. You can write your own to launch any service inside a container.

Place the script in your local `helper-scripts/` directory (shown by `condatainer helper --path`) and run it with:

```bash
condatainer helper my-service
```

## Minimal Example

```bash
#!/bin/bash
#WHATIS: My Custom Service
#NCPUS: 2
#MEM: 8G
#TIME: 4:00:00

my-service \
    --host "${CNT_HELPER_BIND_ADDR:-127.0.0.1}" \
    --port "$CNT_HELPER_PORT" &
PID=$!

condatainer _server_ready --port "$CNT_HELPER_PORT"

wait $PID
```

## Metadata Headers

### Resource defaults

```bash
#NCPUS: 4
#MEM:   16G
#TIME:  12:00:00   # also: 12h, 4d, 1-12:00:00
#GPU:   a100:1     # optional
```

### `#WHATIS:` — description

Shown in `condatainer helper --list` and the dashboard.

```bash
#WHATIS: Jupyter Lab
```

### `#PARAM:` — script-specific flags

Defines a configurable variable resolved from: CLI flag → saved config → default.

```bash
#PARAM: KEY=default --long-flag,-s "Description shown in --help"
```

| Default form | Behavior |
|---|---|
| `KEY=value` | Fixed default, overridable via flag or saved config |
| `KEY=?` | Uses first entry of `#VALUE:` list automatically |
| `KEY=` | Required — prompts interactively if not supplied |

The resolved value is exported as `$KEY` inside the script.

### `#VALUE:` — allowed values list

```bash
#VALUE: VERSION=3.13,3.12,3.11          # comma-separated, sorted newest-first
#VALUE: MODE=light | full | minimal     # pipe-separated, order preserved
```

Used for the interactive prompt selector and `--help` display.

### `#REQUIRED_OVERLAYS:` — SquashFS overlays

Checked and auto-installed before the job starts. `{KEY}` tokens are substituted from resolved `#PARAM:` values. Declaration order is the Apptainer overlay stack order (last = topmost layer).

```bash
#REQUIRED_OVERLAYS: my-app/{VERSION} build-essential
```

### `#IMG_PACKAGES:` — Conda overlay required

Signals that a writable Conda overlay is needed. The value is the default package list for guided creation; omit the value to require an overlay without a package check.

```bash
#IMG_PACKAGES: my-app={VERSION} numpy
#IMG_PACKAGES:                          # overlay required, no package check
```

### `#BIND:` — extra bind mounts

`src:dest` or `src:dest:opts` format. `$VAR` references are expanded at runtime.

```bash
#BIND: $CNT_HELPER_STATE_DIR/config:/app/config
#BIND: /shared/data:/data:ro
```

### `#SINGLETON:` — single instance only

Prevents starting a second instance while one is running.

```bash
#SINGLETON: true
```

## Environment Variables

The Go runner injects these before your script runs:

| Variable | Value |
|----------|-------|
| `CNT_HELPER_NAME` | Script name |
| `CNT_HELPER_STATE_DIR` | NFS directory for runtime files (readable from login node) |
| `CNT_HELPER_SCRIPT_DIR` | Directory containing the helper script |
| `CNT_HELPER_PORT` | Free TCP port on the compute node — bind your service here |
| `CNT_HELPER_BIND_ADDR` | `127.0.0.1` normally; `0.0.0.0` when `helper_bind_all` is set |
| `CNT_HELPER_BIND_ALL` | `1` when `helper_bind_all` is set; unset otherwise |
| `CNT_JOB_TMPDIR` | Node-local scratch dir (for Unix sockets, etc.); cleaned up on exit |
| `SCRATCH` | `$SCRATCH` or `$HOME` if unset |
| `$KEY` | One variable per `#PARAM:` key |

Always use `${CNT_HELPER_BIND_ADDR:-127.0.0.1}` as the bind address so the script works in both normal and `helper_bind_all` modes.

## Signalling Ready

### `condatainer _server_ready`

Call once the service is accepting connections. Writes the state file that the CLI monitor reads.

```bash
# Standard service (port-based)
condatainer _server_ready \
    --port "$CNT_HELPER_PORT" \
    --url-path "?token=$TOKEN"

# External URL (no port forwarding — e.g. tunnel services)
condatainer _server_ready \
    --port 0 \
    --external-url "https://example.com/tunnel/$NAME"
```

### `condatainer _server_message`

Send status lines to the terminal and dashboard while the job is running.

```bash
condatainer _server_message "Service started, loading data..."
condatainer _server_message --level warn "Low disk space"
condatainer _server_message --level error "Failed to start"
```

## Full Example: Versioned Conda App

```bash
#!/bin/bash
#WHATIS: My App (conda-managed)
#NCPUS: 4
#MEM:   16G
#TIME:  12:00:00
#IMG_PACKAGES: my-app={VERSION}
#PARAM: VERSION=? --version,-v "App version"
#VALUE: VERSION=2.1,2.0,1.9

condatainer _server_message "Starting my-app $VERSION..."

my-app \
    --host "${CNT_HELPER_BIND_ADDR:-127.0.0.1}" \
    --port "$CNT_HELPER_PORT" \
    --workdir "$CNT_HELPER_CWD" &
PID=$!

condatainer _server_ready \
    --port "$CNT_HELPER_PORT" \
    --url-path "/"

wait $PID
```

## Saving Runtime State

For params computed at runtime (e.g. auto-generated names), compute on first run and save to config so subsequent runs reuse the value:

```bash
if [ -z "${NAME:-}" ]; then
    NAME="$(whoami)-$(hostname -s | tr '.' '-')"
    condatainer helper "$CNT_HELPER_NAME" config set NAME "$NAME"
fi
```

## Reference

For the full header reference and NFS state file format, see [internal/helper/README.md](https://github.com/Justype/condatainer/blob/main/internal/helper/README.md).
