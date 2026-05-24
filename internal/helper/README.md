# internal/helper

Orchestrates the full lifecycle of a helper service job: resolve params, check/create overlays, submit to scheduler (or run headless), monitor NFS state files, and clean up.

## Files

| File | Purpose |
|------|---------|
| `run.go` | `RunHelper()` — top-level entry point; also server auto-start, scheduler detection, wrapper generation |
| `state.go` | NFS state files (`ready`, `done`, `messages`), JSONL history, `HelperRun` struct |
| `monitor.go` | `PollOnce()` — shared non-blocking read of all three NFS state files |
| `params.go` | `ApplyParamFlags()`, `PromptMissingParams()` — resolve `#PARAM:` headers |
| `config.go` | Per-helper saved config (`~/.local/share/condatainer/helper-config/<name>.conf`) |

---

## Helper Script Headers

Helper scripts are plain bash files with metadata in comment headers. All headers are parsed before the job is submitted.

### Resource headers

Parsed by `scheduler.ReadScriptSpecs()`. Standard scheduler directives also work (`#SBATCH`, `#PBS`, `#BSUB`).

```bash
#NCPUS: 4          # default CPU count
#MEM:   16G        # default memory       — bare number defaults to GB (16 = 16G)
#TIME:  12:00:00   # default walltime     — bare number defaults to hours (12 = 12h)
                   # also accepts: 12h, 2h30m, 4d, 1-12:00:00
#GPU:   a100:1     # optional GPU spec (type:count)
```

### `#WHATIS:` — short description

Shown in `condatainer helper --list` and the server dashboard.

```bash
#WHATIS: Jupyter Lab
```

### `#PARAM:` — helper-specific parameters

Defines a configurable variable. The Go runner resolves each param from (in order): CLI flag → saved config → interactive prompt → default value.

```
#PARAM: KEY=default --long-flag,-s "Description shown in --help and prompts"
```

- `KEY=default` — variable name and default value. Empty default means the user must supply a value.
- `--long-flag,-s` — optional CLI flags (long and short). Omit entirely for prompt-only params.
- `"Description"` — shown in `condatainer helper <name> --help` and the interactive prompt.

The resolved value is exported as `$KEY` inside the job wrapper, so the script sees it as a plain env var.

**First-run save pattern** — for params that must be computed at runtime (e.g. tunnel names), leave the default empty and compute+save inside the script on first run:

```bash
if [ -z "${MACHINE_NAME:-}" ]; then
    MACHINE_NAME="$(whoami)-$(hostname -s)"
    condatainer helper "$CNT_HELPER_NAME" config set MACHINE_NAME "$MACHINE_NAME"
fi
```

Subsequent runs: Go reads the saved config and injects `$MACHINE_NAME` before the script runs.

### `#VALUE:` — allowed values / version list

Provides the ordered list of valid values for a `#PARAM:` key, used for the interactive prompt selector and the web UI dropdown. Updated by `auto.py`.

```bash
#VALUE: CONDA_PYTHON=3.13.0,3.12.7,3.11.9,...
```

Two formats:
- **Comma-separated** `KEY=val1,val2,...` — version list, **sorted descending** automatically (newest first). Order in the script doesn't matter.
- **Pipe-separated** `KEY=opt1 | opt2 | opt3` — ordered option labels, **not sorted**. Order in the script is preserved.

### `#IMG_PACKAGES:` — conda packages for guided overlay creation

Presence signals that a writable conda-env overlay is required. Value is the default package list for guided creation; `{KEY}` tokens are substituted from resolved `#PARAM:` values.

```bash
#IMG_PACKAGES: python={CONDA_PYTHON} jupyterlab
```

- **With overlay (`-e` given)**: Go runs a pre-submission install check (fatal if packages not found).
- **Without overlay**: Go prompts to create one (guided overlay creation flow).

### `#POST_INSTALL_CMD:` — post-install hook

A single shell command run inside the container after guided overlay creation installs the packages. Used to pin package versions or run setup steps.

```bash
#POST_INSTALL_CMD: mm-pin r-base
```

### `#CHECK_PATH:` — pre-submission binary check

One path per line. Checked inside the overlay (fatal) when `#IMG_PACKAGES:` is set, or in the base image (non-fatal warning) otherwise. Optional quoted message shown when the path is missing.

```bash
#CHECK_PATH: /ext3/env/bin/jupyter-lab "Jupyter Lab not found — install with: mm-install jupyterlab"
#CHECK_PATH: /ext3/env/bin/R "Conda R not found — install with: mm-install r-base=<version>"
```

Multiple lines accumulate independently (one path checked per line).

### `#REQUIRED_OVERLAYS:` — named SquashFS overlays

Space-separated overlay names. `{KEY}` tokens are substituted from resolved params. Go checks each overlay exists on disk and installs any that are missing before submitting the job.

```bash
#REQUIRED_OVERLAYS: r{POSIT_R} rstudio-server build-essential
# resolves to: r4.4.3 rstudio-server build-essential
```

### `#BIND:` — extra bind mounts

One bind spec per line, in `src:dest` or `src:dest:opts` form (same as Apptainer `--bind`).
`{KEY}` tokens are substituted from resolved `#PARAM:` values. `$VAR` references are left
verbatim and expanded by the shell when the wrapper runs on the compute node — any of the
[injected env vars](#environment-variables-injected-by-the-wrapper) can be used directly.

```bash
#BIND: $CNT_HELPER_STATE_DIR/vnc:$HOME/.vnc
#BIND: /shared/data:/data:ro
#BIND: /path/to/{KEY_DIR}:/mnt/custom
```

### `#SINGLETON:` — at most one running instance

When `true`, the "start new instance" option is suppressed if an instance is already running. `--new`/`--force` prints an error instead.

```bash
#SINGLETON: true
```

---

## Environment Variables Injected by the Wrapper

The Go runner injects these into every job before the helper script body runs:

| Variable | Value |
|----------|-------|
| `CNT_HELPER_NAME` | Helper script name (e.g. `jupyterlab`) |
| `CNT_HELPER_ID` | Unique run ID: `{name}-{job_id}` or `{name}-{pid}` |
| `CNT_HELPER_PORT` | Free TCP port on the compute node (helper must bind here) |
| `CNT_HELPER_CWD` | Working directory |
| `CNT_HELPER_STATE_DIR` | NFS state directory for this run — write runtime files here |
| `CNT_HELPER_WALLTIME_SECS` | Walltime in seconds |
| `CNT_HELPER_JOB_ID` | Scheduler job ID (empty for headless) |
| `CNT_HELPER_SCRIPT_DIR` | Directory containing the helper script |
| `SCRATCH` | `$SCRATCH` or `$HOME` if unset |
| `$KEY` | One var per resolved `#PARAM:` key |

---

## NFS State Files

Written by the script on the compute node; read by the CLI monitor and server watcher on the login node.

```
~/.local/state/condatainer/helper/{id}/
    ready       — JSON written by _server_ready when service is accessible
    done        — JSON written by the wrapper on exit (always fires, even on crash)
    messages    — append-only JSONL written by _server_message
    job.log     — scheduler stdout/stderr (streamed by server dashboard)
```

### `ready` JSON

```json
{"port": 8888, "node": "cn01", "label": "Jupyter Lab", "timestamp": "...",
 "walltime_secs": 43200, "job_id": "12345",
 "url_path": "?token=abc123", "external_url": ""}
```

`port: 0` + non-empty `external_url` = no proxy tunnel needed (e.g. vscode-tunnel).

### `done` JSON

```json
{"exit_code": 0, "ts": "2024-01-01T12:00:00Z"}
```

### `messages` JSONL

One JSON object per line:

```json
{"level": "info", "text": "Open project: rstudioapi::openProject(...)", "ts": "..."}
{"level": "warn", "text": "Low memory warning", "ts": "..."}
{"level": "error", "text": "Port already in use", "ts": "..."}
```

---

## Hidden Commands Used by Scripts

Scripts call these condatainer subcommands from inside the container:

### `condatainer _server_ready`

Writes the `ready` file and appends to JSONL history.

```bash
condatainer _server_ready \
    --port "$CNT_HELPER_PORT" \
    --label "Jupyter Lab" \
    --url-path "?token=$TOKEN"

# For external-URL services (no port forwarding):
condatainer _server_ready \
    --port 0 \
    --label "VS Code Tunnel: $MACHINE_NAME" \
    --external-url "https://vscode.dev/tunnel/$MACHINE_NAME$CNT_HELPER_CWD"
```

### `condatainer _server_message`

Appends one line to the `messages` file. Shown in the terminal and the server dashboard.

```bash
condatainer _server_message "Open project: rstudioapi::openProject(\"$CNT_HELPER_CWD\")"
condatainer _server_message --level warn "Low disk space"
condatainer _server_message --level error "Failed to start"
```

### `condatainer helper <name> config`

Reads and writes the per-helper saved config. Callable from inside the container (NFS-accessible).

```bash
condatainer helper "$CNT_HELPER_NAME" config set MACHINE_NAME "$MACHINE_NAME"
condatainer helper "$CNT_HELPER_NAME" config get MACHINE_NAME
condatainer helper "$CNT_HELPER_NAME" config show   # all keys
condatainer helper "$CNT_HELPER_NAME" config path   # config file path
```

---

## Minimal Example Script

```bash
#!/bin/bash
#WHATIS: My Custom Service
#NCPUS: 2
#MEM: 8G
#TIME: 4:00:00
#PARAM: MODE=light --mode,-m "Interface mode: light or full"

my-service \
    --mode "$MODE" \
    --bind 127.0.0.1:"$CNT_HELPER_PORT" &
PID=$!

condatainer _server_ready \
    --port "$CNT_HELPER_PORT" \
    --label "My Custom Service ($MODE)"

wait $PID
```

## Example: conda-env service (with overlay)

```bash
#!/bin/bash
#WHATIS: My App (conda)
#IMG_PACKAGES: my-app={VERSION}
#CHECK_PATH: /ext3/env/bin/my-app "my-app not found — install with: mm-install my-app"
#NCPUS: 4
#MEM: 16G
#TIME: 12:00:00
#VALUE: VERSION=2.0,1.9,1.8

my-app --port "$CNT_HELPER_PORT" &
PID=$!

condatainer _server_ready --port "$CNT_HELPER_PORT" --label "My App $VERSION"
wait $PID
```

## Example: external-URL service (no port forwarding)

```bash
#!/bin/bash
#WHATIS: My Tunnel Service
#SINGLETON: true
#NCPUS: 2
#MEM: 4G
#TIME: 24:00:00
#PARAM: NAME= --name,-n "Tunnel name (auto-set on first run)"

if [ -z "${NAME:-}" ]; then
    NAME="$(whoami)-$(hostname -s | tr '.' '-')"
    condatainer helper "$CNT_HELPER_NAME" config set NAME "$NAME"
fi

my-tunnel --name "$NAME" 2>&1 | while IFS= read -r line; do
    echo "$line"
    case "$line" in
        *"device code"*) condatainer _server_message "$line" ;;
        *"tunnel ready"*)
            condatainer _server_ready \
                --port 0 \
                --label "My Tunnel: $NAME" \
                --external-url "https://my-service.example.com/tunnel/$NAME" ;;
    esac
done
```
