# Custom Helper Scripts

A helper script turns any long-running service into a one-command app on HPC.

CondaTainer handles the parts that are tedious to get right: picking a free port on the compute node, submitting the scheduler job, stacking the overlays the service needs, and routing your browser to it through a single SSH tunnel. Your script only has to **start the service and say when it's ready**.

See [Helper Scripts](../helpers/helpers.md) for using the helpers that already ship with CondaTainer. This page is about writing your own.

## Before You Write One

Ask what your service needs to run, because that decides the shape of the script:

| Your service is… | Use | Example |
|---|---|---|
| May keep installing more packages | `#IMG_PACKAGES:` a **writable** env overlay | [JupyterLab](#example-1-jupyterlab-on-a-writable-conda-overlay) |
| A fixed tool stack you don't modify | `#REQUIRED_OVERLAYS:` **read-only** module overlays | [Cytoscape](#example-2-cytoscape-an-app-not-on-conda) |
| Not packaged for Conda at all | A [build script](./custom_app.md) first, then `#REQUIRED_OVERLAYS:` | [Cytoscape](#example-2-cytoscape-an-app-not-on-conda) |

The two can be combined: rstudio-server needs R from image (`#REQUIRED_OVERLAYS:`) and also env overlays (`#IMG_PACKAGES:`) for holding packages.

## Where Scripts Live

```bash
condatainer helper --path      # writable helper-scripts/ directory
```

Drop the file there (executable, no extension) and it becomes a helper immediately:

```bash
condatainer helper --list      # your script now appears, with its #WHATIS: text
condatainer helper my-service
```

## Lifecycle

Understanding what runs where explains most of the rules below:

1. **Login node** — CondaTainer resolves `#PARAM:` values, checks/installs the overlays your script declares, and submits the job.
2. **Compute node** — a wrapper exports the [injected variables](#environment-variables), enters the container with your overlays mounted, redirects all output to `$CNT_HELPER_STATE_DIR/job.log`, and sources your script.
3. Your script binds the service to `$CNT_HELPER_PORT` and calls `condatainer _server_ready`.
4. **Login node** — the CondaTainer server sees the ready file over NFS and starts proxying `http://<id>.localhost:<server-port>/` to the compute node.
5. Your script **must not exit** while the service is up — the job ends when the script does.

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

Three rules this illustrates:

- Bind to `${CNT_HELPER_BIND_ADDR:-127.0.0.1}`, never a hardcoded address — the fallback keeps the script working when `helper_bind_all` is enabled for clusters without inter-node SSH.
- Bind to `$CNT_HELPER_PORT`, never a fixed port — several users share each compute node.
- Keep the script in the foreground (`wait $PID`, or `sleep infinity` for services you can't `wait` on).

## Metadata Headers

### Resource defaults

```bash
#NCPUS: 4
#MEM:   16G        # bare number = GB
#TIME:  12:00:00   # also: 12h, 4d, 1-12:00:00; bare number = hours
#GPU:   a100:1     # optional
```

These are defaults only — users override them with `-c/-m/-t/-g` or in the settings prompt.

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

The resolved value is exported as `$KEY` inside the script, and `{KEY}` tokens in other headers are substituted with it.

### `#VALUE:` — allowed values list

```bash
#VALUE: VERSION=3.13,3.12,3.11          # comma-separated, sorted newest-first
#VALUE: MODE=light | full | minimal     # pipe-separated, order preserved
```

Used for the interactive prompt selector and `--help` display.

### `#REQUIRED_OVERLAYS:` — read-only overlays

Checked and auto-installed before the job starts. Declaration order is the Apptainer overlay stack order — **the last name is the topmost layer** and wins on file conflicts (e.g. `/var/lib/dpkg/status`).

```bash
#REQUIRED_OVERLAYS: my-app/{VERSION} build-essential
```

### `#IMG_PACKAGES:` — writable Conda overlay

Signals that a writable Conda overlay is needed. The value is the default package list for guided creation; omit the value to require an overlay without a package check.

```bash
#IMG_PACKAGES: my-app={VERSION} numpy
#IMG_PACKAGES:                          # overlay required, no package check
```

If the user runs the helper without `-e`, CondaTainer offers to create the overlay and install this list.

### `#POST_INSTALL_CMD:` — post-install hook

One shell command run inside the container after guided overlay creation, for pinning or setup:

```bash
#POST_INSTALL_CMD: mm-pin r-base
```

### `#CHECK_PATH:` — pre-submission binary check

Verifies the binary exists *before* burning queue time. Fatal when `#IMG_PACKAGES:` is set, a warning otherwise. One path per line, with an optional message:

```bash
#CHECK_PATH: /ext3/env/bin/jupyter-lab "Jupyter Lab not found — install with: mm-install jupyterlab"
```

### `#BIND:` — extra bind mounts

`src:dest` or `src:dest:opts` format. `$VAR` references are expanded on the compute node, so any injected variable works.

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
| `CNT_HELPER_ID` | Unique run ID (`{name}-{job_id}`) |
| `CNT_HELPER_STATE_DIR` | NFS directory for runtime files (readable from login node) |
| `CNT_HELPER_SCRIPT_DIR` | Directory containing the helper script |
| `CNT_HELPER_PORT` | Free TCP port on the compute node — bind your service here |
| `CNT_HELPER_BIND_ADDR` | `127.0.0.1` normally; `0.0.0.0` when `helper_bind_all` is set |
| `CNT_HELPER_BIND_ALL` | `1` when `helper_bind_all` is set; unset otherwise |
| `CNT_HELPER_CWD` | Working directory chosen by the user |
| `CNT_HELPER_WALLTIME_SECS` | Walltime in seconds |
| `CNT_JOB_TMPDIR` | Node-local scratch dir (for Unix sockets, etc.); cleaned up on exit |
| `$KEY` | One variable per `#PARAM:` key |

```{tip}
`$CNT_HELPER_STATE_DIR` is on NFS and survives the job — use it for config the login node must read, or for state that should persist across sessions. `$CNT_JOB_TMPDIR` is node-local and faster — use it for sockets, caches, and anything disposable.
```

## Server Hidden Functions

### `condatainer _server_ready`

Call once the service is accepting connections. Writes the state file the CLI monitor and dashboard read; until it lands, the session shows as `[STARTING]`.

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

### `condatainer _pick_port`

Returns another free port, for services that need a second listener (a backend the frontend proxies to). Check it against `$CNT_HELPER_PORT` and re-pick if they collide.

---

## Example 1: JupyterLab on a Writable Conda Overlay

This is the common case for analysis work: you start with Jupyter and a Python version, then keep adding packages — `scanpy` this week, `scvi-tools` next — without rebuilding anything. That is exactly what a **writable env overlay** (`env.img`) is for, and `#IMG_PACKAGES:` is how a helper asks for one.

Source: [helpers/jupyterlab](https://github.com/Justype/cnt-scripts/blob/main/helpers/jupyterlab)

```bash
#!/bin/bash
#WHATIS: Jupyter Lab
#NCPUS: 4
#MEM: 16G
#TIME: 12:00:00

#IMG_PACKAGES: python={CONDA_PYTHON} jupyterlab
#VALUE: CONDA_PYTHON=3.13.13,3.13.12,3.12.13,3.12.12,...
#AUTOUPDATE:CONDA_PYTHON:conda-forge:python>=3.8.0

TOKEN=$(python3 -c "import secrets; print(secrets.token_hex(8))" 2>/dev/null \
    || openssl rand -hex 8 2>/dev/null \
    || echo "$(date +%s)$RANDOM")

jupyter-lab \
    --no-browser \
    --ip="${CNT_HELPER_BIND_ADDR:-127.0.0.1}" \
    --port "$CNT_HELPER_PORT" \
    --IdentityProvider.token "$TOKEN" \
    ${CNT_HELPER_BIND_ALL:+--ServerApp.allow_remote_access=True} &
PID=$!

condatainer _server_ready \
    --port "$CNT_HELPER_PORT" \
    --url-path "?token=$TOKEN"

wait $PID
```

### Why `#IMG_PACKAGES:` and not `#REQUIRED_OVERLAYS:`

`#IMG_PACKAGES: python={CONDA_PYTHON} jupyterlab` does two things:

- **Declares the overlay writable.** Read-only `.sqf` overlays can't take a new package; an `env.img` can. Other packages can be installed later.
- **Drives guided creation.** Launch without `-e` and CondaTainer offers to build an overlay with this package list.

`{CONDA_PYTHON}` is just a placeholder:

- Only used for env overlay creation: it's a **creation-time** choice, changing it won't affect an existing overlay.
- `#VALUE:` gives the available options (the version list offered at creation).
- `#AUTOUPDATE:` is for the CI workflow: it auto-updates the Python version list from conda-forge.

### Details worth copying

**How the variables start the service and route the browser.** The service binds to `--ip="${CNT_HELPER_BIND_ADDR:-127.0.0.1}"` and `--port "$CNT_HELPER_PORT"`, the address and free port CondaTainer picked. `_server_ready` then tells the server where the service is, so your browser lands straight on the running app.

**Add auth when the app supports it.** The service is reachable by anyone on the compute node, so a token is worth having — but this is app-specific.

**`${CNT_HELPER_BIND_ALL:+...}`** adds `--ServerApp.allow_remote_access=True` *only* when `helper_bind_all` is on. Jupyter refuses non-local connections.

### The workflow it enables

```bash
condatainer helper jupyterlab       # first run: offers to create env.img with jupyterlab
```

Then, from a terminal inside JupyterLab (or `condatainer exec -w -o env.img bash`):

```bash
mm-install scanpy leidenalg
```

---

## Example 2: Cytoscape, an App Not on Conda

Cytoscape is a desktop **GUI** app and isn't on Conda. To reach it from a browser, the helper:

1. **Packages the app** into an overlay (custom script).
2. **Runs a KasmVNC server on top of an XFCE4 desktop**, so the GUI has a desktop to draw on and KasmVNC turns that desktop into the web page.
3. **Starts PulseAudio (server + client)** so sound from the desktop reaches the browser.

Source: [helpers/cytoscape](https://github.com/Justype/cnt-scripts/blob/main/helpers/cytoscape)

### 1. Package the app

Cytoscape ships as a vendor tarball on GitHub, so it needs a [custom build script](./custom_app.md#example-2-cytoscape--a-version-template) — see that page for the full walkthrough. The result is a module overlay you can install by name:

```bash
condatainer create cytoscape/3.10.4
```

### 2. Stack the overlays

```bash
#REQUIRED_OVERLAYS: xfce4 openjdk/{CONDA_OPENJDK} cytoscape/{CYTOSCAPE}
#BIND: $CNT_HELPER_STATE_DIR/vnc:$HOME/.vnc

#PARAM: CYTOSCAPE=? --cytoscape,-C "Cytoscape version"
#VALUE: CYTOSCAPE=3.10.4,3.10.3,3.10.2,...

#PARAM: CONDA_OPENJDK=? --openjdk,-j "OpenJDK version (must be 17.x)"
#VALUE: CONDA_OPENJDK=17.0.18,17.0.17,...
```

#### One script, multiple versions

`{KEY}` and `#PARAM: KEY` makes one helper serves many app versions. Take `CYTOSCAPE`:

- `#PARAM: CYTOSCAPE=? ...` accepts a version and exports it as `$CYTOSCAPE`. The `=?` means "default to the first `#VALUE:` entry", so a plain launch just works.
- The value is resolved at submit time and substituted into `cytoscape/{CYTOSCAPE}` in `#REQUIRED_OVERLAYS:`. Pass `-C 3.9.1` and the helper stacks the `cytoscape/3.9.1` overlay instead — no edit to the script.

#### Choose the right directory to bind

The VNC server writes its config, password, and `vncserver.log` to `$HOME/.vnc`. The `#BIND:` line redirects that at `$CNT_HELPER_STATE_DIR/vnc` (on NFS), so each run gets its own copy (two sessions don't fight over one `$HOME/.vnc`) and the log stays readable from the login node.

Match the destination to what the data is:

- **Keep after the job finishes** (config, logs) → `$CNT_HELPER_STATE_DIR`. It's per-run (not shared across sessions) but on NFS, so the login node can read it and it survives after the job ends.
- **Disposable / node-local** (Unix sockets, caches, scratch) → `$CNT_JOB_TMPDIR`. Faster, and cleaned up when the job exits.

### 3. Make the app behave in a container

This part is script specific. Fix things when app is not work as you expected. A few worth showing:

**Derive `JAVA_HOME` from the overlay.** Cytoscape's launcher won't find Java on its own, and the path depends on which OpenJDK overlay got mounted, so resolve it at runtime rather than hardcoding:

```bash
_JAVA_BIN=$(command -v java 2>/dev/null)
if [ -n "$_JAVA_BIN" ]; then
    export JAVA_HOME="$(dirname "$(dirname "$(readlink -f "$_JAVA_BIN")")")"
    condatainer _server_message "Java: JAVA_HOME=$JAVA_HOME"
else
    condatainer _server_message "WARNING: java not found in PATH; Cytoscape may fail to start"
fi
```

**Keep stdin open.** Launched from a desktop entry with `Terminal=false`, Cytoscape's Karaf console reads EOF on stdin and shuts down instantly. Feeding it a pipe that never closes fixes it:

```bash
cat > "$CNT_JOB_TMPDIR/cytoscape-launch.sh" <<'EOF'
#!/bin/bash
exec 0< <(sleep infinity)
exec cytoscape
EOF
chmod +x "$CNT_JOB_TMPDIR/cytoscape-launch.sh"
```

**Add a desktop entry so it shows in the XFCE4 menu.** Write a `.desktop` file into an XDG applications dir under `$CNT_JOB_TMPDIR` (node-local and disposable), pointing `Exec` at the wrapper above. `Terminal=false` is why stdin needed the fix; it's what the app launches with from the menu:

```bash
mkdir -p "$CNT_JOB_TMPDIR/share/applications"
cat > "$CNT_JOB_TMPDIR/share/applications/cytoscape.desktop" <<EOF
[Desktop Entry]
Version=1.0
Type=Application
Name=Cytoscape
Comment=Network Biology Visualization Platform
Exec=$CNT_JOB_TMPDIR/cytoscape-launch.sh
Terminal=false
Categories=Science;Biology;
EOF
```

For XFCE4 to find it, `$CNT_JOB_TMPDIR/share` must be on `XDG_DATA_DIRS`.

**Never assume a port or display is yours.** Nodes are shared, and two launches can race for the same X display, so the script scans for one that works instead of failing on the first collision:

```bash
_port_open() { # Check if the port is available
    (exec 3<>"/dev/tcp/127.0.0.1/$1") 2>/dev/null
}

DISPLAY_NUM=""
for i in $(seq 10 999); do
    [ -f "/tmp/.X${i}-lock" ] && continue
    _port_open $((5900 + i)) && continue
    if vncserver ":${i}" ... >>"$CNT_HELPER_STATE_DIR/vncserver.log" 2>&1; then
        DISPLAY_NUM=$i
        break
    fi
    rm -f "/tmp/.X${i}-lock" "/tmp/.X11-unix/X${i}"
done
[ -z "$DISPLAY_NUM" ] && { condatainer _server_message "ERROR: No free X display"; exit 1; }
export DISPLAY=":${DISPLAY_NUM}"
```

The same applies to the second listener — the VNC web client needs its own port, so it asks for one and re-picks on collision:

```bash
KASM_INTERNAL_PORT=$(condatainer _pick_port)
[ "$KASM_INTERNAL_PORT" = "$CNT_HELPER_PORT" ] && KASM_INTERNAL_PORT=$(condatainer _pick_port)
```

### 4. Hand credentials to the browser via `--url-path`

The whole VNC session is reachable in one click because the ready URL carries the noVNC query string, password included:

```bash
condatainer _server_ready \
    --port "$CNT_HELPER_PORT" \
    --url-path "vnc.html?autoconnect=1&resize=remote&password=${VNC_PASSWD}"

sleep infinity
```

`sleep infinity` rather than `wait` — the desktop session has no single PID worth waiting on. The job ends on walltime or `condatainer helper cytoscape stop`.

---

## Saving Runtime State

For params computed at runtime (e.g. auto-generated names), compute on first run and save to config so subsequent runs reuse the value:

```bash
if [ -z "${NAME:-}" ]; then
    NAME="$(whoami)-$(hostname -s | tr '.' '-')"
    condatainer helper "$CNT_HELPER_NAME" config set NAME "$NAME"
fi
```

## Debugging

Everything your script writes to stdout/stderr goes to `job.log` in the run's state directory:

```bash
ls ~/.local/state/condatainer/helper/           # one dir per run ID
tail -f ~/.local/state/condatainer/helper/<id>/job.log
```

That directory also holds `ready`, `done`, and `messages` — plus whatever your script writes there itself, which is why service logs are worth pointing at `$CNT_HELPER_STATE_DIR` (as the Cytoscape helper does with `vncserver.log`).

Two habits that shorten the loop:

- **Test the payload outside the helper first.** `condatainer exec -w -o env.img -o <module overlays> bash` puts you in the same container the helper gets, without waiting in the queue.
- **Use `_server_message` liberally.** Messages surface live in the terminal and dashboard, so you learn *where* a slow start is stuck rather than only that it never got ready.

## Related

- [Helper Scripts](../helpers/helpers.md) — using the built-in helpers
- [Sharing Your Scripts](../deployment/share_scripts.md) — upstreaming or hosting your own source
- [Custom App Build Scripts](./custom_app.md) — packaging an app that isn't on Conda
- [Writable Environment Overlays](../user_guide/environment_overlays.md) — the `env.img` model
- [internal/helper/README.md](https://github.com/Justype/condatainer/blob/main/internal/helper/README.md) — full header reference and NFS state file format
