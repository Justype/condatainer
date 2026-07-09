# Helper Scripts

Helper scripts launch interactive web applications (VS Code, RStudio, Jupyter, desktop VNC) inside CondaTainer. They run on both HPC clusters (via scheduler job submission) and headless servers (directly).

## How It Works

When you run `condatainer helper ...`, a local **CondaTainer server** is auto-started on the current machine (login node or server). This server:

- Connects to helper services running on compute nodes
- Exposes a **web dashboard** at `http://localhost:<server-port>/` — (start, stop, configure helpers)
- Gives each helper session its own URL: `http://<session-id>.localhost:<server-port>/`

You forward **one port** for the CondaTainer server — no per-helper port juggling.

```
Your Browser → http://<id>.localhost:<server-port>/
                      ↓ SSH port forward (one-time setup)
               CondaTainer server  (login/server node)
                      ↓ (automatic)
               Helper service  (compute node or same server)
```

## Prerequisites

### Install CondaTainer

```bash
curl -fsSL https://get-condatainer.justype.net | bash
```

### Download Helper Scripts

```bash
condatainer helper -u
```

### SSH Port Forwarding

Set up a **single** `LocalForward` for the CondaTainer server port. The server prints its port on first start; you can also check it with `condatainer server status`.

```
# HPC is a shared system, do not use common port like 8080, 8787
Host hpc
    HostName hpc.university.edu
    User your_username
    LocalForward 13182 localhost:13182
```

To use a fixed port, start the server once with `-p` — the choice is saved for all future starts:

```bash
condatainer server start -p 13182
```

Then set your `LocalForward` to that port. Everything else is handled automatically.

````{tip}
Add the following at the end of your local `~/.ssh/config` to prevent idle timeouts while working in the browser:

```
Host *
    ServerAliveInterval 60
    ServerAliveCountMax 5
```
````

````{note}
**The server exits when your SSH session ends.** This is by design: on HPC systems with multiple login nodes, the server must run on the same node your SSH tunnel lands on. Tying the server to the SSH session ensures it is always reachable through your forwarded port.

To resume after reconnecting, simply run any helper command or restart the server explicitly:

```bash
condatainer server start # restart the server
condatainer helper       # or just run a helper — it auto-starts the server
```

If you want the server to survive logout, start it with `--daemon`:

```bash
condatainer server start --daemon
```
````

### HPC Scheduler

On HPC clusters, helpers submit batch jobs. Verify your scheduler is detected:

```bash
condatainer scheduler
```

> Slurm is fully tested. PBS, LSF, and HTCondor are experimentally supported.

### Clusters Without Inter-Node SSH

If your cluster does not allow SSH between nodes, enable `helper_bind_all`:

```bash
condatainer config set helper_bind_all true
```

Helper services will then bind to all interfaces on the compute node. The compute node must be reachable from the login node for this to work.

## Quick Start

```bash
# Update helpers
condatainer helper -u

# Run a helper (submits scheduler job on HPC, runs directly on a server)
# CondaTainer server auto-starts on first run and prints the access URL
condatainer helper vscode-server
```

On a headless server without a scheduler, use `--local`:

```bash
condatainer --local helper -u
condatainer --local helper vscode-server
```

Or set `submit_job: false` in your CondaTainer config to make headless mode permanent.

## Flags

### Resource Flags

These control the scheduler job resources on HPC. They have no effect in headless mode.

| Flag | Description | Default |
|------|-------------|---------|
| `-c, --cpus` | CPUs per job | (helper default) |
| `-m, --mem` | Memory | (helper default) |
| `-t, --time` | Walltime | (helper default) |
| `-g, --gpu` | GPU spec (e.g. `a100:1`, `h100:2`) | none |

Pass `-g ''` to explicitly clear a saved GPU setting.

### Overlay Flags

These apply in all modes (HPC and headless).

| Flag | Description |
|------|-------------|
| `-b, --base` | Override base image |
| `-e, --env` | Writable overlay (default: `env.img`) |
| `-o, --overlay` | Additional read-only overlay (repeatable) |
| `-w, --cwd <path>` | Set working directory (e.g. `-w .` for current directory) |
| `--new` | Skip reuse prompt, force new session |

When `-e` is `env.img` (the default), the helper searches in order:

1. `env.img` (current directory)
2. `overlay/env.img`
3. `src/overlay/env.img`

### Helper-Specific Flags

Each helper may define its own flags (e.g. `--rversion` for `rstudio-server`, `--vnc` for `xfce4`). View them with:

```bash
condatainer helper <name> --help
```

## Configuration

Helpers save your preferences to `~/.config/condatainer/helper/<name>` (or `$XDG_CONFIG_HOME/condatainer/helper/<name>`).

```bash
# View saved config and available keys
condatainer helper vscode-server config

# Set a value
condatainer helper vscode-server config set NCPUS 8

# View raw config file
condatainer helper vscode-server config show
```

Config files use `KEY="VALUE"` format and can be edited directly. The next run picks up any changes.

### Session History

When you run a helper, CondaTainer shows your recent sessions (up to 5, deduplicated by working directory) and lets you pick one:

```
Recent vscode-server sessions:
  1. /scratch/user/project-a   4c 32G        [RUNNING] started 2h ago
     → http://<id>.localhost:13182/...
  2. /scratch/user/project-b   4c 16G        [PENDING] submitted 5m ago
  3. /scratch/user/project-c   8c 64G        last used 3d ago
  n. Start new

[?] Choose [1/2/3/n]:
```

- Active sessions are tagged with their scheduler state: 
  - `[PENDING]` (queued, no URL yet);
  - `[STARTING]` (job started, service coming up)
  - `[RUNNING]` (ready, URL shown). 
- Picking
  - Pick a **running** session prints its access URL and exits; picking a **pending/starting** one waits for the service and prints the URL once it is ready.
  - Pick a **past** session reuses its settings (working directory, resources, overlays) for the new job.
  - Choose `n` to start fresh with config defaults and the current directory.

Pass `--new` to skip the session picker and go straight to this prompt with fresh defaults.

### Settings Confirmation

Before submitting a new job, CondaTainer prints all current settings and lets you edit them inline:

```
[CNT] Settings:
  c: 4     CPUs   m: 16G  Memory   t: 12:00:00  Walltime   g: (none)  GPU
  w: /scratch/user/project-a      Working directory
  e: env.img                      Writable overlay (.img)
  r: 4.5.2                        R version

[Enter to continue, or key value to update]:
```

Type `key value` to update a field before submitting (e.g. `c 8`, `g h100:1`, `r 4.4`, `w /scratch/project`), then press Enter on an empty line to proceed.

Use `-` to clear a field (e.g. `e -` to run without a writable overlay, `g -` to clear GPU, `o -` to clear extra overlays). For `o`, the value replaces the entire overlay list — space-separate multiple names (e.g. `o build-essential extra-deps`).

## Available Applications

| Helper | Description | Guide |
|--------|-------------|-------|
| `vscode-server` | Full VS Code in browser | [VS Code](./vscode.md#vscode-server) |
| `vscode-tunnel` | Full VS Code via Microsoft relay (no port forwarding) | [VS Code](./vscode.md#vscode-tunnel) |
| `code-server` | Open-source VS Code (Open VSX extensions only) | [VS Code](./vscode.md#code-server) |
| `rstudio-server` | RStudio with Posit R image overlays | [RStudio Server](./rstudio-server.md) |
| `rstudio-server-conda` | RStudio with Conda-managed R | [RStudio Server](./rstudio-server.md#rstudio-server-conda-conda-r) |
| `jupyterlab` | Jupyter Lab | — |
| `xfce4` | XFCE desktop via VNC, browser-accessible | [XFCE Desktop](./xfce4.md) |
| `igv` | XFCE desktop with IGV pre-launched | [XFCE Desktop](./xfce4.md#igv) |
| `cytoscape` | XFCE desktop with Cytoscape pre-launched | [XFCE Desktop](./xfce4.md#cytoscape) |
