# Helper Scripts

Helper scripts launch interactive web applications (e.g. RStudio, desktop VNC) inside CondaTainer. They can run on both HPC clusters (via scheduler job submission) and headless servers (directly).

## How It Works

When you run `condatainer helper <name>` (e.g. `rstudio-server`), the helper generates a job script that runs the web application inside CondaTainer and launches it:

- **On HPC**, it submits the job to the scheduler, which dispatches it to a compute node.
- **On a headless server**, it runs the service directly.

A local **CondaTainer server** is auto-started to handle the connection to the helper service. You forward **one port** for it — the same port serves both the dashboard and every launched helper app, so there is no per-helper port juggling.

```
Your Browser → http://<id>.localhost:<server-port>/
                      ↓ SSH port forward (one-time setup)
               CondaTainer server  (login node/normal server)
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
condatainer helper --update # or -u
```

### HPC Scheduler

On HPC clusters, helpers submit batch jobs. Verify your scheduler is detected:

```bash
condatainer scheduler
```

```{note}
*Slurm* is fully tested. *PBS*, *LSF*, and *HTCondor* are experimentally supported.
```

## Quick Start

```bash
# Update helpers
condatainer helper -u

# Run a helper (submits scheduler job on HPC or runs directly on server without scheduler)
condatainer helper vscode-server
```

The command prints an access URL like `http://<id>.localhost:13182/`. To open it in your browser, set up port forwarding as described next.

The printed port will be saved for future reuse.

## Access From Your Browser

The CondaTainer server (and every helper app it launches) is served on a port on the HPC login node. Forward that port over SSH and your local browser can reach the server.

The helper you just ran printed this port in its URL (`13182` in `http://<id>.localhost:13182/`). You can also check it with `condatainer server status`. 

### SSH Port Forwarding (on your local machine)

After setting up port forwarding, any traffic you send to `localhost:<port>` on your computer is forwarded over SSH and delivered to `localhost:<port>` on the HPC side.

There are two ways: 
- Use `ssh -L <port>:localhost:<port>` for a one-off connection;
- Add a `LocalForward` entry to `~/.ssh/config` so it applies every time you connect:

```
Host hpc
    HostName hpc.university.edu
    User your_username
    # Use the port the helper printed (13182 here)
    LocalForward 13182 localhost:13182
```

Then you can use `ssh hpc` for connection.

Now you can open http://localhost:13182 in your browser.

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

To resume after reconnecting, restart the server explicitly:

```bash
condatainer server start # start the server (and use the saved port)
condatainer helper       # or run a helper — it auto-starts the server
```

If you want the server to survive logout, start it with `--daemon`:

```bash
condatainer server start --daemon
```
````

## Managing Helpers

There are two ways to start, monitor, and stop helpers. They act on the same sessions — use whichever is convenient.

- **[From the Web Dashboard](#from-the-web-dashboard)** — point-and-click in your browser.
- **[From the CLI](#from-the-cli)** — run `condatainer helper` in your terminal.

### From the Web Dashboard

Open `http://localhost:<server-port>/` in your browser (after [forwarding the port](#ssh-port-forwarding-on-your-local-machine)). The sidebar has five sections. The three you use to manage helpers — **Start**, **Jobs**, and **History** — are covered below. The other two are general-purpose tools:

- **Overlay** — create and remove overlays (see [Overlays](../user_guide/concepts.md) for the concept).
- **Files** — browse, upload, download, and manage files on the cluster.

#### Web Start

Launch a new helper session on a compute node (or directly on a headless server):

1. Pick a helper from the list;
2. Set its resources (CPUs, memory, walltime, GPU);
3. Choose an env overlay and any module overlays;
4. then click **Start**. 

**Defaults** opens a dialog to change the default settings.

```{image} assets/dashboard-start.png
:alt: CondaTainer dashboard — Start tab
:width: 100%
```

#### Web Jobs

Monitor and open your running (and pending) sessions.

Each card shows the helper name and status, node, working directory, env overlay, and remaining walltime. From a card you can:

- **Open** the access URL, or **copy** the link;
- **Stop** the session.

```{image} assets/dashboard-jobs.png
:alt: CondaTainer dashboard — Jobs tab
:width: 100%
```

Each job is reachable through two URLs:

- `<id>.localhost` — points to that **specific** helper job (the full unique id).
- `<helper-name>.localhost` — the 🔗 button opens this **stable** URL, which always points to the *first running* helper job of that name.

E.g. with these jobs running:

```
# <name>-YYYYMMDD-HHMMSS
myapp-20260722-143015     <- myapp.localhost resolves here (first running "myapp")
another-20260722-150458   <- another.localhost resolves here
myapp-20260722-161032
```

Here `myapp.localhost` points to `myapp-20260722-143015` (the earliest still-running one).

**Why the stable `<helper-name>.localhost` URL matters:** browser-based apps like code-server (VS Code) persist some of their data (settings, layout, ...) in the browser keyed by **origin (the hostname)**. If you connect through the per-id URL, each new session gets a fresh hostname and therefore a fresh workspace. Using the stable `<helper-name>.localhost` hostname keeps the same origin across sessions, so your app state carries over from one run to the next.

#### Web History

Review and rerun past sessions.

The table lists each run's name, status (`running` / `done` / `failed`), node, working directory (CWD), and when it ended. From here you can:

- **Rerun** a session, reusing its settings;
- **Delete** a single row;
- **Clear finished** to prune all completed rows.

```{image} assets/dashboard-history.png
:alt: CondaTainer dashboard — History tab
:width: 100%
```

### From the CLI

Launch or reattach to a helper by name; the command prints the access URL:

```bash
condatainer helper vscode-server      # show picker: launch or show status
condatainer helper                    # show status of all running helpers
condatainer helper vscode-server stop # stop this helper's sessions
condatainer helper stop --all         # stop every running helper
```

#### Resource Flags

These control the scheduler job resources on HPC. They have no effect in headless mode.

| Flag | Description | Default |
|------|-------------|---------|
| `-c, --cpus` | CPUs per job | (helper default) |
| `-m, --mem` | Memory | (helper default) |
| `-t, --time` | Walltime | (helper default) |
| `-g, --gpu` | GPU spec (e.g. `a100:1`, `h100:2`) | none |

Pass `-g ''` to explicitly clear a saved GPU setting.

#### Overlay Flags

These apply in all modes (HPC and headless).

| Flag | Description |
|------|-------------|
| `-b, --base` | Override base image |
| `-e, --env` | Writable overlay (default: `env.img`) |
| `-o, --overlay` | Additional read-only overlay (repeatable) |
| `-w, --cwd <path>` | Set working directory (e.g. `-w .` for current directory) |
| `--new` | Skip reuse prompt, force new session |

When `-e` is unset, the helper searches these directories in order and uses the first overlay found, preferring a per-user `env-$USER.img` over a shared `env.img` in each:

1. current directory
2. `overlay/`
3. `src/overlay/`

Use `-e -` or `-e ''` to disable auto search.

#### Helper-Specific Flags

Each helper may define its own flags (e.g. `--rversion` for `rstudio-server`). View them with:

```bash
condatainer helper <name> --help
```

#### Session History

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
  - Pick a **running** session prints its access URL and exits.
  - Pick a **pending/starting** one waits for the service and prints the URL once it is ready.
  - Pick a **past** session reuses its settings (working directory, resources, overlays) for the new job.
  - Choose `n` to start fresh with config defaults and the current directory.

Pass `--new` to skip the session picker and go straight to this prompt with fresh defaults.

#### Settings Confirmation

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

Use `-` to clear a field (e.g. `e -` to run without a writable overlay, `g -` to clear GPU, `o -` to clear extra overlays).

For `o`, the value replaces the entire overlay list — space-separate multiple names (e.g. `o build-essential extra-deps`).

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

## Available Applications

| Helper | Description | Guide |
|--------|-------------|-------|
| `code-server` | Open-source VS Code | [VS Code](./vscode.md#code-server) |
| `vscode-server` | Full VS Code in browser | [VS Code](./vscode.md#vscode-server) |
| `vscode-tunnel` | Full VS Code via Microsoft relay | [VS Code](./vscode.md#vscode-tunnel) |
| `rstudio-server` | RStudio with Posit R image overlays | [RStudio Server](./rstudio-server.md) |
| `rstudio-server-conda` | RStudio with Conda-managed R | [RStudio Server](./rstudio-server.md#rstudio-server-conda-conda-r) |
| `jupyterlab` | Jupyter Lab | — |
| `xfce4` | XFCE desktop via VNC | [XFCE Desktop](./xfce4.md) |
| `igv` | XFCE with IGV pre-launched | [XFCE Desktop](./xfce4.md#igv) |
| `cytoscape` | XFCE with Cytoscape pre-launched | [XFCE Desktop](./xfce4.md#cytoscape) |

## Potential Issues

### Clusters Without Inter-Node SSH

If your cluster does not allow SSH between nodes, enable `helper_bind_all`:

```bash
condatainer config set helper_bind_all true
```

Helper services will then bind to all interfaces on the compute node. The compute node must be reachable from the login node for this to work.

### Server has scheduler but not functional

If you wants helper explitly on the same machine, set `submit_job: false` in your CondaTainer config to make headless mode permanent.

```bash
condatainer config set submit_job false
```

then run `server` or `helper`
