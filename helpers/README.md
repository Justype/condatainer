# CondaTainer Helper Scripts

Helper scripts for launching interactive services (RStudio Server, VS Code, etc.) inside CondaTainer environments on HPC systems.

## Directory Structure

```
helpers/
├── headless/ # Run directly on the current machine
├── slurm/    # Submit as SLURM batch jobs
├── pbs/      # Submit as PBS batch jobs
└── lsf/      # Submit as LSF batch jobs
```

## Quick Start

These scripts are designed to be executed by the **CondaTainer**. They should not be run directly from the host system.

```bash
condatainer helper -u # Update helper scripts (based on the current scheduler)

condatainer helper rstudio-server -r r4.4.3 # Launch RStudio Server with R 4.4.3
```

## Headless vs Scheduler Modes

If your system has a scheduler, you should not use the headless scripts.

- **headless/** -- Runs the service directly on the current server. Use on headless server which does not have a job scheduler. The script blocks while the service is running.
- **slurm/** -- Submits a SLURM batch job, waits for it to start, then opens an SSH tunnel with port forwarding. Re-running the script while a job is active will reconnect to the existing session.
- **pbs/** -- Same behavior as `slurm/` but submits via PBS (`qsub`/`qstat`/`qdel`).
- **lsf/** -- Same behavior as `slurm/` but submits via LSF (`bsub`/`bjobs`/`bkill`). Interactive helpers in `helpers/lsf/` force single-node placement using `#BSUB -R "span[hosts=1]"` so the service and port-forwarding run on one compute node.

## Desktop / GUI Apps

- **xfce4** -- Launches an XFCE desktop session via VNC / noVNC, accessible through a browser.
- **igv** -- Launches an XFCE desktop with IGV (Integrative Genomics Viewer) via noVNC.

## RStudio Server Variants

- **rstudio-server** -- Uses Posit's R image overlays (`r4.4.3`, `r4.5.2`, etc.). Specify R version with `-r`. If omitted, detects R in the overlay or falls back to the latest available R overlay.
- **rstudio-server-conda** -- Expects R installed inside the conda environment overlay (via `mm-install r-base`). No separate R overlay needed.

## Common Options

| Option | Description |
|--------|-------------|
| `-e <overlay>` | Environment overlay image file (default: `env.img`) |
| `-o <overlay>` | Additional overlay files (repeatable) |
| `-p <port>` | Port number (auto-selected if omitted) |
| `-w` | Use current directory as working directory |
| `-h` | Show help |
| `config` | Show config file path and contents |

Scheduler scripts also accept:

| Option | Description |
|--------|-------------|
| `-c <cpus>` | Number of CPUs (default: 4) |
| `-m <memory>` | Memory allocation (default: 16G or 32G) |
| `-t <time>` | Time limit (default: 12:00:00) |
| `-g <gpus>` | gpu specs (e.g. `type:count`) |
| `-v`        | View mode: less cpu/memory, suitable for light use |

## Configuration

Each script saves its defaults to `$XDG_CONFIG_HOME/condatainer/helper/<script-name>` (defaults to `~/.config/condatainer/helper/<script-name>`) on first run. Subsequent runs load from this file.

```bash
# View current config
./helpers/slurm/rstudio-server config

# Reset to defaults (delete config file, next run recreates it)
rm ~/.config/condatainer/helper/rstudio-server
```

Config files use simple `KEY="VALUE"` format and can be edited directly:

```bash
NCPUS="4"
MEM="32G"
TIME="12:00:00"
OVERLAY="env.img"
R_VERSION=""
PORT=""
BASE_IMAGE=""
OVERLAYS=""
```

## Port Forwarding

Scheduler scripts handle this automatically by SSH-ing to the compute node after the job starts.

## Logs

Scheduler job logs are written to `~/logs/`: 

```
~/logs/<service>-<jobid>.log
```

## Shared Library (.common.sh)

### Config helpers
| Function | Description |
|---|---|
| `config_init` | Create defaults file on first run with KEY=VALUE pairs. |
| `config_load` | Source saved defaults, set `CWD` to pwd, and auto-save `_CONFIG_*` copies for change detection. |
| `config_update` | Update specific keys in defaults without rewriting the file. |
| `config_require` | Validate required variables are set (exits on failure). |
| `config_show` | Print config file path and contents, then exit. |
| `config_path` | Return the config file path for a helper. |
| `state_path` | Return the state file path for a helper. |

### Port helpers
| Function | Description |
|---|---|
| `choose_port` | Pick a random available port in 20000-65000. |
| `validate_port` | Validate a numeric port is in 1024-65535 range. |
| `check_port_available` | Exit if the requested port is already in use. |

### Overlay & filesystem helpers
| Function | Description |
|---|---|
| `check_writable` | Probe that file is writable and not exclusively locked (`flock`). |
| `check_readable` | Probe that file is readable and not locked for writing (`flock`). |
| `require_writable` | Enforce writable availability; **headless** prompts to kill local PIDs, **scheduler** helpers suggest the cluster job query command (e.g. `squeue`/`qstat`/`bjobs`). |
| `check_overlay_in_use` | Convenience wrapper ensuring `.img` overlays are available for writing. |
| `check_overlay_integrity` | Runs `check_overlay_in_use` then `e2fsck -p` to validate/repair overlay (exit on unfixable errors). |
| `check_and_install_overlays` | Install missing named overlays using `condatainer create`. |
| `resolve_overlay_list` | Convert colon-separated overlay lists to absolute file paths. |
| `build_overlays_arg` | Build `-o` arguments from a colon-separated overlay list. |

### Prompt helpers
| Function | Description |
|---|---|
| `confirm_default_yes` | Prompt with `[Y/n]`. Returns 0 on yes (default). |
| `confirm_default_no` | Prompt with `[y/N]`. Returns 0 on yes. |

### Message helpers
| Function | Description |
|---|---|
| `print_msg` | Print `[MSG] message`. |
| `print_info` | Print `[INFO] message`. |
| `print_warn` | Print `[WARN] message`. |
| `print_error` | Print `[ERR] message` to stderr. |
| `print_pass` | Print `[PASS] message`. |

### Misc
| Function | Description |
|---|---|
| `check_condatainer` | Exit if `condatainer` is not in PATH. |

## Scheduler helpers (`helpers/*/.common.sh`)

### Display & reuse
| Function | Description |
|---|---|
| `spec_line` | Print one spec line with column alignment. Shows `(default)` or `(default: X)` by comparing `_ARG_*` and `_CONFIG_*` prefixes. |
| `print_specs` | Print all current job settings. Override in helper scripts for custom fields. |
| `countdown` | Print a countdown with Ctrl+C hint, overwriting the same line each tick. |
| `handle_reuse_mode` | Handle `REUSE_MODE` (always/ask/never) when a previous job is not running. Shows specs, prompts user, sets `REUSE_PREVIOUS_CWD`. |

### Job lifecycle
| Function | Description |
|---|---|
| `read_job_state` | Source state file and detect PENDING/RUNNING job state (returns codes 0-3). |
| `wait_for_job` | Poll `squeue` until job transitions to RUNNING; sets global `NODE` and `JOB_LOG`. |
