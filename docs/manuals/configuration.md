# Configuration

**CondaTainer** uses a layered configuration system with two independent mechanisms.

- [Configuration Priority](#configuration-priority)
- [Data Directory Search Paths](#data-directory-search-paths)

## Configuration Priority

### Environment Variables (admin override)

`CNT_*` environment variables always **replace** the config file value for that key. When set (e.g. via a module file), config files are not consulted for that key. This is the recommended way for sysadmins to enforce cluster-wide settings without editing user configs.

### Config Files (layered)

All three config files are loaded and merged when they exist:

1. **Command-line flags** (highest priority)
2. **Environment variables** (`CNT_*`) — replaces; config files not consulted
3. **User config file** (`~/.config/condatainer/config.yaml`)
4. **Extra-root config** (`$CNT_EXTRA_ROOT/config.yaml`, group/lab layer)
5. **App-root config** (`$CNT_ROOT/config.yaml` or `<install-dir>/config.yaml`)
6. **System config file** (`/etc/condatainer/config.yaml`)
7. **Defaults** (lowest priority)

**Scalar keys** (`apptainer_bin`, `base`, `submit_job`, etc.): the highest-priority config file that sets the key wins.

**`sources`**: **merged** across all config files. Entries from user config appear first (higher search priority), followed by extra-root, app-root, then system. This lets a sysadmin publish shared recipe collections in an app-root or system config without requiring every user to copy them into their own config. See [How Array Settings Merge](#how-array-settings-merge) for a worked example.

**`channels`**: overwrite — the highest-priority config file that sets it wins (not merged), since channel order controls conda package resolution priority.

## Quick Start

Initialize a config file with auto-detected settings:

```bash
condatainer config init
```

View current configuration:

```bash
condatainer config show
```

## Config File Locations

| Type | Path | Use Case |
|------|------|----------|
| User | `~/.config/condatainer/config.yaml` | Personal settings |
| Extra-root | `$CNT_EXTRA_ROOT/config.yaml` | Group/lab layer (requires `$CNT_EXTRA_ROOT`) |
| App-root | `$CNT_ROOT/config.yaml` or `<install-dir>/config.yaml` | Shared cluster/group installation |
| System | `/etc/condatainer/config.yaml` | System-wide defaults |

To create a config file at a specific layer:

```bash
# User config (default for home installations)
condatainer config init -l user

# App-root config (for shared installations; not under home directory)
condatainer config init -l app-root

# Extra-root config (group/lab layer)
CNT_EXTRA_ROOT=/shared/labA/condatainer condatainer config init -l extra-root

# System config (requires appropriate permissions)
condatainer config init -l system
```

### Where Writes Go

Without `-l`, `config set` / `append` / `prepend` / `remove` write to the highest-priority config file that **exists** — your user config if you have one, otherwise the next layer down.

If that file is read-only, the write goes to your user config instead (created if needed), with a warning:

```
[CNT!] app-root config is read-only; saving to your user config instead (applies only to you).
```

With `-l`, a read-only target is an error instead — an explicit layer is never silently retargeted.

## Configuration Options

### Directories

| Key | Default | Description |
|-----|---------|-------------|
| `logs_dir` | `$HOME/logs` | Directory for build job logs |

### Binaries

| Key | Default | Description |
|-----|---------|-------------|
| `apptainer_bin` | Auto-detected | Path to apptainer or singularity binary |
| `scheduler_bin` | Auto-detected | Path to job scheduler binary (sbatch, qsub, bsub, condor_submit, etc.). |

### Recipe Sources

| Key | Default | Description |
|-----|---------|-------------|
| `sources` | the public `cnt` collection | Ordered recipe collections, first match wins |

Each entry is a single-key mapping of a handle to a local directory or a base
URL. Order is priority — like `PATH`, the first collection holding a name wins:

```yaml
sources:
  - lab: /shared/labA/recipes
  - cnt: https://raw.githubusercontent.com/condatainer/recipes/main
```

The public `cnt` collection is appended automatically when nothing else claims
that handle, so recipes resolve on a fresh install. Give an entry the handle
`cnt` to replace it outright.

A collection is a directory (or URL prefix) with `recipes/` and a `source.json`
declaring `schema`, `repository` and `default_base` — the base recipe used when
`base` is unset.

### Options

| Key | Default | Description |
|-----|---------|-------------|
| `submit_job` | `true` | Submit builds as scheduler jobs (disabled if no scheduler found) |
| `parse_module_load` | `false` | Parse `module load` / `ml` lines as dependencies in `check` and `run` |
| `autoload_gpu` | `true` | Pass `--nv` / `--rocm` when the host has the device node. Set `false` if the driver is present but unusable |
| `base` | first source's `default_base` | Base recipe for the container root, e.g. `ubuntu24` → `ubuntu24/base` |
| `scheduler_timeout` | `0` | Seconds to wait for a scheduler command before erroring. `0` disables the timeout |
| `notification` | `web` | Alert when a helper job starts: `web`, `terminal`, `both`, `none` |
| `metadata_cache_ttl` | `7` | Days to cache remote recipe metadata. `0` always fetches |
| `proxy_perjob` | `false` | Auto-start a per-job SOCKS5 proxy inside submitted jobs. See [Proxy](condatainer.md#proxy). |
| `helper_bind_all` | `false` | Bind helper services to `0.0.0.0` for direct TCP instead of an SSH tunnel |

### Build Configuration

| Key | Default | Description |
|-----|---------|-------------|
| `build.ncpus` | `4` | CPUs for build jobs |
| `build.mem` | `8192` | Memory for build jobs (supports units: `8g`, `8192`) |
| `build.time` | `2h` | Time limit for builds |
| `build.compress_args` | Auto-detected | mksquashfs compression arguments (zstd-medium for apptainer≥1.4; lz4 otherwise, including Singularity) |
| `build.block_size` | `128k` | mksquashfs block size for app/env/external overlays (e.g. `128k`, `512k`) |
| `build.data_block_size` | `512k` | mksquashfs block size for data overlays (e.g. `512k`, `1m`) |
| `build.app_tmp_overlay` | `false` | Assemble an **app** build inside a temporary ext3 overlay instead of host directories. Ignored for `data`, `os` and `base` |
| `build.always_submit` | `false` | Always submit builds as scheduler jobs even if the script has no scheduler directives |
| `build.app_tmp_overlay_size` | `20480` | Size of that overlay (supports units: `20g`, `20480`); only used when `app_tmp_overlay` is `true` |
| `channels` | `[conda-forge, bioconda]` | Conda channels passed to micromamba in priority order (first = highest priority) |

> `build.compress_args` also accepts shortcuts: `gzip`, `lz4`, `zstd`, `zstd-fast`, `zstd-medium`, `zstd-high`
>
> `build.block_size` and `build.data_block_size` must be a power of two between `4k` and `1m` (mksquashfs `-b` limit). Larger blocks improve compression ratio but increase random-read latency.

## Managing Configuration

### View Configuration

```bash
# Show all settings
condatainer config show

# Show only config file path
condatainer config show --path

# Get a specific value
condatainer config get apptainer_bin
condatainer config get build.ncpus
```

### Set Configuration Values

Shell completion is available for config keys - press Tab to see available options.

```bash
# Set apptainer binary path
condatainer config set apptainer_bin /usr/bin/apptainer

# Set build resources
condatainer config set build.ncpus 8
condatainer config set build.mem 16g

# Set build time limit (supports Go-style or HPC-style formats)
condatainer config set build.time 4h
condatainer config set build.time 02:00:00

# Disable job submission (run builds locally)
condatainer config set submit_job false

# Set scheduler command timeout (seconds)
condatainer config set scheduler_timeout 10
```

### Manage Array Config Values

Array keys (`sources`, `channels`) use dedicated subcommands:

```bash
# Recipe collections are written as handle=location
condatainer config append sources lab=/shared/lab/recipes

# Prepend (higher priority — checked first)
condatainer config prepend sources labA=/shared/labA/recipes

# Remove
condatainer config remove sources labA=/shared/labA/recipes
```

Shell completion for `remove` offers the current values of the array as candidates.

### Validate Configuration

```bash
condatainer config validate
```

This command checks that key binaries are accessible and build settings are sane.

## Environment Variables

Environment variables **always replace** the corresponding config file value — including merged array keys. This makes env vars suitable for module-file-based admin control: setting a variable in a module file gives a predictable, reproducible environment regardless of what users have in their configs.

Most configuration settings may be overridden by environment variables.
The name is derived automatically from the Viper key by
upper‑casing, replacing `.` with `_`, and prefixing with
`CNT_`.  For example:

* `logs_dir` → `CNT_LOGS_DIR`
* `build.mem` → `CNT_BUILD_MEM`
* `build.app_tmp_overlay_size` → `CNT_BUILD_APP_TMP_OVERLAY_SIZE`

You can list the supported variables with
`condatainer config show` (it prints any that are currently set).

A build works under one of two roots:

| root | where | used by |
|---|---|---|
| **fast** | `$CNT_TMPDIR` → scheduler scratch (`SLURM_TMPDIR`, `PBS_TMPDIR`, `LSF_TMPDIR`, `_CONDOR_SCRATCH_DIR`) → `$TMPDIR` → `/tmp`, plus `cnt-$USER` | `app` builds |
| **stable** | first writable `<data-dir>/tmp` (extra-root → root → scratch → user) | `data` builds, definitions and the base image |

`$CNT_TMPDIR` selects the **fast** root only. It does not redirect the stable
one: collapsing the two would put a large data payload, or a multi-GB `.sif`, on
node-local scratch that the job wipes when it ends. To move the stable root, move
the data directory (`CNT_ROOT` / `CNT_EXTRA_ROOT`).

An external build (`-f`) is the exception: `app` takes the fast root, while
`data` and `.def` builds keep their intermediates beside the target prefix, whose
location you chose.

A few common overrides are shown below for clarity, but the
mapping is consistent for every key handled by the CLI:

| Environment Variable               | Config Key             |
|-----------------------------------|------------------------|
| `CNT_APPTAINER_BIN`        | `apptainer_bin`        |
| `CNT_SUBMIT_JOB`           | `submit_job`           |
| `CNT_AUTOLOAD_GPU`         | `autoload_gpu`         |
| `CNT_BASE`                 | `base`                 |
| `CNT_BUILD_MEM`            | `build.mem`            |
| `CNT_BUILD_ALWAYS_SUBMIT`  | `build.always_submit`  |
| `CNT_BUILD_BLOCK_SIZE`     | `build.block_size`     |
| `CNT_BUILD_DATA_BLOCK_SIZE`| `build.data_block_size`|
| `CNT_ROOT`                 | Cluster/system root dir (loads `config.yaml` + data dirs; replaces bin/ heuristic) |
| `CNT_EXTRA_ROOT`           | Group/lab root dir — single path, loads `config.yaml` + data dirs |
| `CNT_SOURCES`              | `sources` (pipe-separated `handle=location`; replaces the list) |
| `CNT_CHANNELS`             | `channels` (pipe or colon-separated) |
| `CNT_SCHEDULER_TIMEOUT`    | `scheduler_timeout`    |
| `CNT_NOTIFICATION`         | `notification`         |
| `CNT_METADATA_CACHE_TTL`   | `metadata_cache_ttl`   |
| `CNT_PROXY_PERJOB`         | `proxy_perjob`         |
| `CNT_HELPER_BIND_ALL`      | `helper_bind_all`      |
| `CNT_TMPDIR`               | (fast build root; no config key) |

Example:

```bash
# Group/lab root (set in module file)
export CNT_ROOT=/cluster/condatainer      # cluster-level
export CNT_EXTRA_ROOT=/shared/lab/tools   # group-level
```

## Data Directory Search Paths

**CondaTainer** searches multiple directories for images, build scripts, and helper scripts. The search order determines which files are used when duplicates exist. See [Data Layers](../deployment/data_layers.md) for how this interacts with shared installations.

### Search Priority

**Images:**

| # | Directory | Layer |
|---|---|---|
| 1 | `$CNT_EXTRA_ROOT/images/` (group/lab) | `extra-root` |
| 2 | `$CNT_ROOT/images/` or `<install>/images/` | `app-root` |
| 3 | `$SCRATCH/condatainer/images/` | `user` |
| 4 | `~/.local/share/condatainer/images/` | `user` |

Scratch and the XDG data directory are **two directories in one `user` layer** — scratch is preferred when `$SCRATCH` is set, and `-l u` selects both.

Recipes have no directory layer — they come from the `sources` list, which is ordered on its own.

**Helper scripts:** the same layers with `helper-scripts/`.

These layer names are what `condatainer list` tags each directory with, and what `condatainer remove -l` accepts (`u`, `r`, `e`).

### View Search Paths

```bash
condatainer config paths
```

This shows all search paths for:
- **Images**: `.sif` and `.sqf` files
- **Build scripts**: Build recipe files
- **Helper scripts**: Runtime helper scripts

Each path is tagged with its layer — `(extra)`, `(extra-root)`, `(app-root)`, `(user)` — plus whether it exists, is writable, and is the write target.

### Directory Structure

Each base directory follows this structure:

```
<base_dir>/
  images/           # Container images and overlays
  build-scripts/    # Build recipes
  helper-scripts/   # Runtime helpers
  tmp/              # Temporary files during builds
```

## Example Config File

```yaml
# ~/.config/condatainer/config.yaml

# Log directory for build jobs
logs_dir: /home/user/logs

# Binary paths (scheduler type is auto-detected from binary)
apptainer_bin: /usr/bin/apptainer
scheduler_bin: /usr/bin/sbatch

# Submit builds as scheduler jobs
submit_job: true

# Recipe collections, in priority order (first match wins).
# The public cnt collection is appended automatically unless redefined here.
sources:
  - lab: /shared/labA/recipes
  - cnt: https://raw.githubusercontent.com/condatainer/recipes/main

# Parse "module load" / "ml" lines as dependencies in 'check' and 'run' (default: false)
parse_module_load: false

# Pass --nv / --rocm when the host has the matching device node (default: true)
# Set false on a node whose driver is installed but unusable
autoload_gpu: true

# Base recipe for the container root (default: the first source's default_base)
base: ubuntu24

# Maximum seconds to wait for scheduler CLI commands (default: 0 = disabled)
scheduler_timeout: 0

# Days to cache remote recipe metadata (default: 7 = 1 week, 0 = disabled)
metadata_cache_ttl: 7

# Notification when a helper job starts running (default: web)
# Values: web | terminal | both | none (or empty)
# notification: web

# Extra base directories (standard layout: images/, helper-scripts/)
# For a group/lab root with standard layout, set in module file:
# export CNT_EXTRA_ROOT=/project/shared/condatainer

# Build configuration
build:
  ncpus: 4
  mem: 8g
  time: 2h
  compress_args: -comp zstd -Xcompression-level 8
  block_size: 128k       # SquashFS block size for app/env/external overlays
  data_block_size: 512k  # SquashFS block size for data overlays
  app_tmp_overlay: false   # Assemble an app build inside an ext3 overlay (app only)
  always_submit: false    # Always submit as scheduler jobs even without directives
  app_tmp_overlay_size: 20g  # Only used when app_tmp_overlay is true

# proxy_perjob: true   # auto-start per-job proxy inside submitted jobs

# Default: conda-forge then bioconda
channels:
  - conda-forge
  - bioconda
```

## Time Duration Format

The `build.time` setting accepts two formats:

**Go-style durations:**
- `2h` - 2 hours
- `30m` - 30 minutes
- `1h30m` - 1 hour 30 minutes
- `90s` - 90 seconds

**HPC-style durations:**
- `02:00:00` - 2 hours (HH:MM:SS)
- `2:30:00` - 2 hours 30 minutes
- `1:30` - 1 hour 30 minutes (HH:MM)

## Standalone Installations

For shared group installations, CondaTainer supports a standalone layout where the config and data live alongside the executable:

```
/project/group/condatainer/
  bin/
    condatainer         # Executable
  config.yaml           # Root config (loaded alongside user, extra-root, and system configs)
  images/               # Shared images
  build-scripts/        # Shared build scripts
```

All config files (user, extra-root, app-root, system) are loaded simultaneously. For `sources`, entries from all configs are **merged** — so a group admin can publish shared recipe collections in the app-root config and every user automatically searches them, even if they also have a personal config.

For scalar keys like `apptainer_bin`, the user config takes priority; users can override app-root/system defaults in their own config without affecting other users.

**Explicit root via `$CNT_ROOT`:** Instead of relying on the `bin/` layout detection, set `$CNT_ROOT` to point directly to the installation directory. This is useful when the binary is installed to a standard location (e.g. `/usr/local/bin`) but the data lives elsewhere:

```bash
# In a module file:
export CNT_ROOT=/shared/cluster/condatainer
```

`$CNT_ROOT` takes priority over the executable-location heuristic. The directory does not need a `bin/` subdirectory.

To set up an app-root config for a shared installation:

```bash
condatainer config init -l app-root
```

## Multi-Tier Setup (System → Group → User)

On HPC systems, configuration is typically layered across three scopes. CondaTainer supports this natively — all config files are loaded and merged simultaneously. For the deployment guides that set these up, see [Data Layers](../deployment/data_layers.md).

| Tier | Data layer | Scope | Sets |
|---|---|---|---|
| System / cluster | `app-root` | Sysadmin | `apptainer_bin`, `scheduler_bin`, shared images, `channels` |
| Group / lab | `extra-root` | Lab admin | Lab-specific images, build scripts, helper scripts |
| User | `user` | Individual | Personal overrides, personal scratch dirs |

Priority: **user > group > system > defaults**

---

**Filesystem layout:**
```
/cluster/condatainer/          ← system tier (CNT_ROOT or bin/ detection)
  bin/condatainer
  config.yaml                  ← apptainer_bin, scheduler_bin, channels
  images/                      ← cluster-wide base images

/shared/labA/condatainer/      ← group tier (CNT_EXTRA_ROOT)
  config.yaml                  ← sources
  images/                      ← lab-specific images
  recipes/                     ← lab-specific recipes (as a `sources` entry)

~/.config/condatainer/config.yaml   ← user tier (auto-loaded)
```

**System config** (`/cluster/condatainer/config.yaml`):
```yaml
apptainer_bin: /usr/local/bin/apptainer
scheduler_bin: /usr/bin/sbatch
channels:
  - conda-forge
  - bioconda
```

**Group config** (`/shared/labA/condatainer/config.yaml`):
```yaml
sources:
  - labA: /shared/labA/condatainer/recipes
```

**User config** (`~/.config/condatainer/config.yaml`):
```yaml
logs_dir: /scratch/myuser/logs
build:
  ncpus: 8
  mem: 16g
```

Initialize the group config:
```bash
CNT_EXTRA_ROOT=/shared/labA/condatainer condatainer config init -l extra-root
```

> **Environment variables are always the highest priority.** Any `CNT_*` variable set in the shell overrides the corresponding key from all config files. This is useful for one-off overrides or admin control via module files — e.g. `export CNT_BUILD_NCPUS=16` overrides `build.ncpus` from every config layer.

### How Array Settings Merge

Scalar keys (like `build.ncpus`) are *overridden* — the highest tier that sets one wins. `sources` instead **merges: its entries are concatenated across every tier, user entries first**. A user adds to the merged list; they never replace what an admin published. (`channels` is the exception — an array that *overwrites*, since channel order decides package resolution.)

Say each tier contributes one recipe collection:

```yaml
# System config  (/cluster/condatainer/config.yaml)
sources:
  - cluster: /cluster/condatainer/recipes

# Group config   (/shared/labA/condatainer/config.yaml)
sources:
  - labA: /shared/labA/condatainer/recipes
```

The user prepends their own, so it is consulted first:

```bash
condatainer config prepend sources mine=/scratch/myuser/recipes
```
```yaml
# User config    (~/.config/condatainer/config.yaml)
sources:
  - mine: /scratch/myuser/recipes
```

The effective `sources` every lookup sees is all three, in **user → group → system** order:

```
1. mine    → /scratch/myuser/recipes        (user)
2. labA    → /shared/labA/condatainer/recipes (group)
3. cluster → /cluster/condatainer/recipes   (system)
```

First match wins, like `PATH`. The user was able to *prepend* their collection but cannot remove the group's or system's, so the shared collections stay in every user's resolution order. Confirm it any time with `condatainer config paths`.

Data directories do not merge this way — they come from the fixed tier list (`CNT_EXTRA_ROOT` → app-root → scratch → user), and `condatainer create` writes to the first one it can write.

---

## Compression Settings

CondaTainer auto-detects the best compression based on your runtime:

- **Apptainer >= 1.4**: Uses zstd compression (`-comp zstd -Xcompression-level 8`)
- **Apptainer < 1.4**: Uses lz4 compression (`-comp lz4`)
- **Singularity**: Uses lz4 compression (`-comp lz4`)

lz4 is the floor. Singularity and Apptainer below 1.4 cannot mount a
zstd-compressed SquashFS, so only Apptainer 1.4+ is moved up — an image the
runtime cannot read is worse than one that compresses less.

`gzip` is still available if you ask for it; nothing selects it automatically.

To override:

```bash
# explicit mksquashfs options
condatainer config set build.compress_args "-comp gzip -Xcompression-level 9"

# shorthand names (completion will offer these)
condatainer config set build.compress_args gzip
condatainer config set build.compress_args zstd-fast
```

## Block Size Settings

The SquashFS block size controls how data is chunked during compression. Two separate defaults are used based on overlay type:

| Overlay type | Config key | Default | CLI flag |
|---|---|---|---|
| App / Env / External | `build.block_size` | `128k` | `--block-size` |
| Data / Reference | `build.data_block_size` | `512k` | `--data-block-size` |

- Valid values: power of two between `4k` and `1m` (mksquashfs constraint). Common choices: `4k`, `8k`, `16k`, `32k`, `64k`, `128k`, `256k`, `512k`, `1m`.
- Smaller blocks (`128k`) give better random-read performance.
- Larger blocks (`1m`) give better compression ratios and sequential read.

To override per-build via CLI:

```bash
condatainer create samtools/1.22 --block-size 256k
condatainer create grch38/gtf-gencode/49 --data-block-size 1m
```

To set persistent defaults:

```bash
condatainer config set build.block_size 256k
condatainer config set build.data_block_size 1m
```

## Troubleshooting

### Config file not found

Run `condatainer config init` to create a config file with auto-detected settings.

### Apptainer/Singularity not found

Run `condatainer config init` — it will automatically search environment modules via `module avail` if no binary is found in `$PATH`. If detection still fails, verify the module name with `module avail apptainer` and set the path manually:

```bash
condatainer config set apptainer_bin /path/to/apptainer
```

### Scheduler not detected

If your HPC uses a non-standard scheduler path, set the binary and the type will be auto-detected:

```bash
condatainer config set scheduler_bin /custom/path/sbatch  # or qsub, bsub, condor_submit
```

### Scheduler commands timing out

If CondaTainer reports a scheduler timeout error, the scheduler daemon may be slow to respond. Increase the timeout or disable it entirely:

```bash
condatainer config set scheduler_timeout 30  # increase to 30 seconds
condatainer config set scheduler_timeout 0   # disable timeout entirely
```

### Remote metadata unavailable or stale

CondaTainer caches remote build script and helper script metadata for 1 week by default. If the cache is expired and the network is unavailable, the stale cache is used with a warning. To force a refresh:

```bash
condatainer update           # re-fetch build + helper metadata (default)
condatainer update --build   # build metadata only
condatainer update --helper  # helper metadata only
```

To disable caching entirely (always fetch live):

```bash
condatainer config set metadata_cache_ttl 0
```

To extend the cache lifetime (e.g. 2 weeks):

```bash
condatainer config set metadata_cache_ttl 14
```

### Using multiple recipe collections

To consult an institutional or personal collection before the public one:

```bash
# handle=location; prepend puts it ahead of everything already configured
condatainer config prepend sources myorg=https://raw.githubusercontent.com/MyOrg/recipes/main

# Remove when no longer needed
condatainer config remove sources myorg=https://raw.githubusercontent.com/MyOrg/recipes/main
```

Or via environment variable (pipe-separated), which replaces the whole list:

```bash
export CNT_SOURCES="myorg=https://raw.githubusercontent.com/MyOrg/recipes/main|cnt=https://raw.githubusercontent.com/condatainer/recipes/main"
```

Earlier entries win. Each remote collection caches its index separately, and
`condatainer update` refreshes them.

### Disable job submission

For systems without a scheduler or for local builds:

```bash
condatainer config set submit_job false
```
