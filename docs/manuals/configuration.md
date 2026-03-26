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

**Scalar keys** (`apptainer_bin`, `default_distro`, `submit_job`, etc.): the highest-priority config file that sets the key wins.

**Directory and source array keys** (`extra_image_dirs`, `extra_build_dirs`, `extra_helper_dirs`, `extra_scripts_links`): **merged** across all config files. Entries from user config appear first (higher search priority), followed by extra-root, app-root, then system. This lets a sysadmin publish shared directories in an app-root or system config without requiring every user to copy them into their own config.

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
| Extra-root | `$CNT_EXTRA_ROOT/config.yaml` | Group/lab layer (requires `CNT_EXTRA_ROOT`) |
| App-root | `$CNT_ROOT/config.yaml` or `<install-dir>/config.yaml` | Shared cluster/group installation |
| System | `/etc/condatainer/config.yaml` | System-wide defaults |

To create a config file at a specific location:

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

## Configuration Options

### Directories

| Key | Default | Description |
|-----|---------|-------------|
| `logs_dir` | `$HOME/logs` | Directory for build job logs |
| `extra_image_dirs` | `[]` | Explicit image directories (direct paths). Entries support `:ro` (search-only) or `:rw` (writable, default) markers. |
| `extra_build_dirs` | `[]` | Explicit build-scripts directories (direct paths). |
| `extra_helper_dirs` | `[]` | Explicit helper-scripts target directories (direct paths). Entries support `:ro` (search-only) or `:rw` (writable, default) markers. |
| `extra_image_dirs` entries with `:ro` | — | Search-only image dirs (no writes) |

### Binaries

| Key | Default | Description |
|-----|---------|-------------|
| `apptainer_bin` | Auto-detected | Path to apptainer or singularity binary |
| `scheduler_bin` | Auto-detected | Path to job scheduler binary (sbatch, qsub, bsub, condor_submit, etc.). |

### Remote Sources

| Key | Default | Description |
|-----|---------|-------------|
| `scripts_link` | `https://raw.githubusercontent.com/Justype/cnt-scripts/main` | Base URL for remote build and helper scripts (lowest-priority remote) |
| `extra_scripts_links` | `[]` | Additional remote URLs prepended before `scripts_link` (first entry = highest priority); array, set via `config append/prepend` or `CNT_EXTRA_SCRIPTS_LINKS` |
| `prebuilt_link` | `https://github.com/Justype/cnt-scripts/releases/download` | Base URL for downloading prebuilt overlays |
| `prefer_remote` | `false` | Remote build scripts take precedence over local |

### Options

| Key | Default | Description |
|-----|---------|-------------|
| `submit_job` | `true` | Submit builds as scheduler jobs (disabled if no scheduler found) |
| `parse_module_load` | `false` | Parse `module load` / `ml` lines as dependencies in `check` and `run` |
| `scheduler_timeout` | `5` | Seconds to wait for a scheduler command (sbatch, qsub, etc.) before returning an error. Set to `0` to disable the timeout. |
| `metadata_cache_ttl` | `7` | Days to keep the cached remote build script metadata (default: 7 days = 1 week). Set to `0` to disable caching and always fetch from the network. |
| `default_distro` | `ubuntu24` | Base OS distro for the base image and bare-name expansion. Accepted values: `ubuntu20`, `ubuntu22`, `ubuntu24`. Determines the base image filename (e.g. `ubuntu24--base_image.sif`) and the distro prefix added to bare package names (e.g. `igv` → `ubuntu24/igv`). |

### Build Configuration

| Key | Default | Description |
|-----|---------|-------------|
| `build.ncpus` | `4` | CPUs for build jobs |
| `build.mem` | `8192` | Memory for build jobs (supports units: `8g`, `8192`) |
| `build.time` | `2h` | Time limit for builds |
| `build.compress_args` | Auto-detected | mksquashfs compression arguments (gzip for singularity; zstd-medium for apptainer≥1.4; lz4 otherwise) |
| `build.block_size` | `128k` | mksquashfs block size for app/env/external overlays (e.g. `128k`, `512k`) |
| `build.data_block_size` | `1m` | mksquashfs block size for data overlays (e.g. `512k`, `1m`) |
| `build.use_tmp_overlay` | `false` | Use a temporary ext3 overlay during builds instead of host directories |
| `build.tmp_overlay_size` | `20480` | Temporary ext3 overlay size (supports units: `20g`, `20480`); only used when `use_tmp_overlay` is `true` |
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

Array keys (`extra_image_dirs`, `extra_build_dirs`, `extra_helper_dirs`, `extra_scripts_links`, `channels`) use dedicated subcommands:

```bash
# Explicit image/scripts directories
condatainer config append extra_image_dirs /shared/lab/images:ro   # search-only
condatainer config append extra_image_dirs /fast/scratch/images    # writable
condatainer config append extra_build_dirs /shared/lab/scripts

# Prepend (higher priority — checked first)
condatainer config prepend extra_image_dirs /fast/images
condatainer config prepend extra_scripts_links https://raw.githubusercontent.com/MyOrg/my-scripts/main

# Remove
condatainer config remove extra_image_dirs /shared/lab/images:ro
```

Shell completion for `remove` offers the current values of the array as candidates.

### Validate Configuration

```bash
condatainer config validate
```

This command checks that key binaries are accessible, build settings are sane, and the `default_distro` value is one of the supported distro names.

## Environment Variables

Environment variables **always replace** the corresponding config file value — including merged array keys. This makes env vars suitable for module-file-based admin control: setting a variable in a module file gives a predictable, reproducible environment regardless of what users have in their configs.

Most configuration settings may be overridden by environment variables.
The name is derived automatically from the Viper key by
upper‑casing, replacing `.` with `_`, and prefixing with
`CNT_`.  For example:

* `logs_dir` → `CNT_LOGS_DIR`
* `scheduler.time` → `CNT_SCHEDULER_TIME`
* `build.tmp_overlay_size` → `CNT_BUILD_TMP_OVERLAY_SIZE`

You can list the supported variables with
`condatainer config show` (it prints any that are currently set).

`CNT_TMPDIR` is a special override used by the build system to select
temporary directories across build types. If set, it takes precedence over
scheduler-provided scratch, writable tmp auto-detection, and
`TMPDIR`/`TEMP`/`TMP`. CondaTainer will append `cnt-$USER` to the path to avoid
user collisions, including external `.sh`/`.bash` builds.

A few common overrides are shown below for clarity, but the
mapping is consistent for every key handled by the CLI:

| Environment Variable               | Config Key             |
|-----------------------------------|------------------------|
| `CNT_APPTAINER_BIN`        | `apptainer_bin`        |
| `CNT_SUBMIT_JOB`           | `submit_job`           |
| `CNT_SCRIPTS_LINK`         | `scripts_link`         |
| `CNT_PREBUILT_LINK`        | `prebuilt_link`        |
| `CNT_PREFER_REMOTE`        | `prefer_remote`        |
| `CNT_SCHEDULER_NODES`      | `scheduler.nodes`      |
| `CNT_BUILD_MEM`            | `build.mem`            |
| `CNT_BUILD_BLOCK_SIZE`     | `build.block_size`     |
| `CNT_BUILD_DATA_BLOCK_SIZE`| `build.data_block_size`|
| `CNT_ROOT`                 | Cluster/system root dir (loads `config.yaml` + data dirs; replaces bin/ heuristic) |
| `CNT_EXTRA_ROOT`           | Group/lab root dir — single path, loads `config.yaml` + data dirs |
| `CNT_EXTRA_IMAGE_DIRS`     | `extra_image_dirs` (pipe-separated; entries support `:ro`/`:rw`) |
| `CNT_EXTRA_BUILD_DIRS`     | `extra_build_dirs` (pipe or colon-separated; plain paths) |
| `CNT_EXTRA_HELPER_DIRS`    | `extra_helper_dirs` (pipe-separated; entries support `:ro`/`:rw`) |
| `CNT_EXTRA_SCRIPTS_LINKS`  | `extra_scripts_links` (pipe-separated) |
| `CNT_CHANNELS`             | `channels` (pipe or colon-separated) |
| `CNT_SCHEDULER_TIMEOUT`    | `scheduler_timeout`    |
| `CNT_METADATA_CACHE_TTL`   | `metadata_cache_ttl`   |
| `CNT_TMPDIR`               | (special override)     |

Example:

```bash
# Group/lab root (set in module file)
export CNT_ROOT=/cluster/condatainer      # cluster-level
export CNT_EXTRA_ROOT=/shared/lab/tools   # group-level

# Explicit image directories (pipe-separated; supports :ro/:rw markers)
export CNT_EXTRA_IMAGE_DIRS="/shared/lab/images:ro|/fast/scratch/images"
```

## Data Directory Search Paths

**CondaTainer** searches multiple directories for images, build scripts, and helper scripts. The search order determines which files are used when duplicates exist.

### Search Priority

**Images:**
1. `extra_image_dirs` — explicit image directories (`:ro` entries skipped for writes)
2. **Extra-root** → `$CNT_EXTRA_ROOT/images/` (group/lab layer)
3. **App-root** → `$CNT_ROOT/images/` or `<install>/images/`
4. **Scratch** → `$SCRATCH/condatainer/images/`
5. **User** → `~/.local/share/condatainer/images/`

**Build / Helper Scripts:**
**Build scripts:**
1. `extra_build_dirs` — explicit build-scripts directories
2. **Extra-root**, **App-root**, **Scratch**, **User** (same pattern)

**Helper scripts:**
1. `extra_helper_dirs` — explicit helper-scripts directories
2. **Extra-root**, **App-root**, **Scratch**, **User** (same pattern)

### View Search Paths

```bash
condatainer config paths
```

This shows all search paths for:
- **Images**: `.sif` and `.sqf` files
- **Build scripts**: Build recipe files
- **Helper scripts**: Runtime helper scripts

### Directory Structure

Each base directory follows this structure:

```
<base_dir>/
  images/           # Container images and overlays
  build-scripts/    # Build recipes
  helper-scripts/   # Runtime helpers
  cache/            # Cached remote metadata (remote-scripts-<hash>.json, helper-scripts-<hash>.json)
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

# Base URL for remote build scripts and helper scripts (includes branch)
# Change to use a private or institutional scripts repo
scripts_link: https://raw.githubusercontent.com/Justype/cnt-scripts/main

# Additional remote sources (higher priority than scripts_link; first entry wins on conflict)
# extra_scripts_links:
#   - https://raw.githubusercontent.com/MyOrg/my-scripts/main

# Base URL for downloading prebuilt images and overlays
prebuilt_link: https://github.com/Justype/cnt-scripts/releases/download

# Remote build scripts take precedence over local
prefer_remote: false

# Parse "module load" / "ml" lines as dependencies in 'check' and 'run' (default: false)
parse_module_load: false

# Maximum seconds to wait for scheduler CLI commands (default: 5, 0 = disabled)
scheduler_timeout: 5

# Days to cache remote build script metadata (default: 7 = 1 week, 0 = disabled)
metadata_cache_ttl: 7

# Base OS distro: ubuntu20, ubuntu22, or ubuntu24 (default: ubuntu24)
# Sets the base image (e.g. ubuntu24--base_image.sif) and prefix for bare package names
default_distro: ubuntu24

# Explicit image directories (direct paths; :ro = search-only, :rw = writable default)
extra_image_dirs:
  - /shared/lab/images:ro        # shared read-only store
  - /fast/scratch/images         # writable personal store

# Explicit build-scripts directories
# extra_build_dirs:
#   - /shared/lab/scripts

# Explicit helper-scripts directories
# extra_helper_dirs:
#   - /shared/lab/helpers

# Extra base directories (standard layout: images/, build-scripts/, etc.)
# For a group/lab root with standard layout, set in module file:
# export CNT_EXTRA_ROOT=/project/shared/condatainer

# Build configuration
build:
  ncpus: 4
  mem: 8g
  time: 2h
  compress_args: -comp zstd -Xcompression-level 8
  block_size: 128k       # SquashFS block size for app/env/external overlays
  data_block_size: 1m    # SquashFS block size for data overlays
  use_tmp_overlay: false  # Use ext3 tmp overlay instead of host directories
  tmp_overlay_size: 20g  # Only used when use_tmp_overlay is true

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

All config files (user, extra-root, app-root, system) are loaded simultaneously. For `extra_image_dirs` and other directory keys, entries from all configs are **merged** — so a group admin can add shared directories to the app-root config and every user automatically searches those directories, even if they also have a personal config.

For scalar keys like `apptainer_bin`, the user config takes priority; users can override app-root/system defaults in their own config without affecting other users.

**Explicit root via `CNT_ROOT`:** Instead of relying on the `bin/` layout detection, set `CNT_ROOT` to point directly to the installation directory. This is useful when the binary is installed to a standard location (e.g. `/usr/local/bin`) but the data lives elsewhere:

```bash
# In a module file:
export CNT_ROOT=/shared/cluster/condatainer
```

`CNT_ROOT` takes priority over the executable-location heuristic. The directory does not need a `bin/` subdirectory.

To set up an app-root config for a shared installation:

```bash
condatainer config init -l app-root
```

## Multi-Tier Setup (System → Group → User)

On HPC systems, configuration is typically layered across three scopes. CondaTainer supports this natively — all config files are loaded and merged simultaneously.

| Tier | Scope | Sets |
|---|---|---|
| System / cluster | Sysadmin | `apptainer_bin`, `scheduler_bin`, shared images, `channels` |
| Group / lab | Lab admin | Lab-specific images, build scripts, helper scripts |
| User | Individual | Personal overrides, personal scratch dirs |

Priority: **user > group > system > defaults**

---

**Filesystem layout:**
```
/cluster/condatainer/          ← system tier (CNT_ROOT or bin/ detection)
  bin/condatainer
  config.yaml                  ← apptainer_bin, scheduler_bin, channels
  images/                      ← cluster-wide base images

/shared/labA/condatainer/      ← group tier (CNT_EXTRA_ROOT)
  config.yaml                  ← extra_image_dirs, extra_build_dirs
  images/                      ← lab-specific images
  build-scripts/               ← lab-specific build recipes

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
extra_image_dirs:
  - /shared/labA/condatainer/images
extra_build_dirs:
  - /shared/labA/condatainer/build-scripts
scripts_link: https://raw.githubusercontent.com/LabA/cnt-scripts/main
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

---

## Compression Settings

CondaTainer auto-detects the best compression based on your runtime:

- **Singularity**: Uses gzip compression (`-comp gzip`)
- **Apptainer >= 1.4**: Uses zstd compression (`-comp zstd -Xcompression-level 8`)
- **Apptainer < 1.4**: Uses lz4 compression (`-comp lz4`)

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
| Data / Reference | `build.data_block_size` | `1m` | `--data-block-size` |

- Valid values: power of two between `4k` and `1m` (mksquashfs constraint). Common choices: `4k`, `8k`, `16k`, `32k`, `64k`, `128k`, `256k`, `512k`, `1m`.
- Smaller blocks (`128k`) give better random-read performance — ideal for executables loaded at runtime.
- Larger blocks (`1m`) give better compression ratios — ideal for large reference files that are read sequentially.

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

Run `condatainer config init` — it will automatically search environment modules via `module avail` if no binary is found in `PATH`. If detection still fails, verify the module name with `module avail apptainer` and set the path manually:

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

### Using multiple remote script sources

To add institutional or personal script repositories alongside the default:

```bash
# Add via CLI (takes priority over the default scripts_link)
condatainer config prepend extra_scripts_links https://raw.githubusercontent.com/MyOrg/my-scripts/main

# Remove when no longer needed
condatainer config remove extra_scripts_links https://raw.githubusercontent.com/MyOrg/my-scripts/main
```

Or via environment variable (pipe-separated):

```bash
export CNT_EXTRA_SCRIPTS_LINKS=https://raw.githubusercontent.com/MyOrg/my-scripts/main
```

Earlier entries take priority over `scripts_link`. Each remote gets its own cache file (`remote-scripts-<hash>.json`, `helper-scripts-<hash>.json`). Cache files for removed remotes are cleaned up automatically on `condatainer update`.

### Disable job submission

For systems without a scheduler or for local builds:

```bash
condatainer config set submit_job false
```
