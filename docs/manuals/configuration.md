# Configuration

**CondaTainer** uses a layered configuration system with multiple sources. Configuration values are applied in priority order, with higher-priority sources overriding lower ones.

- [Configuration Priority](#configuration-priority)
- [Data Directory Search Paths](#data-directory-search-paths)

## Configuration Priority

1. **Command-line flags** (highest priority)
2. **Environment variables** (`CNT_*`)
3. **User config file** (`~/.config/condatainer/config.yaml`)
4. **Portable config** (`<install-dir>/config.yaml`)
5. **System config file** (`/etc/condatainer/config.yaml`)
6. **Defaults** (lowest priority)

## Quick Start

Initialize a config file with auto-detected settings:

```bash
condatainer config init
```

View current configuration:

```bash
condatainer config show
```

````{note}
If you have Apptainer modules instead of system-wide binaries, load the module first:

```bash
ml apptainer
condatainer config init
```
````

## Config File Locations

| Type | Path | Use Case |
|------|------|----------|
| User | `~/.config/condatainer/config.yaml` | Personal settings |
| Portable | `<install-dir>/config.yaml` | Shared group installation |
| System | `/etc/condatainer/config.yaml` | System-wide defaults |

To create a config file at a specific location:

```bash
# User config (default for home installations)
condatainer config init -l user

# Portable config (for shared installations; not under home directory)
condatainer config init -l portable

# System config (requires appropriate permissions)
condatainer config init -l system
```

## Configuration Options

### Directories

| Key | Default | Description |
|-----|---------|-------------|
| `logs_dir` | `$HOME/logs` | Directory for build job logs |
| `extra_base_dirs` | `[]` | Additional directories to search for images and scripts |

### Binaries

| Key | Default | Description |
|-----|---------|-------------|
| `apptainer_bin` | Auto-detected | Path to apptainer or singularity binary |
| `scheduler_bin` | Auto-detected | Path to job scheduler binary (sbatch, qsub, bsub, condor_submit, etc.). |

### Remote Sources

| Key | Default | Description |
|-----|---------|-------------|
| `scripts_link` | `https://raw.githubusercontent.com/Justype/cnt-scripts/main` | Base URL for remote build and helper scripts |
| `prebuilt_link` | `https://github.com/Justype/cnt-scripts/releases/download` | Base URL for downloading prebuilt overlays |
| `prefer_remote` | `false` | Remote build scripts take precedence over local |

### Options

| Key | Default | Description |
|-----|---------|-------------|
| `submit_job` | `true` | Submit builds as scheduler jobs (disabled if no scheduler found) |
| `parse_module_load` | `false` | Parse `module load` / `ml` lines as dependencies in `check` and `run` |
| `scheduler_timeout` | `5` | Seconds to wait for a scheduler command (sbatch, qsub, etc.) before returning an error. Set to `0` to disable the timeout. |
| `default_distro` | `ubuntu24` | Base OS distro for the base image and bare-name expansion. Accepted values: `ubuntu20`, `ubuntu22`, `ubuntu24`. Determines the base image filename (e.g. `ubuntu24--base_image.sif`) and the distro prefix added to bare package names (e.g. `igv` → `ubuntu24/igv`). |

### Build Configuration

| Key | Default | Description |
|-----|---------|-------------|
| `build.ncpus` | `4` | CPUs for build jobs |
| `build.mem_mb` | `8192` | Memory in MB (8GB) |
| `build.time` | `2h` | Time limit for builds |
| `build.compress_args` | Auto-detected | mksquashfs compression arguments (gzip for singularity; zstd-medium for apptainer≥1.4; lz4 otherwise) |
| `build.overlay_type` | `ext3` | Overlay filesystem type: `ext3` or `squashfs` |
| `build.block_size` | `128k` | mksquashfs block size for app/env/external overlays (e.g. `128k`, `512k`) |
| `build.data_block_size` | `1m` | mksquashfs block size for data overlays (e.g. `512k`, `1m`) |
| `build.use_tmp_overlay` | `false` | Use a temporary ext3 overlay during builds instead of host directories |
| `build.tmp_overlay_size_mb` | `20480` | Temporary ext3 overlay size in MB (20GB); only used when `use_tmp_overlay` is `true` |

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
condatainer config set build.mem_mb 16384

# Set build time limit (supports Go-style or HPC-style formats)
condatainer config set build.time 4h
condatainer config set build.time 02:00:00

# Disable job submission (run builds locally)
condatainer config set submit_job false

# Set scheduler command timeout (seconds)
condatainer config set scheduler_timeout 10
```

### Edit Config File Directly

```bash
condatainer config edit
```

This opens the config file in your default editor (`$EDITOR`).

### Validate Configuration

```bash
condatainer config validate
```

This command checks that key binaries are accessible, build settings are sane, and the `default_distro` value is one of the supported distro names.

## Environment Variables

Most configuration settings may be overridden by environment variables.
The name is derived automatically from the Viper key by
upper‑casing, replacing `.` with `_`, and prefixing with
`CNT_`.  For example:

* `logs_dir` → `CNT_LOGS_DIR`
* `scheduler.time` → `CNT_SCHEDULER_TIME`
* `build.tmp_overlay_size_mb` → `CNT_BUILD_TMP_OVERLAY_SIZE_MB`

You can list the supported variables with
`condatainer config show` (it prints any that are currently set).

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
| `CNT_BUILD_MEM_MB`         | `build.mem_mb`         |
| `CNT_BUILD_BLOCK_SIZE`     | `build.block_size`     |
| `CNT_BUILD_DATA_BLOCK_SIZE`| `build.data_block_size`|
| `CNT_EXTRA_BASE_DIRS`      | `extra_base_dirs` (colon-separated) |
| `CNT_SCHEDULER_TIMEOUT`    | `scheduler_timeout`    |

Example:

```bash
# Add extra search directories (highest priority)
export CNT_EXTRA_BASE_DIRS=/shared/tools:/project/common
```

## Data Directory Search Paths

**CondaTainer** searches multiple directories for images, build scripts, and helper scripts. The search order determines which files are used when duplicates exist.

### Search Priority

1. **Extra base directories** (from `extra_base_dirs` config or `CNT_EXTRA_BASE_DIRS`)
2. **Portable** (group/shared directory next to the executable; if not under `$HOME`)
3. **Scratch** (`$SCRATCH/condatainer` on HPC systems)
4. **User** (`~/.local/share/condatainer`)

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

# Base URL for downloading prebuilt images and overlays
prebuilt_link: https://github.com/Justype/cnt-scripts/releases/download

# Remote build scripts take precedence over local
prefer_remote: false

# Parse "module load" / "ml" lines as dependencies in 'check' and 'run' (default: false)
parse_module_load: false

# Maximum seconds to wait for scheduler CLI commands (default: 5, 0 = disabled)
scheduler_timeout: 5

# Base OS distro: ubuntu20, ubuntu22, or ubuntu24 (default: ubuntu24)
# Sets the base image (e.g. ubuntu24--base_image.sif) and prefix for bare package names
default_distro: ubuntu24

# Additional search directories
extra_base_dirs:
  - /project/shared/condatainer
  - /apps/bioinformatics/condatainer

# Build configuration
build:
  ncpus: 4
  mem_mb: 8192
  time: 2h
  compress_args: -comp zstd -Xcompression-level 8
  overlay_type: ext3
  block_size: 128k       # SquashFS block size for app/env/external overlays
  data_block_size: 1m    # SquashFS block size for data overlays
  use_tmp_overlay: false  # Use ext3 tmp overlay instead of host directories
  tmp_overlay_size_mb: 20480  # Only used when use_tmp_overlay is true
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

## Portable Installations

For shared group installations, CondaTainer supports "portable" mode where the config and data live alongside the executable:

```
/project/group/condatainer/
  bin/
    condatainer         # Executable
  config.yaml           # Portable config
  images/               # Shared images
  build-scripts/        # Shared build scripts
```

Users can still have personal configs (`~/.config/condatainer/config.yaml`) that override the shared settings.

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

Load the apptainer (or singularity) module first, then reinitialize:

```bash
module load apptainer   # or: module load singularity
condatainer config init
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

### Disable job submission

For systems without a scheduler or for local builds:

```bash
condatainer config set submit_job false
```
