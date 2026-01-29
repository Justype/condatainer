# CondaTainer Manual

📦 **CondaTainer** is a CLI tool designed to streamline the management of Apptainer (Singularity) containers backed by Conda environments and SquashFS overlays.

## Table of Contents

- [Naming Convention](#naming-convention)
- [Mount Points](#mount-points)
- [Create](#create)
- [Overlay](#overlay)
- [Container Management (Avail, List, Remove)](#container-management-avail-list-remove)
- [Exec](#exec)
- [Runtime (Check, Run)](#runtime-check-run)
- [Info](#info)
- [Helper](#helper)
- [Config](#config)
- [Scheduler](#scheduler)
- [Update](#update)
- [Completion](#completion)

## Overall Command Structure

```
CondaTainer: Use Apptainer/Conda/Overlays/SquashFS to manage tools/data/env for HPC users.

Usage:
  condatainer [command]

Available Commands:
  avail       Check available build scripts (local and remote)
  check       Check if the dependencies of a script are installed
  completion  Generate shell completion script
  config      Manage condatainer configuration
  create      Create a new SquashFS overlay
  exec        Execute a command using overlays
  helper      Manage and run helper scripts
  info        Show information about a specific overlay
  list        List installed overlays matching search terms
  o           Shortcut for 'overlay create'
  overlay     Manage persistent overlay images (create, resize, check, info)
  remove      Remove installed overlays matching search terms
  run         Run a script and auto-solve the dependencies by #DEP tags
  scheduler   Display scheduler information
  self-update Update condatainer to the latest version from GitHub

Flags:
      --debug     Enable debug mode with verbose output
  -h, --help      help for condatainer
      --local     Disable job submission (run locally)
  -q, --quiet     Suppress verbose overlay creation messages (errors/warnings still shown)
  -v, --version   version for condatainer

Use "condatainer [command] --help" for more information about a command.
```

## Naming Convention

**CondaTainer** classifies overlays into three distinct categories based on their naming structure.

### System Applications (Sys)

Used for system-level applications that do not follow the standard application or reference naming conventions.

Like text editors, IDEs, build-essential tools, etc.

* **Format:** `name`
* **Constraints:**
  * Must **not** contain slashes (`/`) or double-dash sequences (`--`).
* **Example:** `rstudio-server`, `texlive`, `code-server`

You can also create these system apps using the `-n, --name` flag with the `create` command.

### Custom Environments (Env)

Used when creating a custom environment with a user-defined name (`create -p`, or `overlay` action).

* **Format:** `custom_name`
* **Constraints:**
  * Must **not** contain the double-dash sequence (`--`).
  * Should **not** conflict with common application names (to avoid ambiguity with app overlays).
* **Example:** `my_analysis_env`, `project_x_utils`

> [!NOTE]
> Use `-n` to install system apps only (e.g., `nvim`). To create an environment `.sqf`, use `-p` to specify a prefix/path.

### Application Overlays (App)

Used for standard software packages and tools managed by **CondaTainer** build scripts.

* **Format:** `name/version`
* **Structure:**
  * **name**: The software package or tool name (e.g., `bcftools`).
  * **version**: The specific version of the software (e.g., `1.22`).
* **Example:** `cellranger/9.0.1`

### Reference Overlays (Ref)

Used for reference datasets, genome assemblies, or indices.

* **Format:** `assembly/data-type/version`
* **Structure:**
  * **assembly**: The genome assembly or project (e.g., `grch38`).
  * **data-type**: The type of data (e.g., `gtf-gencode`).
  * **version**: The release or build version (e.g., `47`).
* **Example:** `grch38/gtf-gencode/47`

> [!NOTE]
> `--` is reserved for separating components in overlays. For example, `bcftools/1.22` overlay will be located at `images/bcftools--1.22.sqf`.

### Rename or not to Rename?

TL;DR: Do not rename `sqf` files. `img` files can be renamed.

**CondaTainer** relies on the overlay names to determine mount points and environment variable settings. Renaming `.sqf` files can lead to inconsistencies and runtime errors.

e.g.

- `bcftools/1.22` is mounted at `/cnt/bcftools/1.22`
- `grch38/gtf-gencode/47` is mounted at `/cnt/grch38/gtf-gencode/47`
- `my_analysis_env.sqf` is mounted at `/cnt/my_analysis_env`

For `.img` files, the mount point is always `/ext3/env`, so renaming is allowed.

## Mount Points

**CondaTainer** uses clear terminology for overlay types and purposes. Read-only overlays use the `.sqf` extension (SquashFS) and writable overlays use the `.img` extension (ext3).

**Overlay Terminology**

| Term | Extension | R/W | Content Path | Purpose |
|------|-----------|-----|--------------|---------|
| OS | `.sqf` | R/O | `/bin`, `/lib`, `/usr` | System Foundation — minimal system root that can run standalone or serve as a base for Modules/Bundles. |
| Module | `.sqf` | R/O | `/cnt/<name>/<version>` | Individual tool overlay (single package) mounted under `/cnt` and layered on top of an OS or Bundle. |
| Bundle | `.sqf` | R/O | `/cnt/<env_name>` | Frozen Conda environment (prebuilt collection of packages) mounted under `/cnt` as a named environment. |
| Workspace | `.img` | R/W | `/ext3/env` | Writable Conda environment (ext3 image) for interactive work and runtime package changes. |

Read-only `.sqf` overlays are ideal for distributing immutable software and reference data. Writable `.img` overlays are for live development or when packages must be changed at runtime.

Mount points and examples

- Read-only overlays (`.sqf`) are mounted under `/cnt` following their naming convention.
  - Module example: `bcftools/1.22` → mounted at `/cnt/bcftools/1.22`.
  - Bundle example: `my_project_env.sqf` → mounted at `/cnt/my_project_env`.
  - OS overlays present system directories (e.g., binaries under `/bin`) inside the container filesystem exposed by the overlay.
- Writable images (`.img`) mount at `/ext3/env` and expose a full ext3 filesystem where packages may be installed or modified.

Writability: `.sqf` overlays are read-only. To enable write access for a `.img` overlay use the `-w` / `--writable-img` flag with `exec` or `run`.

## Overlay

Manage overlay writable ext3 images:

- **create**: Create an image with a conda environment inside.
- **info**: Display disk usage and filesystem statistics.
- **check**: Verify filesystem integrity.
- **chown**: Change ownership of files inside the image.
- **resize**: Resize the image file.

### Overlay Create

Create an ext3 `.img` with a conda environment inside. This is useful for creating writable conda environments.

**Usage:**

```
condatainer overlay create [OPTIONS] [NAME]
```

**Options:**

* `-s`, `--size [SIZE]`: Size of the overlay image (default: 10G). Supports units like `20G`, `2048M`.
* `-t`, `--type [TYPE]`: Overlay profile: `small`, `balanced`, or `large` (default: balanced).
* `--fs [FS]`: Filesystem type: `ext3` or `ext4` (default: ext3).
* `-f`, `--file [FILE]`: Initialize with a Conda environment file (.yml or .yaml).
* `--fakeroot`: Create image compatible with fakeroot (owned by root, must use with `--fakeroot` later).
* `--sparse`: Create a sparse image file.
* NAME: Name of the overlay image (`env.img` by default if not specified).

**Examples:**

```bash
# Create a 20GB overlay img with a custom name
condatainer overlay create -s 20G my_analysis_env

# Create with a specific profile (small/balanced/large)
condatainer overlay create -t large data_env.img

# o is a shortcut for 'overlay create'
condatainer o -f environment.yml project_env.img

# Create a fakeroot-compatible sparse overlay
condatainer o --fakeroot --sparse
```

Then you can use the `mm-<operation>` helper commands to create/install/update/remove conda packages inside the writable container.

```bash
# install more packages
mm-install numpy pandas

# pin the package version
mm-pin numpy # after installation
mm-pin -d numpy # to unpin
mm-pin -l # list pinned packages

# update packages
mm-update

# remove packages
mm-remove pandas

# list installed packages
mm-list

# clean all conda cache
mm-clean -ay

# export the env
mm-export --no-builds > my_env.yaml
```

### Overlay Info

Display disk usage and filesystem statistics for an ext3 `.img` overlay.

**Usage:**

```
condatainer overlay info [image_path]
```

**Output includes:**

* File info (path, size, type)
* Filesystem info (format, state, UUID, creation time)
* Ownership (root or current UID)
* Disk usage
* Inode usage

**Example:**

```bash
condatainer overlay info env.img
```

### Overlay Check

Verify filesystem integrity of an ext3 `.img` overlay.

**Usage:**

```
condatainer overlay check [image_path]
```

**Options:**

* `-f`, `--force`: Force check even if filesystem appears clean.

**Example:**

```bash
condatainer overlay check env.img
condatainer overlay check -f env.img
```

### Overlay Chown

Change the ownership of files inside an ext3 `.img` overlay. This is useful when sharing a writable overlay with other users: the recipient can reset ownership to their UID/GID so files are writable with their account.

**Usage:**

```
condatainer overlay chown [OPTIONS] [image]
```

**Options:**

* `-u`, `--uid UID`    : Set the owner UID (default: current user's UID).
* `-g`, `--gid GID`    : Set the group GID (default: current user's GID).
* `--root`             : Set UID/GID to 0 (root).
* `-p`, `--path PATH`  : Path inside the overlay to change (default: /ext3).
* `image`              : Path to the overlay `img` file.

**Examples:**

```bash
# Set ownership inside env.img to the current user
condatainer overlay chown env.img

# Explicitly set UID and GID
condatainer overlay chown -u 1001 -g 1001 env.img

# Set ownership to root
condatainer overlay chown --root env.img
```

### Overlay Resize

Resize an ext3 `.img` overlay to a new size. Supports both expanding and shrinking.

**Usage:**

```
condatainer overlay resize -s SIZE [image]
```

**Options:**

* `-s`, `--size SIZE`  : New size for the overlay image (e.g., `20G`, `2048M`). **Required.**
* `image`              : Path to the overlay `img` file.

**Examples:**

```bash
# Resize env.img to 20GB
condatainer overlay resize -s 20G env.img
```

## Create

Initialize and build a new **CondaTainer** SquashFS overlay. You can build from existing recipes (local/remote), a Conda environment file, or a remote container source.

**Usage:**

```
condatainer create [OPTIONS] [packages...]
```

**Aliases:** `install`, `i`

**Options:**

* `-n`, `--name [NAME]`: Custom name for the resulting overlay file. If used, all specified packages are bundled into one overlay.
* `-p`, `--prefix [PATH]`: Custom prefix path for the overlay file.
* `-f`, `--file [FILE]`: Path to definition file (.yaml, .sh, .def).
* `-b`, `--base-image [PATH]`: Base image to use instead of default.
* `-s`, `--source [URI]`: Remote source URI (e.g., `docker://ubuntu:22.04`).
* `--temp-size [SIZE]`: Size of temporary overlay (default: 20G).
* `packages`: List of packages to install (e.g., `bcftools/1.22` or `samtools=1.10` or `grch38/genome/gencode`).

**Compression Options:**

* `--lz4`: Use LZ4 compression (default).
* `--zstd-fast`: Use Zstandard compression level 3.
* `--zstd-medium`: Use Zstandard compression level 8.
* `--zstd`: Use Zstandard compression level 14.
* `--zstd-high`: Use Zstandard compression level 19.
* `--gzip`: Use Gzip compression.

**Build Modes:**

* **Default:** Each package gets its own `.sqf` via the build system.
* **`--name`:** Create a single `.sqf` with multiple packages bundled together.
* **`--prefix` + `--file`:** Create `.sqf` from external source file (.sh, .def, .yml).
* **`--source`:** Create `.sqf` from a remote container source URI.

### System level Examples

```bash
# Create from a specific recipe (App Overlay)
condatainer create bcftools/1.22

# Create multiple overlays at once
condatainer create samtools/1.16 bcftools/1.15

# Create a reference overlay (Ref Overlay)
condatainer create grch38/gtf-gencode/47

# Create from a yaml file with a custom name (Custom Env)
# NOTE: `-f/--file` must be used with `-p/--prefix`.
condatainer create -p my_analysis_env -f environment.yml

# Create a custom env with multiple packages bundled together
condatainer create -n tools samtools=1.16 bcftools=1.15

# Create from a remote container source
condatainer create --source docker://ubuntu:22.04 -n myubuntu
```

**Features**:

- Automatic Fetching: If a build script is not found locally, **CondaTainer** attempts to fetch it from the remote repository.
- Conda Fallback: If no build script exists, **CondaTainer** attempts to create the module by installing the package with the requested name and version from conda-forge or bioconda.
- Metadata Parsing: Parses `#ENV` and `#ENVNOTE` tags from build scripts to inject environment variables and help text into the generated modulefile.

### Project level Examples

Create a read-only overlay with a Conda environment using a Conda YAML file.

```bash
condatainer create -p my_project_env -f environment.yml
```

Create a read-only overlay using a custom Apptainer definition file. See [Custom OS Overlays](../advanced_usage/custom_os.md) for more details.

```bash
condatainer create -p my_project_env -f custom_def.def
```

Create a read-only overlay using a shell script that installs packages. See [Custom Build Script using Build Scripts](../advanced_usage/custom_bundle.md) for more details.

```bash
condatainer create -p my_project_env -f install_packages.sh
```

These methods will generate `my_project_env.sqf` in the current directory.

## Container Management (Avail, List, Remove)

Manage your local library of built containers and available recipes.

### Avail

Search for available build scripts (both local scripts in build-scripts/ and remote metadata). Search terms are combined using **AND** logic (all terms must be present).

**Aliases:** `av`

```
condatainer avail [search_terms...] [flags]
```

**Options:**

* `--local`: Show only local build scripts.
* `--remote`: Show only remote build scripts from GitHub.
* `-i`, `--install`: Install any found packages that are not currently installed.
* `-a`, `--add`: Alias for `--install`.

**Examples:**

```bash
# Search for all resources matching 'grcm' and 'M33'
$ condatainer avail grcm M33
grcm39/gtf-gencode/M33
grcm39/salmon/1.10.2/gencodeM33
grcm39/star-2.7.11b/gencodeM33-101
grcm39/star-2.7.11b/gencodeM33-151
grcm39/transcript-gencode/M33

# Search for salmon using specific reference tags
$ condatainer avail grcm M33 salmon
grcm39/salmon/1.10.2/gencodeM33

# Directly install from search results (will have prompt)
$ condatainer avail grcm M33 salmon -i
grcm39/gtf-gencode/M33
[CondaTainer] Do you want to install the above overlays? [y/N]:
```

### List

List installed overlays stored in the `images/` directory (app and reference overlays now share that location).

**Aliases:** `ls`

```
condatainer list [search_terms...] [flags]
```

**Options:**

* `-d`, `--delete`: Prompt to delete listed overlays after displaying them.
* `-r`, `--remove`: Alias for `--delete`.

**Features:**

* Lists both app overlays (with versions) and reference overlays.
* Searches across all image directories.
* Marks system apps and env overlays.
* Uses AND logic for multiple search terms.

### Remove

Deletes specific overlays and their associated .env files.

**Aliases:** `rm`, `delete`, `uninstall`

```
condatainer remove [search_terms...]
```

**Features:**

* Exact match mode for single overlay names.
* Search mode with AND logic for multiple terms.
* Prompts for confirmation before deletion.
* Only removes from writable directories.

## Exec

Execute a command inside a containerized environment composed of one or more overlays.

**Usage:**

```
condatainer exec [flags] [overlays...] [--] [command...]
```

**Aliases:** `e`

**Options:**

* `-o`, `--overlay [OVERLAY]`: Overlay file to mount (can be used multiple times).
* `-w`, `--writable-img`: Mount `.img` overlays as writable (default: read-only).
* `-r`, `--read-only`: Mount `.img` overlays as read-only (only for `e` shortcut).
* `-k`, `--keep`: Disable automatic command parsing (don't convert commands to overlay names).
* `-n`, `--no-autoload`: Disable autoloading `env.img` from current directory (only for `e` shortcut).
* `-b`, `--base-image [PATH]`: Base image to use instead of default.
* `--env [KEY=VALUE]`: Set environment variable inside the container (can be used multiple times).
* `--bind [HOST:CONTAINER]`: Bind mount path into the container (can be used multiple times).
* `--fakeroot`: Run container with fakeroot privileges.
* COMMAND: The command to execute inside the container. If not provided, opens an interactive bash shell.

**Special Features:**

* **Auto-detection:** Overlay names are automatically converted from command names (e.g., `samtools` → `samtools/1.22` if installed).
* **Shortcut `e`:** The `e` alias provides intelligent parsing with automatic `env.img` loading from the current directory.

**Environment Variables (inside container):**

* `IN_CONDATAINER=1`: Set inside the container.
* `CNT_CONDA_PREFIX`: Path to the `.img` default conda path, if `-w`, `--writable-img` is used.

**Examples:**

```bash
# Run a command using a specific overlay
condatainer exec -o bcftools/1.22 bcftools --version

# Run with overlay specified as positional argument
condatainer exec samtools/1.22

# Run with writable overlay and custom environment variable
condatainer exec -w env.img --env MYVAR=value bash

# Use the 'e' shortcut (auto-loads env.img if present)
e samtools/1.22

# Disable autoload of env.img
e --no-autoload samtools/1.22 bash

# Run with fakeroot privileges
condatainer exec --fakeroot env.img bash
```

### Autoload and Shell Completion

- **Autoload local env images:** When running `condatainer e` in a directory that contains `env.img`, **CondaTainer** will automatically mount it into the container and open an interactive shell. Use `-n` / `--no-autoload` to disable this behavior.
- **Enable shell completion:** Generate and install the completion script for your shell. See the [Completion](#completion) section for details.

## Runtime (Check, Run)

Utilities for running scripts with automatic dependency handling via  `#DEP:` tags and `module load` or `ml` commands.

### Check

Parses a script for `#DEP:` tags and `module load` or `ml` commands and checks if the required overlays are installed.

```
condatainer check [SCRIPT] [-a]
```

* `-a`, `--auto-install`: Automatically attempt to build/install missing dependencies found in the script.
* `-i`, `--install`: Alias for `--auto-install`.

### Run

Executes a script inside the **CondaTainer** environment, mounting dependencies defined in the script. Autosolves dependencies based on `#DEP:` tags and `module load` or `ml` commands within the script.

```
condatainer run [SCRIPT] [SCRIPT_ARGS...]
```

**Options:**

* `-w`, `--writable-img`: Make `.img` overlays writable (default: read-only).
* `-b`, `--base-image [PATH]`: Use custom base image.
* `-a`, `--auto-install`: Automatically install missing dependencies.
* `-i`, `--install`: Alias for `--auto-install`.

**Script Tags:**

Scripts can use special comment tags to declare dependencies and configure the container:

| Tag | Description |
|-----|-------------|
| `#DEP: package/version` | Declare a dependency overlay |
| `#CNT [args]` | Additional arguments passed to condatainer |

**Available `#CNT` arguments:**

* `-w`, `--writable-img`: Make `.img` overlays writable.
* `-b`, `--base-image PATH`: Use custom base image.
* `--env KEY=VALUE`: Set environment variable.
* `--bind HOST:CONTAINER`: Bind mount path.
* `--fakeroot`: Run with fakeroot privileges.

**Script Example:**

```bash
#!/bin/bash
#DEP: bcftools/1.22
#CNT --writable-img
#CNT --env MYVAR=value
bcftools view input.vcf | head
```

### CondaTainer is compatible with module systems

**CondaTainer** will scan your script for `module load` or `ml` commands and mount the corresponding overlays automatically.

**Example:**

```bash
#!/bin/bash
module load bcftools/1.22
bcftools --version
```

You can run `check` or `run` commands to automatically handle the dependencies.

```bash
# Install missing dependencies
condatainer check my_script.sh -a

# Run the script with automatic dependency resolution
condatainer run my_script.sh
```

## Info

Displays metadata regarding a specific overlay.

**Usage:**

```
condatainer info [OVERLAY]
```

**Output includes:**

* Writable status.
* File size.
* Compression type (for `.sqf` files).
* Timestamps (status change, modification, access).
* Mount path (e.g., `/cnt/bcftools/1.22` or `/ext3/env`).
* File ownership (for `.img` files).
* Environment variables defined within the overlay's `.env` file.

**Examples:**

```bash
condatainer info samtools/1.22
condatainer info env.img
```

## Helper

Download and manage small helper scripts stored in the `helper-scripts/` folder inside the CondaTainer repository.

**Usage:**

```
condatainer helper [-u|--update] [-p|--path] [SCRIPT_NAME] [SCRIPT_ARGS...]
```

Options:

* `-u`, `--update`: Update helper scripts from remote metadata.
* `-p`, `--path`: Print the absolute path of the helper folder or a specific helper script and exit.
* `SCRIPT_NAME`: Name of the helper script to run or update (optional).
* `SCRIPT_ARGS...`: Remaining arguments are passed directly to the helper script when running it.

**Examples**

```bash
# Print helper folder path
condatainer helper --path

# Download/Update all helper scripts (auto selects sbatch or headless based on availability)
condatainer helper --update
# Forcely Update headless scripts version
condatainer --local helper --update

# Run a helper script with arguments
condatainer helper code-server -p 18230
```

## Config

Manage **CondaTainer** configuration settings.

### Config Show

Display current configuration including file paths, settings, and environment variable overrides.

```
condatainer config show [--path]
```

* `--path`: Show only the config file path.

### Config Get

Get a specific configuration value.

```
condatainer config get <key>
```

**Examples:**

```bash
condatainer config get apptainer_bin
condatainer config get build.default_cpus
condatainer config get submit_job
```

### Config Set

Set a configuration value.

```
condatainer config set <key> <value>
```

**Examples:**

```bash
condatainer config set apptainer_bin /usr/bin/apptainer
condatainer config set build.default_cpus 8
condatainer config set build.default_time 4h
condatainer config set submit_job false
```

**Time formats:** `2h`, `30m`, `1h30m`, `90s`, `02:00:00`, `HH:MM:SS`

### Config Init

Create a config file with auto-detected defaults.

```
condatainer config init [-l|--location user|portable|system]
```

* `-l`, `--location`: Config location (default: auto-detect).

**Auto-detects:**
* Apptainer binary
* Scheduler binary
* Compression support (zstd vs lz4)

### Config Edit

Edit configuration file in your default editor.

```
condatainer config edit
```

Opens the config file in `$EDITOR` (falls back to `vi`).

### Config Paths

Show data search paths for images, build scripts, and helper scripts.

```
condatainer config paths
```

### Config Validate

Validate the current configuration.

```
condatainer config validate
```

**Checks:**
* Apptainer binary accessibility
* Scheduler binary accessibility
* Build configuration (CPUs > 0, Memory > 0)

### Configuration Priority

Configuration is loaded in the following order (highest to lowest priority):

1. Command-line flags
2. Environment variables (`CONDATAINER_*`)
3. User config file (`~/.config/condatainer/config.yaml`)
4. Portable config (`<install-dir>/config.yaml`)
5. System config (`/etc/condatainer/config.yaml`)
6. Built-in defaults

### Configuration File Example

```yaml
base_dir: /path/to/data

# Binary paths
apptainer_bin: /usr/bin/apptainer
scheduler_bin: /usr/bin/srun
scheduler_type: slurm

# Submission settings
submit_job: true

# Build configuration
build:
  default_cpus: 8
  default_mem_mb: 16384
  default_time: 4h
  tmp_size_mb: 20480
  compress_args: "-comp zstd -Xcompression-level 8"
  overlay_type: squashfs

# Additional data directories
extra_base_dirs:
  - /path/to/shared/data
  - /path/to/other/location
```

### Environment Variables

| Variable | Description |
|----------|-------------|
| `CONDATAINER_BASE_DIR` | Override base directory |
| `CONDATAINER_LOGS_DIR` | Override logs directory |
| `CONDATAINER_APPTAINER_BIN` | Override apptainer binary path |
| `CONDATAINER_SCHEDULER_BIN` | Override scheduler binary path |
| `CONDATAINER_SUBMIT_JOB` | Override job submission setting (true/false) |
| `CONDATAINER_BUILD_DEFAULT_CPUS` | Override default CPUs for builds |
| `CONDATAINER_BUILD_DEFAULT_MEM_MB` | Override default memory for builds |
| `CONDATAINER_BUILD_DEFAULT_TIME` | Override default build time |
| `CONDATAINER_EXTRA_BASE_DIRS` | Colon-separated list of extra base directories |

## Scheduler

Display information about the configured job scheduler.

**Usage:**

```
condatainer scheduler
```

**Output includes:**

* Scheduler type (SLURM, PBS, etc.)
* Binary path
* Version
* Availability status
* Available GPUs (if any)
* Resource limits (CPUs, memory, GPUs, time)

**Example:**

```bash
$ condatainer scheduler
Scheduler: SLURM
Binary: /usr/bin/srun
Version: 23.02.4
Status: Available
GPUs: nvidia_a100 (4), nvidia_v100 (8)
```

## Update

Updates the **CondaTainer** binary to the latest version from GitHub releases.

**Aliases:** `update`

**Usage:**

```
condatainer self-update [-y]
```

* `-y`, `--yes`: Skip confirmation prompt and auto-update.

**Features:**

* Downloads latest binary from GitHub releases.
* Detects current OS and architecture.
* Compares versions before updating.
* Creates temporary backup before replacing.
* Supports symlink resolution.

## Completion

Generates shell completion scripts for CondaTainer.

**Usage:**

```
condatainer completion [SHELL]
```

* SHELL: Specify the shell type (`bash`, `zsh`, or `fish`). If not provided, auto-detected from `$SHELL`.

**Bash:**

Add the following to your `~/.bashrc`:

```bash
source <(condatainer completion bash)
```

**Zsh:**

Add the following to your `~/.zshrc`:

```zsh
source <(condatainer completion zsh)
# Or generate to fpath:
# condatainer completion zsh > "${fpath[1]}/_condatainer"
```

**Fish:**

Run the following command:

```fish
condatainer completion fish | source
# Or save to config:
# condatainer completion fish > ~/.config/fish/completions/condatainer.fish
```
