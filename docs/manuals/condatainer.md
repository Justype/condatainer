# CondaTainer Manual

📦 **CondaTainer** is a CLI tool designed to streamline the management of Apptainer (Singularity) containers backed by Conda environments and SquashFS overlays.

## Table of Contents

- [Naming Convention](#naming-convention)
- [Mount Points](#mount-points)
- [Create](#create)
- [Overlay](#overlay)
- [Container Management (Avail, List, Remove)](#container-management-avail-list-remove)
- [Exec](#exec)
- [E (Quick Exec)](#e-quick-exec)
- [Runtime (Check, Run)](#runtime-check-run)
- [Instance](#instance)
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
  e           Quick shortcut for executing commands with overlays
  exec        Execute a command using overlays (explicit -o flag)
  helper      Manage and run helper scripts
  info        Show information about a specific overlay
  instance    Manage Apptainer instances
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
  -q, --quiet     Suppress messages (warnings/errors are still shown)
  -v, --version   version for condatainer
  -y, --yes       Automatically answer yes to all prompts

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

Writability: `.sqf` overlays are read-only. To enable write access for a `.img` overlay use the `-w` / `--writable` flag with `exec` or `run`.

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
* `--root`             : Set UID/GID to 0 (root); overrides `-u` and `-g`.
* `-p`, `--path PATH`  : Path inside the overlay to change (can specify multiple, default: `/ext3:/opt`).
* `image`              : Path to the overlay `img` file.

**Examples:**

```bash
# Set ownership inside env.img to the current user (default paths /ext3 and /opt)
condatainer overlay chown env.img

# Explicitly set UID and GID
condatainer overlay chown -u 1001 -g 1001 env.img

# Set ownership to root
condatainer overlay chown --root env.img

# Set entire image to root
condatainer overlay chown --root -p / env.img

# Chown only /ext3 path
condatainer overlay chown -p /ext3 env.img

# Chown multiple specific paths
condatainer overlay chown -p /ext3 -p /data env.img
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
* `--remote`: Remote build scripts take precedence over local.
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

* `--remote`: Remote build scripts take precedence over local (on duplicates).
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
* `-e`, `--exact`: Require exact match instead of substring match.

**Features:**

* Lists both app overlays (with versions) and reference overlays.
* Searches across all image directories.
* Marks system apps and env overlays.
* Uses AND logic for multiple search terms (substring mode).
* Exits with code 1 if no overlays match the search terms.

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

Execute a command inside a containerized environment using explicit overlay specifications.

**Usage:**

```
condatainer exec [flags] [command...]
```

**Options:**

* `-o`, `--overlay [OVERLAY]`: Overlay file to mount (can be used multiple times).
* `-w`, `--writable`: Mount `.img` overlays as writable (default: read-only).
* `-b`, `--base-image [PATH]`: Base image to use instead of default.
* `-f`, `--fakeroot`: Run container with fakeroot privileges.
* `--env [KEY=VALUE]`: Set environment variable inside the container (can be used multiple times).
* `--bind [HOST:CONTAINER]`: Bind mount path into the container (can be used multiple times).

**Features:**

* Use `-o/--overlay` to explicitly specify overlays.
* All positional arguments are treated as commands.
* Read-only by default for `.img` overlays (use `-w` for writable).
* Defaults to bash if no command specified.

**Environment Variables (inside container):**

* `IN_CONDATAINER=1`: Set inside the container.
* `CNT_CONDA_PREFIX`: Path to the `.img` default conda path, if `-w/--writable` is used.

**Examples:**

```bash
# Run bash with samtools overlay
condatainer exec -o samtools/1.22

# Run samtools command with overlay
condatainer exec -o samtools/1.22 samtools view file.bam

# Use writable .img overlay
condatainer exec -w -o env.img bash

# Multiple overlays
condatainer exec -o samtools/1.22 -o bcftools/1.20 bash

# Set environment variables
condatainer exec --env MYVAR=value -o samtools/1.22 bash

# Bind mount directories
condatainer exec --bind /data:/mnt -o env.img bash

# Run with fakeroot privileges (short flag)
condatainer exec -f -o env.img bash

# Pass apptainer flags (use --flag=value format)
condatainer exec --nv --home=/custom -o samtools/1.22 python gpu_script.py
```

## E (Quick Exec)

Quick shortcut for executing commands with overlays using simplified syntax.

**Usage:**

```
condatainer e [flags] [overlays...] [--] [command...]
```

**Options:**

* `-r`, `--read-only`: Mount `.img` overlays as read-only (default: writable).
* `-n`, `--no-autoload`: Disable autoloading `env.img` from current directory.
* `-b`, `--base-image [PATH]`: Base image to use instead of default.
* `-f`, `--fakeroot`: Run container with fakeroot privileges.
* `--env [KEY=VALUE]`: Set environment variable inside the container (can be used multiple times).
* `--bind [HOST:CONTAINER]`: Bind mount path into the container (can be used multiple times).

**Key Differences from `exec`:**

* Overlays are positional arguments before `--`.
* Commands go after `--`.
* Writable by default (use `-r` for read-only).
* Auto-loads `env.img` unless `-n` is specified.
* Defaults to bash if no command specified.

**Environment Variables (inside container):**

* `IN_CONDATAINER=1`: Set inside the container.
* `CNT_CONDA_PREFIX`: Path to the `.img` default conda path (writable by default).

**Examples:**

```bash
# Auto-load env.img if present, run bash
condatainer e

# Multiple overlays
condatainer e samtools/1.22 bcftools/1.20

# Run specific command
condatainer e samtools/1.22 -- samtools view file.bam

# Read-only .img overlay
condatainer e -r env.img

# Disable env.img auto-loading
condatainer e -n samtools/1.22

# Set environment variables
condatainer e --env MYVAR=value samtools/1.22

# Bind mount directories
condatainer e --bind /data:/mnt env.img

# Run with fakeroot privileges (short flag)
condatainer e -f env.img

# Pass apptainer flags (use --flag=value format)
condatainer e --home=/custom samtools/1.22
```

### Autoload and Shell Completion

- **Autoload local env images:** When running `condatainer e` in a directory that contains `env.img`, **CondaTainer** will automatically mount it into the container and open an interactive shell. Use `-n` / `--no-autoload` to disable this behavior.
- **Enable shell completion:** Generate and install the completion script for your shell. See the [Completion](#completion) section for details.

## Instance

Manage persistent Apptainer instances that run in the background. Instances are useful for long-running services or when you need to execute multiple commands in the same containerized environment without restarting the container each time.

**Key Features:**

* Persistent background containers
* State preservation (environment variables, overlays, bind paths)
* Multiple instances with unique names
* Pattern-based operations (wildcards supported)

### Instance Start

Start a named instance with specified overlays.

**Usage:**

```
condatainer instance start [flags] <name>
```

**Options:**

* `-o`, `--overlay [OVERLAY]`: Overlay file to mount (can be used multiple times)
* `-w`, `--writable`: Mount `.img` overlays as writable (default: read-only)
* `-b`, `--base-image [PATH]`: Base image to use instead of default
* `-f`, `--fakeroot`: Run instance with fakeroot privileges
* `--env [KEY=VALUE]`: Set environment variable (can be used multiple times)
* `--bind [HOST:CONTAINER]`: Bind mount path (can be used multiple times)

**Examples:**

```bash
# Start instance with multiple overlays
condatainer instance start -o samtools/1.22 -o bcftools/1.22 myinstance

# Start with writable .img overlay
condatainer instance start -w -o env.img data_analysis

# Pass apptainer flags
condatainer instance start --home=/custom -o samtools/1.22 pipeline1
```

**State Preservation:**

When you start an instance, CondaTainer saves the instance configuration (overlays, environment variables, bind paths) to a state file. This allows `instance exec` to automatically apply the same environment when executing commands in the instance.

### Instance Stop

Stop one or more running instances.

**Usage:**

```
condatainer instance stop [flags] [name]
```

**Options:**

* `-a`, `--all`: Stop all user's instances
* `-F`, `--force`: Force kill instance (may corrupt data)
* `-s`, `--signal [SIGNAL]`: Signal to send to instance (e.g., SIGTERM, TERM, 15)
* `-t`, `--timeout [SECONDS]`: Timeout before force kill (default: 10)

**Examples:**

```bash
# Stop a specific instance
condatainer instance stop myinstance

# Stop all instances matching pattern
condatainer instance stop mysql*

# Stop all instances
condatainer instance stop --all

# Force stop with custom signal
condatainer instance stop --force --signal TERM myinstance

# Custom timeout before force kill
condatainer instance stop --timeout 30 myinstance
```

**Wildcard Support:**

You can use shell wildcards (`*`, `?`, `[]`) to stop multiple instances at once. State files are automatically cleaned up for stopped instances.

### Instance List

List all running instances.

**Usage:**

```
condatainer instance list
```

Shows all active instances with their names and process information.

### Instance Stats

Display resource usage statistics for a running instance.

**Usage:**

```
condatainer instance stats <name>
```

**Example:**

```bash
condatainer instance stats myinstance
```

Shows CPU, memory, and I/O statistics for the specified instance.

### Instance Exec

Execute a command in a running instance.

**Usage:**

```
condatainer instance exec [flags] <name> <command> [args...]
```

**Options:**

* `--env [KEY=VALUE]`: Set additional environment variables (can be used multiple times)

**Features:**

* Commands run in the same environment as when the instance was started
* Automatically applies saved environment variables from instance state
* Overlays and bind paths from instance start are preserved
* Defaults to bash if no command specified
* Unknown apptainer flags can be passed using `--flag=value` format

**Examples:**

```bash
# Run a command in an instance
condatainer instance exec myinstance samtools view file.bam

# Open interactive shell (default)
condatainer instance exec myinstance

# Set additional environment variables
condatainer instance exec --env VAR1=val1 --env VAR2=val2 myinstance bash

# Pass apptainer flags
condatainer instance exec --home=/custom myinstance bash
```

**State Restoration:**

When you use `instance exec`, CondaTainer automatically loads the saved state from when the instance was started. This includes:

* Environment variables (`#ENV` tags from overlays)
* Overlay mount points
* Bind paths
* Base image configuration

Environment variables set during `instance exec` are added on top of the saved state, allowing you to override or add variables as needed.

**Note:** Apptainer's native `apptainer exec instance://name` does not preserve environment variables from instance start. CondaTainer solves this by storing instance state and reapplying it during exec.

### Instance Workflow Example

Here's a complete workflow showing how to use instances:

```bash
# 1. Start an instance with your desired overlays and writable env
condatainer instance start -o lxde -o igv -o env.img -w desktop

# 2. Run multiple commands in the same instance
condatainer instance exec desktop websockify --web /usr/share/novnc ...
condatainer instance exec desktop vncserver :1 -geometry 1920x1080 -depth 24 ...
condatainer instance exec desktop igv

# 3. Check instance status
condatainer instance stats desktop

# 4. Stop the instance when done
condatainer instance stop desktop
```

## Runtime (Check, Run)

Utilities for running scripts with automatic dependency handling via `#DEP:` tags. Parsing of `module load` / `ml` lines is disabled by default and can be enabled with `--module` or `parse_module_load: true` in config.

### Check

Parses a script for `#DEP:` tags and checks if the required overlays are installed.

```
condatainer check [SCRIPT] [-a] [--module]
```

* `-a`, `--auto-install`: Automatically attempt to build/install missing dependencies found in the script.
* `-i`, `--install`: Alias for `--auto-install`.
* `--module`: Also parse `module load` / `ml` lines as dependencies.

### Run

Executes a script inside the **CondaTainer** environment, mounting dependencies defined in the script. Autosolves dependencies based on `#DEP:` tags within the script.

```
condatainer run [SCRIPT] [SCRIPT_ARGS...]
```

**Options:**

* `-w`, `--writable`: Make `.img` overlays writable (default: read-only).
* `-b`, `--base-image [PATH]`: Use custom base image.
* `-a`, `--auto-install`: Automatically install missing dependencies.
* `-i`, `--install`: Alias for `--auto-install`.
* `--module`: Also parse `module load` / `ml` lines as dependencies.

**Script Tags:**

Scripts can use special comment tags to declare dependencies and configure the container:

| Tag | Description |
|-----|-------------|
| `#DEP: package/version` | Declare a dependency overlay |
| `#DEP: path.img` | Declare a external overlay (.sqf and .img) |
| `#CNT [args]` | Additional arguments passed to condatainer |
| `#SBATCH [args]` | SLURM scheduler directives (auto-submit as job) |
| `#PBS [args]` | PBS scheduler directives (auto-submit as job) |
| `#BSUB [args]` | LSF scheduler directives (auto-submit as job) |

**Available `#CNT` arguments:**

* `-w`, `--writable`: Make `.img` overlays writable.
* `-b`, `--base-image PATH`: Use custom base image.
* `--env KEY=VALUE`: Set environment variable.
* `--bind HOST:CONTAINER`: Bind mount path.
* `-f`, `--fakeroot`: Run with fakeroot privileges.

**Script Example:**

```bash
#!/bin/bash
#DEP: bcftools/1.22
#CNT --writable
#CNT --env MYVAR=value
#CNT -f
bcftools view input.vcf | head
```

### Scheduler Integration

If your script contains scheduler directives (`#SBATCH`, `#PBS`, or `#BSUB`), `condatainer run` will automatically submit it as a scheduler job instead of running it locally. HTCondor uses native `.sub` submit files instead of in-script directives.

**Behavior:**

1. **Script has scheduler specs + scheduler available**: Submits the script as a job
2. **Script has scheduler specs + already inside a job**: Runs locally (avoids nested job submission)
3. **Script has scheduler specs + `--local` flag**: Runs locally
4. **Script without scheduler specs**: Runs locally (as before)

**Scheduler Script Example:**

```bash
#!/bin/bash
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --time=2:00:00
#DEP: samtools/1.22
#DEP: bcftools/1.22

samtools view -@ 4 input.bam | bcftools call -mv -o output.vcf
```

When you run this script:

```bash
# Shows missing deps, then submits as SLURM job
condatainer run analysis.sh

# Auto-install missing deps first, then submit as job
condatainer run analysis.sh -a
```

**Log Files:**

* If a log output path is specified in the script (e.g., `#SBATCH --output=...`), logs go to that path
* Otherwise, logs are written to the global logs directory (`~/logs` by default)
* Job scripts are created alongside the log files

### Exit Codes (script and job-submission behavior) ⚠️

CondaTainer uses specific exit codes so automation and downstream tooling can detect special states:

- `0` — Success (all requested builds completed locally or nothing to do)
- `1` — Generic error (invalid arguments, build failures, or other fatal errors)
- `3` — **Jobs submitted to scheduler** — overlays will be created asynchronously by scheduler jobs

Commands that may return exit code `3` when scheduler jobs were submitted include:

- `condatainer create ...`
- `condatainer check -a ...` (auto-install missing deps)
- `condatainer avail -i ...` (install from search results)

Quick example for shell scripts that detect the job-submitted state:

```bash
condatainer create samtools/1.22
if [ $? -eq 3 ]; then
  echo "Jobs submitted to scheduler — overlays will be created asynchronously"
  # Optionally: exit 0 or wait/monitor jobs here
fi
```

Note: When exit code `3` is returned, CondaTainer prints a message showing the number of jobs submitted (and can be extended to emit JSON or write job metadata for automation).

**Disabling Job Submission:**

```bash
# Run locally even with scheduler specs
condatainer --local run analysis.sh

# Or set in config
condatainer config set submit_job false
```

### MPI Auto-Detection

When a scheduler script requests more than one task (`--ntasks-per-node`, `--ntasks`, or PBS/LSF equivalents), `condatainer run` automatically detects the host MPI and wraps the job command with `mpiexec`:

**Detection order:**

1. `mpiexec` already in `PATH` → use it directly
2. `module avail -t openmpi` → pick the highest available version, load it with `module purge && module load`
3. Neither found → warn and run without MPI wrapper

**Generated job command:**

```bash
# Via direct PATH:
mpiexec condatainer run script.sh

# Via module system:
module purge && module load openmpi/4.1.5 && mpiexec condatainer run script.sh
```

**MPI Script Example:**

```bash
#!/bin/bash
#SBATCH --nodes=2
#SBATCH --ntasks-per-node=3
#SBATCH --mem=1G
#SBATCH --time=00:30:00

#DEP: mpi.img

python my_mpi_script.py
```

Run with:

```bash
condatainer run mpi_job.sh
```

CondaTainer detects `ntasks = 6`, finds `mpiexec`, and submits:

```bash
mpiexec condatainer run mpi_job.sh
```

Each MPI rank launches its own container, all sharing the same MPI communicator via SLURM's process management interface.

````{important}
You need to have the same major and minor version of OpenMPI installed inside the container as on the host.

```bash
ml av openmpi
# openmpi/4.1.5
condatainer e mpi.img -- mm-install mpi4py openmpi=4.1 -y
```
````

### CondaTainer is compatible with module systems

**CondaTainer** can scan your script for `module load` or `ml` commands and mount the corresponding overlays automatically. This is disabled by default; enable it with `--module` or by setting `parse_module_load: true` in your config.

**Example:**

```bash
#!/bin/bash
module load bcftools/1.22
bcftools --version
```

Use `--module` with `check` or `run` to handle these dependencies:

```bash
# Install missing dependencies (including module load lines)
condatainer check my_script.sh --module -a

# Run the script with automatic dependency resolution
condatainer run my_script.sh --module
```

To enable it permanently:

```bash
condatainer config set parse_module_load true
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

```{note}
Helper commands are not available inside a container or a scheduler job (except `--path` and `--list`).
```

**Usage:**

```
condatainer helper [FLAGS] [SCRIPT_NAME] [SCRIPT_ARGS...]
```

Options:

* `-u`, `--update`: Update helper scripts from remote metadata.
* `-p`, `--path`: Show all helper script search paths and the writable directory. If a `SCRIPT_NAME` is given, print the absolute path of that specific helper script and exit.
* `-l`, `--list`: List available helper scripts with their descriptions (from `#WHATIS` tags).
* `SCRIPT_NAME`: Name of the helper script to run (optional).
* `SCRIPT_ARGS...`: Remaining arguments are passed directly to the helper script when running it.

**Examples**

```bash
# List all available helper scripts with descriptions
condatainer helper --list

# Print helper script search paths
condatainer helper --path

# Print path to a specific helper script
condatainer helper --path code-server

# Download/Update all helper scripts (auto selects sbatch or headless based on availability)
condatainer helper --update
# Force update headless scripts version
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
condatainer config get build.ncpus
condatainer config get scheduler.ncpus_per_task
condatainer config get submit_job
```

### Config Set

Set a configuration value and save to the active config file. For a full list of supported keys, see the [Configuration manual](configuration.md).

```
condatainer config set <key> <value>
```

**Examples:**

```bash
condatainer config set apptainer_bin /usr/bin/apptainer
condatainer config set submit_job false
condatainer config set scheduler.ncpus_per_task 8
condatainer config set build.ncpus 8
condatainer config set build.time 4h
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

See the [Configuration manual](configuration.md) for a full reference of all available keys.

```yaml
# Binary paths (scheduler type is auto-detected from binary)
apptainer_bin: /usr/bin/apptainer
scheduler_bin: /usr/bin/sbatch

# Submission settings
submit_job: true

# Default scheduler specs (used when scripts lack explicit directives)
scheduler:
  ncpus_per_task: 4
  mem_mb_per_node: 8192
  time: 4h

# Build configuration
build:
  ncpus: 8
  mem_mb: 16384
  time: 4h
  tmp_size_mb: 20480
  # compress_args options (gzip, lz4, zstd, zstd-fast, zstd-medium, zstd-high)
  # Or explicitly set mksquashfs arguments:
  compress_args: "-comp zstd -Xcompression-level 8"
  overlay_type: ext3

# Additional data directories
extra_base_dirs:
  - /path/to/shared/data
  - /path/to/other/location
```

## Scheduler

Display information about the configured job scheduler.

**Usage:**

```
condatainer scheduler [FLAGS]
```

**Options:**

* `-p`, `--partitions`: Show per-partition resource limits.
* `-Q`, `--queue`: Show per-queue resource limits (alias for `-p`).
* `--cpu`: Show only CPU-only partitions (no GPUs); automatically enables `-p`.
* `--gpu`: Show only GPU partitions; automatically enables `-p`.

**Output includes:**

* Scheduler type (SLURM, PBS, LSF, HTCondor)
* Binary path
* Version
* Availability status
* Max resource limits across all partitions (CPUs, memory, time)
* Available GPUs (type, total, available)

**Examples:**

```bash
# Show scheduler information
condatainer scheduler

# Show per-partition/queue resource limits
condatainer scheduler -p

# Show only GPU partitions
condatainer scheduler --gpu

# Show per-partition limits for CPU-only partitions
condatainer scheduler -p --cpu
```

## Update

Updates the **CondaTainer** binary to the latest version from GitHub releases.

**Aliases:** `update`

**Usage:**

```
condatainer self-update [FLAGS]
```

**Options:**

* `-y`, `--yes`: Skip confirmation prompt and auto-update.
* `-f`, `--force`: Force update even if already on the latest version.
* `--dev`: Include pre-release versions (also enabled automatically when the config `branch` is set to `dev`).

**Features:**

* Downloads latest binary from GitHub releases.
* Detects current OS and architecture.
* Compares versions before updating.
* Updates base image when minor or major version changes (not for patch updates).
* Warns about major version upgrades and suggests rebuilding def-built containers.
* Supports symlink resolution.

**Examples:**

```bash
# Update to latest stable version
condatainer self-update

# Update without confirmation
condatainer self-update --yes

# Force update even if already on latest version
condatainer self-update -f

# Include pre-release versions
condatainer self-update --dev
```

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
