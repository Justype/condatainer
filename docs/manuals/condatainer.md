# CondaTainer Manual

📦 **CondaTainer** is a wrapper script designed to streamline the management of Apptainer (Singularity) containers backed by Conda environments and SquashFS overlays.

## Table of Contents

- [Naming Convention](#naming-convention)
- [Mount Points](#mount-points)
- [Create](#create)
- [Overlay](#overlay)
- [Container Management (Avail, List, Remove)](#container-management-avail-list-remove)
- [Exec](#exec)
- [Runtime (Check, Run)](#runtime-check-run)
- [Info](#info)
- [Update](#update)
- [Completion](#completion)
- [Helper](#helper)

## Overall Command Structure

```
usage: condatainer [-h] [-v] [--debug] [--local] COMMAND ...

CondaTainer: Use Apptainer/Conda/SquashFS to manage tools for HPC users.

positional arguments:
  COMMAND              Available actions
    overlay            Manage overlay images (create, chown, resize)
    o                  Shortcut for 'overlay create'
    create (install, i)
                       Create a new SquashFS overlay using available build scripts or Conda
    helper             Run or update helper scripts from helper-scripts/
    avail (av)         Check available local and remote build scripts
    list (ls)          List installed overlays matching search terms
    remove (rm, delete)
                       Remove installed overlays matching search terms
    exec               Execute a command using overlays
    e                  Run bash using writable overlays
    check              Check if the dependencies of a script are installed
    run                Run a script and auto-solve the dependencies by #DEP tags
    info               Show information about a specific overlay
    self-update        Update CondaTainer to the latest version
    completion         Generate shell completion script 'source <(condatainer completion)'

options:
  -h, --help           show this help message and exit
  -v, --version        Show program's version number and exit
  --debug              Enable debug mode with verbose output
  --local              Disable job submission
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

- create: create an image with a conda environment inside.
- chown: change ownership of files inside the image.
- resize: resize the image file.

### Overlay Create

Create an ext3 `.img` with a conda environment inside. This is useful for creating writable conda environments.

**Usage:**

```
condatainer overlay create [OPTIONS] [NAME]
```

**Options:**

* `-s`, `--size [SIZE]`: Size (MB) of the overlay image (default: 10G).
* `-f`, `--file [FILE]`: Path to a Conda environment file (.yml or .yaml).
* `--fakeroot`: Create image compatible with fakeroot (must use with `--fakeroot` later).
* `--sparse`: Create a sparse image file.
* NAME: Name of the overlay image (`env.img` by default if not specified).

**Examples:**

```bash
# Create a 20GB overlay img with a custom name
condatainer overlay create -s 20G my_analysis_env

# o is a shortcut for 'overlay create'
condatainer o -f environment.yml project_env.img
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

### Overlay Chown

Change the ownership of files inside an ext3 `.img` overlay. This is useful when sharing a writable overlay with other users: the recipient can reset ownership to their UID/GID so files are writable with their account.

**Usage:**

```
condatainer overlay chown [-h] [-u UID] [-g GID] [-p PATH] image
```

**Options:**

* `-u`, `--uid UID`    : Set the owner UID (default: current user's UID).
* `-g`, `--gid GID`    : Set the group GID (default: current user's GID).
* `-p`, `--path PATH`  : Path inside the overlay to change (default: /ext3).
* `image`              : Path to the overlay `img` file.

**Examples:**

```bash
# Set ownership inside env.img to the current user
condatainer overlay chown env.img

# Explicitly set UID and GID
condatainer overlay chown -u 1001 -g 1001 env.img
```

### Overlay Resize

Resize an ext3 `.img` overlay to a new size.

```
condatainer overlay resize [-h] -s SIZE image
```

**Options:**

* `-s`, `--size SIZE`  : New size (MB) for the overlay image.
* `image`              : Path to the overlay `img` file.

**Examples:**

```bash
# Resize env.img to 10GB
condatainer overlay resize -s 10G env.img
```

## Create

Initialize and build a new **CondaTainer** SquashFS overlay. You can build from existing recipes (local/remote) or a Conda environment file.

**Usage:**

```
condatainer create [OPTIONS] [NAME_VERSIONS...]
```

**Options:**

* `-n`, `--name [NAME]`: Custom name for the overlay file. If used, all specified packages are bundled into one overlay.
* `-p`, `--prefix [PATH]`: Custom prefix path for the overlay file.
* `-f`, `--file [FILE]`: Path to a Conda environment file (.yml/.yaml/.def/.sh).
* NAME_VERSIONS: List of packages to install (e.g., `bcftools/1.22` or `samtools=1.10` or `grch38/genome/gencode`).

**Compression Options:**

* `--zstd` / `--zstd-medium` / `--zstd-fast`: Use Zstandard compression (levels 14, 8, or 3).
* `--gzip`: Use Gzip compression.
* `--lz4`: Use LZ4 compression.

### System level Examples

```bash
# Create from a specific recipe (App Overlay)
condatainer create bcftools/1.22

# Create a reference overlay (Ref Overlay)
condatainer create grch38/gtf-gencode/47

# Create from a yaml file with a custom name (Custom Env)
# NOTE: `-f/--file` must be used with `-p/--prefix` (the script requires `--prefix` when using `--file`).
condatainer create -p my_analysis_env -f environment.yml

# Create a custom env with multiple packages
condatainer create -n tools.sqf samtools=1.16 bcftools=1.15
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

```
condatainer avail [search_terms...] [-i|--install]
```

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

```
condatainer list [-d] [search_terms...]
```

* `-d`, `--delete`: Prompt to delete listed overlays after displaying them.
* `-r`, `--remove`: Alias for `--delete`.

### Remove

Deletes specific overlays and their associated .env files.

```
condatainer remove [search_terms...]
```

* Requires at least one search term.
* Prompts for confirmation before deletion.

## Exec

Execute a command inside a containerized environment composed of one or more overlays.

**Usage:**

```
condatainer exec [OPTIONS] [COMMAND...]
```

**Options:**

* `-o`, `--overlay [OVERLAY]`: Specify specific overlays to mount (can be used multiple times).
* `-w`, `--writable-img`: Make a `.img` overlay writable (default: read-only). Only one `.img` overlay may be used at a time; use this flag to enable writability.
* `-k`, `--keep`: Do not attempt to parse the command itself as an overlay name.
* `--env [ENV_VAR=VALUE]`: Set additional environment variables inside the container (can be used multiple times).
* `--bind [HOST_PATH:CONTAINER_PATH]`: Bind additional host paths into the container (can be used multiple times).
* COMMAND: The command to execute inside the container. If not provided, opens an interactive bash shell.

**Environment Variables:**

* `IN_CONDATAINER=1`: Set inside the container.
* `CNT_CONDA_PREFIX`: Paths to the `.img` default conda path, if `-w`, `--writable-img` is used.

**Example:**

```bash
# Run a command using a specific overlay
condatainer exec -o bcftools/1.22 bcftools --version
```

### Autoload and Bash completion

- **Autoload local env images:** When running `condatainer e` in a directory that contains `env.img`, **CondaTainer** will automatically mount it into the container and open an interactive shell.
- **Enable Bash completion on the host:** Generate and install the completion script for your shell. Example for Bash (per-user):

Add the following to your `~/.bashrc`:

```bash
source <(condatainer completion bash)
```

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

* **Dependency Injection:** Reads `#DEP:`and `module load` or `ml` commands lines in the script to determine which overlays to mount.
* **Argument Injection:** Reads `#CNT` lines in the script to inject arguments into the condatainer command itself.
* `-w`, `--writable-img`: Enable writable mounting for .img overlays.
* `-a`, `--auto-install`: Automatically attempt to build/install missing dependencies found in the script.
* `-i`, `--install`: Alias for `--auto-install`.

**Script Example:**

```bash
#!/bin/bash
#DEP: bcftools/1.22
#CNT --writable-img
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

* File size.
* Potential mount path (e.g., `/cnt/bcftools/1.22` or `/ext3/env`).
* Environment variables defined within the overlay's `.env` file.

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

## Update

Updates the **CondaTainer** script itself to the latest version from the GitHub repository.

**Usage:**

```
condatainer self-update [-y]
```

* `-y`, `--yes`: Automatically confirm the download and replacement of the current script.

## Completion

Generates shell completion scripts for CondaTainer. (Bash, Zsh)

**Usage:**

```
condatainer completion [SHELL]
```

* SHELL: Specify the shell type (e.g., `bash`, `zsh`). If not provided, auto-detected.

For Bash, you can enable completion by adding the following to your `~/.bashrc`:

```bash
source <(condatainer completion bash)
```

For Zsh, add the following to your `~/.zshrc`: (not tested)

```zsh
source <(condatainer completion zsh)
```
