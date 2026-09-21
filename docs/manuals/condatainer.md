# CondaTainer Manual

📦 **CondaTainer** is a CLI tool designed to streamline the management of Apptainer (Singularity) containers backed by Conda environments and SquashFS overlays.

## Table of Contents

- [Naming Convention](#naming-convention)
- [Mount Points](#mount-points)
- [Overlay](#overlay)
- [Create](#create)
- [Container Management (Avail, List, Remove, Search)](#container-management-avail-list-remove-search)
- [Exec](#exec)
- [E (Quick Exec)](#e-quick-exec)
- [Runtime (Check, Run)](#runtime-check-run)
- [Info](#info)
- [Export](#export)
- [Helper](#helper)
- [Config](#config)
- [Registry](#registry)
- [Store](#store)
- [Project](#project)
- [Scheduler](#scheduler)
- [Update](#update)
- [Proxy](#proxy)
- [Completion](#completion)

## Overall Command Structure

```
Single-file Conda env, tools, and data for HPC — plus app helpers to run RStudio and more.

Usage:
  condatainer [command]

Available Commands:
  avail             List available build scripts
  check             Check if the dependencies of script(s) are installed
  completion        Generate shell completion script
  config            Manage condatainer configuration
  create            Create a new SquashFS overlay
  e                 Shortcut for exec with overlays, writable by default
  env               Manage the Conda environment of the mounted overlay
  exec              Execute a command with overlays
  helper            Run web apps like RStudio on HPC
  info              Show details about an overlay
  list              List installed overlays
  o                 Shortcut for 'overlay create'
  overlay           Manage ext3 overlays (create, resize, check, info)
  project           Pin, restore and publish a project's artifacts
  proxy             Manage proxy tunnels for compute nodes
  registry          Publish and fetch artifacts through an OCI registry
  remove            Remove installed overlays matching search terms
  run               Run a script and auto-solve the dependencies by #DEP tags
  scheduler         Display scheduler information
  search            Search conda packages via anaconda.org
  self-update       Update condatainer to the latest version
  server            Manage the condatainer dashboard server
  store             Inspect overlays kept under an exact build identity
  update            Update script metadata caches or the base image

Flags:
      --debug     Enable debug mode with verbose output
  -h, --help      help for condatainer
  -q, --quiet     Suppress messages (warnings/errors are still shown)
  -v, --version   version for condatainer
  -y, --yes       Automatically answer yes to all prompts

Use "condatainer [command] --help" for more information about a command.
```

## Naming Convention

**CondaTainer** classifies overlays into three distinct categories based on their naming structure.

### System Applications (OS)

Used for system-level applications that do not follow the standard application or reference naming conventions.

Like text editors, IDEs, build-essential tools, etc.

* **Format:** `<distro>/<name>`
* **Example:** `ubuntu22/rstudio-server`, `ubuntu24/code-server`
* **Shortcut:** when `<distro>` is the configured `default_distro` (`ubuntu24`), it can be omitted — `code-server` resolves to `ubuntu24/code-server`. Standing in a project, its own pinned root takes over instead — see [Project Select-Distro](#project-select-distro).
* **Version:** `exec -o`, `e`, `create` and helper `#REQUIRED_OVERLAYS:` also accept a name without a version or with a partial one — `samtools` or `samtools/1.22` picks the newest matching version already installed, and only then the newest the recipes offer.

### Custom Environments (Bundle/Environment)

Used when creating a custom environment with a user-defined name (`create -p`, or `overlay`).

* **Format:** `custom_name`
* **Constraints:** Must **not** contain the double-dash sequence (`--`).
* **Example:** `my_analysis_env`, `project_x_utils`

### Application Overlays (Module)

Used for standard software packages managed by **CondaTainer** build scripts or from conda.

* **Format:** `name/version`
* **Structure:**
  * **name**: The software package or tool name (e.g., `bcftools`).
  * **version**: The specific version of the software (e.g., `1.22`).
* **Example:** `cellranger/9.0.1`

If you want a fixed combination of conda packages for everyday use (e.g. an editor), use `-n <name> <conda packages>` to bundle them into a single overlay:

```bash
condatainer create -n nvim nvim nodejs  # → nvim.sqf (Neovim + Node.js for plugins)
condatainer exec -o nvim nvim           # run it
```

### Reference Overlays (Module)

Used for reference datasets, genome assemblies, or indices.

* **Format:** `<assembly|project>/<datatype>/<version>`
* **Structure:**
  * **assembly**: The genome assembly or project (e.g., `grch38`).
  * **datatype**: The type of data (e.g., `gtf-gencode`).
  * **version**: The release or build version (e.g., `47`).
* **Example:** `grch38/gtf-gencode/47`

**Version delimiters:** the following are accepted when specifying a version — `/`, `--`, `=`, `@` (e.g. `bcftools/1.22`, `bcftools--1.22`, `bcftools=1.22`, `bcftools@1.22` all resolve to the same overlay).

```{note}
`--` is reserved for separating components in overlays. For example, `bcftools/1.22` overlay will be located at `images/bcftools--1.22.sqf`. Because `--` serves as a delimiter, it must not be used within names.
```

### Rename or not to Rename?

TL;DR: Do not rename `sqf` files. `img` files can be renamed.

**CondaTainer** relies on the overlay names to determine mount points and environment variable settings. Renaming `.sqf` files can lead to inconsistencies and runtime errors.

e.g.

- `bcftools/1.22` is mounted at `/cnt/bcftools/1.22`
- `grch38/gtf-gencode/47` is mounted at `/cnt/grch38/gtf-gencode/47`
- `my_analysis_env.sqf` is mounted at `/cnt/my_analysis_env`

For `.img` files, the mount point is always `/cnt_env`, so renaming is allowed.

## Mount Points

**CondaTainer** uses clear terminology for overlay types and purposes. Read-only overlays use the `.sqf` extension (SquashFS) and writable overlays use the `.img` extension (ext3).

**Overlay Terminology**

| Term | Extension | R/W | Content Path | Purpose |
|------|-----------|-----|--------------|---------|
| OS | `.sqf` | R/O | `/bin`, `/lib`, `/usr` | System Foundation — minimal system root that can run standalone or serve as a base for Modules/Bundles. |
| Module | `.sqf` | R/O | `/cnt/<name>/<version>` | Individual tool overlay (single package) mounted under `/cnt` and layered on top of an OS or Bundle. |
| Bundle | `.sqf` | R/O | `/cnt/<env_name>` | Frozen Conda environment (prebuilt collection of packages) mounted under `/cnt` as a named environment. |
| Environment | `.img` | R/W | `/cnt_env` | Writable Conda environment (ext3 image) for interactive work and runtime package changes. |

Read-only `.sqf` overlays are ideal for distributing immutable software and reference data. Writable `.img` overlays are for live development or when packages must be changed at runtime.

Mount points and examples

- Read-only overlays (`.sqf`) are mounted under `/cnt` following their naming convention.
  - Module example: `bcftools/1.22` → mounted at `/cnt/bcftools/1.22`.
  - Bundle example: `my_project_env.sqf` → mounted at `/cnt/my_project_env`.
  - OS overlays present system directories (e.g., binaries under `/bin`) inside the container filesystem exposed by the overlay.
- Writable images (`.img`) expose the Conda environment at `/cnt_env`, where packages may be installed or modified.

Writability: `.sqf` overlays are read-only. To enable write access for a `.img` overlay use the `-w` / `--writable` flag with `exec` or `run`.

## Overlay

Manage overlay writable ext3 overlay images:

- **create**: Create an overlay with a conda environment inside.
- **info**: Display disk usage and filesystem statistics.
- **check**: Verify filesystem integrity.
- **chown**: Change ownership of files inside the overlay.
- **resize**: Resize the overlay file.
- **export**: Export the recipe that produced an overlay (definition, build script, or conda env spec). See [Export](#export).

### Overlay Create

Create an ext3 `.img` with a conda environment inside. This is useful for creating writable conda environments.

**Usage:**

```
condatainer overlay create [OPTIONS] [path]
```

**Options:**

* `-s`, `--size [SIZE]`: Size of the overlay image (default: 10G). Supports units like `20G`, `2048M`.
* `-p`, `--profile [PROFILE]`: Overlay profile: `small`, `balanced`, or `large` (default: balanced). Aliases: `conda`/`python` = `small`; `data`/`genome` = `large`.
* `-f`, `--file [FILE]`: Initialize with a Conda environment file (.yml or .yaml).
* `--fakeroot`: Create image compatible with fakeroot (owned by root, must use with `--fakeroot` later).
* `-S`, `--sparse`: Create a sparse image file (no pre-allocation).
* `--no-tmp`: Create the overlay directly at the target path instead of staging at a local tmp directory first.
* `-- [packages...] [-c channel...]`: Initialize with inline conda packages instead of a YAML file. Mutually exclusive with `-f`. (`conda-forge` will be appended if missing)
* NAME: Name of the overlay image (`env.img` if not specified).

**Examples:**

```bash
# Create a 20GB overlay img with a custom name
condatainer overlay create -s 20g my_analysis_env

# Create with a specific profile (small/balanced/large)
condatainer overlay create -p large data_env.img

# o is a shortcut for 'overlay create'
condatainer o -f environment.yml project_env.img

# Create a fakeroot-compatible sparse overlay
condatainer o --fakeroot --sparse

# Initialize with inline conda packages (faster than writing a YAML)
condatainer o -- python=3.11 numpy pandas
condatainer o myenv.img -- python=3.11

# Specify channels with -c; conda-forge is always appended if not listed
condatainer o -- pytorch torchvision torchaudio pytorch-cuda=12.4 -c pytorch -c nvidia
```

When neither packages nor `--file` are supplied, overlay creation skips conda initialization. The
first `mm install` inside the writable container creates the project conda environment.

Inside the writable container, use `mm` to install, update, and remove conda packages. `mm` is a
bash shell function available in a bash session; in every other shell, run `condatainer env`
directly instead — same subcommands, same arguments.

```bash
# install more packages (uses channels saved during creation)
mm install numpy pandas
mm install pytorch-cuda=12.4 -c pytorch -c nvidia  # extra channels merged with saved ones

# manage channels (/cnt_env/.condarc)
mm channels list              # show configured channels
mm channels prepend pytorch  # move/add channel to highest priority
mm channels append bioconda  # move/add channel to lowest priority
mm channels remove nvidia    # remove a channel

# search for packages (uses saved channels)
mm search pytorch-cuda

# pin the package version
mm pin add numpy # after installation
mm pin remove numpy # to unpin
mm pin list # list pinned packages

# update packages
mm update --all

# remove packages
mm remove pandas

# list installed packages
mm list

# clean all conda cache
mm clean -ay

# export the env
mm export --no-builds > my_env.yaml
```

`mm install`/`update`/`remove` automatically re-run any `activate.d`/`deactivate.d` scripts a
package sets up (e.g. `JAVA_HOME`), so a bash session picks up the change without exiting and
re-entering the container. Run this manually after a plain `condatainer env install`/`update`/
`remove`, or after a `mm` call that failed partway through:

```bash
eval "$(mm reactivate)"               # bash
eval "$(condatainer env reactivate)"  # zsh
condatainer env reactivate --shell fish | source  # fish
```

`--shell` defaults to auto-detecting the shell you're running in; pass it explicitly if that's
wrong for your session. Fish only re-runs a hook shipped as `*_activate.fish`/
`*_deactivate.fish` — same as real conda's own fish support — so a package whose hook is
`*_activate.sh` (the common case) sets nothing under fish either way; exit and re-enter the
container to pick that up.

### Overlay Info

Alias for [`condatainer info`](#info) scoped to ext3 `.img` overlays. Produces identical output.

```
condatainer overlay info <overlay>
```

```bash
condatainer overlay info env.img
# same as:
condatainer info env.img
```

### Overlay Check

Verify filesystem integrity of an ext3 `.img` overlay.

**Usage:**

```
condatainer overlay check <overlay>
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
condatainer overlay chown [OPTIONS] <overlay>
```

**Options:**

* `-u`, `--uid UID`    : Set the owner UID (default: current user's UID).
* `-g`, `--gid GID`    : Set the group GID (default: current user's GID).
* `--root`             : Set UID/GID to 0 (root); overrides `-u` and `-g`.
* `-p`, `--path PATH`  : Path inside the overlay to change (repeatable, default: `/`).
* `path`               : Path to the ext3 overlay file.

**Examples:**

```bash
# Set ownership inside env.img to the current user (entire overlay by default)
condatainer overlay chown env.img

# Explicitly set UID and GID
condatainer overlay chown -u 1001 -g 1001 env.img

# Set ownership to root
condatainer overlay chown --root env.img

# Chown only /ext3 path
condatainer overlay chown -p /ext3 env.img

# Chown multiple specific paths
condatainer overlay chown -p /ext3 -p /data env.img
```

### Overlay Resize

Resize an ext3 `.img` overlay to a new size. Supports both expanding and shrinking.

**Usage:**

```
condatainer overlay resize -s SIZE <overlay>
```

**Options:**

* `-s`, `--size SIZE`  : New size for the overlay (e.g., `20g`, `2048M`). **Required.** Case insensitive
* `-S`, `--sparse`     : Leave the image sparse instead of pre-allocating its blocks.
* `path`               : Path to the ext3 overlay file.

Blocks are pre-allocated by default, as in `overlay create`, so a sparse image
cannot promise space the host filesystem is unable to supply. Allocation uses
`fallocate`; where that is unsupported (NFSv3, Lustre before 2.14) it warns and
leaves the image sparse.

**Examples:**

```bash
# Resize env.img to 20GB
condatainer overlay resize -s 20g env.img

# Grow without reserving blocks
condatainer overlay resize -s 20g env.img --sparse

# Pre-allocate an existing sparse image, size unchanged
condatainer overlay resize -s 20g env.img   # env.img is already 20g
```

### Overlay Export

Export the Conda environment in a writable `.img` overlay. See [Export](#export).

```bash
condatainer overlay export env.img > environment.yml
```

### Overlay Freeze

Pack a writable `.img` overlay into an immutable `.sqf` artifact — a snapshot. The environment is kept exactly as it is, not rebuilt: packages installed by hand and files deleted from the base are all preserved.

A frozen environment can be pinned, pushed and restored; a writable overlay cannot. It has no recipe behind it, so the `.sqf` is the only copy — back it up, or publish it with [`project push`](#project-push).

Loading a writable `.img` looks beside it for a matching frozen snapshot and mounts it automatically underneath — a personal `<name>-<user>.sqf` line is checked before a shared `<name>.sqf` line, and whichever is found becomes the read-only base the `.img` builds on.

With no destination, `overlay freeze` is the routine snapshot loop: the target is whichever of those two lines the `.img` was actually loaded against, and freezing replaces it — the source `.img` is then removed, so the next `overlay create` of the same name starts fresh and thin, autoloading the new snapshot again. `--keep`, or giving an explicit destination, keeps the `.img` instead — that's for a deliberate artifact meant to be found and used elsewhere, not the routine loop.

**Usage:**

```
condatainer overlay freeze [OPTIONS] <overlay.img> [artifact.sqf]
```

**Options:**

* `-d`, `--description TEXT` : Description recorded in the artifact.
* `--block-size SIZE`        : SquashFS block size (default: `build.block_size`).
* `--use-tmp`                : Copy the payload to the temp directory and pack from there. Faster; needs payload-sized free space.
* `--keep`                   : Keep the source `.img` after a bare freeze (default: remove it).
* `--<compression>`          : Any `create` compression flag, e.g. `--zstd-high` (default: `build.compress_args`).
* `overlay.img`              : The writable overlay to pack.
* `artifact.sqf`             : Where to write it. A directory puts `<overlay>.sqf` inside it. Omitted, the target is resolved as described above.

The artifact's type is `environment`; `condatainer info` reads back what it holds. Deletions are carried into the artifact and stay deleted when it is mounted, and freeze reports how many it translated.

**Refused:**

* an explicit destination that already exists (a bare freeze may replace its resolved target instead — see above);
* a bare freeze whose resolved target is occupied by something that isn't a snapshot — move it aside or give an explicit destination;
* a target inside an images directory — a frozen environment is addressed by path, not filed by name;
* an overlay with nothing written into it;
* an overlay a writable session is using, or a target another session is actively reading. A read-only (`chmod a-w`) overlay, or a `.sqf` protected the same way, still freezes and is never replaced.

**Examples:**

```bash
# Pack env.img into env.sqf beside it, then remove env.img
condatainer overlay freeze env.img

# Same pack, but keep env.img
condatainer overlay freeze env.img --keep

# Into a project's overlay directory, with a description (source kept)
condatainer overlay freeze env.img ./overlays/ -d "paper revision 2"

# Pack harder than the build default
condatainer overlay freeze env.img --zstd-high
```

To install with `apt` before freezing, give the overlay to root, install under `--fakeroot`, then take it back:

```bash
condatainer overlay chown env.img --root
condatainer exec --fakeroot -o env.img -- apt install <package>
condatainer overlay freeze env.img
condatainer overlay chown env.img
```

Freeze itself does not need `--fakeroot`.

### Overlay Unfreeze

Rebuild a writable `.img` overlay from a frozen environment, to keep developing in it. The payload comes back at the paths it had, deletions included.

The result has no identity: it cannot be pinned, locked or published until `overlay freeze` packs it again. The `.sqf` is left as it is, and a project whose lock names it still resolves that artifact.

**Usage:**

```
condatainer overlay unfreeze [OPTIONS] <frozen.sqf> [overlay.img]
```

**Options:**

* `-s`, `--size SIZE` : Size of the resulting overlay (e.g. `40g`). Default: the payload plus 5 GB.
* `-S`, `--sparse`    : Leave the image sparse instead of pre-allocating its blocks.
* `frozen.sqf`        : The artifact `overlay freeze` produced.
* `overlay.img`       : Where to write the overlay. A directory puts `<artifact>.img` inside it. Omitted, it lands beside the artifact with the extension changed.

An `<overlay>.img.env` sidecar records which artifact the overlay came from.

**Refused:**

* a target that already exists;
* a `.sqf` that is not a frozen environment;
* a `--size` smaller than the payload, before any work is done.

**Examples:**

```bash
# Rebuild env.img beside the artifact
condatainer overlay unfreeze env.sqf

# To a separate path, with room to install more
condatainer overlay unfreeze env.sqf dev.img -s 40G
```

## Create

Initialize and build a new **CondaTainer** SquashFS overlay. You can build from existing recipes (local/remote), a Conda environment file, or a remote container source.

**Usage:**

```
condatainer create [OPTIONS] [packages...]
```

**Aliases:** `install`, `i`

Named overlays are written to the first writable [data layer](../deployment/data_layers.md), which `create` reports before building:

```
[CNT◇] Installing to /opt/condatainer/images (app-root)
```

Check the layer if it matters who can see the result — `(app-root)` or `(extra-root)` is shared, `(user)` is yours alone. A shared directory you cannot write to is skipped silently, so an install meant for everyone can land in your own directory instead. Use `-p` to write to an exact path.

**Arguments:**

* `packages`: List of packages to install (e.g., `bcftools/1.22` or `samtools=1.10` or `grch38/genome/gencode`). Supports conda channel annotations: `bioconda::star=2.7.11b` (version required in default mode; optional with `-n`/`-p`).

**Flags:**

* `-n`, `--name [NAME]`: Custom name for the resulting overlay file. If used, all specified packages are bundled into one overlay.
* `-p`, `--prefix [PATH]`: Custom prefix path for the overlay file. When `-f` is used, this can be omitted — the prefix is inferred from the file name.
* `-f`, `--file [FILE]`: Path to definition file (.yaml, .sh, .def), an already-built `.sif`, or an Apptainer sandbox directory.
* `--from [URI]`: Build from an external image URI (e.g., `docker://ubuntu:22.04`).
* `-c`, `--channel [CHANNEL]`: Conda channel to use, overriding config channels. Repeatable: `-c conda-forge -c bioconda`.
* `-u`, `--update`: Rebuild overlays even if they already exist (atomic `.new` swap). Useful for refreshing a package to the latest version.
* `-l`, `--layer [LAYER]`: Build into a chosen [data layer](../deployment/data_layers.md) — `u`/`user`, `r`/`app-root`, `e`/`extra-root` — instead of the first writable one. Cannot be combined with `-p`, which already names where the overlay goes.

**Build Flags:**

* `--block-size [SIZE]`: SquashFS block size for app/env/external overlays (e.g. `128k`, `512k`; default: `128k`). Must be a power of two between `4k` and `1m`.
* `--data-block-size [SIZE]`: SquashFS block size for data/reference overlays (e.g. `512k`, `1m`; default: `512k`). Must be a power of two between `4k` and `1m`.
* `--no-prebuilt`: Build from the recipe even when a registry publishes a prebuilt artifact for it. Without it, a matching prebuilt is pulled instead of built.
* `--store`: Build into the [store](#store), filed under this build's identity instead of taking the plain name.

  It never skips and never replaces. Ordinarily a build of an already-installed name is skipped, and `-u` swaps the installed one out; `--store` is how you get a second build of that name installed alongside the first. It behaves the same when nothing holds the name yet — the build is filed by identity either way, so it stays out of `condatainer list` and out of `exec -o <name>` until you promote it with [`store use`](#store). That is the point of the flag: a build that changes nothing anyone else resolves.

  Before building, it works out the identity the build *would* produce and stops if that exact artifact is already installed — a recipe build knows its identity from the recipe and its dependencies, and a Conda build learns it from a solve (`--dry-run`), without creating the environment. So re-running `--store` after a successful one costs a solve, not a build.

  Works with any build that lands in an images directory — a `name/version` recipe, `-n` with packages, `-n -f environment.yml`, `-n --from docker://…`. The one conflict is `-p`/`--prefix`, which names an exact output file while the store generates its filename from the artifact's keys. A bare `-f environment.yml` with no `-n` derives a prefix from the file name, so it conflicts too; give it a `-n`.
* `--always-submit-data`: Submit `data` builds as scheduler jobs even when the recipe has no scheduler directives. Other builds run here unless their recipe carries directives. Conda builds, definitions and image imports never go to the scheduler.
* `--no-submit`: Disable job submission; build locally even if the build script has scheduler directives.
* `-s`, `--source`: Use only this configured recipe collection; repeat to set the lookup order.

When a selected recipe source declares `oci.pull` endpoints, `create` tries them
in order before building locally. A candidate is accepted only when its
equivalence key matches the selected recipe. Missing artifacts, missing host
architecture, or an unavailable registry fall back to the local build;
authentication, compatibility, and integrity failures are reported.

**Compression Options:**

* `--lz4`: Use LZ4 compression (default).
* `--zstd-fast`: Use Zstandard compression level 3.
* `--zstd-medium`: Use Zstandard compression level 8.
* `--zstd`: Use Zstandard compression level 14.
* `--zstd-high`: Use Zstandard compression level 19.
* `--gzip`: Use Gzip compression.

**Build Modes:**

* **Default:** Each package gets its own `.sqf` via the build system. Build script lookup is skipped when a package uses channel annotation (`bioconda::star=2.7.11b`).
* **`--name`:** Create a single `.sqf` with multiple packages bundled together.
* **`--name` + `--file`:** Create `.sqf` from a source file under that name, in the managed images directory.
* **`--prefix` + packages:** Create a conda env `.sqf` at a custom path, like `conda create -p`.
* **`--file` only:** Create `.sqf` from external source file; prefix inferred from file name (e.g. `condatainer create -f r-collect.sh` → `r-collect.sqf`). The name is inferred only when neither `--name` nor `--prefix` is given.
* **`--prefix` + `--file`:** Create `.sqf` from external source file at a custom path.
* **`--from`:** Create `.sqf` from an external container image URI.
* **`--file` with a `.sif` or sandbox directory:** Import an already-built Apptainer/Singularity
  root into a real `.sqf` — no rebuild. Only a root originally built from `docker://`, `oras://`, or
  `library://` can be imported; anything else is refused with a clear message. A `.sif` built for a
  different architecture than this host is refused rather than packed unusable.

### System level Examples

**Conda environment:**

```bash
# Single package module overlay (conda fallback)
condatainer create bcftools/1.22

# Multiple package module overlays at once
condatainer create samtools/1.16 bcftools/1.15

# Create a conda env sqf with packages at a custom path (like conda create -p)
condatainer create python=3.11 numpy -p /scratch/myenv

# Create from a yaml file (prefix inferred from filename)
condatainer create -f environment.yml
# Create from a yaml file with a custom prefix
condatainer create -p my_analysis_env -f environment.yml

# Rebuild an existing conda overlay (force update)
condatainer create -n nvim nvim nodejs -u

# Override channels for this build only
condatainer create python=3.11 -n myenv -c conda-forge -c bioconda

# Channel-annotated package: skips build script lookup, goes direct to conda
# (version required in default mode)
condatainer create bioconda::star=2.7.11b

# Channel annotation in bundled env
condatainer create -p rnaseq bioconda::star bioconda::salmon=1.10.0
```

**Build script:**

```bash
# App overlay via build script
condatainer create cellranger/9.0.1

# Bare package name is expanded using the configured default_distro (e.g.
# ubuntu24/build-essential), or a standing project's pinned root instead
condatainer create build-essential

# Data overlay
condatainer create grch38/gtf-gencode/47

# From a remote container source
condatainer create --from docker://ubuntu:22.04 -n myubuntu

# Template script — prompts for each placeholder interactively
condatainer create grch38/salmon-gencode
# [CNT] Placeholder template: grch38/salmon-gencode
# [CNT] Salmon GRCh38 GENCODE{gencode_version} index for transcript quantification
#   Target: grch38/salmon/{salmon_version}/gencode{gencode_version}
#   salmon_version [1.0.0-1.11.4] (default: 1.11.4): 1.10.2
#   gencode_version [23-49] (default: 49):
#   → Creating grch38/salmon/1.10.2/gencode49

# Template script — skip prompts by specifying the resolved target name directly
condatainer create grch38/salmon/1.10.2/gencode49
```

**Features**:

- Automatic Fetching: If a build script is not found locally, **CondaTainer** attempts to fetch it from the remote repository.
- Conda Fallback: If no build script exists, **CondaTainer** attempts to create the module by installing the package with the requested name and version from conda-forge or bioconda.
- Metadata Parsing: Parses `#ENV:` tags (with their inline `## ` notes) from recipes to inject environment variables and help text into the generated modulefile.
- Template Resolution: If the name matches a template script (`#PH:` / `#TARGET:`), **CondaTainer** prompts for each placeholder interactively, then builds the resolved concrete overlay. You can also bypass prompts by specifying the resolved target name directly.

### Project level Examples

Create a read-only overlay with a Conda environment using a Conda YAML file. The prefix can be omitted — it is inferred from the file name.

```bash
# Prefix inferred: generates environment.sqf in the current directory
condatainer create -f environment.yml

# Custom prefix: generates my_project_env.sqf
condatainer create -p my_project_env -f environment.yml
```

Create a read-only overlay using a custom Apptainer definition file. See [Custom OS Overlays](../advanced_usage/custom_os.md) for more details.

```bash
condatainer create -f custom_def.def
```

Create a read-only overlay using a shell script that installs packages. See [Custom Build Script using Build Scripts](../advanced_usage/custom_bundle.md) for more details.

```bash
condatainer create -f install_packages.sh
```

Import an already-built `.sif` or Apptainer sandbox directory into a real `.sqf`, instead of
rebuilding it from source:

```bash
condatainer create -f alpine.sif -p alpine
condatainer create -f ./my-sandbox-dir -p my_base
```

### Exit Codes (script and job-submission behavior)

CondaTainer uses specific exit codes so automation and downstream tooling can detect special states:

- `0` — Success (all requested builds completed locally or nothing to do)
- `1` — Generic error (build failures, runtime errors, invalid arguments, or other fatal errors)
- `3` — **Jobs submitted to scheduler** — overlays will be created asynchronously by scheduler jobs

Commands that may return exit code `3` when scheduler jobs were submitted include:

- `condatainer create ...`
- `condatainer check -a ...` (auto-install missing deps)

Quick example for shell scripts that detect the job-submitted state:

```bash
condatainer create grch38/salmon/1.10.2/gencode49
if [ $? -eq 3 ]; then
  echo "Jobs submitted to scheduler — overlays will be created asynchronously"
  # Optionally: exit 0 or wait/monitor jobs here
fi
```

Note: When exit code `3` is returned, CondaTainer prints a message showing the number of jobs submitted.

**Disabling Job Submission:**

```bash
# Run locally with scheduler specs
condatainer run --no-submit analysis.sh

# Or set in config
condatainer config set submit_job false
```

## Container Management (Avail, List, Remove, Search)

Manage your local library of built containers and available recipes.

### Avail

Search the recipes of every configured [recipe collection](configuration.md#recipe-sources).

**Aliases:** `av`

```
condatainer avail [search_terms...] [flags]
```

**Options:**

* `-s`, `--source`: Search only this configured collection; repeat to set the lookup order.
* `-e`, `--expand`: Expand template groups to show individual concrete entries instead of the collapsed template header.
* `--description`: Show the description (`#DESC:`) for each entry.

**Search rules:**

| Input | Mode |
|---|---|
| Single plain string | Substring match |
| Single term with `*` / `?` | Wildcard (e.g. `cell*`) |
| Single term with `^` `$` `(` `[` `+` `{` `\|` | Regex (e.g. `^cell.*9\.0`) |
| Multiple terms, first is an exact name | Each term matched exactly (OR) |
| Multiple terms, first not found | All terms AND substring |

**Template display:**

Template scripts (`#PH:` / `#TARGET:`) are shown collapsed by default as a group header with variant count and placeholder value summaries. Use `-e` to expand all concrete combinations instead.

```
grcm39/salmon-gencode  [594 variants]
  Salmon GRCm39 GENCODEM{gencode_version} index for transcript quantification
  → grcm39/salmon/{salmon_version}/gencodeM{gencode_version}
  - salmon_version:   1.0.0-1.11.4  (18 values)
  - gencode_version:  6-38  (33 values)
```

Placeholder lists of more than 5 concrete values collapse to a `first-last  (n values)` range; shorter lists are shown in full. `*` marks an open-ended placeholder that accepts any value.

**What a search term matches:**

Terms match entry **names**, plus **descriptions** while those are shown. You can search whatever is displayed:

| Invocation | Descriptions shown | Searched |
|---|---|---|
| `avail <term>` | yes | name + description |
| `avail -e <term>` | no | name only |
| `avail -e --description <term>` | yes | name + description |
| `avail --description=false <term>` | no | name only |

Without `-e`, only templates and plain entries are searched, so templates stay collapsed no matter how many variants they have:

```bash
$ condatainer avail grch star gencode
grch38/star-gencode  [1176 variants]
  → grch38/star/{star_version}/gencode{gencode_version}-{read_length}
  - star_version:     2.7.0b-2.7.11b  (21 values)
  ...
```

Read the values you need from the header, then pass the assembled name to `install`.

Use `-e` to match variant names, including version strings:

```bash
$ condatainer avail gencode47-101
[CNT!] No matching build scripts found.
[CNT◇] Templates are listed collapsed — use -e to search individual variants.

$ condatainer avail -e gencode47-101
grch38/star/2.7.0b/gencode47-101
grch38/star/2.7.0d/gencode47-101
...
```

Under `-e`, a template matched by name contributes all of its variants; otherwise each variant must match on its own.

**Installing what you find:**

Pass any name from the listing to `condatainer install`:

- **Concrete entry** (e.g. `grch38/genome/gencode`): built directly.
- **Template entry** (e.g. `grch38/star-gencode`): prompts for each placeholder, then builds the resolved name. Defaults prefer the latest already-installed version of each dependency, falling back to the latest available.
- **Filled-in template name** (e.g. `grch38/star/2.7.11b/gencode47-101`): already concrete, built with no prompting.

Values for open-ended (`*`) placeholders cannot be enumerated, so they never appear in `avail -e` output — but `install` still accepts them, e.g. `condatainer install grch38/star/2.7.11b/gencode47-75`.

**Examples:**

```bash
# Substring search
$ condatainer avail cellranger

# Wildcard
$ condatainer avail 'cell*'

# AND search with multiple terms
$ condatainer avail cellranger 9

# Install one of the results
$ condatainer install cellranger/9.0.1

# Search variants for all resources matching 'grcm' and 'M33' (-e required:
# 'M33' spans the literal 'M' and the placeholder value '33')
$ condatainer avail -e grcm M33
grcm39/gtf-gencode/M33
grcm39/salmon/1.0.0/gencodeM33
grcm39/salmon/1.1.0/gencodeM33
...

# Search descriptions too (shown by default, so this needs no flag)
$ condatainer avail java

# Hide descriptions (and stop matching against them)
$ condatainer avail star --description=false

# Expand all template combinations
$ condatainer avail star -e
```

### List

List installed overlays stored in the `images/` directory.

**Aliases:** `ls`

```
condatainer list [search_terms...] [flags]
```

**Options:**

* `-d`, `--delete`: Prompt to delete listed overlays after displaying them (requires search terms).
* `-r`, `--remove`: Alias for `--delete`.
* `-l`, `--layer`: Limit to one data layer — `u`/`user`, `r`/`app-root`, `e`/`extra-root`.
* `-D`, `--dir`: Limit listing (and deletion) to specific image directories.

Matching rules: 
* no leading `/` → substring match (`scratch` matches `/scratch/user/images`)
* leading `/` → exact match
* leading `/` with `*`/`?` → wildcard `*` (`/scratch/*` matches `/scratch/user/images`).

**Search rules:**

Same as `avail` (single term: substring/wildcard/regex; multiple terms: exact-first or AND substring). Distro-prefix aliases are recognised for exact multi-term matching (e.g. `rstudio-server` matches `ubuntu24/rstudio-server`).

**Features:**

* Output is grouped by image directory with a full-width header per directory, tagged with its data layer — `(user)`, `(extra-root)`, or `(app-root)`. Directories appear in read order, nearest first. The same overlay may appear under several directories; the nearest copy is the one commands resolve by name, so a personal rebuild wins over a shared one.
* Lists OS overlays, app overlays, and data overlays.
* An overlay is listed under the name recorded inside it. If that is not the name its filename spells, the listing says so: the recorded name is what the image is, and the filename is what `exec -o` resolves it by. Rename the file or reinstall to make them agree.
* Missing directories are skipped; existing but empty directories show a `(no overlays)` line.
* Exits with code `1` if search terms are given but no overlays match.

**Delete mode (`-d`/`-r`):**

* Requires search terms — `list -r` alone only lists.
* If the same overlay name exists in multiple directories, `--layer` (or `--dir`) is required to avoid ambiguity.
* Checks file lock and write permission before each deletion.

**Examples:**

```bash
$ condatainer list                            # all overlays, grouped by dir
$ condatainer list cellranger                 # substring match
$ condatainer list 'cell*'                    # wildcard
$ condatainer list cellranger 9               # AND search
$ condatainer list --dir scratch              # only dirs with "scratch" in path
$ condatainer list --dir /scratch/user/images # exact dir match
$ condatainer list --dir '/scratch/*'         # wildcard: all dirs under /scratch
$ condatainer list cellranger/9.0.1 -d        # show then prompt to delete
$ condatainer list cellranger -r --dir /scratch  # delete from /scratch only
```

### Remove

Deletes specific overlays and their associated `.env` files.

**Aliases:** `rm`, `delete`, `uninstall`

```
condatainer remove [search_terms... | file_paths...]
```

**Options:**

* `-l`, `--layer`: Limit to one data layer — `u`/`user`, `r`/`app-root`, `e`/`extra-root`.
* `-D`, `--dir`: Limit to specific image directories (search mode only). No leading `/` → substring match; leading `/` → exact match; leading `/` with `*`/`?` → wildcard.

Both are optional and compose. When an overlay of the same name exists in more than one layer, `remove` refuses and lists the copies with their layers; pick one with `-l` (or `-D`):

```
[CNT!] The following overlays exist in multiple layers. Use --layer (or --dir) to pick one:
  samtools/1.22
    /opt/condatainer/images (app-root)
    /scratch/me/condatainer/images (user)

$ condatainer remove -l u samtools/1.22     # removes your copy only
```

**Modes:**

**Search mode** — args have no overlay extension:

| Input | Mode |
|---|---|
| Single plain string | Exact match (including distro-prefix alias) |
| Single term with `*` / `?` | Wildcard |
| Single term with `^` `$` `(` `[` `+` `{` `\|` | Regex |
| Multiple terms, first is an exact name | Each term matched exactly (OR), warns if not found |
| Multiple terms, first not found | All terms AND substring |

**File mode** — args end with `.img`, `.sqf`, `.sqsh`, or `.squashfs`:

Removes the specified files directly, bypassing the name-based search. Useful for overlays that are not in any configured image directory.

**Features:**

* Overlays to be removed are displayed grouped by directory before confirmation.
* If the same overlay name exists in multiple directories, `--layer` (or `--dir`) is required to avoid ambiguity.
* Checks file lock and write permission before each deletion, and says which one refused: an image
  in use by a running container, an image that is write-protected, or one that is missing.
* Also removes the associated `.env` file if present.

**Pinning an image.** Clearing the write bit protects an image from `remove` and from
`build --update`, including for its owner, and it stays readable and mountable while pinned:

```bash
$ chmod a-w ~/.local/share/condatainer/images/grch38--genome--gencode.sqf
$ condatainer rm grch38/genome/gencode
[CNT!] Cannot remove grch38/genome/gencode: ... is write-protected; chmod +w to allow changes
```

**Examples:**

```bash
$ condatainer rm cellranger/9.0.1                    # exact version
$ condatainer rm cellranger                          # all cellranger versions (exact)
$ condatainer rm 'cell*'                             # wildcard
$ condatainer rm cellranger/9.0.1 cellranger/8.0.1   # multiple exact versions
$ condatainer rm cellranger 9                        # AND search
$ condatainer rm cellranger --dir scratch            # only from dirs with "scratch" in path
$ condatainer rm cellranger --dir '/scratch/*'       # wildcard: all dirs under /scratch
$ condatainer rm /path/to/cellranger--9.0.1.sqf      # direct file path
$ condatainer rm *.img                               # all .img files in current dir
# path and name/version cannot be mixed
```

### Search

Search conda packages via the anaconda.org API.

```
condatainer search <package> [flags]
```

**Options:**

* `--json`: Output results in JSON format.
* `-f`, `--fuzzy`: Search package via api.anaconda.org.
* `-l`, `--limit N`: Number of fuzzy search results before filter (default: 100).
* `-c`, `--channel [CHANNEL]`: Channel to search, overriding config channels. Repeatable: `-c bioconda -c conda-forge`.

**Notes:**

* Uses the anaconda.org REST API — no base image or micromamba required.
* **Exact match** (default): queries each configured channel's package API in order and returns the first hit, matching install behaviour.
* **Fuzzy match** (`-f`): queries the anaconda.org search API with two requests — one for the current platform and one for `noarch` — then merges and sorts results alphabetically. Results are capped at `--limit`; a warning is shown if the limit was hit.
* Uses channels from the config (`build.channels`, default: `conda-forge`, `bioconda`). `-c` fully replaces the config channel list for this invocation.

**Examples:**

```bash
$ condatainer search samtools
$ condatainer search samtools --json
$ condatainer search -f samtool
$ condatainer search -f samtool -l 200
$ condatainer search star -c bioconda
```

## Exec

Execute a command inside a containerized environment using explicit overlay specifications.

**Usage:**

```
condatainer exec [flags] [command...]
```

**Options:**

* `-o`, `--overlay [OVERLAY]`: Overlay file to mount (repeatable).
* `-w`, `--writable`: Mount `.img` overlays as writable (default: read-only).
* `-f`, `--fakeroot`: Run container with fakeroot privileges.
* `--env [KEY=VALUE]`: Set environment variable inside the container (repeatable).
* `--bind [HOST:CONTAINER]`: Bind mount path into the container (repeatable).
* `--gpu`: Force GPU flags (`--nv`/`--rocm`) even if `autoload_gpu` is disabled.
* `--activation [all|env|none]`: Which activate.d scripts to source before the
  command (default `all`). `all` sources every mounted app overlay's own plus the
  mounted conda environment's; `env` sources only the conda environment's; `none`
  sources neither.
* `--project [DIR]`: Act as if standing in `DIR` instead of the current directory —
  every relative argument, and project lookup, resolves against it.
* `--no-project`: Ignore any project found above the current directory, even one
  standing inside it would otherwise resolve. Cannot be combined with `--project`.

**Features:**

* Use `-o/--overlay` to explicitly specify overlays.
* All positional arguments are treated as commands.
* Read-only by default for `.img` overlays (use `-w` for writable).
* Defaults to bash if no command specified.
* There is no separate base-image flag: the container root is a `-o` overlay
  that is a `.sif` if one is given (only one `.sif` is allowed per command),
  otherwise the first `-o` overlay that is itself an OS layer, if any;
  otherwise, inside a project, the root that project has pinned; otherwise
  the configured default (built automatically if missing).
* Inside a project — a directory holding `cnt-lock/`, or any directory below
  one — every `-o` name resolves through that project's lock instead of by
  installed name, and so does the root when nothing else supplied one. See
  below.

**Inside a project:**

When the current directory holds `cnt-lock/`, or sits anywhere below a
directory that does, `-o` resolves against that project's lock, so a run
mounts the exact artifact the project pinned rather than whatever currently
answers to the name. `e`, `run` and `check` behave the same way, and all four
print `Project: <root>` once, the first time anything in the run actually
resolves through it — including when nothing was passed with `-o` at all and
only the container root came from the project.

* A name the lock does not pin is an error naming the remedy, not a
  fallback to the installed copy.
* An artifact that is pinned but not present here is an error naming
  `condatainer project restore`. Nothing is fetched or built to satisfy a mount.
* A version constraint such as `-o star/2.7.11b>=2.7.0` is refused: the lock says
  which version, so a range is not a request.
* A project-relative `.sqf` is checked against the lock before it is mounted.
* A writable `.img` has no identity to pin and is mounted as written, relative to
  the project root.
* **The root is no different.** When nothing you passed is itself an OS layer,
  the project's pinned root is used — not the configured default_distro — so
  `exec`, `e` and `run` all put you in the same container regardless of how
  your own machine is configured. If that root is not installed here, it is
  the same error naming `condatainer project restore`, not a silent fall back
  to your own configuration.

Standing anywhere outside a project's tree is the ordinary opt-out — no flag
required. `--project DIR` runs the command as if you had `cd`'d to `DIR`
first, and `--no-project` disables project lookup outright even when standing
inside one, for the rare case where you need the ambient, by-installed-name
behavior anyway. These are distinct from `condatainer project`'s own
`--project DIR`, which names a project to administer (lock, pin, restore,
push) rather than a directory to relocate into.

**Environment Variables (inside container):**

* `IN_CONDATAINER=1`: Set inside the container.
* `CNT_CONDA_ROOT`: Path to the mounted project conda environment.
* `CNT_CONDA_WRITABLE`: `1` when the `.img` is writable, otherwise `0`.

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
* `-f`, `--fakeroot`: Run container with fakeroot privileges.
* `--env [KEY=VALUE]`: Set environment variable inside the container (repeatable).
* `--bind [HOST:CONTAINER]`: Bind mount path into the container (repeatable).
* `--activation [all|env|none]`: same as [`exec`](#exec).
* `--project [DIR]`, `--no-project`: same as [`exec`](#exec).

**Key Differences from `exec`:**

* Overlays are positional arguments before `--`.
* Commands go after `--`.
* Writable by default (use `-r` for read-only).
* Auto-loads `env.img` unless `-n` is specified.
  * looks in the current directory only; `env-$USER.img` takes priority over `env.img`
  * if neither `.img` exists but its frozen snapshot does (`env-$USER.sqf` or `env.sqf`), that is loaded instead, read-only
* Defaults to bash if no command specified.

**Environment Variables (inside container):**

* `IN_CONDATAINER=1`: Set inside the container.
* `CNT_CONDA_ROOT`: Path to the mounted project conda environment.
* `CNT_CONDA_WRITABLE`: `1` when the `.img` is writable, otherwise `0`.

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

## Runtime (Check, Run)

Utilities for running scripts with automatic dependency handling via `#DEP:` tags.

### Check

Parses scripts for `#DEP:` tags and checks if the required overlays are installed.

**Usage:**

```
condatainer check <script|dir|name> [script|dir|name ...] [-a]
```

Each argument can be:
- A local `.sh` file path
- A directory — all `.sh` files **directly inside** it are checked (non-recursive)
- A package name (e.g., `grch38/salmon/1.10.2/gencode49`) resolved from the configured recipe collections

**Options:**

* `-a`, `--auto-install`: Automatically attempt to build/install missing dependencies.
* `-i`, `--install`: Alias for `--auto-install`.
* `--no-submit`: Disable job submission; build missing dependencies locally.
* `--project [DIR]`, `--no-project`: same as [`exec`](#exec).

**Inside a project:**

Run from a directory holding `cnt-lock/`, or anywhere below one, `check <script>` 
resolves that script's declarations through the lock instead of by
name, and reports whether the script can run right now. It goes through the
same code `run` does, so the two cannot disagree — including opening a
project-relative `.sqf` to confirm it still carries the pinned identity. `-a`
is refused there: installing by name is what a lock exists to prevent, so the
answer is [`condatainer project restore`](#project-restore).

**Output:**

Dependencies are grouped into two sections:

```
Module Overlays:
  ✓ samtools/1.22
  ✗ bcftools/1.21

External Overlays:
  ✓ src/overlay/r-env.sqf
  ✗ src/overlay/r-collect.sqf  (.yml found)
  ✗ src/overlay/custom.sqf     (.sh found, but has #DEP - create manually)
```

**Auto-install (`-a`):**

- **Module overlays**: resolved via build scripts or conda, same as `condatainer create`.
- **External overlays**: auto-created from the sibling source file if found and has no `#DEP:` tags.
  - `.yml` / `.yaml` → conda environment
  - `.def` → Apptainer build and extract the squashfs
  - `.sh` (no `#DEP:`) → shell build (cannot have `#DEP:` tags)

### Run

Executes a script inside the **CondaTainer** environment, mounting dependencies defined in the script. Autosolves dependencies based on `#DEP:` tags within the script.

The `#DEP:` tags read here are **your script's** — they name the overlays to mount for this run. An overlay does not carry dependencies of its own: a recipe's `#DEP:` is a build-time edge and is never re-expanded at run time.

A `#DEP:` may be a full `name/version`, a partial version (`samtools/1.22`) or a bare name (`samtools`); a bare or partial name uses the newest matching version already installed. A version range (`samtools>=1.20`) is for build recipes only and aborts the run.

```
condatainer run [OPTIONS] SCRIPT [SCRIPT_ARGS...]
```

```{note}
All options (`-a`, `-o`, `--afterok`, etc.) must appear **before** `SCRIPT`. Arguments after the script name are forwarded to the script.
```

**Inside a project:** run from a directory holding `cnt-lock/`, or anywhere below
one, `run` resolves the script's `#DEP:` and its container root through that
project's lock, the same way `exec`/`e` do — see [Exec, Inside a project](#exec)
for the full rule and what an unresolved pin does.

**Container Flags:**

* `-w`, `--writable`, `--writable-img`: Make `.img` overlays writable (default: read-only).
* `-f`, `--fakeroot`: Run with fakeroot privileges.
* `--project [DIR]`, `--no-project`: same as [`exec`](#exec).
* `--bind HOST:CONTAINER`: Bind mount a path into the container (repeatable).
* `--env KEY=VALUE`: Set an environment variable inside the container (repeatable).

**Resource Override Flags:**

These override the script's scheduler directives (`#SBATCH`, `#PBS`, `#BSUB`) for this run.

* `-c`, `--cpu INT`: Override CPUs per task (e.g. `4`).
* `-m`, `--mem STRING`: Override memory per task (e.g. `4G`, `8192M`).
* `-t`, `--time STRING`: Override walltime (e.g. `4d12h`, `2h30m`, `01:30:00`).
* `-g`, `--gpu SPEC`: Override GPUs per node. Formats: `N` (any type), `TYPE:N`, or `TYPE` (count=1). E.g. `1`, `a100:2`, `a100`. When running (whether submitted or local), a GPU request here also forces `--nv`/`--rocm` on the container even if `autoload_gpu` is disabled.

**Job Flags:**

* `-n`, `--name NAME`: Override the job name shown in the scheduler queue (e.g. `squeue`). Takes priority over `#SBATCH --job-name` and similar directives.
* `-o`, `--output PATH`: Override the job stdout path (creates parent directory if needed). Takes priority over scheduler stdout settings.
* `-e`, `--error PATH`: Override the job stderr path. Takes priority over scheduler stderr settings.
* `-A`, `--account STRING`: Override the billing/allocation account. Falls back to the script's own directive, then `scheduler.account` in config, when unset.
* `-p`, `--partition STRING`: Override the partition/queue. Falls back to the script's own directive, then `scheduler.partition` in config, when unset.
* `--afterok IDS`: Submit job that runs only if all listed jobs **succeed**. Colon-separated IDs: `123:456:789`.
* `--afternotok IDS`: Submit job that runs only if any listed job **fails**. Colon-separated IDs.
* `--afterany IDS`: Submit job that runs after all listed jobs finish **regardless of outcome**. Colon-separated IDs.
* `--array FILE`: Input file for an array job — one subjob per line, tokens become positional args.
* `--array-limit N`: Max concurrently running subjobs (0 = unlimited).
* `--dry-run`: Preview what would be submitted without executing anything, including the base, dependencies, and how nested running would get apptainer (`nested_run`).
* `--no-submit`: Disable job submission; run the script locally even if it has scheduler directives.

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

| Condition | Behavior |
|-----------|----------|
| Already inside a running job or container | Always runs locally (no nested submission) |
| `--no-submit` flag or `submit_job: false` in config | Always runs locally |
| Script has no scheduler specs | Runs locally (prints a note) |
| Script has `#SBATCH`/`#PBS`/`#BSUB` + scheduler available | Submits as a scheduler job |
| Script has scheduler specs but scheduler not found/available | Runs locally (prints a note) |

**Scheduler Script Example:**

```bash
#!/bin/bash
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --time=2:00:00
#DEP: samtools/1.22
#DEP: bcftools/1.22

samtools view -@ $NCPUS input.bam | bcftools call -mv -o output.vcf
```

Install missing dependencies first:

```bash
condatainer check -a salmon_quant.sh
```

After all dependencies are installed, submit as a job:

```bash
condatainer run salmon_quant.sh
```

**Log Files:**

* `-o`/`-e` CLI flags take highest priority — they override any `#SBATCH --output`/`--error` in the script
* If no CLI flag, the path from the script directive (e.g., `#SBATCH --output=...`) is used
* Otherwise, logs are written to the global logs directory (`~/logs` by default)
* If only `-o` is set (no `-e`), stderr is merged into the same file
* Job scripts are created in the same directory as the log file

**Resource Override Examples:**

```bash
# Override CPU/memory/time on top of script's #SBATCH directives
condatainer run -c 8 -m 32G -t 4h analysis.sh

# Override GPU spec (any GPU type, count 2)
condatainer run -g 2 gpu_job.sh

# Override with specific GPU type and count
condatainer run -g a100:4 gpu_job.sh

# GPU type only (count defaults to 1)
condatainer run -g h100 gpu_job.sh
```

### Script Arguments

Anything after the script name is passed into the script as `$1`, `$2`, ...

```bash
#!/bin/bash
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --time=1:00:00
#DEP: samtools/1.22

# $1 = input BAM, $2 = output directory
mkdir -p "$2"
samtools sort -@ $NCPUS "$1" -o "$2/sorted.bam"
```

Submit with:

```bash
condatainer run sort.sh /data/raw/sample1.bam results/sample1/
# inside sort.sh: $1=/data/raw/sample1.bam  $2=results/sample1/
```

### Array Jobs

Run the same script over a list of inputs with `--array`. Each non-empty line in the input file becomes one subjob; its space-separated tokens arrive as positional arguments (`$1`, `$2`, …) inside the script.

```bash
# samples.txt — one entry per line; tokens become $1, $2, ...
sample1 condition_A
sample2 condition_B
sample3 condition_C
```

```bash
# quant.sh — $1=sample name, $2=condition
salmon quant -i $SALMON_INDEX_DIR -p $NCPUS -l A \
  -r reads/${1}.fq -o quants/${1}_${2}/
```

```bash
condatainer run --array samples.txt quant.sh
condatainer run --array samples.txt --array-limit 4 quant.sh  # max 4 concurrent
condatainer run --dry-run --array samples.txt quant.sh        # preview args
```

**Mixing array args with extra CLI args:**

When extra args are provided after the script name alongside `--array`, array line tokens are **prepended** before the CLI args inside the script:

```bash
condatainer run --array samples.txt quant.sh genome_version
# line "sample1 treated" → $1=sample1  $2=treated  $3=genome_version
```

**Rules:**
- All non-empty lines must have the same number of shell-split tokens (quoted strings count as one).
- Blank lines are flagged on `--dry-run` and block submission — remove them before submitting.
- Per-subjob logs go to `{logsDir}/{jobname}_{idx}_{tag}.log` where `{tag}` is all args joined, with non-alphanumeric characters replaced by `_`, capped at 20 characters.
- Array jobs always require scheduler submission (cannot run locally).

### Job Chaining

All `[CNT]` messages go to stderr, only job ID is printed to stdout, so you can capture it for downstream job submission.

```bash
TRIM=$(condatainer run trim.sh sample1)
ALIGN=$(condatainer run --afterok "$TRIM" align.sh sample1)
condatainer run --afterok "$ALIGN" quant.sh sample1
```

Three dependency flags are available:

| Flag | Runs when upstream job… |
|---|---|
| `--afterok IDS` | Succeeds |
| `--afternotok IDS` | Fails |
| `--afterany IDS` | Finishes (any outcome) |

Multiple job IDs can be passed as a colon-separated list (e.g. `--afterok 123:456:789`). All three flags can be combined in a single submission.

### Array Jobs + Chaining

`--array` and `--afterok` can be combined: each stage submits an array job and waits for the previous stage to **succeed** before starting. A final single job can collect results after all subjobs complete.

```bash
# Stage 1: trim all samples (no dependency)
JOB=$(condatainer run --array samples.txt --array-limit 10 trim.sh)

# Stage 2: align — waits for ALL trim subjobs to finish
JOB=$(condatainer run --array samples.txt --array-limit 10 --afterok "$JOB" align.sh)

# Stage 3: quant — waits for ALL align subjobs to finish
JOB=$(condatainer run --array samples.txt --array-limit 10 --afterok "$JOB" quant.sh)

# Final: single job collecting results — waits for ALL quant subjobs
condatainer run --afterok "$JOB" collect_results.sh samples.txt
```

Each array stage fans out across all samples in parallel (up to the concurrency limit), and `--afterok` ensures the next stage only starts once every subjob in the previous stage has succeeded. Use `--afterany` instead to proceed even if some subjobs failed.

If you don't like this all-or-nothing behavior, you can chain individual subjobs.

```bash
declare -a FINAL_JOBS

while IFS= read -r line; do
  if [[ -z "$line" ]]; then continue; fi  # skip blank lines
  TRIM=$(condatainer run trim.sh $line)
  ALIGN=$(condatainer run --afterok "$TRIM" align.sh $line)
  QUANT=$(condatainer run --afterok "$ALIGN" quant.sh $line)
  FINAL_JOBS+=("$QUANT")
done < samples.txt
JOB_IDS=$(IFS=:; echo "${FINAL_JOBS[*]}")
condatainer run --afterok "$JOB_IDS" collect_results.sh samples.txt
```

Do not quote the line variable in the loop if you want the tokens to be split into separate arguments.

### Usable ENV

When running scripts with scheduler directives, the following environment variables are automatically available inside the container:

| Variable | Description | Example |
|----------|-------------|---------|
| `NNODES` | Number of compute nodes | `2` |
| `NTASKS_PER_NODE` | Number of MPI tasks per node | `4` |
| `NTASKS` | Total number of MPI tasks | `8` |
| `NCPUS` | CPUs per task | `4` |
| `MEM` | Memory per task in MB | `8192` |
| `MEM_GB` | Memory per task in GB | `8` |

**Priority Order (Scheduler ENV > Script > Default):**

Environment variables are resolved with the following priority:

1. **Scheduler ENV** (highest priority) — Environment variables from the scheduler (e.g., `$SLURM_CPUS_PER_TASK`)
2. **Script** — Scheduler directives in the script (`#SBATCH`, `#PBS`, `#BSUB`)
3. **Default** (lowest priority) — Built-in default values in config (`scheduler.*`)

### MPI Auto-Detection

When a scheduler script requests more than one task (`--ntasks-per-node`, `--ntasks`, or PBS/LSF equivalents), `condatainer run` automatically detects the host MPI and wraps the job command with `mpiexec`:

**Detection:** `mpiexec` must be available in `$PATH` at submission time. If it is not found, `condatainer run` exits with an error.

**Generated job command:**

```bash
mpiexec condatainer run script.sh
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
# ml openmpi # or other MPI module to ensure mpiexec is in PATH
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
condatainer e mpi.img -- mm install mpi4py openmpi=4.1 -y
```
````

### Module commands are not dependencies

**CondaTainer** reads `#DEP:` and nothing else. A `module load` or `ml` line names your site's own
module tree — a separate tool, with its own builds, that only happens to spell names the same way —
so it is left alone and never mounts an overlay. Declare what the script needs:

```bash
#!/bin/bash
#DEP: bcftools/1.22
bcftools --version
```

## Info

Display detailed metadata about an installed overlay, the base image, or an external file. Accepts an installed overlay name (`name/version`), the configured base recipe name, or a direct file path (`.sqf` / `.img`, or a `.sif` from elsewhere).

**Usage:**

```
condatainer info [--verify] <overlay>
```

**Examples:**

```bash
condatainer info samtools/1.22
condatainer info ubuntu24/base          # the base image
condatainer info env.img
condatainer info ./ubuntu--22.04.sqf
condatainer info nvim --verify          # also check its files against the payload key
```

`--verify` reads a `.sqf` overlay's files and checks them against the payload key its manifest records,
shown as `Payload`. It reads every byte, so it takes as long as reading the overlay, and it needs
`squashfuse` and `unshare`. Only script builds record a payload key; for a Conda or definition build it
says so.

CondaTainer builds every image as a `.sqf`, the base included. A `.sif` from
elsewhere reads the same way: its payload is a SquashFS partition starting
partway into the file, so the archive reads take that offset. The `Type` line is
the type from the image's embedded metadata — `app`, `os`, `data`, or
`environment` for one produced by `overlay freeze`. An image built before that
metadata existed reports `unknown`.

### Image output

| Section | Fields |
|---------|--------|
| **File** | Name, Path, file Size, Type (`app` / `os` / `data` / `environment` / `unknown`, Read-Only), and Created timestamp |
| **SquashFS** | Compression algorithm (with level if set), Block Size, Inode count, Fragment count, Deduplication flag |
| **Payload** | `Prefix` (`/cnt/<name>/<version>`, for `app` and `data`); `Deletions` for a frozen environment. Read from the image's metadata, never from its filename |
| **Environment** | Variables from the image's embedded metadata, with their notes |


### ext3 (`.img`) output

| Section | Fields |
|---------|--------|
| **File** | Name, Path, file Size, Type (always `environment`, Writable; sparse images also show actual on-disk size) |
| **Filesystem** | Format, State, Block Size, Created, Modified, Last Mounted |
| **Ownership** | Inner UID/GID of files inside the image (or `root` for fakeroot-compatible images) |
| **Disk Usage** | Used / Total (%), Reserved blocks, Free |
| **Inode Usage** | Used / Total (%), Free |
| **Payload** | `Prefix` — always `/cnt_env` for a writable overlay |
| **Conda Env** | Channels and explicitly-installed packages, from `conda-meta/history` |
| **Environment** | Variables from the `.env` sidecar file (`KEY=value ## note`) |

When the `.img` autoloads a paired snapshot (a personal `<name>-<user>.sqf` line, else a shared
`<name>.sqf` line — see [Overlay Freeze](#overlay-freeze)), both the **Conda Env** and **Environment**
sections show the merged view — the snapshot's own channels, packages, and variables, with the `.img`'s
own edits on top — not the thin `.img` alone.

## Export

Export the Conda environment in a writable `.img` overlay, to stdout (or a file
with `-p`). Runs `micromamba env export` against `/cnt_env`, so it captures the
environment as it is now — including anything installed since the overlay was
created. That is why it is limited to `.img`: an installed image is immutable,
carries its own metadata, and is reproduced by rebuilding from its recipe rather
than by recovering one from the image.

Exporting an installed image reports what to do instead:

```
export needs a writable .img overlay; …/samtools--1.22.sqf is an installed image.
Rebuild it from its recipe instead: condatainer info samtools--1.22.sqf shows what it is
```

Definition-built images keep their definition regardless: it is recorded at
`/.cnt/recipe`, and Apptainer's own copy travels in the payload at
`/.singularity.d/Singularity`.

**Usage:**

```
condatainer overlay export [OPTIONS] <overlay.img>
```

**Options:**

* `-p`, `--prefix <path>`: Write to `<path>.<ext>` instead of stdout.

* `-e`, `--explicit`: Explicit (URL-pinned) format — exact, no re-solve. Round-trips via `create -f spec.txt`.
* `--no-build`, `--no-builds`: Strip build strings from the spec.
* `--channel-subdir`: Prefix each dependency with its `channel/subdir`.
* `--no-md5`: Disable MD5 checksums in explicit output.
* `--from-history`: Reconstruct spec from install history only.

**Channel order:** `micromamba env export` sorts the `channels:` block alphabetically, which can break re-solve. When an overlay carries a `.condarc`, its channel priority is used to reorder the block (matching `conda env export`); channels not listed there are prepended.

**In-use overlays:** exporting an ext3 `.img` while a writable session holds it fails with a clear "is open for writing by another process" error. (If you are in an ext3 overlay, use `mm export`)

**Examples:**

```bash
# Conda environment -> environment.yml on stdout
condatainer overlay export env.img > environment.yml

# Write to a file; extension chosen by format (here .yml)
condatainer overlay export env.img -p ./env

# Explicit, fully reproducible (round-trips via create -f)
condatainer export env.img -e -p ./spec        # writes ./spec.txt
condatainer create -f ./spec.txt -p rebuilt

# From history (only user-requested packages)
condatainer export env.img --from-history > minimal.yml
```

## Helper

Download and manage small helper scripts stored in the `helper-scripts/` folder inside the CondaTainer repository.

```{note}
Helper commands are not available inside a container or a scheduler job.
```

**Usage:**

```
condatainer helper [FLAGS] [SCRIPT_NAME] [SCRIPT_ARGS...]
```

Options:

* `-u`, `--update`: Update helper scripts from remote metadata.
* `--path`: Show all helper script search paths and the writable directory. If a `SCRIPT_NAME` is given, print the absolute path of that specific helper script and exit.
* `-l`, `--list`: List available helper scripts with their descriptions (from `#DESCRIPTION` tags).
* `--no-submit`: Disable job submission for this run — run headless on this node instead.
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

# Download/Update all helper scripts
condatainer helper --update

# Run a helper script with arguments (e.g. request 4 CPUs)
condatainer helper code-server -c 4
```

To run a single helper without submitting a job, pass `--no-submit`
(`condatainer helper --no-submit code-server`) — useful for a quick inspection or debug session,
or to keep working when the scheduler itself is unavailable. To disable submission for every run,
set `submit_job: false` in the config (`condatainer config set submit_job false`) instead. The
dashboard offers the same override as a "Run headless" checkbox next to Start, shown only when a
scheduler is actually available to opt out of.

Pass `-A`/`--account` and `-p`/`--partition` to submit under a billing account or a specific
partition/queue for one run (`condatainer helper code-server -A myproj -p gpu`), or set
`scheduler.account`/`scheduler.partition` in the config for every run — see
[Configuration](configuration.md). Both can also be typed at the interactive settings prompt
(`a`/`p`), and the dashboard exposes matching Account/Partition fields next to Start.

## Config

Manage **CondaTainer** configuration settings.

### Config Show

Display current configuration including file paths, settings, and environment variable overrides. (including available distros)

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
condatainer config get build.system_apptainer
condatainer config get build.ncpus
condatainer config get submit_job
```

### Config Set

Set a scalar configuration value and save to the active config file. For a full list of supported keys, see the [Configuration manual](configuration.md).

```
condatainer config set <key> <value>
```

**Examples:**

```bash
condatainer config set build.system_apptainer /usr/bin/apptainer
condatainer config set submit_job false
condatainer config set build.ncpus 8
condatainer config set build.time 4h
```

**Time formats:** `2h`, `30m`, `1h30m`, `90s`, `02:00:00`, `HH:MM:SS`

**Where it writes:** the highest-priority config file that **exists** — your user config if you have one, otherwise the next layer down (extra-root, app-root, system). If that file is read-only, the write goes to your user config instead (created if needed) and says so:

```
[CNT!] app-root config is read-only; saving to your user config instead (applies only to you).
```

Pass `-l/--layer` to target a layer explicitly (`user`, `app-root`, `extra-root`, `system`). A read-only target is then an error, not a fallback.

### Config Append / Prepend / Remove

Manage array config keys (`sources`, `channels`) from the CLI.

```
condatainer config append  <key> <value>
condatainer config prepend <key> <value>
condatainer config remove  <key> <value>
```

* `append` — add a value to the **end** of the array (lower priority).
* `prepend` — add a value to the **beginning** of the array (higher priority).
* `remove` — delete all occurrences of a value from the array.

**Examples:**

```bash
# Add an institutional scripts source (takes priority over the default)
condatainer config prepend sources myorg=https://raw.githubusercontent.com/MyOrg/recipes/main

# Add a lab collection at lower priority
condatainer config append sources lab=/shared/lab/recipes

# Remove entries
condatainer config remove sources lab=/shared/lab/recipes

# Show result
condatainer config show
condatainer config get sources
```

### Config Init

Create a config file with auto-detected defaults.

```
condatainer config init [-l|--layer user|app-root|extra-root|system]
```

* `-l`, `--layer`: Config layer (default: auto-detect).

**Auto-detects:**
* Apptainer binary
* Scheduler binary
* Compression support (zstd vs lz4)

### Config Paths

Show the recipe sources, then the data search paths for images and helper scripts, in read order (nearest first). Writes go the opposite way, to the first writable directory starting from the furthest-out layer.

```
condatainer config paths
```

Each entry is tagged with its data layer — `(user)`, `(extra-root)`, `(app-root)` — and its status, including which directory receives writes:

```
Images:
  1. /shared/labA/condatainer/images (extra-root) (writable, target)
  2. /opt/condatainer/images (app-root) (read-only)
  3. /scratch/me/condatainer/images (user) (not found)
```

The layer names are the values `-l`/`--layer` accepts (`u`, `r`, `e`).

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
2. Environment variables (`CNT_*`)
3. User config file (`~/.config/condatainer/config.yaml`)
4. Extra-root config (`$CNT_EXTRA_ROOT/config.yaml`, group/lab layer)
5. App-root config (`$CNT_ROOT/config.yaml` or `<install-dir>/config.yaml`)
6. System config (`/etc/condatainer/config.yaml`)
7. Built-in defaults

### Configuration File Example

See the [Configuration manual](configuration.md) for a full reference of all available keys.

```yaml
build:
  system_apptainer: /usr/bin/apptainer  # scheduler type is auto-detected from scheduler.bin
scheduler:
  bin: /usr/bin/sbatch

# Submission settings
submit_job: true

# Default distro for the container root (default: the first source's default_distro)
default_distro: ubuntu24

# Recipe collections, in priority order (first match wins).
# The public cnt collection is appended automatically unless redefined here.
sources:
  - cnt: https://raw.githubusercontent.com/condatainer/recipes/main

# Days to cache remote metadata (default: 7, set 0 to always fetch live)
metadata_cache_ttl: 7

# Build configuration
build:
  ncpus: 8
  mem: 16g
  time: 4h
  # compress_args options (gzip, lz4, zstd, zstd-fast, zstd-medium, zstd-high)
  compress_args: "-comp zstd -Xcompression-level 8"
# For a group/lab root with standard layout, set in module file:
# export CNT_EXTRA_ROOT=/project/shared/condatainer
```

## Registry

Publish and fetch completed `.sqf` images — overlays and bases alike — through an
OCI registry. Writable `.img` overlays are never distributable.

```
condatainer registry push|pull|tags|resolve|login|logout
```

Registry bases include the host, owner, and optional prefix, for example
`ghcr.io/my-lab/condatainer`. Push takes either an installed `name/version` or a
path to an artifact.

Only an installed `name/version` infers its destination, from the artifact's
recorded recipe source and that configured source's descriptor; `--registry` is
then an explicit override, and is required when the provenance matches no single
configured source. **A path always requires `--registry`** — a file you point at
is typically a one-off build, a mirror, or a project's own image, so its
destination is stated rather than guessed. Pull, tags, and resolve require it
too, because they do not start from a local artifact.

```bash
# Save a token for a whole registry
printf '%s\n' "$TOKEN" | condatainer registry login ghcr.io \
  --username "$USER" --password-stdin

# Save a different token for one repository
condatainer registry login ghcr.io/my-lab/rnaseq --username "$USER"

# Save a token in a shared layer, where everyone who reads that layer finds it
condatainer registry login ghcr.io -l extra-root --username group-bot

# Show the saved credentials and their layers (never the secrets)
condatainer registry list

# Publish by installed name, inferring the selected source's oci.push endpoint.
condatainer registry push grch38/genome/gencode49

# A path states its endpoint. Public is the safe default; say restricted only
# when a known set of people are the only ones who can pull.
condatainer registry push ./licensed-app.sqf \
  --registry registry.lab.example/cnt --audience restricted

# Inspect a repository, or resolve one exact platform manifest
condatainer registry tags grch38/genome \
  --registry ghcr.io/my-lab/condatainer
condatainer registry resolve grch38/genome:gencode49 \
  --registry ghcr.io/my-lab/condatainer

# Install an exact tag/digest into the managed images directory
condatainer registry pull grch38/genome:gencode49 \
  --registry ghcr.io/my-lab/condatainer
condatainer registry pull ubuntu24/base@sha256:<digest> \
  --registry ghcr.io/my-lab/condatainer

# Override the managed name, or write to an exact external path
condatainer registry pull grch38/genome:gencode49 --name grch38/genome/gencode49 \
  --registry ghcr.io/my-lab/condatainer
condatainer registry pull grch38/genome:gencode49 --prefix /project/images/gencode49 \
  --registry ghcr.io/my-lab/condatainer
```

`pull` requires an exact address. A bare name belongs to `create`, which owns
version selection and may build when no published artifact exists. Pull never
falls back to a build. Placement precedence is `--prefix`, `--name`, the exact
address when it contains a complete name, then the published OCI title.

Credentials come from `GITHUB_TOKEN` for `ghcr.io`, then from those saved by
`registry login`, then anonymous access. A saved credential is keyed by a host
(`ghcr.io`) or by one repository (`ghcr.io/my-lab/rnaseq`). A push or pull uses the
most specific key that covers its repository, and for the same key the nearest
config layer: user, extra-root, app-root, system. So a group's read-only token on
the host and your write token on one repository can be saved side by side.

`registry login -l <layer>` chooses the layer; the default is `user`. Credentials
are saved in `registry-auth.json` beside that layer's `config.yaml`, readable only
by you (mode 0600). Saving into a shared layer warns that anyone who can read the
file can use the token, asks for confirmation (`-y` answers yes), and leaves the
file's permissions to you.

`registry list` shows who can read each file: `you`, `group <name>`, or `everyone`,
from the file's and its directory's permissions. `login` and `list` also warn when:

- the `user` layer's file can be read by anyone else (`chmod 600` fixes it);
- anyone other than you can replace a credential, because its file or its
  directory is writable by others. In a shared layer, group members replacing it
  is only a note.

Access control lists are not checked.

Versioned tags are immutable unless `push --force` is used; version-less OS
artifacts publish a `YYYYMMDD` tag and `latest`.

### Large pushes

A large artifact is split into layers, and `push` prints the plan it settled on
before sending anything:

```
[CNT] upload plan size=26.00 GB layer-size=2.00 GB (floor) layers=13 requests=~39 registry=ghcr.io
```

There is no setting for the layer size. It is derived from the artifact, because
the limit that decides it is a count of requests rather than of bytes: a fixed
size would make the layer count grow with the artifact and eventually exceed what
a registry accepts in one run.

Registries refuse a client that writes too often, so `push` **pauses and
resumes** rather than failing:

```
[CNT] upload rate limited layer=4/11 retry=1/4 wait=1m0s
[CNT] upload resumed layer=4/11
```

That is normal for a multi-hour transfer and needs no action. Four such pauses
end the push. A genuine permission failure is reported immediately instead, and
never waited on. Interrupting during a pause exits at once.

`pull` waits out the same limits. A rate-limited pull never silently falls back
to building locally — it is a wait, not a missing artifact.

### Tag schemes

Two different naming schemes publish to a registry, and which one you are looking
at is decided by the command that pushed, not by the registry:

| | **catalog scheme** (`registry push`) | **project scheme** (`project push`) |
|---|---|---|
| repository | the artifact name: `grch38/genome` | one per project: `my-lab/rnaseq-2026/cnt` |
| tag | the version segment: `gencode49` | the whole name, `/` → `--`: `grch38--genome--gencode49` |
| version-less | `YYYYMMDD` plus `latest` | no special case — every artifact has an identity |
| extra tag | none | `<name>__<12 hex>`, the artifact's identity |
| packages | one per artifact name | one per project |
| mutability | a versioned tag is immutable without `--force` | the plain tag moves; the identity tag never does |

So `grch38--genome--gencode49__a31f902c12ab` in a project package and
`grch38/genome:gencode49` in a collection package can be the same bytes; they are
addressed differently because they are retained differently. See
[Project Push](#project-push) for why a project needs the identity tag.

## Store

Most overlays answer to a plain name: one `star/2.7.11b`, in one file, found by
every command that takes a name. The **store** is where the extra ones go when
that is not enough — a second build of a name, filed under its exact identity so
both can be installed at once.

```
condatainer store add | list | path | use | validate | rm | gc
```

An identity names one exact build. Any command taking `--identity` accepts the
complete `scheme@sha256:…` key or any unambiguous prefix of its digest;
`condatainer info` prints both keys for an overlay.

### Getting something into the store

Two ways, depending on whether it exists yet:

```bash
# Build a second one, keeping the installed build of that name
condatainer create --store star/2.7.11b

# Install an overlay file you already have
condatainer store add ./star-2.7.11b.sqf
condatainer store add /scratch/builds/star.sqf --layer user
```

`store add` copies — never moves, never links — so the source file stays exactly
where it is. The name and identity come from the file's own metadata, never from
its filename, and inside the store the filename is generated from the keys. Use
`-l`/`--layer` to choose *which* store; to put an overlay at a path of your own
choosing use `create -p` or `registry pull -p` instead.

If the identity is already installed, nothing is copied and the existing path is
reported instead. With `--layer` that check is limited to the layer you named, so
asking for a build in a particular directory puts it there even when another
directory already has one — which is how you lift a shadow before `store use`.

### Choosing which build a name means

A store entry is reachable only by its identity: `condatainer list` and
`exec -o star/2.7.11b` see the plain-name build. `store use` swaps them.

```bash
condatainer store list star/2.7.11b            # what else is installed
condatainer store use star/2.7.11b --identity 9f2c1ab
```

The chosen build takes the plain name and the one that held it moves into the
same directory's store. Both stay installed and both stay resolvable by
identity — **nothing is deleted, and no project lock breaks**, because a lock
pins an identity rather than a path.

Two rules follow from how overlays are found. Reads resolve nearest-first across
[data layers](../deployment/data_layers.md), so:

* a build in a directory that a nearer one already shadows is **refused** — promoting it there would change nothing. Copy it into the nearer directory first with `store add --layer`, then run `store use`.
* a build in a *nearer* directory than the current holder simply wins, and the farther copy is left alone and reported as shadowed.

`store use` only ever renames, and only inside one directory. That is what keeps
it safe for other people: nothing leaves the directory, so no identity vanishes
from anyone's view whatever their own layer configuration. In a shared directory
it does change what the name means for everyone reading it, so it asks first
unless `-y` is given.

### Inspecting and reclaiming

```bash
condatainer store list                                 # name, identity, size, path
condatainer store list star/2.7.11b --equiv 9f2c1ab    # what could stand in for it
condatainer store path star/2.7.11b --identity 9f2c1ab # one exact path, for scripts
condatainer store validate                             # re-check every entry's keys
condatainer store validate --payload                   # also read each entry's files (slow)
condatainer store rm star/2.7.11b --identity 9f2c1ab   # delete one entry
condatainer store gc                                   # what could be reclaimed
condatainer store gc --layer user --apply              # reclaim it
```

`store rm` deletes store entries only — an overlay under a plain name belongs to
`condatainer remove`. Neither touches an overlay whose write bit is clear
(`chmod a-w` is how an artifact is pinned) or one a running container is reading.

`store validate --payload` also reads each entry's files and checks them against the payload key its
manifest records (shown by `condatainer info` as `Payload`). It reads every byte, so it takes as long
as reading the store, and it needs `squashfuse` and `unshare`. Only script builds record a payload key;
Conda and definition builds are checked by their other keys.

`gc` reports by default and needs `--dir` or `--layer` before `--apply`, so a
directory shared with people who are not at the keyboard is never collected by
omission.

## Project

Pin the exact artifacts a project mounts, and carry that pin in Git.

```
condatainer project lock | pin | unpin | select-distro | select-match | list | validate | restore | registry [set|unset] | push
```

A project is any directory holding `cnt-lock/`, which `project lock` creates.
Every subcommand acts on the current directory unless `--project DIR` names
another. `run`, `check`, `exec -o` and `e -o` instead act *in* whatever project
you are standing in, with no flag either way.

Commit `cnt-lock/` alongside the code. It holds the pins plus the
manifest and rebuild sources of every pinned artifact, and no payload, no
absolute path and nothing machine-local, so a checkout reproduces the same
artifacts on a machine that has never run CondaTainer:

```text
project/
  analysis.sh
  cnt-lock/
    lock.json
    provenance/
      star--2.7.11b@a31f902c12ab/
        manifest.json
        recipe
```

### Identity and equivalence

Every artifact has two keys, both computed from how it was built and never from
the files inside it, so a rebuild from the same inputs gets the same keys on any
machine.

- **Identity** is which exact build this is. A pin records it.
- **Equivalence** is whether this build can stand in for another. It moves only
  when what you asked for changes, not when the build environment does.

**What changes them**

| change | identity | equivalence |
|---|---|---|
| edit the recipe's commands | changes | changes |
| choose another `#PH:` or `#ENV:` value | changes | changes |
| an upstream `#SOURCE:` file changes | changes | same |
| a `data` dependency is rebuilt as an equivalent, not identical, build | changes | same |
| an `app` or `os` dependency named in the artifact's own name changes version | changes | changes |
| the same version of that dependency is a different build | changes | same |
| a dependency used only to build it (a mounted `build-essential`) changes | changes | same |
| reword a comment or `#DESC:`, edit a scheduler directive, extend a `#PH:` menu | same | same |

Edits that cannot change what is built leave both keys alone, so they do not create a
new artifact.

An `os` definition is the same, with its upstream image in place of `#SOURCE:`: when
the `From:` image is updated, rebuilding gives a new identity and the same
equivalence.

A Conda environment has two keys from its two exports. Identity is `explicit.txt`
(the exact package URLs). Equivalence is `environment.yml` (the versions asked for).
A new solve that lands on the same versions with different build strings is
equivalent, not identical. A frozen environment has one key for both, so nothing
substitutes for it.

**What you see**

| where | what it shows |
|---|---|
| `condatainer info <name>` | the `Identity:` and `Equivalence:` lines |
| `project restore` | `(equivalent, not sha256:…; differs: src:gtf)` when a substitute was used or built, with the inputs that differ |
| `project restore --dry-run` | the same note on a step a substitute already answers |
| `run`, `exec` | `is mounted from an equivalent artifact, not <identity>` |
| `create` | `this exact build is already installed, skipping` when the identity is present |

Under an [`identity` project](#project-select-match) each of these is an error instead
of a note.

### Project Lock

Rescan every script for `#DEP:` declarations, pin what they name, and reconcile
the result with the lock. Pins nothing declares any more are dropped, and so are
pins whose artifact no longer verifies.

It creates `cnt-lock/` in the current directory when there is none, so this is
also how a project starts. The other subcommands act on pins that must
already exist, and refuse rather than create one.

```bash
# Start a project, or rescan an existing one
condatainer project lock
```

Scan rules:

- `.sh` and `.bash` files only;
- `cnt-lock/` and dot directories are skipped;
- build recipes are skipped;
- symlinks are not followed.

A script that expands `$CNT_PREFIX` is a build recipe and is not scanned. Its
`#DEP:` are the build dependencies of the artifact it produces, already recorded
in that artifact's provenance, so reading them would make the project pin what
it never mounts. Comments are stripped before the check, so mentioning the
variable in a comment does not skip a script.

Declaring the built overlay is unaffected: `#DEP: overlays/tool.sqf` in an
analysis script is pinned as usual.

A declaration whose overlay is not installed fails the lock. Other installed
builds of the same name are listed rather than pinned; `project pin` takes one
of those instead.

Every run also pins the project's **root** — which distro's `.../base` a
restore builds inside — from the configured `default_distro`, whether or not
anything in the project actually needs one. See
[Project Select-Distro](#project-select-distro).

### Project Pin

Resolve one declaration to one exact local artifact and vendor its sources.

```bash
# Pick one exact artifact by identity prefix, digest, or identity reference
condatainer project pin star/2.7.11b a31f902c12ab

# A project-relative .sqf takes no identity: it already names the file it means
condatainer project pin overlays/combined.sqf
```

The request is written the way the `#DEP:` writes it — a name, or a
project-relative path, classified by extension exactly as a declaration is.

A **file is never accepted as the identity.** Restore and push both find a
pinned artifact through the image roots and the store, so pointing at a file
elsewhere would record an identity only a rebuild could satisfy. To pin a `.sqf`
sitting outside the project, install it into an images root first and pin it by
identity.

Every copy in every readable image root is a candidate, not just the nearest one,
so a pin can name a copy that ordinary name resolution would hide.

A pin also vendors the artifact's source closure into `cnt-lock/provenance/`,
and — best effort, over the network — records any collection endpoint that
already publishes that exact identity (see [Two kinds of remote](#two-kinds-of-remote)).

A dependency nothing can pin is always a scan finding: `project validate` fails
on it, and `project lock` warns and publishes what it could. Nothing written in
a script silences one.

A writable `.img` is the usual case, and the finding names the remedy:

```
run.sh:2: env.img is writable, so it has no identity to pin; freeze it into the
project with `condatainer overlay freeze env.img overlays/<name>.sqf` and
declare that instead
```

For a `.sqf` outside the project, copy it under the project and declare that
path — restore only writes to paths the project owns.

#### Manual pins

A pin is **manual** when nothing in the project declares it. Pinning is checked
against a scan once, when the pin is made, and the answer is recorded in the lock:

```json
"path:overlays/env.sqf": { "artifact": "provenance/env@0a398f4644c7", "manual": true }
```

Artifacts reach a project without a `#DEP:` routinely: a helper script names its
service overlays in `#REQUIRED_OVERLAYS:`, and a frozen environment is named by
nobody. Pin those by hand:

```bash
condatainer project pin overlays/env.sqf
condatainer project pin ubuntu24/build-essential
```

`project lock` sweeps every pin its scan no longer produces, and a manual pin is
exempt from that sweep. That is the whole of what the flag does, so removing one
takes `project unpin`.

### Project Unpin

```bash
condatainer project unpin overlays/env.sqf
condatainer project unpin ubuntu24/build-essential
```

Removes a manual pin, and with it every vendored artifact directory and remote
that nothing else reaches once it is gone.

**Only a manual pin.** A pin a declaration produced is refused here: the next
`project lock` would re-pin the same declaration, possibly at a different
identity. Delete the `#DEP:` line and run `project lock`, which drops every pin
nothing asks for any more.

**No image is deleted.** A `path:` pin keeps its file where it is; a named pin
keeps its overlay in the images root. Only the lock's record of it goes.

### Project Select-Distro

```bash
condatainer project select-distro rocky9
condatainer project select-distro --auto
```

Overrides which distro's `.../base` `project restore` builds inside, instead
of the one `project lock` derives from the configured `default_distro`. There
is one root per project, so re-running **replaces** the choice rather than
adding to it. `--auto` clears the override and re-derives immediately, rather
than waiting for the next `project lock`.

Every collaborator restoring the project gets this root regardless of their
own `default_distro` — that is the whole point. An artifact reaching outside
`/cnt_env` (a frozen environment, or an `os` overlay layered above the root)
still couples to whichever release it actually ran against; recording a root
does not check that, it only removes "my machine is configured differently"
as a way two checkouts of the same project can disagree.

Standing in a project, this is also what every bare `<distro>/<name>` shortcut
expands against — `avail`, `list`, `remove`, `info`, `overlay`, shell
completion, and `create`'s own bare-name expansion all follow the selected
root instead of the configured `default_distro`.

### Project Select-Match

```bash
condatainer project select-match identity
condatainer project select-match equivalence
```

Records in the lock which copy of a pinned artifact the project accepts.
`equivalence`, the default, takes any equivalent build: one that can stand in for the pinned one.
`identity` takes only the exact pinned build.

For example, a Conda environment is pinned and later rebuilt. The package
versions are the same, but the build strings differ:

| mode | the rebuilt copy |
|---|---|
| `equivalence` | accepted; the result line says it is equivalent, not identical |
| `identity` | refused; only a copy with the pinned build strings is used |

A data artifact behaves the same way. If a `#SOURCE:` file was re-released
upstream under the same name, the rebuilt artifact is equivalent but not
identical. See [Identity and equivalence](#identity-and-equivalence) for what
each key covers.

The mode applies wherever the pins are used: `project restore` refuses a
substitute, `project validate --installed` reports one, and `run` and `exec`
stop when only a substitute is installed. Commit the change and every
collaborator is held to it.

### Project List

```bash
condatainer project list
condatainer project list --json
```

Lists every pin and the artifact it holds:

```text
  base (ubuntu24/base)         sha256:8f1c02de41ab
  star/2.7.11b                 sha256:a31f902c12ab
  path:overlays/env.sqf (env)  sha256:0a398f4644c7  manual
[CNT] 3 pin(s), 1 manual
```

`base` is the reserved root pin every project carries — see
[Project Select-Distro](#project-select-distro) — labelled by its role rather
than by its raw lock key.

Reads `cnt-lock/` and nothing else — no scan, no network, no image opened — so
it answers on a fresh clone that has restored nothing.

The first column is the **pin key**, which is what `project pin` and
`project unpin` take. A path key addresses a file and says nothing about what is
in it, so the artifact's name follows in parentheses; a name key already is the
name and is not repeated. Rows go to stdout and the summary to stderr, so
`project list | grep manual` gets the pins and nothing else.

A pin whose artifact is missing or does not verify is listed as `unreadable`.
Why it is unreadable is `project validate`'s answer; whether an artifact is
actually *available* is `project restore --dry-run`'s.

### Project Validate

```bash
condatainer project validate
condatainer project validate --json
condatainer project validate --installed
condatainer project validate --payload
```

By default it checks the lock alone:

- every declaration in the scripts has a pin
- the project's root is pinned
- every pin points at a vendored artifact
- every vendored artifact regenerates the keys it records
- every dependency edge resolves to another vendored artifact

By default it reads `cnt-lock/` and the project's scripts, and nothing else — no overlay, no
store, no catalog, no configuration, no network — so it answers the same on a
fresh clone that has restored nothing. Every key is recomputed from the vendored
sources rather than trusted, and every dependency edge is followed by name *and*
identity.

By default it never asks whether an artifact is installed here.

| flag | effect |
|---|---|
| `--installed` | also fail for every pin that is not installed here as an artifact the project's [match mode](#project-select-match) accepts |
| `--payload` | also read each installed artifact's files and check them against its recorded payload key; implies `--installed` |

`--installed` uses the lookup [`project restore --dry-run`](#project-restore) plans with, and a build
dependency that only a rebuild would need is not a problem. `--payload` reads every byte, so it takes as
long as reading the overlays, and it needs `squashfuse` and `unshare`. Only script builds record a
payload key; a Conda or definition build is skipped.

### Project Restore

```bash
condatainer project restore
condatainer project restore --dry-run
condatainer project restore --no-prebuilt
condatainer project restore --keep-build-deps
```

Reuses, fetches, or rebuilds each locked artifact until the project can run.
Nothing is mounted and `cnt-lock/` is never modified.

| flag | effect |
|---|---|
| `--dry-run` | report what would happen and acquire nothing |
| `--no-prebuilt` | build missing artifacts from source where possible instead of downloading a recorded one |
| `--keep-build-deps` | install newly produced build dependencies instead of discarding them |
| `--replace` | overwrite a project path holding something the lock does not name |
| `--json` | print JSON |

**A published overlay is downloaded; anything else is rebuilt.** An artifact with a
recorded location is fetched by content, so what arrives is exactly the recorded
build and nothing upstream is consulted. One with no recorded location is rebuilt
from the recipe in `cnt-lock/`, which downloads its `#SOURCE:` files again. If an
upstream file has changed since, the rebuild is a different build.

Under the default `equivalence` match it is accepted, and the result line names what differs:

```
  built    grch38/gtf-gencode/49 → …/49.sqf (equivalent, not sha256:41ab1c2d3e4f; differs: src:gtf)
```

`--json` carries the digests of each input that differs. Under an `identity`
project the restore fails instead. An artifact addressed by name is kept in the store under
its own identity, and `condatainer project pin` can pin it.

A rebuild can also match the locked identity and still produce different files, because the
recipe is not deterministic. The restore keeps the result and warns:

```
  built    grch38/gtf-gencode/49 → …/49.sqf
  warning: grch38/gtf-gencode/49 rebuilt with the locked identity but different files: its recipe does not produce the same output twice
```

To keep an analysis reproducible, publish its overlays with
[`project push`](#project-push).

**A project path is the only file a restore can destroy.** An artifact addressed
by name goes to a store name no other identity holds, so nothing there is ever
clobbered. A `path:` pin is materialized at exactly the path it declares, and
whatever is there is renamed over.

So a path already holding something the lock does not name — a `.sqf` copied in
from elsewhere, a stale build — is **refused during planning**, before anything
is fetched or built:

```
[CNT✗] overlays/combined.sqf already holds sha256:7c4e11ab27f0, which is not
       what combined pins; pass --replace to overwrite it
```

`--dry-run` reports the same, and shows what each step would replace. `--replace`
is the deliberate act that proceeds.

A named pin lands where the store's destination rule puts it — the flat name when
it is free, `store/` when that name is already held at a different identity. An
artifact nothing pins is a build dependency: it exists only so its dependent can
be built, so it is produced in a temporary directory and removed when the restore
ends.

`--no-prebuilt` is not an offline mode: every build needs the network, which is
what the [proxy](#proxy) is for on a compute node without egress.

**A frozen environment cannot be rebuilt.** It was captured from a writable
overlay rather than built from a recipe, so it vendors no sources and a registry
copy is the only thing that produces it. A restore that cannot find one is
refused during planning:

```
[CNT✗] env is a frozen environment and cannot be rebuilt; it is only obtainable
       from a registry, so publish it with `condatainer project push` from a
       checkout that has it
```

`--no-prebuilt` chooses building over downloading, so it does not apply here:
a published frozen environment is still fetched under it, because there is no
source to build from instead. A project carrying one must publish it for any
other checkout to restore.

A rebuild whose recipe carries scheduler directives is submitted rather than run
here, and so is anything waiting on it; dependencies become `afterok` edges. The
command then exits with the jobs-submitted code, having made nothing available
yet. Re-running the restore is how it resumes. Restore is atomic per artifact,
not across the project: a later failure leaves earlier results in place, and
re-running adopts them.

### Two kinds of remote

The lock's `remotes` map says where each artifact can be **fetched** from, as a
repository coordinate plus a platform manifest digest — never a mutable tag.
Restore tries those locations before it builds. Entries get written two ways, and
a project normally uses both:

| | **upstream remote** | **project remote** |
|---|---|---|
| where | the collection's `oci.pull` endpoints, from its `source.json` | the endpoint `project registry set` recorded |
| written by | `project pin`, best effort | `project push` |
| naming | the catalog scheme | the project scheme |
| costs | nothing — the bytes are already there | the upload |
| covers | artifacts whose pinned identity is *literally* what the collection publishes | everything else |
| survives | as long as the collection keeps it | as long as the project's own package lives |

An upstream remote is recordable only when the pinned artifact's **complete
identity** matches what the endpoint advertises, which in practice means the
artifact was pulled from there in the first place. Order in the list is retry
priority, and the free location ends up first because a pin records it
before any push runs.

Remotes are written by machines, never typed: one is recorded only after
something confirmed the artifact is actually there. There is no flag for entering
a coordinate by hand.

### Project Registry

Show, set, or clear the publish destination.

```bash
# Show what is recorded
condatainer project registry

# Record the destination — one repository for the whole project
condatainer project registry set ghcr.io/my-lab/rnaseq-2026/cnt \
  --audience restricted

# Forget it; already-recorded fetch locations are left alone
condatainer project registry unset
```

The endpoint is a complete repository coordinate with no tag or digest. It is
tracked in `lock.json`, not machine-local: every collaborator publishes to the
same package, or the recorded fetch locations become a set of places *some* of
the artifacts are.

`--audience` sets what `push` may publish to this destination — see
[Publishing rules](#publishing-rules). It defaults to `public`, the stricter
rule; a `restricted` destination takes anything, so say it only if the registry
really is limited to people who may receive everything. It does **not** control
who can pull, and it is not a GitHub package's visibility, which is a separate
per-package setting CondaTainer never reads or changes.

`--source` records the project's code repository, defaulting to its GitHub
origin, and becomes the published packages' source link.

Setting a destination writes to the lock and contacts nothing. Changing it does
not invalidate locations already recorded — those digests still resolve where
they were published.

### Project Push

```bash
# See the plan and its cost without uploading
condatainer project push --dry-run

# Publish the pins a collection does not already serve
condatainer project push

# Publish everything, build dependencies included
condatainer project push --all --closure
```

| flag | effect |
|---|---|
| `--all` | publish every pin, including ones a collection already serves |
| `--closure` | also publish build dependencies |
| `--registry` | publish to this repository instead of the recorded one |
| `--dry-run` | report what would be published and upload nothing |
| `--json` | emit machine-readable output |

The project must pass `condatainer project validate` first: a lock inconsistent
with what it vendors, a declaration nothing can pin, or a declaration with no pin
each stop the upload, and no flag overrides that. A published project is one
another checkout can restore, so a gap in the lock is a gap in what was
published. `--dry-run` reports the same problems and still shows the plan.

Every artifact goes into the one repository `project registry set` recorded, with
its name carried in the tag, so the whole project is one registry package rather
than one per artifact. That is the reason for the flat scheme: on GitHub a
distinct repository path is a distinct package, and each package's visibility is
its own setting that follows neither the linked code repository nor anything
CondaTainer records — so one package is one setting to manage by hand instead of
a dozen. See [Distributing Artifacts](../deployment/distribution.md).

Each artifact also gets an **identity tag**, `<name>__<12 hex>`. Nothing fetches
by tag — a `remotes` entry records the platform manifest digest — so what a tag
has to do is keep a manifest *referenced*, because what no tag reaches is the
registry's to reclaim. A plain name tag can only keep one identity alive, so
after a re-pin moved it, an older commit's recorded digest would resolve to
nothing. The plain tag is written only for the current pin, as a human
handle that says which identity the project uses now.

Two pins can share a manifest name — a `path:` pin of a locally built copy, or
two frozen environments, since every one of them is named `env`. Neither takes
the plain tag then; both keep their identity tags, and the push proceeds with a
warning naming what lost it.

The default push set is the pins **no collection already serves**: restore
tries the upstream location first, so a second copy of those buys nothing but the
bytes. `--all` publishes them anyway, which is worth it when the project must
outlive the collection's retention — an archived paper, a collection you do not
control, or a registry your compute nodes can reach when the collection's is not.

Push publishes what is already here and **never builds**. Anything missing at its
locked identity is reported with the restore that would produce it, and the same
goes for `--closure`, which will not acquire a build dependency that a
`restore --keep-build-deps` never installed. The project must also validate
first: what is published has to be what the checkout already describes.

`--dry-run` reports four dispositions per artifact:

```text
[CNT] Publishing to ghcr.io/my-lab/rnaseq-2026/cnt (restricted)
[CNT]   packages link to https://github.com/my-lab/rnaseq-2026
[CNT]   upload   star/2.7.11b  star--2.7.11b__a31f902c12ab star--2.7.11b
[CNT]   present  samtools/1.23.1  already published
[CNT]   upstream grch38/genome/gencode49  ghcr.io/cnt-recipes/cnt/grch38/genome
[CNT]   refused  cellranger/9.0.1  app artifacts are not published to a public registry…
[CNT] 3 to upload
```

`upstream` and `present` are decided by asking the network, and both only ever
*remove* work — an unreachable registry leaves the plan as computed rather than
failing a push nobody could complete offline. A `refused` step stops the push
before anything uploads.

There is no `--force`. An identity tag carries the content key, so it can only be
re-pushed with the same content, and the plain tag is a moving pointer replaced
by design. A content-keyed tag holding *different* content is refused outright,
because that is a corrupted package rather than something to overwrite.

Each successful upload records its location in the lock as its own transaction,
so an interrupted push leaves every artifact that did land recorded and
re-running skips them.

### Publishing rules

What may be pushed — by either `registry push` or `project push` — depends on who
can pull it and on what the artifact says about itself. Everything is read from
the embedded manifest, never guessed from a filename, and push **refuses** rather
than warns.

A `restricted` endpoint takes anything. At a `public` one:

| the recipe declared | may publish publicly |
|---|---|
| `#REDISTRIBUTE: no` | **never**, whatever the type, and no flag overrides it |
| `#REDISTRIBUTE: yes` | yes, whatever the type |
| nothing, and it is `os` | yes — a container root, and packages from a public distribution |
| nothing, and it is `data` | yes — the type asserts public reference data and the indexes built from it |
| nothing, and it is an `app` | **no** — someone else's software with unstated terms |
| nothing, and it is a Conda build | yes — it embeds no recipe, so it could never carry the declaration |
| nothing, and it is a frozen environment | yes — same reason: it has no recipe to declare in |

The declaration lives in the recipe rather than on the command line, so it is
authored once, reviewed in a commit, and travels with every artifact built from
it. Nothing verifies it. `#LICENSE:` is published verbatim as
`org.opencontainers.image.licenses` and never gates anything.

A Conda build's channels are **reported** at push time rather than judged — a
private or vendor channel is a question for a person, not for an allowlist:

```text
[CNT]   upload   rnaseq/1.0  rnaseq--1.0__6f21c0b81a4d
[CNT]            channels: conda-forge, bioconda
```

Writable `.img` overlays are never published to any endpoint. That is structural
rather than policy: a writable overlay has no identity, so there is nothing to
publish it *as*.

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
* GPUs (type, total)

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

Refreshes the cached recipe and helper indexes of every remote collection, or refreshes the
self-provisioned toolchain (mksquashfs, squashfuse, apptainer).

**Usage:**

```
condatainer update [FLAGS]
```

**Options:**

* `--build`: Refresh the build script metadata cache.
* `--helper`: Refresh the helper script metadata cache.
* `--libexec [package...]`: Update the self-provisioned toolchain. It always holds `micromamba`. Named
  packages are installed with it: `apptainer`, `squashfs-tools`, `squashfuse`.

By default (no flags), both `--build` and `--helper` are enabled.

**Features:**

* Prints each remote URL as it fetches metadata.
* Downloads and caches metadata locally per remote URL (default TTL: 7 days).
* Cached metadata is reused by `avail` and `create` without a network round-trip.
* Supports multiple recipe collections (`sources`); each remote gets its own cache file.
* Removes cache files for remotes no longer configured (orphan cleanup).
* `--libexec` refuses rather than waits if any condatainer session is currently using the
  toolchain — stop those sessions first, then retry.
* `--libexec` with no packages checks what the system already has (a usable apptainer, `mksquashfs`
  and `unsquashfs`, `squashfuse`), installs only what is missing (together with `micromamba` when there
  is no toolchain yet), and updates what is installed.
* `--libexec` prints the version of each installed tool after a successful update.

**Examples:**

```bash
# Refresh build + helper metadata (default)
condatainer update

# Build script metadata only
condatainer update --build

# Helper script metadata only
condatainer update --helper

# Update the installed toolchain
condatainer update --libexec

# Install packages into it
condatainer update --libexec apptainer squashfs-tools squashfuse
```

The default root image — `<distro>/base` (e.g. `ubuntu24/base`), where `<distro>` is `default_distro`
in config or a source's `default_distro` — is otherwise never
updated on its own: it is built the first time something needs it, the same as any other
named artifact, and reused until you rebuild it with `condatainer create --update <name>/<version>`.
Every other command treats it as a prerequisite — `create` builds it alongside the images
that need it, and `exec`/`e`/`run` build it before starting a container.

## Self-Update

Updates the **CondaTainer** binary to the latest version from GitHub releases.

**Usage:**

```
condatainer self-update [FLAGS]
```

**Options:**

* `-y`, `--yes`: Skip confirmation prompt and auto-update.
* `-f`, `--force`: Force update even if already on the latest version.
* `--dev`: Include pre-release versions.

**Features:**

* Downloads latest binary from GitHub releases.
* Detects current OS and architecture.
* Compares versions before updating.
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

## Proxy

SSH tunnel + dual-protocol (HTTP CONNECT + SOCKS5) proxy so compute node jobs can reach the internet through the login node.

The daemon automatically selects the best available SSH tunnel method:
1. **Go SSH library** — pure Go, uses hostbased auth (no user keys needed), SSH agent, or key files
2. **`ssh -D` Unix socket** — delegates to the system `ssh` binary
3. **`ssh -D` TCP port** — fallback for older OpenSSH (< 6.7)

```
condatainer proxy start|stop|status|show
```

### Two modes

| | Shared | Per-job |
|---|---|---|
| **Started on** | Login node | Compute node (inside job) |
| **Bind** | `0.0.0.0:PORT` | `127.0.0.1:PORT` |
| **Scope** | All compute nodes | This node only |
| **PID file** | `~/.local/share/condatainer/proxy.pid` (NFS) | `$SLURM_TMPDIR/cnt-$USER/proxy.pid` (local) |
| **Lifetime** | Until stopped / logout | Until job ends (tmpdir cleaned by scheduler) |
| **Start** | `condatainer proxy start` | `condatainer proxy start --via <login-node>` |

### proxy start

```
condatainer proxy start [--host HOST] [--via VIA] [--port PORT]
```

| Flag | Description |
|---|---|
| `--host HOST` | Run the daemon on `HOST` instead of the current node (login node only; delegates via SSH). |
| `--via VIA` | SSH server to tunnel through. Default: current node (shared mode); required for per-job mode inside a job. |
| `--port PORT` | Listen port. Default: OS-assigned. |

Shared mode attempts `loginctl enable-linger` so the daemon survives logout. A warning is printed if it fails.

```bash
condatainer proxy start                 # shared, OS picks port
condatainer proxy start --port 1080     # shared, fixed port
condatainer proxy start --host login02  # shared, delegate to login02
condatainer proxy start --via gateway   # shared, custom SSH gateway
condatainer proxy start --via login01   # per-job (inside a job)
```

### proxy stop / status / show / export

```bash
condatainer proxy stop      # kill shared proxy (SSH to other login node if needed)
condatainer proxy status    # check per-job then shared; clean stale PID files
condatainer proxy show      # print active URL (per-job first); empty if none
condatainer proxy export    # print export statements for all proxy env vars

# Export all proxy env vars (http_proxy, https_proxy, all_proxy, no_proxy)
# into the current shell; no-op when no proxy is active:
source <(condatainer proxy export)
```

### Automatic injection (inside jobs)

When running inside a scheduler job, condatainer auto-injects the proxy URL into:
- **container env** — `http_proxy`, `https_proxy`, `HTTP_PROXY`, `HTTPS_PROXY`, `all_proxy`, `ALL_PROXY`
- **condatainer's HTTP clients** — script fetching and metadata downloads use the proxy transport directly

Lookup order: per-job proxy → shared proxy.

The proxy port speaks both HTTP CONNECT (`http_proxy`) and SOCKS5 (`all_proxy`) — compatible with all tools (curl, wget, pip, micromamba, etc.).

### `proxy_perjob` config option

When `proxy_perjob: true`, condatainer injects `condatainer proxy start --via <login-node>` at the top of every generated scheduler script (build, run, helper jobs). The per-job proxy starts automatically on the compute node before the job body runs.

```yaml
proxy_perjob: true   # default: false
```

`CNT_PROXY_PERJOB=1` enables it for a single invocation without changing the config file.

## Completion

Generates shell completion scripts for CondaTainer.

**Usage:**

```
condatainer completion [SHELL]
```

* SHELL: Specify the shell type (`bash`, `zsh`, or `fish`). If not provided, the current shell is auto-detected. An unsupported shell is rejected.

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
