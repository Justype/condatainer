# CondaTainer Manual

**CondaTainer** is a wrapper script designed to streamline the management of Apptainer (Singularity) containers backed by Conda environments and SquashFS overlays.

## Table of Contents

- [Naming Convention](#naming-convention)
- [Mount Points](#mount-points)
- [Create](#create)
- [Container Management (Avail, List, Remove)](#container-management-avail-list-remove)
- [Exec](#exec)
- [Runtime (Check, Run)](#runtime-check-run)
- [Info](#info)
- [Apptainer](#apptainer)
- [Update](#update)

## Overall Command Structure

```
usage: condatainer [-h] [-v] [--debug] {create,avail,list,remove,exec,check,run,info,apptainer,update,help} ...

CondaTainer: Use apptainer/conda/squashFS to manage tools for HPC users.

positional arguments:
  {create,avail,list,remove,exec,check,run,info,apptainer,update,help}
                        Available actions
    create              Create a new SquashFS overlay using conda or available build scripts
    avail               Check available local and remote build scripts
    list                List installed overlays based on search terms
    remove              Remove installed overlays based on search terms
    exec                Execute a command using a group of overlays
    check               Check if the dependencies of a script are installed
    run                 Run a script and auto-solve the dependencies by #DEP tags
    info                Show information about a specific overlay
    apptainer           Get latest Apptainer executable from conda-forge
    update              Update CondaTainer to the latest version
    help                Show help information about CondaTainer

options:
  -h, --help            show this help message and exit
  -v, --version         Show the version of CondaTainer
  --debug               Enable debug mode with verbose output
```

## Naming Convention

CondaTainer classifies overlays into three distinct categories based on their naming structure.

### Custom Environments (Env)

Used when creating a custom environment with a user-defined name (e.g., via \-n or \-p).

* **Format:** `custom_name`
* **Constraints:**
  * Must **not** contain the double-dash sequence (`--`).
  * Should **not** conflict with common application names (to avoid ambiguity with app overlays).
* **Example:** `my_analysis_env`, `project_x_utils`

### Application Overlays (App)

Used for standard software packages and tools managed by CondaTainer build scripts.

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

## Mount Points

CondaTainer supports two types of overlay files, each with specific purposes and mount locations inside the container.

### SquashFS Overlays (.sqf)

* **Type:** Read-only, highly compressed.
* **Usage:** The default format for all Applications, References, and most Environments.
* **Mount Point:** `/cnt/[name]/[version]`
  * Files are accessible at the path defined by their naming convention.
  * *Example:* `bcftools/1.22` is mounted at `/cnt/bcftools/1.22`.

### Writable Images (.img)

* **Type:** Writable (ext3 filesystem), uncompressed.
* **Usage:** **Strictly for Conda environments** that require runtime modifications (e.g., installing new packages on the fly).
* **Mount Point:** `/ext3/[name]`
  * *Example:* `my_env.img` is mounted at `/ext3/my_env`.
* **Writability:** These are read-only by default. To enable writing, you must use the `-w` / `--writable-img` flag during exec or run.

## Create

Initialize and build a new CondaTainer SquashFS overlay. You can build from existing recipes (local/remote) or a Conda environment file.

**Usage:**

```
condatainer create [OPTIONS] [NAME_VERSIONS...]
```

**Options:**

* `-n`, `--name [NAME]`: Custom name for the overlay file. If used, all specified packages are bundled into one overlay.
* `-p`, `--prefix [PATH]`: Custom prefix path for the overlay file.
* `-f`, `--file [FILE]`: Path to a Conda environment file (.yml or .yaml).
* NAME_VERSIONS: List of packages to install (e.g., `bcftools/1.22` or `samtools=1.10` or `grch38/genome/gencode`).

**Compression Options:**

* `--zstd` / `--zstd-medium` / `--zstd-fast`: Use Zstandard compression (levels 14, 8, or 3).
* `--gzip`: Use Gzip compression.
* `--lz4`: Use LZ4 compression.

**Examples:**

```bash
# Create from a specific recipe (App Overlay)
condatainer create bcftools/1.22

# Create a reference overlay (Ref Overlay)
condatainer create grch38/gtf-gencode/47

# Create from a yaml file with a custom name (Custom Env)
condatainer create -f environment.yml -n my_analysis_env

# Create a custom env with multiple packages
condatainer create -n tools.sqf samtools=1.16 bcftools=1.15
```

**Features**:

- Automatic Fetching: If a build script is not found locally, ModGen attempts to fetch it from the remote repository.
- Conda Fallback: If no build script exists, ModGen attempts to create the module by installing the package named name with version version from conda-forge or bioconda.
- Metadata Parsing: Parses `#ENV` and `#ENVNOTE` tags from build scripts to inject environment variables and help text into the generated modulefile.

## Container Management (Avail, List, Remove)

Manage your local library of built containers and available recipes.

### Avail

Search for available build scripts (both local scripts in build-scripts/ and remote metadata). Search terms are combined using **AND** logic (all terms must be present).

```
condatainer avail [search_terms...] [-i|--install]
```

* `-i`, `--install`: Automatically install any found packages that are not currently installed.

**Examples:**

```bash
# Search for all resources matching 'grcm' and 'M33'
$ condatainer avail grcm M33
grcm39/gtf-gencode/M33
grcm39/salmon-1.10.2/gencodeM33
grcm39/star-2.7.11b/gencodeM33-101
grcm39/star-2.7.11b/gencodeM33-151
grcm39/transcript-gencode/M33

# Search for salmon using specific reference tags
$ condatainer avail grcm M33 salmon
grcm39/salmon-1.10.2/gencodeM33

# Directly install from search results (will have prompt)
$ condatainer avail grcm M33 salmon -i
grcm39/gtf-gencode/M33
[CondaTainer] Do you want to install the above overlays? [y/N]:
```

### List

Lists installed overlays found in the images/ and ref-images/ directories.

```
condatainer list [-d] [search_terms...]
```

* `-d`, `--delete`: Prompt to delete listed overlays after displaying them.

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
* `-w`, `--writable-img`: Mount .img overlays as writable (default is read-only). Only the last specified .img is made writable.
* `-k`, `--keep`: Do not attempt to parse the command itself as an overlay name.

**Environment Variables:**

* `IN_CONDATINER=1`: Set inside the container.
* `WRITABLE_PATH`: Paths to the writable directory if `-w`, `--writable-img` is used.

**Example:**

```bash
# Run a command using a specific overlay
condatainer exec -o bcftools/1.22 bcftools --version
```

## Runtime (Check, Run)

Utilities for running scripts with automatic dependency handling via \#DEP: tags.

### Check

Parses a script for `#DEP:` tags and `module load` or `ml` commands and checks if the required overlays are installed.

```
condatainer check [SCRIPT] [-a]
```

* `-a`, `--auto-install`: Automatically attempt to build/install missing dependencies found in the script.

### Run

Executes a script inside the CondaTainer environment, mounting dependencies defined in the script.

```
condatainer run [SCRIPT] [SCRIPT_ARGS...]
```

* **Dependency Injection:** Reads `#DEP:` lines in the script to determine which overlays to mount.
* **Argument Injection:** Reads `#CNT` lines in the script to inject arguments into the condatainer command itself.
* `-w`, `--writable-img`: Enable writable mounting for .img overlays.

**Script Example:**

```bash
#!/bin/bash
#DEP: bcftools/1.22
#CNT --writable-img
bcftools view input.vcf | head
```

## Info

Displays metadata regarding a specific overlay.

**Usage:**

```
condatainer info [OVERLAY]
```

**Output includes:**

* File size.
* Potential mount path (e.g., `/cnt/bcftools/1.22` or `/ext3/my_env`).
* Environment variables defined within the overlay's `.env` file.

## Apptainer

Manages the local installation of Apptainer. If Apptainer is not found on the host system, CondaTainer can install it via Micromamba into a local directory. But this is not recommended, because it cannot build singularity `sif` images.

**Usage:**

```
condatainer apptainer [-y] [-f]
```

**Options:**

* `-y`, `--yes`: Automatically confirm installation.
* `-f`, `--force`: Force re-installation even if it already exists.

## Update

Updates the CondaTainer script itself to the latest version from the GitHub repository.

**Usage:**

```
condatainer update [-y]
```

* `-y`, `--yes`: Automatically confirm the download and replacement of the current script.
