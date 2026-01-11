# ModGen Manual

🏭 **ModGen** is a wrapper script designed to streamline the management of Environment Modules or Lmod modules backed by Conda environments.

## Table of Contents

- [Initialization](#initialization)
- [Naming Convention](#naming-convention)
- [Create](#create)
- [Module Management (Avail, List, Remove)](#module-management-avail-list-remove)
- [Exec](#exec)
- [Runtime (Check, Run)](#runtime-check-run)
- [Update](#update)

## Overall Command Structure

```
usage: modgen [-h] [-v] [--debug] [-n] COMMAND ...

ModGen: Use conda and build scripts to create environment-modules or Lmod modules.

positional arguments:
  COMMAND               Available actions
    init                Initialize ModGen environment
    create (install, i)
                        Create a modules using conda or available build scripts
    avail (av)          Check available local and remote build scripts
    list (ls)           List installed modules matching search terms
    remove (delete)     Remove installed modules matching search terms
    exec (e)            Execute a command with modules loaded
    check (c)           Check if the dependencies of a script are installed
    run (r)             Run a script and auto-solve the dependencies by #DEP tags
    self-update         Update ModGen to the latest version
    condatainer         Install CondaTainer

optional arguments:
  -h, --help            show this help message and exit
  -v, --version         Show the version of ModGen
  --debug               Enable debug mode with verbose output
  -n, --no-sbatch, --local
                        Run all operations locally without using sbatch
```

## Initialization

Before using ModGen, you must initialize the shell environment. This ensures that the module system (Lmod or Environment Modules) knows where to look for the modules created by ModGen.

**Usage:**

```
modgen init [-m] [-c]
```

**Options**:

- `-m`, `--modules`: Add module use commands to your shell profile (e.g., `~/.bashrc`). This allows module load to find your ModGen libraries.
- `-c`, `--conda`: Initialize the internal Micromamba environment used to build packages. (not required for just module management)

## Naming Convention

**ModGen** organizes modules into two specific categories. This structure dictates where files are installed and how modulefiles are generated.

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

## Create

Create a new module. **ModGen** can build modules from existing build scripts (local/remote) or directly from Conda packages. It automatically generates the corresponding modulefile (Lua or Tcl) with appropriate PATH and environment variable settings.

**Usage:**

```
modgen create [OPTIONS] [NAME_VERSIONS...]
```

**Arguments**:

- `NAME_VERSIONS`: List of modules to create (e.g., `bcftools/1.22` or `grch38/genome/gencode`).

**Features**:

- Automatic Fetching: If a build script is not found locally, **ModGen** attempts to fetch it from the remote repository.
- Conda Fallback: If no build script exists, **ModGen** attempts to create the module by installing the package named name with version version from conda-forge or bioconda.
- Metadata Parsing: Parses `#ENV` and `#ENVNOTE` tags from build scripts to inject environment variables and help text into the generated modulefile.

## Module Management (Avail, List, Remove)

Manage the library of installed modules and available recipes.

### Avail

Search for available build scripts (both local scripts in build-scripts/ and remote metadata). Search terms are combined using **AND** logic (all terms must be present).

```
modgen avail [search_terms...] [-i|--install]
```

* `-i`, `--install`: Automatically install any found packages that are not currently installed.
* `-a`, `--add`: Alias for `--install`.

**Examples:**

```bash
# Search for all resources matching 'grcm' and 'M33'
$ modgen avail grcm M33
grcm39/gtf-gencode/M33
grcm39/salmon/1.10.2/gencodeM33
grcm39/star-2.7.11b/gencodeM33-101
grcm39/star-2.7.11b/gencodeM33-151
grcm39/transcript-gencode/M33

# Search for salmon using specific reference tags
$ modgen avail grcm M33 salmon
grcm39/salmon/1.10.2/gencodeM33

# Directly install from search results (will have prompt)
$ modgen avail grcm M33 salmon -i
grcm39/gtf-gencode/M33
[ModGen] Do you want to install the above modules? [y/N]:
```

### List

List currently installed modules generated by ModGen. (`apps`|`app-modules` and `refs`|`ref-modules` directories)

```
modgen list [-d] [search_terms...]
```

* `-d`, `--delete`: Delete the listed modules after confirmation.
* `-r`, `--remove`: Alias for `--delete`.

> [!TIP]
> You can also use Lmod to check available modules with `module avail` or `ml av`.

### Remove

Uninstall specific modules and remove their modulefiles. (`apps`|`app-modules` and `refs`|`ref-modules`)

```
modgen remove [search_terms...]
```

* Requires at least one search term.
* Prompts for confirmation before deletion.

## Exec

Execute a command in a clean environment with specific modules loaded. This is useful for one-off commands without modifying your current shell environment.

**Usage:**

```
modgen exec [OPTIONS] [COMMAND...]
```

**Options:**

* `-m`, `--module [MODULE]`: Specify modules to load (can be used multiple times).
* `-k`, `--keep`: Do not attempt to parse the command itself as a module name.

**Environment Variables:**

* `IN_CONDATINER=1`: Set when running the command.

**Behavior**: **ModGen** performs a module purge followed by module load for the requested modules before running the command.

**Example:**

```bash
# Run a command using a specific module
modgen exec -m bcftools/1.22 bcftools --version
```

> [!TIP]
> You can use module load directly in your shell for persistent module usage.

## Runtime (Check, Run)

Utilities for running scripts with automatic dependency handling via  `#DEP:` tags and `module load` or `ml` commands.

### Check

Parses a script for `#DEP:` tags and `module load` or `ml` commands and checks if the required overlays are installed.

```
modgen check [SCRIPT] [-a]
```

* `-a`, `--auto-install`: Automatically attempt to build/install missing dependencies found in the script.
* `-i`, `--install`: Alias for `--auto-install`.

### Run

Run a script with automatic dependency resolution based on `#DEP:` tags and `module load` or `ml` commands within the script.

```
modgen run [SCRIPT] [SCRIPT_ARGS...]
```

* **Dependency Injection:** Reads `#DEP:`and `module load` or `ml` commands lines in the script to determine which overlays to mount.
* **Argument Injection:** Reads `#CNT` lines in the script to inject arguments into the condatainer command itself.
* `-w`, `--writable-img`: Placeholder for CondaTainer compatibility (no effect).
* `-a`, `--auto-install`: Automatically attempt to build/install missing dependencies found in the script.
* `-i`, `--install`: Alias for `--auto-install`.

**Script Example:**

```bash
#!/bin/bash
#DEP: bcftools/1.22
#DEP: samtools/1.16
bcftools view input.vcf | head
```

## Update

Updates the **ModGen** script itself to the latest version from the GitHub repository.

**Usage:**

```
modgen self-update [-y]
```

* `-y`, `--yes`: Automatically confirm the download and replacement of the current script.
