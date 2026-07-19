# Build Script Manual

This document gives instructions on how to create your own build scripts for **CondaTainer**.

## Table of Contents

- [Naming Conventions](#naming-conventions)
- [Available Variables](#available-variables)
- [Headers](#headers)
  - [WhatIs and URL](#whatis-and-url)
  - [Set Dependencies](#set-dependencies)
  - [Scheduler Parameters](#scheduler-parameters)
  - [Type Tag](#type-tag)
  - [Template Tags](#template-tags)
  - [Auto-Update Tag](#auto-update-tag)
  - [Environment Variables](#environment-variables) and [ENV Naming Guidelines](#env-naming-guidelines)
  - [Interactive Tag](#interactive-tag)
- [Apps](#apps)
- [Data](#data)
- [OS](#os)

## Naming Conventions

The file path must follow the naming convention below to be recognized by CondaTainer:

`build_scripts/<name_conversion>` (should not include the .sh suffix)

Where `<name_conversion>` is defined as:

### OS

Apptainer definition files for distro-level system tools.

* **Format:** `<distro>/<name>`
* **Structure:**
  * **distro**: The base OS distribution (e.g., `ubuntu24`).
  * **name**: The tool name (e.g., `igv`, `r4.4.3`).
* **Example:** `ubuntu24/igv`

### Apps

Apps not available as conda packages, or specific versions not in conda.

**Single version per file:**

* **Format:** `<name>/<version>`
* **Example:** `cellranger/9.0.1`

**Template (multiple versions, one file):**

* **Format:** `<name>` (a single file with `#PL:` and `#TARGET:` headers)
* **Example:** `cytoscape` → expands to `cytoscape/3.10.3`, `cytoscape/3.10.4`, etc.
* Use this when the install logic is identical across versions and only the download URL changes.

### Data

Any data, including genome reference indexes.

* **Format:** `<assembly|project>/<datatype>/<version>`
* **Structure:**
  * **assembly/project**: The genome assembly or project name (e.g., `grch38`).
  * **datatype**: The type of data (e.g., `gtf-gencode`).
  * **version**: The release or build version (e.g., `47`).
* **Example:** `grch38/genome/gencode`

**Template (multiple versions, one file):**

* **Format:** `<assembly|project>/<datatype>` (a single file with `#PL:` and `#TARGET:` headers)
* **Example:** `grch38/star-gencode` → expands to `grch38/star/2.7.11b/gencode47-101`, `grch38/star/2.7.11b/gencode50-151`, etc.
* Use this when the install/create logic is identical across versions.

### Example

- [cellranger/9.0.1](https://github.com/Justype/cnt-scripts/blob/main/build-scripts/cellranger/9.0.1) (App)
- [cytoscape](https://github.com/Justype/cnt-scripts/blob/main/build-scripts/cytoscape) (App template)
- [grch38/cellranger/2024-A](https://github.com/Justype/cnt-scripts/blob/main/build-scripts/grch38/cellranger/2024-A) (Data)
- [grch38/star-gencode](https://github.com/Justype/cnt-scripts/blob/main/build-scripts/grch38/star-gencode) (Data templdate)

## Available Variables

| Variable | Description |
| -------- | ----------- |
| `NCPUS` | Number of CPUs (from script directives, or `build.ncpus` config setting) |
| `MEM` | Memory per task in MB (from script directives, or `build.mem` config setting) |
| `MEM_GB` | Memory per task in GB (integer) |
| `app_name` | app: name; data/ref: assembly/datatype |
| `version` | apps: app version; ref: data version |
| `target_dir` | Target installation directory (managed by **CondaTainer**) |
| `tmp_dir` | Temporary working directory (managed by **CondaTainer**) |

> `tmp_dir` paths:
>
> - For OS `.def`: intermediate `.sif` will created to condatainer root tmp, or next to the target dir.
> - For app module: use  scheduler tmp/`$TMPDIR`/`/tmp`
> - For data module: use condatainer root tmp
> - For external script: determined by `#TYPE:` tag, if app(default) use tmp, if data use target dir
> - if `CNT_TMPDIR` is set, it overrides all tmp behaviors.

| Function | Description |
| -------- | ----------- |
| `print_stderr` | Print message to stderr with current time |
| `pigz_or_gunzip` | Decompress gz file using pigz or gunzip |
| `tar_xf_pigz` | Extract tar.gz using pigz if available |
| `pigz_or_gunzip_pipe` | Decompress gz file and pipe to stdout using pigz or gunzip |

## Headers

Headers are special comments at the beginning of build scripts that provide metadata and instructions for CondaTainer.

**Example Header**: [grch38/bowtie2/ucsc_no_alt](https://github.com/Justype/cnt-scripts/blob/main/build-scripts/grch38/bowtie2/ucsc_no_alt)

```bash
#!/usr/bin/bash
#DEP:bowtie2/2.5.5>=2.3
#AUTOUPDATE:bowtie2:bioconda:bowtie2
#DEP:grch38/genome/ucsc_no_alt

#WHATIS:bowtie2 index for GRCh38 UCSC no alt reference genome
#URL:https://bowtie-bio.sourceforge.net/bowtie2/index.shtml

#ENV:BOWTIE2_PREFIX=$app_root/GRCh38_no_alt_analysis_set
#ENVNOTE:GRCh38 reference genome

#SBATCH --cpus-per-task=12
#SBATCH --mem=12G
#SBATCH --time=2:00:00
#SBATCH --job-name=bowtie2-build
#SBATCH --output=%x-%j.log

install() {
    ...
}
```

### WhatIs and URL

`#WHATIS:` is a one-line description of what the overlay provides; `#URL:` points to its homepage or documentation.

**CondaTainer** displays `#WHATIS:` in:

- `condatainer avail -w` (or `--whatis`) — shows the description beside each build script.
- `condatainer info <overlay>` — shows the description of a built overlay. This is read from the overlay's `.env` sidecar, which is written from `#WHATIS:` at build time.
- `condatainer create` — shown during interactive template resolution.

Both lines also replace the `{WHATIS}` and `{HELP}` placeholders in modulefiles generated by [ModGen](https://github.com/Justype/condatainer/blob/main/assets/modgen/manual.md).

### Set Dependencies

`#DEP:` lines specify dependencies that must be installed before building the current overlay.

When **CondaTainer** processes the build script, it will ensure that all specified dependencies are available and load them in the same order as listed.

**Basic (exact version):**

```bash
#DEP:samtools/1.21
```

Requires exactly `samtools/1.21`. If not installed, builds it.

**With version constraint:**

```bash
#DEP:samtools/1.22.1>=1.10
```

- **Preferred version**: `1.22.1` — built if no compatible version is installed.
- **Minimum version** (`>=1.10`): any installed version in the range `[1.10, 1.22.1]` is accepted.
- **Upper bound**: the preferred version acts as an implicit upper bound. A version higher than `1.22.1` (e.g., `2.0`) will not be used and the preferred version will be built instead.

`>` is also supported for a strict lower bound:

```bash
#DEP:samtools/1.22.1>1.10
```

**Version format**: partial versions are accepted — `1`, `1.10`, `1.22.1`. Missing components are treated as `0` (so `1.10` matches `1.10.0`, `1.10.3`, etc.).

When multiple installed versions satisfy the constraint, the latest one is used.

### Scheduler Parameters

Scheduler directive lines (`#SBATCH`, `#PBS`, or `#BSUB`) allow you to specify job parameters for the build process. (HTCondor uses native `.sub` submit files, so its directives is not supported.)

If a supported scheduler is available, **CondaTainer** will submit the build job with the specified parameters.

**Slurm Example:**

```bash
#SBATCH --cpus-per-task=16
#SBATCH --mem=42G
#SBATCH --time=2:00:00
#SBATCH --job-name=star-index
#SBATCH --output=%x-%j.log
```

- `--cpus-per-task`, `--mem`, and `--time` should be set according to the expected resource requirements.
- `--nodes`, `--ntasks`: **must not be set** (always single-task).
- `--output`: will always be overwritten to point to the `logs` directory.

**PBS Example:**

```bash
#PBS -l select=1:ncpus=16:mem=42gb
#PBS -l walltime=2:00:00
#PBS -N star-index
```

**LSF Example:**

```bash
#BSUB -n 16
#BSUB -R "span[hosts=1] rusage[mem=2688MB]"
#BSUB -W 2:00
#BSUB -J star-index
```

- `span[hosts=1]` is **required** — it keeps all 16 CPUs on a single node. Without it, `-n 16` means 16 MPI tasks and will be randomly distributed.
- `rusage[mem=]` is **per slot** (here, per CPU): `2688MB` × 16 = 42 GB total.

### Type Tag

`#TYPE:` only controls temporary build path behavior for **external** `.sh`/`.bash` builds.

- `#TYPE:app` (default): build in scratch tmp (`utils.GetTmpDir()`, affected by `CNT_TMPDIR`).
- `#TYPE:data`: build alongside the target prefix directory.

If `CNT_TMPDIR` is set, it overrides both `#TYPE:app` and `#TYPE:data` behavior and forces scratch tmp resolution under `CNT_TMPDIR/cnt-$USER`.

Accepted aliases (case-insensitive):

- App aliases: `app`, `env`, `tool`, `conda`, `small`
- Data aliases: `data`, `ref`, `large`

Examples:

```bash
#TYPE:app
#TYPE:data
#TYPE:ref
```

### Template Tags

`#PL:` and `#TARGET:` turn a single script into a parameterized **template**. Running `condatainer create <template-path>` prompts the user for each placeholder value and then creates the concrete overlay.

- `#PL:<name>:<values>` — declares a placeholder. The separator sets the ordering:
  - `,` — comma list, **sorted descending** (use for versions): `3.10,3.9,3.8`
  - `|` — pipe list, **author's order preserved** (use for labels): `kasm | turbo`

  In either form:
  - Integer ranges (`a-b`) expand to every integer between `a` and `b`, inclusive (`22-25` → `22,23,24,25`). Only pure digits are a range, so `2024-A` stays literal.
  - A `*` entry makes the list open-ended (any value accepted); it is always listed last.
  - Duplicates are dropped and whitespace around entries is ignored.
- `#TARGET:<pattern>` — the **module path** of the resulting overlay, once the placeholders are filled.
  Use `{name}` tokens matching the `#PL:` names.
  This is the name users refer to from then on — what they pass to `condatainer create`/`build`, what `#DEP:` lines point at, and the modulefile path [ModGen](https://github.com/Justype/condatainer/blob/main/assets/modgen/manual.md) generates. It is independent of the script's own file name.

Every `#PL:` name must appear as a `{name}` token in `#TARGET:`, and every token must have a matching `#PL:`. Otherwise **CondaTainer** warns and skips the expansion, because unused placeholders would silently collapse every value onto the same target.

**Example** — STAR index template [grch38/star-gencode](https://github.com/Justype/cnt-scripts/blob/main/build-scripts/grch38/star-gencode):

```bash
#!/usr/bin/bash
#PL:star_version:2.7.0b,...,2.7.11a,2.7.11b
#PL:gencode_version:22-49
#PL:read_length:101,151,*
#TARGET:grch38/star/{star_version}/gencode{gencode_version}-{read_length}

#DEP:grch38/genome/gencode
#DEP:grch38/gtf-gencode/{gencode_version}
#DEP:star/{star_version}

#WHATIS:STAR GRCh38 GENCODE{gencode_version} index for read length {read_length}
#URL:https://github.com/alexdobin/STAR/blob/master/doc/STARmanual.pdf

#ENV:STAR_INDEX_DIR=$app_root
#ENVNOTE:STAR index for GRCh38 GENCODE v{gencode_version} with read length {read_length}

#SBATCH --cpus-per-task=16
#SBATCH --mem=42G
#SBATCH --time=2:00:00
#SBATCH --job-name=star-index

install() {
    ...
}
```

When the user runs `condatainer create grch38/star-gencode`, **CondaTainer** shows the `#WHATIS:` description, the `#TARGET:` pattern (with `{placeholder}` tokens highlighted), then prompts for each placeholder in declaration order:

```
[CNT] Placeholder template: grch38/star-gencode
[CNT] STAR GRCh38 GENCODE{gencode_version} index for read length {read_length}
  Target: grch38/star/{star_version}/gencode{gencode_version}-{read_length}
  star_version [2.7.0b-2.7.11b] (default: 2.7.11b):
  gencode_version [22-49] (default: 49): 47
  read_length [suggested: 151, 101, or any value] (default: 151): 
  → Creating grch38/star/2.7.11b/gencode47-151
```

If the user already has a compatible dependency installed (e.g. `star/2.7.10` is installed), the default for `star_version` will be `2.7.10` instead of the latest available `2.7.11b`.

See [grch38/star-gencode](https://github.com/Justype/cnt-scripts/blob/main/build-scripts/grch38/star-gencode) for a real template example.

User can directly use the target to fill the placeholders:

```bash
condatainer create grch38/star/2.7.11b/gencode47-101
# Will directly set:
#   star_version=2.7.11b
#   gencode_version=47
#   read_length=101
```

### Auto-Update Tag

`#AUTOUPDATE:` opts a script into automatic version maintenance. A CI workflow runs twice every month, fetches the latest versions from the specified source, and rewrites the version list in place.

```
#AUTOUPDATE:{key}:{source}:{identifier}[>={min}][<{max}|<={max}]
```

`{key}` must match an existing `#PL:`, `#DEP:`, or (in helpers) `#VALUE:` header in the same file. The target type is detected automatically:

| Header matched | Behavior |
|---|---|
| `#PL:key:` | Rewrites the full version list (all versions ≥ min) |
| `#DEP:key/` | Rewrites the pinned version to latest only; preserves `>=constraint` |

**Supported sources:**

| Source | Format | Example |
|---|---|---|
| `github` | `github:{org}/{repo}` | `github:cytoscape/cytoscape>=3.9.0` |
| `bioconda` | `bioconda:{package}` | `bioconda:star>=2.7.0b` |
| `conda-forge` | `conda-forge:{package}` | `conda-forge:r-base>=4.0.0` |
| `docker` | `docker:{image}:{tag_regex}` | `docker:posit/r-base:^(\d+\.\d+\.\d+)-noble(?:-[^-]+)?$>=4.0.0` |

The `docker` source requires a full capture-group regex as the tag pattern — the first capture group is extracted as the version string.

**Examples:**

```bash
# PL template — full version list, all 3.9+
#PL:cytoscape_version:3.9.0,3.9.1,3.10.0,3.10.3,3.10.4
#AUTOUPDATE:cytoscape_version:github:cytoscape/cytoscape>=3.9.0
#TARGET:cytoscape/{cytoscape_version}

# DEP pin — latest samtools; min constraint preserved from the #DEP: line
#DEP:samtools/1.23.1>=1.10
#AUTOUPDATE:samtools:bioconda:samtools

# DEP pin with upper bound — stay on openjdk 17.x, never upgrade to 18+
#DEP:openjdk/17.0.12>=17
#AUTOUPDATE:openjdk:bioconda:openjdk>=17<18
```

Place `#AUTOUPDATE:` immediately after the `#PL:` or `#DEP:` line it manages.

### Environment Variables

- `#ENV:` lines define environment variables to be set when the overlay is loaded.
- `#ENVNOTE:` lines provide descriptions for the environment variables, which will be included in the modulefile help text and **CondaTainer** `.env` file.
  - `#ENVNOTE:` must directly follow its corresponding `#ENV:` line.
  - Only one `#ENVNOTE:` line is allowed per `#ENV:` variable.

`$app_root` is a special placeholder that will be replaced with the actual installation path of the overlay when loaded.

**Example:**

```bash
#ENV:CELLRANGER_REF_DIR=$app_root
#ENVNOTE:cellranger reference dir
#ENV:GENOME_FASTA=$app_root/fasta/genome.fa
#ENVNOTE:genome fasta
#ENV:ANNOTATION_GTF_GZ=$app_root/genes/genes.gtf.gz
#ENVNOTE:10X modified gtf
```

#### ENV Naming Guidelines

For common data: genome fasta, gtf, etc., use standard variable names like `GENOME_FASTA`, `ANNOTATION_GTF_GZ`.

- If the file is compressed, add `_GZ` suffix. e.g. `ANNOTATION_GTF` for uncompressed gtf, `ANNOTATION_GTF_GZ` for gzipped gtf.

For tool-specific references, use the tool name as a prefix.

- If the index is a directory, use `_DIR` suffix.
- If the index is a file prefix, use `_PREFIX` suffix.
- If the index is a specific file, use appropriate suffix based on file type.

**Examples:**

- `CELLRANGER_REF_DIR` for Cellranger references.
- `STAR_INDEX_DIR` for STAR indices.
- `BOWTIE2_PREFIX` for Bowtie2 indices.
- `BWA_MEM2_FASTA` for BWA-MEM2 genome fasta with `bwa-mem2` indices.

### Interactive Tag

- `#INTERACTIVE:<Prompt>` tag indicates that the build script requires input from the user during execution.
- It is common for apps that need license agreement acceptance or custom configuration.
- When **CondaTainer** encounters this tag, it will prompt the user with the specified `<Prompt>` message before building. It will take the user input and pass it to the build script during execution.
- You can use `\n` to add new lines in the prompt message.

**Example:**

```bash
#!/usr/bin/bash
#WHATIS:10X Genomics Single Cell Software Suite
#URL:https://www.10xgenomics.com/support/software/cell-ranger/downloads/previous-versions
#INTERACTIVE:⚠️ 10X links only valid for one day. Please go to the link below and get tar.gz link.\nhttps://www.10xgenomics.com/support/software/cell-ranger/downloads/previous-versions
```

Example: [cellranger/9.0.1](https://github.com/Justype/cnt-scripts/blob/main/build-scripts/cellranger/9.0.1)

## Apps

- Do not try to manually download apps that are already available via conda-forge or bioconda.
- Also, I don't recommend compiling apps from source unless absolutely necessary.
  - HPC systems often lack required build tools or dependencies unless you load specific modules.
  - To maximize compatibility (**CondaTainer**), it's better to rely on pre-compiled packages.

Template: [build-template-apps](https://github.com/Justype/cnt-scripts/blob/main/assets/build-template-apps)

### Tips

You can use `tar_xf_pigz` and `pigz_or_gunzip` functions to speed up decompression of large files if `pigz` is available on your system.

If the app requires specific environment variables to function properly, make sure to add them using `#ENV:` and `#ENVNOTE:` tags. e.g. [orad/2.7.0](https://github.com/Justype/cnt-scripts/blob/main/build-scripts/orad/2.7.0)

### Examples

- [cellranger/9.0.1](https://github.com/Justype/cnt-scripts/blob/main/build-scripts/cellranger/9.0.1)
- [orad/2.7.0](https://github.com/Justype/cnt-scripts/blob/main/build-scripts/orad/2.7.0)

## Data

- Data often require downloading large files from external sources.
- Indices may need to be built using specific versions of software.
  - If indices are version dependent, ensure the app version is included in the name. e.g. `grch38/star/2.7.11b/gencode47-101`
  - If indices require building, ensure you have the scheduler parameters (e.g. `#SBATCH`) set appropriately to allocate sufficient resources.
- Always add environment variables using `#ENV:` and `#ENVNOTE:` to help users locate and understand the reference data.

Template: [build-template-ref](https://github.com/Justype/cnt-scripts/blob/main/assets/build-template-ref)

### Tips

You can use `tar_xf_pigz` and `pigz_or_gunzip` functions to speed up decompression of large files if `pigz` is available on your system.

### Examples

- [grch38/genome/ucsc_no_alt](https://github.com/Justype/cnt-scripts/blob/main/build-scripts/grch38/genome/ucsc_no_alt)
- [grch38/transcript-gencode](https://github.com/Justype/cnt-scripts/blob/main/build-scripts/grch38/transcript-gencode) (template)
- [grch38/star-gencode](https://github.com/Justype/cnt-scripts/blob/main/build-scripts/grch38/star-gencode)
- [grcm39/salmon-gencode](https://github.com/Justype/cnt-scripts/blob/main/build-scripts/grcm39/salmon-gencode)

## OS

OS scripts are Apptainer definition files (`.def`) for distro-level system tools, built into a `.sif` image. The `.def` suffix is not part of the overlay name — `ubuntu24/igv.def` is referred to as `ubuntu24/igv`.

- Use these for tools that are not available as conda packages and need a full distro environment.
- Prefer an existing upstream container image (via `Bootstrap: docker`) over building from source.

### OS Templates

`.def` files support [Template Tags](#template-tags) too. Since a `.def` has no `install()` function, placeholders are substituted **throughout the whole file** at build time — including the `Bootstrap`/`From` header — so one file can cover every upstream image tag.

[ubuntu24/posit-r.def](https://github.com/Justype/cnt-scripts/blob/main/build-scripts/ubuntu24/posit-r.def) builds every R version from a single definition:

```
#PL:version:4.4.3,4.5.0,4.5.1
#AUTOUPDATE:version:docker:posit/r-base:^(\d+\.\d+\.\d+)-noble(?:-[^-]+)?$>=3.1.3

#TARGET:ubuntu24/r{version}
#WHATIS: Ubuntu noble R {version}

Bootstrap: docker
From: posit/r-base:{version}-noble
```

- `#TARGET:ubuntu24/r{version}` expands to `ubuntu24/r4.4.3`, `ubuntu24/r4.5.0`, …
- `{version}` in `From:` is substituted at build time, so each expansion pulls its own upstream tag.
- `#AUTOUPDATE:` can track a Docker tag to keep the `#PL:` list current. See [Auto-Update Tag](#auto-update-tag).

### Examples

- [ubuntu24/code-server.def](https://github.com/Justype/cnt-scripts/blob/main/build-scripts/ubuntu24/code-server.def)
- [ubuntu24/posit-r.def](https://github.com/Justype/cnt-scripts/blob/main/build-scripts/ubuntu24/posit-r.def) (template)
- [ubuntu24/xfce4.def](https://github.com/Justype/cnt-scripts/blob/main/build-scripts/ubuntu24/xfce4.def)
