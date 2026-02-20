# Build Script Manual

This document gives instructions on how to create your own build scripts for **CondaTainer**.

## Table of Contents

- [Naming Conventions](#naming-conventions)
- [Available Variables](#available-variables)
- [Headers](#headers)
  - [WhatIs and URL](#whatis-and-url)
  - [Set Dependencies](#set-dependencies)
  - [Scheduler Parameters](#scheduler-parameters)
  - [Environment Variables](#environment-variables) and [ENV Naming Guidelines](#env-naming-guidelines)
  - [Interactive Tag](#interactive-tag)
- [Apps](#apps)
- [References](#references)

## Naming Conventions

The file path must follow the naming convention below to be recognized by CondaTainer:

`build_scripts/<name_conversion>` (should not include the .sh suffix)

Where `<name_conversion>` is defined as:

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

### Example

- [cellranger/9.0.1](https://github.com/Justype/condatainer/blob/main/build-scripts/cellranger/9.0.1) (App)
- [grch38/cellranger/2024-A](https://github.com/Justype/condatainer/blob/main/build-scripts/grch38/cellranger/2024-A) (Ref)
- [grch38/star/2.7.11b/gencode47-101](https://github.com/Justype/condatainer/blob/main/build-scripts/grch38/star/2.7.11b/gencode47-101) (Ref)

## Available Variables

| Variable | Description |
| -------- | ----------- |
| `NCPUS` | Number of CPUs (from script directives, or `Build.DefaultCPUs`) |
| `app_name` | apps: app name; ref: assembly/data-type |
| `version` | apps: app version; ref: data version |
| `target_dir` | Target installation directory (managed by **CondaTainer**) |
| `tmp_dir` | Temporary working directory (managed by **CondaTainer**) |

| Function | Description |
| -------- | ----------- |
| `print_stderr` | Print message to stderr with current time |
| `pigz_or_gunzip` | Decompress gz file using pigz or gunzip |
| `tar_xf_pigz` | Extract tar.gz using pigz if available |
| `pigz_or_gunzip_pipe` | Decompress gz file and pipe to stdout using pigz or gunzip |

## Headers

Headers are special comments at the beginning of build scripts that provide metadata and instructions for CondaTainer.

**Example Header**: [star/2.7.11b/gencode47-101](https://github.com/Justype/condatainer/blob/main/build-scripts/grch38/star/2.7.11b/gencode47-101)

```bash
#!/usr/bin/bash
#DEP:grch38/genome/gencode
#DEP:grch38/gtf-gencode/22
#DEP:star/2.7.11b
#WHATIS:STAR grch38 GENCODE22 index for read length 101
#URL:https://github.com/alexdobin/STAR/blob/master/doc/STARmanual.pdf
#ENV:STAR_INDEX_DIR=$app_root
#ENVNOTE:STAR index for grch38 GENCODE v22 with read length 101

#SBATCH --cpus-per-task=16
#SBATCH --mem=42G
#SBATCH --time=2:00:00
#SBATCH --job-name=star-index
#SBATCH --output=%x-%j.log

install() {
    ...
}
```

### WhatIs and URL

`#WHATIS:` and `#URL:` lines are used to replace modulefile's `{WHATIS}` and `{HELP}` placeholders.

It will only be used by [ModGen](./modgen.md) when generating modulefiles. (**CondaTainer** does not support these tags yet.)

### Set Dependencies

`#DEP:` lines specify dependencies that must be installed before building the current overlay.

When **CondaTainer** processes the build script, it will ensure that all specified dependencies are available and load them in the same order as listed.

### Scheduler Parameters

Scheduler directive lines (`#SBATCH`, `#PBS`, or `#BSUB`) allow you to specify job parameters for the build process. HTCondor uses native `.sub` submit files instead of in-script directives.

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
- `--nodes`, `--ntasks`: should not be set (build jobs do not support MPI).
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
#BSUB -M 43008
#BSUB -W 2:00
#BSUB -J star-index
```

For LSF, **CondaTainer** will add `-R "span[hosts=1]"` to ensure all CPUs are allocated on the same node.

> **Note**: Certain SLURM flags are not supported and will cause the build to fail: `--topology-plugin`, `--switches`, `--gpus-per-socket`, `--sockets-per-node`, `--cores-per-socket`, `--threads-per-core`, `--ntasks-per-socket`, `--ntasks-per-core`, `--distribution`. Remove these from your build script.

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

Example: [cellranger/9.0.1](https://github.com/Justype/condatainer/blob/main/build-scripts/cellranger/9.0.1)

## Apps

- Do not try to manually download apps that are already available via conda-forge or bioconda.
- Also, I don't recommend compiling apps from source unless absolutely necessary.
  - HPC systems often lack required build tools or dependencies unless you load specific modules.
  - To maximize compatibility (**CondaTainer**), it's better to rely on pre-compiled packages.

Template: [build-template-apps](https://github.com/Justype/condatainer/blob/main/assets/build-template-apps)

### Tips

You can use `tar_xf_pigz` and `pigz_or_gunzip` functions to speed up decompression of large files if `pigz` is available on your system.

If the app requires specific environment variables to function properly, make sure to add them using `#ENV:` and `#ENVNOTE:` tags. e.g. [orad/2.7.0](https://github.com/Justype/condatainer/blob/main/build-scripts/orad/2.7.0)

### Examples

- [cellranger/9.0.1](https://github.com/Justype/condatainer/blob/main/build-scripts/cellranger/9.0.1)
- [orad/2.7.0](https://github.com/Justype/condatainer/blob/main/build-scripts/orad/2.7.0)

## References

- References often require downloading large files from external sources.
- indices may need to be built using specific versions of software.
  - If indices are version dependent, ensure the app version is included in the name. e.g. [grch38/star/2.7.11b/gencode47-101](https://github.com/Justype/condatainer/blob/main/build-scripts/grch38/star/2.7.11b/gencode47-101)
  - If indices require building, ensure you have the scheduler parameters (`#SBATCH`, `#PBS`, or `#BSUB`) set appropriately to allocate sufficient resources.
- Always add environment variables using `#ENV:` and `#ENVNOTE:` to help users locate the reference data.

Template: [build-template-ref](https://github.com/Justype/condatainer/blob/main/assets/build-template-ref)

### Tips

You can use `tar_xf_pigz` and `pigz_or_gunzip` functions to speed up decompression of large files if `pigz` is available on your system.

### Examples

- [grch38/genome/ucsc_no_alt](https://github.com/Justype/condatainer/blob/main/build-scripts/grch38/genome/ucsc_no_alt)
- [grch38/transcript-gencode/47](https://github.com/Justype/condatainer/blob/main/build-scripts/grch38/transcript-gencode/47)
- [grch38/star/2.7.11b/gencode47-101](https://github.com/Justype/condatainer/blob/main/build-scripts/grch38/star/2.7.11b/gencode47-101)
- [grcm39/salmon/1.10.2/gencodeM36](https://github.com/Justype/condatainer/blob/main/build-scripts/grcm39/salmon/1.10.2/gencodeM36)
