# CondaTainer & ModGen

<img src="assets/logos_gemini_generated.png" width="250" alt="CondaTainer and ModGen Logo">

An (HPC) toolkit for managing softwares and genome references efficiently.

This repository contains two complementary tools designed to simplify package management on HPC systems where users often face challenges with file quotas (inodes), complex dependency trees, and environment isolation.

- **CondaTainer**: Manages software using [Apptainer](https://apptainer.org/) ([Singularity](https://sylabs.io/docs/)), [Micromamba](https://mamba.readthedocs.io/en/latest/user_guide/micromamba.html) and [SquashFS](https://github.com/plougher/squashfs-tools) overlays. Ideal for saving file quota and creating isolated, reproducible containers.
- **ModGen**: Automates the creation of [Environment Modules](https://modules.readthedocs.io/en/latest/) (Tcl) or [Lmod](https://lmod.readthedocs.io/en/latest/) (Lua) modules. Ideal for users who prefer `module load`.

## 🛠️ Installation

By default, The script will be installed to `$SCRATCH/condatainer/` and add the installation path to your `.bashrc` or `.zshrc`.

If `$SCRATCH` is not defined, it will be installed to `$HOME/condatainer/`.

```bash
curl -sL https://raw.githubusercontent.com/Justype/condatainer/main/assets/install.sh | bash
```

## 📦 CondaTainer

**CondaTainer** solves the "too many files" problem inherent to Conda environments. Instead of creating thousands of small files, it packs Conda environments into single, highly compressed SquashFS files (`.sqf`) and mounts them inside an Apptainer container.

**Key Features**

- **Quota Saver**: Converts heavy Conda environments (1k+ files) into 1 SquashFS file.
- **Auto-Dependency Resolution**: Detects `#DEP:` tags and `module load` or `ml` commands in scripts and auto-mounts required overlays.
- **Writable Mode**: Supports writable ext3 images (.img) for environments that need runtime modification.
- **Smart Environment & Info**: Supports loading genome data and defining common environments variables. Helpful descriptions are automatically displayed when loading the module in an interactive shell (TTY).

**Quick Usage**

```bash
# 1. Create an overlay for samtools
condatainer create samtools/1.16

# 2. Run a command inside the container
condatainer exec -o samtools/1.16 samtools --version

# 3. List available recipes (local & remote)
condatainer avail

# 4. Create a custom environment from a YAML file
condatainer create -f environment.yml -p my_analysis

# 5. Cellranger reference example
bin/condatainer exec -o grch38--cellranger--2024-A bash
# [CondaTainer] Overlay envs:
#   CELLRANGER_REF_DIR: cellranger reference dir
#   GENOME_FASTA      : genome fasta
#   ANNOTATION_GTF_GZ : 10X modified gtf
#   STAR_INDEX_DIR    : STAR index dir
```

[Read the full CondaTainer Manual](assets/CNT_MANUAL.md)

## 🧩 ModGen

**ModGen** bridges the gap between Conda and HPC module systems. It installs packages via Conda into a centralized directory and automatically generates modulefiles (Lua or Tcl), allowing users to load them via standard module load commands.

**Key Features**

- **Native HPC Feel**: Uses module load and module avail.
- **Auto-Generation**: No need to write modulefiles manually; ModGen parses metadata and builds them for you.
- **Dependency Checking**: Scans scripts for dependencies and installs missing modules on the fly.
- **Smart Environment & Info**: Supports loading genome data and defining common environments variables. Helpful descriptions are automatically displayed when loading the module not in a SLURM job.

```bash
# 1. Initialize shell hooks (only needed once)
modgen init

# 2. Install bcftools and generate the modulefile
modgen create bcftools/1.16

# 3. Load the module (standard HPC way)
module load bcftools/1.16
bcftools --version

# 4. Search for packages to install
modgen avail salmon grcm M36

# 5. 10X reference example
ml grch38/cellranger/2024-A
# Available environment variables:
#   CELLRANGER_REF_DIR (cellranger reference dir)
#   GENOME_FASTA       (genome fasta)
#   ANNOTATION_GTF_GZ  (10X modified gtf)
#   STAR_INDEX_DIR     (STAR index dir)
```

[Read the full ModGen Manual](assets/MG_MANUAL.md)

## ⚖️ Which one should I use?

Feature | CondaTainer 🐳 | ModGen 🧩
-- | -- | --
Technology | Apptainer + SquashFS | Conda + Lmod/EnvModules
File Count (Inodes) | Very Low (1 file per env) | High (Standard Conda install)
Performance | High (LZ4 Compressed) | High
Isolation | Perfect (Containerized) | Shared (Host environment)
Best For | Saving quota | Rapid prototyping

## 🚀 Script Runtime Automation

Both tools support a "Headless" mode for running scripts. You can define dependencies directly in your bash scripts using metadata tags.

Example Script (`analysis.sh`):

```bash
#!/bin/bash
#DEP: salmon/1.10.3
#DEP: grcm39/salmon-1.10.3/gencodeM33

salmon quant -i $SALMON_INDEX_DIR ...
```

In shell or other scripts

```bash
# Run with CondaTainer
condatainer run analysis.sh

# Run with ModGen
modgen run analysis.sh
```

Integration with SLURM:

```bash
#!/bin/bash
#SBATCH --time=2:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=1GB
#DEP:samtools/1.22.1

if [ -z "$IN_CONDATINER" ]; then
    condatainer run "$0" "$@"
    exit $?
fi

samtools --version
```

## 📂 Naming Convention

Both tools utilize a standardized naming schema to organize software and reference data.

- Apps: `name/version` (e.g., `bcftools/1.16`)
- References: `assembly/datatype/version` (e.g., `grcm39/salmon-1.10.3/gencodeM33` salmon index for GRCm39 using gencode M33 transcripts)

## 🔗 Links & Resources

- [NYU HPC Using Containers on HPC](https://services.rt.nyu.edu/docs/hpc/containers/containers/)
- [NYU HPC Singularity with Conda](https://services.rt.nyu.edu/docs/hpc/containers/singularity_with_conda/)
- [NYU HPC SquashFS and Singularity](https://services.rt.nyu.edu/docs/hpc/containers/squash_file_system_and_singularity/)
- [Compression Method Benchmarks](https://github.com/inikep/lzbench)
- [AllianceCAN Multi-Instance GPU](https://docs.alliancecan.ca/wiki/Multi-Instance_GPU)
- [AllianceCAN CPU RAM GPU ratio](https://docs.alliancecan.ca/wiki/Allocations_and_compute_scheduling#Ratios_in_bundles)
