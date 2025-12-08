# CondaTainer & ModGen

<img src="assets/logos_gemini_generated.png" width="250" alt="CondaTainer and ModGen Logo">

A toolkit for efficient resource management on HPC systems, specifically targeting inode limits and environment isolation.

- **CondaTainer**: Wraps environments in Apptainer overlays using [Apptainer](https://apptainer.org/) and [SquashFS](https://github.com/plougher/squashfs-tools).
- **ModGen**: Automates [Lmod](https://lmod.readthedocs.io/en/latest/)/[Environment Modules](https://modules.readthedocs.io/en/latest/) modules creation.

If you want to reduce inode usage, use **CondaTainer**. If you prefer `Lmod` and using system available modules, use **ModGen**.

## 🛠️ Installation

```bash
curl -sL https://raw.githubusercontent.com/Justype/condatainer/main/assets/install.sh | bash
```

The installation script is **interactive**. You will be prompted to confirm the installation path (defaulting to `$SCRATCH/condatainer/` or `$HOME/condatainer/`). The script will also edit shell config.

## 📦 CondaTainer

**CondaTainer** solves the "too many files" problem inherent to Conda environments. Instead of creating thousands of small files. It packs Conda environments into single, highly compressed SquashFS files (`.sqf`) and mounts them inside an Apptainer container.

- **Inode Efficiency**: Compresses heavy conda environments (10k+ files) into 1 file.
- **Auto Generation**: Builds overlays from recipes or Conda environment files.
- **Smart Dependencies**: Scans scripts to detect and install missing modules on the fly. (by `#DEP:` tags and `module load`)
- **Hybrid Modes**: Supports read-only SquashFS (`.sqf`) for production and writable ext3 (`.img`) for development.
- **Context-Aware**: Manages genome references/indexes, and displays helpful info when in interactive sessions (TTY).

```bash
# 1. Create an overlay for samtools
condatainer create samtools/1.16

# 2. Run a command inside the container
condatainer exec -o samtools/1.16 samtools --version

# 3. List available recipes (local & remote)
condatainer avail

# 4. Create a custom environment from a YAML file
condatainer create -f environment.yml -p my_analysis

# 5. auto install the dependencies
condatainer check analysis.sh -a

# 6. Helpful info example: Cellranger reference
condatainer exec -o grch38--cellranger--2024-A bash
# [CondaTainer] Overlay envs:
#   CELLRANGER_REF_DIR: cellranger reference dir
#   GENOME_FASTA      : genome fasta
#   ANNOTATION_GTF_GZ : 10X modified gtf
#   STAR_INDEX_DIR    : STAR index dir
```

[Read the full CondaTainer Manual](assets/CNT_MANUAL.md)

## 🏭 ModGen

**ModGen** streamlines HPC software management by automatically converting apps and genome references into modules. It allows users to access Conda and non-Conda software (e.g., Cell Ranger, ORAD) but also genome references via standard `module` load commands without manual configuration.

- **Auto-Generation**: Installs packages and builds Lua/Tcl modulefiles automatically.
- **Smart Dependencies**: Scans scripts to detect and install missing modules on the fly. (by `#DEP:` tags and `module load`)
- **Native Integration**: seamless `module avail` and `module load` experience.
- **Context-Aware**: Manages genome references/indexes, and displays helpful info when not running within a SLURM job.

```bash
# 1. Initialize shell hooks (only needed once)
modgen init

# 2. Install samtools and generate the modulefile
modgen create samtools/1.16

# 3. Load the module (standard HPC way)
module load samtools/1.16
samtools --version

# 4. Search for packages to install
modgen avail salmon grcm M36

# 5. Auto install dependencies from a script
modgen check analysis.sh -a

# 6. Helpful info example: Cellranger reference
ml grch38/cellranger/2024-A
# Available environment variables:
#   CELLRANGER_REF_DIR (cellranger reference dir)
#   GENOME_FASTA       (genome fasta)
#   ANNOTATION_GTF_GZ  (10X modified gtf)
#   STAR_INDEX_DIR     (STAR index dir)
```

[Read the full ModGen Manual](assets/MG_MANUAL.md)

## ⚖️ Which one should I use?

Feature | CondaTainer 📦 | ModGen 🏭
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
#DEP: salmon/1.10.2
#DEP: grcm39/salmon-1.10.2/gencodeM33

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
- References: `assembly/datatype/version` (e.g., `grcm39/salmon-1.10.2/gencodeM33` salmon index for GRCm39 using gencode M33 transcripts)

## 🔗 Links & Resources

- [Compression Method Benchmarks](https://github.com/inikep/lzbench)
- [NYU HPC Using Containers on HPC](https://services.rt.nyu.edu/docs/hpc/containers/containers/)
- [NYU HPC Singularity with Conda](https://services.rt.nyu.edu/docs/hpc/containers/singularity_with_conda/)
- [NYU HPC SquashFS and Singularity](https://services.rt.nyu.edu/docs/hpc/containers/squash_file_system_and_singularity/)
- [AllianceCAN Multi-Instance GPU](https://docs.alliancecan.ca/wiki/Multi-Instance_GPU)
- [AllianceCAN CPU RAM GPU ratio](https://docs.alliancecan.ca/wiki/Allocations_and_compute_scheduling#Ratios_in_bundles)

## Acknowledgements

This project used computational resources provided by **McMaster University** and the **Digital Research Alliance of Canada**.

<img src="assets/logo_mcmaster.svg" height="50" alt="McMaster University Logo" style="vertical-align: top;">
<img src="assets/logo_alliance.svg" height="40" alt="Digital Research Alliance of Canada Logo" style="margin-left:20px; vertical-align: top;">
