# CondaTainer & ModGen

<img src="assets/logos_gemini_generated.png" width="250" alt="CondaTainer and ModGen Logo">

A toolkit for efficient resource management on HPC systems, targeting inode limits and environment isolation.

- **CondaTainer**: Wraps environments in overlays using [Apptainer](https://apptainer.org/) and [SquashFS](https://github.com/plougher/squashfs-tools).
- **ModGen**: Automates [Lmod](https://lmod.readthedocs.io/en/latest/)/[Environment Modules](https://modules.readthedocs.io/en/latest/) modules creation.

Which one to choose?

- If you want to reduce inode usage, use **CondaTainer**.
- If you prefer `Lmod` and using system available modules, use [**ModGen**](MG_README.md).

## 🛠️ Installation

```bash
curl -sL https://raw.githubusercontent.com/Justype/condatainer/main/assets/install.sh | bash
```

The installation script is **interactive**. You will be prompted to confirm the installation path (defaulting to `$SCRATCH/condatainer/` or `$HOME/condatainer/`). The script will also edit shell config.

## 📦 CondaTainer

**CondaTainer** solves the "too many files" problem inherent to Conda environments. Instead of creating thousands of small files. It packs Conda environments into single, highly compressed SquashFS files (`.sqf`) and mounts them inside an Apptainer container.

- **Inode Efficiency**: Compresses heavy conda environments (10k+ files) into 1 file.
- **Smart Dependencies**: Scans scripts to detect and install missing modules. (`#DEP:` and `module load`)
- **Module Compatibility**: Recognizes `module load` commands and loads corresponding overlays.
- **Hybrid Modes**: Supports read-only SquashFS (`.sqf`) for production and writable ext3 (`.img`) for development.
- **Context-Aware**: Manages references/indexes, and displays helpful info when in interactive sessions.
- **SLURM Integration**: Submits index generation jobs automatically when needed.

(Tips for writable images: see [here](assets/TIPS.md#create-writable-images))

```bash
# 1. Create an overlay for samtools
condatainer create samtools/1.16

# 2. Run a command inside the container
condatainer exec -o samtools/1.16 samtools --version

# 3. List available recipes (local & remote)
condatainer avail

# 4. Create a custom environment from a YAML file
condatainer create -f environment.yml -p my_analysis

# 5. Create a writable image
condatainer overlay -f environment.yml --size 10240 my_analysis_dev.img

# 6. auto install the dependencies
condatainer check analysis.sh -a

# 7. Helpful info example: Cellranger reference
condatainer exec -o grch38--cellranger--2024-A bash
# [CondaTainer] Overlay envs:
#   CELLRANGER_REF_DIR: cellranger reference dir
#   GENOME_FASTA      : genome fasta
#   ANNOTATION_GTF_GZ : 10X modified gtf
#   STAR_INDEX_DIR    : STAR index dir
```

[Read the full CondaTainer Manual](assets/CNT_MANUAL.md)

## 📂 Naming Convention

Both tools utilize a standardized naming schema to organize software and reference data.

- Apps: `name/version` (e.g., `bcftools/1.16`)
- References: `assembly/datatype/version`
  - `grcm39/genome/gencode`: mouse genome GRCm39 with Gencode style naming
  - `grcm39/salmon/1.10.2/gencodeM33`: salmon index for mouse genome GRCm39 with Gencode M33 annotation

> [!NOTE]
> CondaTainer uses overlay names to identify environments. So names must be unique and should not be changed after creation.

## 🚀 Script Runtime Automation

CondaTainer supports a "Headless" mode. You can define dependencies in scripts using tags.

Example Script (`analysis.sh`):

```bash
#!/bin/bash
#DEP: salmon/1.10.2
#DEP: grcm39/salmon/1.10.2/gencodeM33

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

## 🔗 Links & Resources

- [Compression Method Benchmarks](https://github.com/inikep/lzbench)
- Container related: [Using Containers](https://services.rt.nyu.edu/docs/hpc/containers/containers/), [Singularity with Conda](https://services.rt.nyu.edu/docs/hpc/containers/singularity_with_conda/) and [Singularity with SquashFS](https://services.rt.nyu.edu/docs/hpc/containers/squash_file_system_and_singularity/)
- Allocation related: [Multi-Instance GPU](https://docs.alliancecan.ca/wiki/Multi-Instance_GPU), [RAM GPU ratio](https://docs.alliancecan.ca/wiki/Allocations_and_compute_scheduling#Ratios_in_bundles)

## Acknowledgements

This project used computational resources provided by

<div>
<a href="https://www.mcmaster.ca" target="_blank" rel="noopener noreferrer"><img src="assets/logo_mcmaster.svg" height="45" alt="McMaster University Logo"></a>
&nbsp;&nbsp;&nbsp;
<a href="https://alliancecan.ca" target="_blank" rel="noopener noreferrer"><img src="assets/logo_alliance.svg" height="40" alt="Digital Research Alliance of Canada Logo"></a>
</div>
