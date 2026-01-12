# CondaTainer & ModGen

[![Read the Docs](https://readthedocs.org/projects/condatainer/badge/?version=latest)](https://condatainer.readthedocs.io/en/latest/) [![GitHub Release](https://img.shields.io/github/v/release/Justype/condatainer)](https://github.com/Justype/condatainer/releases)

<div>
<a href="#-condatainer"><img src="assets/logo_cnt.png" height="100" alt="CondaTainer Logo"></a>
&nbsp;&nbsp;
<a href="./README_MG.md"><img src="assets/logo_mg.png" height="100" alt="ModGen Logo"></a>
</div>

Both tools aim to automate environment management on HPC, but with different approaches:

- **CondaTainer**: Wraps environments in overlays using [Apptainer](https://apptainer.org/) and [SquashFS](https://github.com/plougher/squashfs-tools).
- **ModGen**: Automates [Lmod](https://lmod.readthedocs.io/en/latest/)/[Environment Modules](https://modules.readthedocs.io/en/latest/) modules creation.

Which one to choose?

- If you care about inode usage and project-level isolation, use [**CondaTainer**](#-condatainer).
- If you want rstudio-server, igv, and other tools on HPC, use [**CondaTainer**](#-condatainer).
- If you prefer Lmod and using system available modules, use [**ModGen**](README_MG.md).

Please go to [Read the Docs](https://condatainer.readthedocs.io/en/latest/) for the full documentation, including installation guides, user manuals, and tutorials.

## 🛠️ Installation

```bash
curl -sL https://raw.githubusercontent.com/Justype/condatainer/main/assets/install.sh | bash
```

You will be prompted to confirm the installation path (defaulting to `$SCRATCH/condatainer/` or `$HOME/condatainer/`). The script will also edit shell config.

## 📦 CondaTainer

**CondaTainer** solves the "too many files" problem inherent to Conda environments: it packs Conda environments into single, highly-compressed SquashFS files and mounts them inside an Apptainer container.

- **Inode Efficiency**: Compresses heavy conda environments (30k+ files) into 1 file.
- **Smart Dependencies**: Scans scripts to detect and install missing modules. (`#DEP:` and `module load`)
- **Hybrid Modes**: Supports read-only (`.sqf`) for production and writable (`.img`) for development.
- **SLURM Integration**: Submits index generation jobs automatically when needed.

```bash
# 1. Create an overlay for samtools
condatainer create samtools/1.16

# 2. Run a command inside the container
condatainer exec -o samtools/1.16 samtools --version

# 3. List available recipes (local & remote)
condatainer avail

# 4. Create a read-only environment from a YAML file
condatainer create -f environment.yml -p my_analysis

# 5. Create a writable image for development
condatainer overlay create -s 10240 dev.img

# 6. Automatically install the dependencies
condatainer check analysis.sh -a

# 7. Helpful info example: Cellranger reference
condatainer exec -o grch38--cellranger--2024-A bash
# [CondaTainer] Overlay envs:
#   CELLRANGER_REF_DIR: cellranger reference dir
#   GENOME_FASTA      : genome fasta
#   ANNOTATION_GTF_GZ : 10X modified gtf
#   STAR_INDEX_DIR    : STAR index dir
```

- 📜 [Read the full CondaTainer Manual](./docs/manuals/condatainer.md)
- 📁 [Naming Conventions](./docs/user_guide/concepts.md#-naming-convention)

## 🚀 Automation

**CondaTainer** supports inline dependency declaration, allowing you to define requirements directly within your scripts using tags `#DEP:`.

Example Script (`analysis.sh`):

```bash
#!/bin/bash
#SBATCH --time=2:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=1GB
#DEP:samtools/1.22.1

if [ -z "$IN_CONDATINER" ] && command -v condatainer >/dev/null 2>&1; then
    condatainer run "$0" "$@"
    exit $?
fi

samtools --version
```

First check and install missing dependencies:

```bash
condatainer check analysis.sh -a
```

Then use sbatch as usual:

```bash
sbatch analysis.sh
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
