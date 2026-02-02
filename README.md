# CondaTainer

<div>
<a href="#-condatainer"><img src="docs/_static/logo_cnt.svg" height="80" alt="CondaTainer Logo"></a>
</div>

[![Read the Docs](https://readthedocs.org/projects/condatainer/badge/?version=latest)](https://condatainer.readthedocs.io/en/latest/) [![GitHub Release](https://img.shields.io/github/v/release/Justype/condatainer)](https://github.com/Justype/condatainer/releases)

**CondaTainer** is an HPC-oriented CLI that streamlines environment and data management by utilizing Apptainer, SquashFS, OverlayFS, and Micromamba.

* **Inode Saver:** Packing 30k+ Conda files into a single image to satisfy inode quotas.
* **Web-App Ready:** Out-of-the-box support for running *RStudio Server*, *code-server* and more on HPC.
* **Environment Isolation:** Supports read-only (`.sqf`) for production and writable (`.img`) for development.
* **Workload Manager Integration:** Native compatibility with Slurm for batch job submission.

## 🛠️ Installation

```bash
curl -fsSL https://get-condatainer.justype.net | bash
```

You will be prompted to confirm the installation path (defaulting to `$SCRATCH/condatainer/` or `$HOME/condatainer/`). The script will also edit shell config.

### ⚙️ Configuration

After installation, run the following command to create a default configuration file:

```bash
source ~/.bashrc # or source your shell config
module load apptainer # if apptainer is provided as a module
condatainer config init
```

This step will let **CondaTainer** save the apptainer path for future use.

## 👀 Quick Look

```bash
condatainer avail # List available recipes
condatainer list  # List installed

# Install and Run
condatainer create samtools/1.16
condatainer exec -o samtools/1.16 samtools --version

# Create a read-only project overlay from a YAML file
condatainer create -f environment.yml -p my_analysis
condatainer exec -o my_analysis.sqf bash

# Automatically install the dependencies
condatainer check analysis.sh -a

# Print helpful info when running with overlays
condatainer exec -o grch38/cellranger/2024-A bash
# [CondaTainer] Overlay envs:
#   CELLRANGER_REF_DIR: cellranger reference dir
#   GENOME_FASTA      : genome fasta
#   ANNOTATION_GTF_GZ : 10X modified gtf
#   STAR_INDEX_DIR    : STAR index dir
```

- 📜 [Read the full CondaTainer Manual](./docs/manuals/condatainer.md)
- 📁 [Naming Conventions](./docs/user_guide/concepts.md#-naming-convention)

## 📦 Writable Overlay

**CondaTainer** allows users to create writable overlays for development purposes.

```bash
# Create a writable overlay of size 5G
condatainer overlay create -s 5G env.img
# Resize an existing overlay to 10G
condatainer overlay resize -s 10G env.img
# Change the file UID/GID to current user in an overlay
# This is useful when the overlay was created by another user
condatainer overlay chown env.img

# Launch a bash shell inside the writable overlay
condatainer exec -w -o env.img bash
condatainer e # shortcut for above

# Run in read-only mode
condatainer exec -o env.img bash
```

**CondaTainer** will set `PATH` and other environment variables for you automatically.

When in writable mode, you can install packages using `mm-*` Micromamba wrappers. (conda-forge and bioconda only)

```bash
mm-install r-base=4.4 r-tidyverse
mm-pin r-base
mm-pin -r r-base # unpin
mm-list
mm-search r-ggplot2
mm-remove r-tidyverse
mm-update
mm-clean -a
mm-export
```

## 🚀 Automation

**CondaTainer** supports inline dependency declaration and automatic job submission. Define requirements with `#DEP:` tags and scheduler directives with `#SBATCH` or `#PBS`.

Example Script (`analysis.sh`):

```bash
#!/bin/bash
#SBATCH --time=4:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=30G
#DEP:salmon/1.10.2
#DEP:grch38/salmon/1.10.2/gencode47

salmon quant \
  -i $SALMON_INDEX_DIR \
  -p $SLURM_CPUS_PER_TASK \
  -l A -r reads.fq -o quants/
```

Auto install dependencies and submit the job with:

```bash
condatainer run analysis.sh -a
```

If no scheduler directives are found or job submission is disabled, the script will run immediately in the current shell.

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
