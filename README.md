# CondaTainer

<div>
<a href="#-condatainer"><img src="docs/_static/logo_cnt.svg" height="80" alt="CondaTainer Logo"></a>
</div>

[![Read the Docs](https://readthedocs.org/projects/condatainer/badge/?version=latest)](https://condatainer.readthedocs.io/en/latest/) [![Go Reference](https://pkg.go.dev/badge/github.com/Justype/condatainer.svg)](https://pkg.go.dev/github.com/Justype/condatainer) [![GitHub Release](https://img.shields.io/github/v/release/Justype/condatainer)](https://github.com/Justype/condatainer/releases)

**CondaTainer** is a rootless CLI designed to manage tools, data, and project environments, and seamlessly launch interactive apps on HPC clusters — perfect for individuals and small teams using institutional or regional compute.

* **Web-App Ready:** Launch *RStudio*, *VS Code*, *noVNC* and more with one command.
* **Unified Management:** Centralize group resources and easily isolate project environments.
* **Inode Saver:** Packing 30k+ Conda files into a single portable image.
* **Scheduler Native:** Out-of-the-box integration with *Slurm*, *PBS*, *LSF*, and *HTCondor*.

> [!NOTE]
> Slurm is the primary, tested scheduler. Others are experimental, bug reports are welcome!

## 🛠️ Installation

```bash
curl -fsSL https://get-condatainer.justype.net | bash
```

You will be prompted to confirm the installation path (defaults to `$SCRATCH/condatainer/` or `$HOME/condatainer/`). The script will also update shell config.

### ⚙️ Configuration

After installation, reload your shell and initialize the default configuration:

```bash
source ~/.bashrc # or source your shell config
ml apptainer # or ml singularity # if using module system
condatainer config init
```

This step ensures **CondaTainer** locates and saves Apptainer path for future use.

## 🧰 Tools/Data Management

```bash
condatainer avail # List available recipes
condatainer list  # List installed tools/data

# Create and run a specific tool
condatainer create samtools/1.16
condatainer exec -o samtools/1.16 samtools --version

# Print helpful info when running with overlays
condatainer exec -o grch38/cellranger/2024-A bash
# [CNT] Overlay envs:
#   CELLRANGER_REF_DIR: cellranger reference dir
```

**CondaTainer** will set `PATH` and other environment variables for you.

## 📦 Project Environment Management

**CondaTainer** supports both read-only (`.sqf`) and writable (`.img`) overlays, ensuring production stability alongside development flexibility.

### Read-Only Environments (Production)

Create an immutable bundle overlay directly from a Conda YAML file:

```bash
condatainer create -f environment.yml -p my_analysis
condatainer exec -o my_analysis.sqf bash
```

### Writable Environments (Development)

Manage a writable workspace overlay for development and testing:

```bash
condatainer overlay create -s 5G env.img   # Create a 5G overlay
condatainer overlay resize -s 10G env.img  # Resize an overlay to 10G
condatainer overlay chown --root env.img   # Make it compatible with --fakeroot

condatainer exec -w -o env.img bash  # Launch a shell in writable mode
condatainer e                        # Quick shortcut for the above command
```

Inside a writable mode, use `mm-*` Micromamba wrappers to manage packages

```bash
mm-install r-base=4.4 r-tidyverse  # Install packages
mm-pin r-base         # Pin a package version
mm-pin -r r-base      # Unpin a package
mm-remove r-tidyverse # Remove a package
mm-update             # Update packages
mm-export             # Export environment to YAML
```

## 🐕‍🦺 Web Apps & GUI Helpers

**CondaTainer** includes built-in helper scripts to launch apps (like RStudio, XFCE4 noVNC) directly on HPC compute nodes.

```bash
condatainer helper --update       # Fetch the latest helper scripts
condatainer helper --list         # View all available apps

condatainer helper vscode-tunnel  # Start a VS Code tunnel
condatainer helper igv -p 10212   # Start IGV via noVNC on port 10212
```

Please check out [ReadTheDocs - Helper Scripts](https://condatainer.readthedocs.io/en/dev/tutorials/helpers_on_HPC.html) for more details and examples.

## 🚀 Automation

**CondaTainer** can parse scripts to resolve dependencies (`#DEP:`) and scheduler directives (`#SBATCH`, `#PBS`, or `#BSUB`).

Example Script (`salmon_quant.sh`):

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
condatainer run -a salmon_quant.sh
```

### Job Chaining

All `[CNT]` messages go to stderr. Only job ID is printed to stdout, so you can capture it for downstream.

```bash
while read sample; do
  JOB=$(condatainer run -o log/trim_${sample}.out trim.sh $sample)
  JOB=$(condatainer run -o log/align_${sample}.out --afterok "$JOB" align.sh $sample)
  condatainer run -o log/quant_${sample}.out --afterok "$JOB" quant.sh $sample
done < samples.txt
```

## 🔗 Links & Resources

- [CondaTainer ReadTheDocs](https://condatainer.readthedocs.io/en/latest/)
- [cnt-scripts Repository](https://github.com/Justype/cnt-scripts)

Related tools and resources:

- [Compression Method Benchmarks](https://github.com/inikep/lzbench)
- Used tools: [Apptainer](https://github.com/apptainer/apptainer), [Micromamba](https://github.com/mamba-org/micromamba-releases), [squashfs-tools](https://github.com/plougher/squashfs-tools) [e2fsprogs](https://github.com/tytso/e2fsprogs)
- Container related: [Using Containers](https://services.rt.nyu.edu/docs/hpc/containers/containers/), [Singularity with Conda](https://services.rt.nyu.edu/docs/hpc/containers/singularity_with_conda/) and [Singularity with SquashFS](https://services.rt.nyu.edu/docs/hpc/containers/squash_file_system_and_singularity/)
- Allocation related: [Multi-Instance GPU](https://docs.alliancecan.ca/wiki/Multi-Instance_GPU), [RAM GPU ratio](https://docs.alliancecan.ca/wiki/Allocations_and_compute_scheduling#Ratios_in_bundles)

## Acknowledgements

This project used computational resources provided by

<div>
<a href="https://www.mcmaster.ca" target="_blank" rel="noopener noreferrer"><img src="assets/logo_mcmaster.svg" height="45" alt="McMaster University Logo"></a>
&nbsp;&nbsp;&nbsp;
<a href="https://alliancecan.ca" target="_blank" rel="noopener noreferrer"><img src="assets/logo_alliance.svg" height="40" alt="Digital Research Alliance of Canada Logo"></a>
</div>
